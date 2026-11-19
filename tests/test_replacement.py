"""Original-file replacement (换件) workflow tests.

Two kinds of coverage:

* API/end-to-end tests through ``/api/batches/{id}/replacements`` -- the
  background validation thread is real, so these poll for terminal states.
* Pure logic tests for the group-rebuild rules, which are deterministic and
  easier to state exhaustively.
"""
from __future__ import annotations

import json
import time
from pathlib import Path

from tests.conftest import DEMO_DIR
from tests.test_api import upload_files, wait_for_terminal


def wait_for_replacements(client, batch_id, timeout=10.0):
    deadline = time.time() + timeout
    while time.time() < deadline:
        data = client.get(f"/api/batches/{batch_id}").json()
        if not data["replace_in_progress"]:
            return data
        time.sleep(0.05)
    raise AssertionError("replacement did not reach a terminal state")


def submit_replacement(client, batch_id, target_file_id, filename):
    handle = open(DEMO_DIR / filename, "rb")
    try:
        return client.post(
            f"/api/batches/{batch_id}/replacements",
            data={"target_file_id": target_file_id},
            files={"file": (filename, handle, "text/plain")},
        )
    finally:
        handle.close()


def file_by_name(data, name):
    return next(f for f in data["files"] if f["original_name"] == name)


def group_by_sample(data, sample):
    return next(g for g in data["dup_groups"] if g["sample_name"] == sample)


def choose(client, batch_id, group, file_name_part):
    file_id = next(
        c["file_id"] for c in group["choices"]
        if file_name_part in c["file_name"]
    )
    resp = client.put(
        f"/api/batches/{batch_id}/selection/{group['id']}",
        data={"file_id": file_id},
    )
    assert resp.status_code == 200, resp.text
    return file_id


# --------------------------------------------------------------------------- #
# End-to-end
# --------------------------------------------------------------------------- #
def test_replace_commits_and_updates_file_status_and_groups(client):
    data = upload_files(client, ["batch_a.vcf", "batch_b.vcf"]).json()
    batch_id = data["id"]
    target = file_by_name(data, "batch_b.vcf")

    # Select S2 from batch_a (an unchanged original) before the swap.
    s2 = group_by_sample(data, "S2")
    a_choice = choose(client, batch_id, s2, "batch_a")

    resp = submit_replacement(client, batch_id, target["id"], "replacement_b.vcf")
    assert resp.status_code == 200, resp.text
    final = wait_for_replacements(client, batch_id)

    replacement = final["replacements"][0]
    assert replacement["status"] == "committed", replacement["error"]
    assert "replacement_b.vcf" == replacement["original_name"]

    # Original updated in place: same file id, new name and contents.
    replaced = next(f for f in final["files"] if f["id"] == target["id"])
    assert replaced["original_name"] == "replacement_b.vcf"
    assert replaced["valid"] is True
    assert replaced["samples"] == ["S2", "S5"]
    assert replaced["record_count"] == 3

    # S3 vanished from the cohort; S5 is new and single-source (auto-selected);
    # S2 was multi-source before and its choice (batch_a, not replaced) stands.
    samples = {g["sample_name"]: g for g in final["dup_groups"]}
    assert "S3" not in samples
    assert samples["S5"]["chosen_file_id"] is not None
    assert samples["S5"]["is_duplicate"] is False
    assert samples["S2"]["chosen_file_id"] == a_choice
    assert final["unresolved_samples"] == []
    assert final["ready_to_generate"] is True
    assert final["batch_valid"] is True


def test_replace_clears_choice_when_selected_original_was_replaced(client):
    data = upload_files(client, ["batch_a.vcf", "batch_b.vcf"]).json()
    batch_id = data["id"]
    target = file_by_name(data, "batch_b.vcf")
    # S2 is selected from batch_b, the file being replaced.
    client.put(
        f"/api/batches/{batch_id}/selection/{group_by_sample(data, 'S2')['id']}",
        data={"file_id": target["id"]},
    )
    submit_replacement(client, batch_id, target["id"], "replacement_b.vcf")
    final = wait_for_replacements(client, batch_id)

    assert final["replacements"][0]["status"] == "committed"
    s2 = group_by_sample(final, "S2")
    assert s2["is_duplicate"] is True
    assert s2["chosen_file_id"] is None
    assert "S2" in final["unresolved_samples"]
    assert final["ready_to_generate"] is False


def test_replace_invalidates_completed_result(client):
    data = upload_files(client, ["batch_a.vcf", "batch_b.vcf"]).json()
    batch_id = data["id"]
    target = file_by_name(data, "batch_b.vcf")
    choose(client, batch_id, group_by_sample(data, "S2"), "batch_a")

    client.post(f"/api/batches/{batch_id}/generate")
    done = wait_for_terminal(client, batch_id)
    assert done["latest_generation"]["status"] == "done"
    old_output = Path(done["latest_generation"]["output_path"])
    assert old_output.is_file()
    version_before = done["selection_version"]

    submit_replacement(client, batch_id, target["id"], "replacement_b.vcf")
    final = wait_for_replacements(client, batch_id)

    assert final["replacements"][0]["status"] == "committed"
    assert final["latest_generation"]["status"] == "superseded"
    assert final["selection_version"] == version_before + 1
    assert not old_output.exists()
    assert client.get(f"/api/batches/{batch_id}/download").status_code == 409

    # A fresh generation under the replacement works.
    client.post(f"/api/batches/{batch_id}/generate")
    final2 = wait_for_terminal(client, batch_id)
    assert final2["latest_generation"]["status"] == "done"
    assert client.get(f"/api/batches/{batch_id}/download").status_code == 200


def test_failed_replacement_shows_reason_and_keeps_batch_usable(client):
    data = upload_files(client, ["batch_a.vcf", "batch_b.vcf"]).json()
    batch_id = data["id"]
    target = file_by_name(data, "batch_b.vcf")
    choose(client, batch_id, group_by_sample(data, "S2"), "batch_a")

    submit_replacement(client, batch_id, target["id"], "bad_replacement_b.vcf")
    final = wait_for_replacements(client, batch_id)

    rep = final["replacements"][0]
    assert rep["status"] == "failed"
    assert "换件后整批校验失败" in rep["error"]
    assert "REF 冲突" in rep["error"]

    # Old original untouched: file name did not change, its old samples and
    # validation status remain, and the batch is still usable.
    replaced = next(f for f in final["files"] if f["id"] == target["id"])
    assert replaced["original_name"] == "batch_b.vcf"
    assert replaced["samples"] == ["S2", "S3"]
    assert final["batch_valid"] is True
    assert final["selection_version"] == 2  # only the earlier manual choice
    assert final["ready_to_generate"] is True
    client.post(f"/api/batches/{batch_id}/generate")
    assert wait_for_terminal(client, batch_id)["latest_generation"]["status"] == "done"


def test_replacement_of_badly_formed_file_fails_per_file_validation(client):
    data = upload_files(client, ["batch_a.vcf", "batch_b.vcf"]).json()
    batch_id = data["id"]
    target = file_by_name(data, "batch_b.vcf")

    bad = ("broken.vcf", b"##fileformat=VCFv4.2\ngarbage\n", "text/plain")
    resp = client.post(
        f"/api/batches/{batch_id}/replacements",
        data={"target_file_id": target["id"]},
        files={"file": bad},
    )
    assert resp.status_code == 200, resp.text
    final = wait_for_replacements(client, batch_id)
    rep = final["replacements"][0]
    assert rep["status"] == "failed"
    assert "新原件未通过单文件校验" in rep["error"]


def test_replace_refused_when_another_original_is_invalid(client):
    from tests.conftest import DEMO_DIR

    # Three originals, one of which fails per-file validation.
    def upload_with_broken():
        with open(DEMO_DIR / "batch_a.vcf", "rb") as a, \
             open(DEMO_DIR / "batch_b.vcf", "rb") as b:
            return client.post(
                "/api/batches",
                files=[
                    ("files", ("batch_a.vcf", a, "text/plain")),
                    ("files", ("batch_b.vcf", b, "text/plain")),
                    ("files", ("broken.vcf", b"nonsense", "text/plain")),
                ],
            )

    # Replacing the broken original itself is allowed: every *other* original
    # already passed validation. The replacement is consistent, so it commits.
    data = upload_with_broken().json()
    batch_id = data["id"]
    broken = file_by_name(data, "broken.vcf")
    assert broken["valid"] is False
    resp = submit_replacement(client, batch_id, broken["id"], "replacement_b.vcf")
    assert resp.status_code == 200, resp.text
    final = wait_for_replacements(client, batch_id)
    assert final["replacements"][0]["status"] == "committed", final["replacements"][0]

    # Replacing a valid original while another original is still invalid is
    # refused synchronously, before any validation work starts.
    data2 = upload_with_broken().json()
    valid_target = file_by_name(data2, "batch_a.vcf")
    refused = submit_replacement(
        client, data2["id"], valid_target["id"], "replacement_b.vcf"
    )
    assert refused.status_code == 400, refused.text
    assert "其余原件尚未通过单文件校验" in refused.json()["detail"]
    # No replacement record was created for the refused attempt.
    assert client.get(f"/api/batches/{data2['id']}").json()["replacements"] == []


def test_selection_can_change_while_replacement_validates(client, app_env):
    from app import config

    data = upload_files(client, ["batch_a.vcf", "batch_b.vcf"]).json()
    batch_id = data["id"]
    target = file_by_name(data, "batch_b.vcf")
    s2 = group_by_sample(data, "S2")

    config.REPLACE_VALIDATE_DELAY_SECONDS = 0.6
    try:
        resp = submit_replacement(
            client, batch_id, target["id"], "replacement_b.vcf"
        )
        assert resp.status_code == 200, resp.text
        assert resp.json()["replace_in_progress"] is True
        # Generation is blocked while validating.
        assert client.post(
            f"/api/batches/{batch_id}/generate"
        ).status_code == 409
        # Selection is still accepted during validation.
        a_file = file_by_name(data, "batch_a.vcf")["id"]
        changed = client.put(
            f"/api/batches/{batch_id}/selection/{s2['id']}",
            data={"file_id": a_file},
        )
        assert changed.status_code == 200, changed.text
        final = wait_for_replacements(client, batch_id)
    finally:
        config.REPLACE_VALIDATE_DELAY_SECONDS = 0.0

    assert final["replacements"][0]["status"] == "committed"
    # The selection made during validation is honoured.
    assert group_by_sample(final, "S2")["chosen_file_id"] == a_file


def test_replacement_committing_first_rejects_the_earlier_unfinished_one(
    client, app_env, monkeypatch
):
    """Requirement: if another replacement commits while this one is still
    validating, this one must be rejected with a refresh hint."""
    import threading

    from app import db, service

    data = upload_files(client, ["batch_a.vcf", "batch_b.vcf", "batch_c.vcf"]).json()
    batch_id = data["id"]
    b_file = file_by_name(data, "batch_b.vcf")["id"]
    c_file = file_by_name(data, "batch_c.vcf")["id"]

    # Gate the *first* replacement's validation so the second one (submitted
    # afterwards) can validate and commit first.
    release = threading.Event()
    first_call = threading.Lock()
    real_validate = service.validate_vcf
    gate = {"used": False}

    def gated_validate(path):
        with first_call:
            should_block = not gate["used"]
            gate["used"] = True
        if should_block:
            assert release.wait(timeout=10)
        return real_validate(path)

    monkeypatch.setattr(service, "validate_vcf", gated_validate)

    first = submit_replacement(client, batch_id, c_file, "replacement_b.vcf")
    assert first.status_code == 200, first.text
    first_id = first.json()["replacement_id"]

    second = submit_replacement(client, batch_id, b_file, "replacement_b.vcf")
    assert second.status_code == 200, second.text
    second_id = second.json()["replacement_id"]

    def wait_status(replacement_id, terminal, timeout=10.0):
        deadline = time.time() + timeout
        row = None
        while time.time() < deadline:
            with db.get_db() as conn:
                row = db.get_replacement(conn, replacement_id)
            if row["status"] in terminal:
                return row
            time.sleep(0.02)
        raise AssertionError(f"replacement stuck in {row['status']}")

    from app import config

    committed = wait_status(second_id, config.REPLACE_TERMINAL_STATUSES)
    assert committed["status"] == "committed", committed["error"]

    # Let the earlier attempt finish: it now sees a bumped replace_seq.
    release.set()
    rejected = wait_status(first_id, config.REPLACE_TERMINAL_STATUSES)
    assert rejected["status"] == "rejected", rejected["error"]
    assert "刷新" in rejected["error"]

    # Only the winner's swap is in effect.
    final = client.get(f"/api/batches/{batch_id}").json()
    assert file_by_name(final, "replacement_b.vcf")["id"] == b_file
    assert file_by_name(final, "batch_c.vcf")["id"] == c_file
    assert final["replace_seq"] == 1


def test_restart_marks_validating_replacement_failed(client, app_env):
    from app import db as db_mod

    data = upload_files(client, ["batch_a.vcf", "batch_b.vcf"]).json()
    batch_id = data["id"]
    target = file_by_name(data, "batch_b.vcf")["id"]
    with db_mod.get_db() as conn:
        repl_id = db_mod.create_replacement(
            conn, batch_id, target, "x.vcf", "/nonexistent", 0
        )
    db_mod.init_db()
    with db_mod.get_db() as conn:
        row = db_mod.get_replacement(conn, repl_id)
    assert row["status"] == "failed"
    assert "重启" in row["error"]


# --------------------------------------------------------------------------- #
# Group rebuild rules in isolation
# --------------------------------------------------------------------------- #
def _add_files(conn, batch_id, specs):
    """specs: list of (name, valid, samples)."""
    from app import db as db_mod

    ids = []
    for name, valid, samples in specs:
        file_id = db_mod.create_file(conn, batch_id, name, f"/tmp/{name}")
        db_mod.update_file_validation(
            conn, file_id, valid=valid, error=None,
            reference_version="ref", sample_count=len(samples),
            record_count=1, samples=samples, contigs=["chr1"],
        )
        ids.append(file_id)
    return ids


def test_group_plan_rules(app_env):
    db = app_env["db"]
    with db.get_db() as conn:
        batch_id = db.create_batch(conn)
        a, b, c = _add_files(conn, batch_id, [
            ("a.vcf", True, ["X", "S"]),
            ("b.vcf", True, ["S"]),
            ("c.vcf", True, ["Y"]),
        ])
        # Before swap: S is multi-source with b chosen; X/Y single-source.
        s_group = db.create_dup_group(conn, batch_id, "S", [a, b])
        db.set_dup_group_choice(conn, s_group, b)
        x_group = db.create_dup_group(conn, batch_id, "X", [a])
        db.set_dup_group_choice(conn, x_group, a)
        y_group = db.create_dup_group(conn, batch_id, "Y", [c])
        db.set_dup_group_choice(conn, y_group, c)

        # Replace A: its samples become X,S. Now S sources are A(now), B.
        # The previously chosen B is untouched -> kept. X single -> auto.
        # If instead the chosen source was the replaced file, clear it.
        old_groups = db.get_dup_groups(conn, batch_id)
        plan = db.build_group_plan(
            conn, batch_id, old_groups=old_groups, target_file_id=a
        )
        by_sample = {sample: (ids, choice) for sample, ids, choice in plan}
        assert by_sample["S"][1] == b
        assert by_sample["X"][1] == a
        assert by_sample["Y"][1] == c


def test_group_plan_clears_when_replaced_file_was_selected(app_env):
    db = app_env["db"]
    with db.get_db() as conn:
        batch_id = db.create_batch(conn)
        a, b = _add_files(conn, batch_id, [
            ("a.vcf", True, ["S"]),
            ("b.vcf", True, ["S"]),
        ])
        s_group = db.create_dup_group(conn, batch_id, "S", [a, b])
        db.set_dup_group_choice(conn, s_group, a)
        old_groups = db.get_dup_groups(conn, batch_id)
        plan = db.build_group_plan(
            conn, batch_id, old_groups=old_groups, target_file_id=a
        )
        assert plan[0][2] is None


def test_group_plan_single_source_becoming_duplicate_needs_confirmation(app_env):
    """A sample that was previously single-source must not inherit that
    implicit source when the swap introduces a second one."""
    db = app_env["db"]
    with db.get_db() as conn:
        batch_id = db.create_batch(conn)
        a, b = _add_files(conn, batch_id, [
            ("a.vcf", True, ["S"]),
            ("b.vcf", True, []),
        ])
        group = db.create_dup_group(conn, batch_id, "S", [a])
        db.set_dup_group_choice(conn, group, a)
        # After the swap b also carries S: S is now multi-source. Simulate by
        # adding S to b's sample list before planning.
        conn.execute(
            "UPDATE vcf_files SET samples_json=? WHERE id=?",
            (json.dumps(["S"]), b),
        )
        old_groups = db.get_dup_groups(conn, batch_id)
        plan = db.build_group_plan(
            conn, batch_id, old_groups=old_groups, target_file_id=b
        )
        assert plan[0][2] is None


def test_group_plan_single_source_auto_selected(app_env):
    db = app_env["db"]
    with db.get_db() as conn:
        batch_id = db.create_batch(conn)
        (a,) = _add_files(conn, batch_id, [("a.vcf", True, ["ONLY"])])
        plan = db.build_group_plan(conn, batch_id, old_groups=[], target_file_id=a)
        assert plan == [("ONLY", [a], a)]
