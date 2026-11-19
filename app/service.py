"""Workflow orchestration: uploads, validation, selections and merge workers.

Selection versioning is the core integrity mechanism:

* every batch has a monotonically increasing ``selection_version``;
* every generation records the version it was started with;
* any choice change bumps the version and immediately marks every existing
* generation -- completed or still running -- as ``superseded``;
* a worker checks the live version between sites and aborts, never publishing
  an output belonging to an outdated choice, even if the user changes back.

Original replacement (lab reissues one corrected VCF) follows the same idea:
an attempt is validated in the background and commits atomically only if no
other replacement committed in the meantime. A commit swaps the original file,
re-validates the whole batch, rebuilds the sample groups and bumps
``selection_version`` so old outputs and in-flight merges are invalidated.
"""
from __future__ import annotations

import json
import logging
import os
import shutil
import threading
import time
from pathlib import Path

from . import config, cross_check, db
from .merge import MergeAborted, SampleSource, merge_vcfs
from .validation import validate_vcf

logger = logging.getLogger(__name__)


class SelectionError(ValueError):
    """Raised when a generation is requested but prerequisites are unmet."""


class ReplacementError(ValueError):
    """Raised when a replacement is rejected before validation starts."""


# --------------------------------------------------------------------------- #
# Upload + validation
# --------------------------------------------------------------------------- #
def _safe_stem(filename: str) -> tuple[str, str]:
    lower = filename.lower()
    if lower.endswith(".vcf.gz"):
        suffix = ".vcf.gz"
    else:
        suffix = Path(filename).suffix.lower()
    stem = Path(filename).name[: -len(suffix)] if suffix else Path(filename).stem
    safe = "".join(ch if ch.isalnum() or ch in "-_." else "_" for ch in stem)
    return safe, suffix


def save_upload(batch_id: str, index: int, filename: str, fileobj) -> tuple[str, str]:
    safe_stem, suffix = _safe_stem(filename)
    stored_name = f"{batch_id}_{index}_{safe_stem}{suffix}"
    stored_path = config.UPLOAD_DIR / stored_name
    config.UPLOAD_DIR.mkdir(parents=True, exist_ok=True)
    with open(stored_path, "wb") as handle:
        shutil.copyfileobj(fileobj, handle)
    return filename, str(stored_path)


def process_uploaded_batch(batch_id: str, saved: list[tuple[str, str]]) -> None:
    """Validate all files, build duplicate groups and run cross-file checks."""
    with db.get_db() as conn:
        file_rows: list = []
        for original_name, stored_path in saved:
            file_id = db.create_file(conn, batch_id, original_name, stored_path)
            _validate_and_store_file(conn, file_id, stored_path)
            file_rows.append(db.get_file(conn, file_id))

        valid_rows = [row for row in file_rows if row["valid"]]

        # Duplicate sample groups are derived from valid files only. Groups
        # with a single source are auto-selected: that original file is the
        # only possible choice. Multi-file groups start unresolved and block
        # generation until the operator picks one.
        plan = _group_plan(valid_rows)
        for sample, file_ids in plan:
            group_id = db.create_dup_group(conn, batch_id, sample, file_ids)
            if len(file_ids) == 1:
                db.set_dup_group_choice(conn, group_id, file_ids[0])

        reference_version, batch_errors = _cross_check_rows(file_rows)
        if reference_version:
            db.set_batch_reference(conn, batch_id, reference_version)
        if batch_errors:
            db.add_batch_errors(conn, batch_id, batch_errors)


def _validate_and_store_file(conn, file_id: str, stored_path: str) -> None:
    result = validate_vcf(stored_path)
    db.update_file_validation(
        conn,
        file_id,
        valid=result.valid,
        error=result.error,
        reference_version=result.reference_version,
        sample_count=len(result.samples),
        record_count=result.record_count,
        samples=result.samples,
        contigs=result.contigs,
    )


def _group_plan(valid_rows) -> list[tuple[str, list[str]]]:
    """Return ``(sample, file_ids)`` pairs sorted by sample name."""
    sample_to_files: dict[str, list[str]] = {}
    for row in valid_rows:
        for sample in json.loads(row["samples_json"]):
            sample_to_files.setdefault(sample, []).append(row["id"])
    return sorted(sample_to_files.items())


def _cross_check_rows(file_rows) -> tuple[str | None, list[str]]:
    """Batch-level checks over the valid files; returns (ref, errors)."""
    valid_rows = [row for row in file_rows if row["valid"]]
    batch_errors: list[str] = []
    reference_version: str | None = None
    if len(valid_rows) >= 2:
        check = cross_check.cross_check(
            [(row["original_name"], row["stored_path"]) for row in valid_rows]
        )
        reference_version = check.reference_version
        batch_errors.extend(check.errors)
    else:
        invalid_names = [
            row["original_name"] for row in file_rows if not row["valid"]
        ]
        batch_errors.append(
            "可参与合并的有效文件不足 2 份"
            + (f"，失效文件：{', '.join(invalid_names)}" if invalid_names else "")
        )
    return reference_version, batch_errors


# --------------------------------------------------------------------------- #
# Selection helpers
# --------------------------------------------------------------------------- #
def selection_summary(conn, batch_id: str) -> dict:
    groups = db.get_dup_groups(conn, batch_id)
    unresolved = [g["sample_name"] for g in groups if not g["chosen_file_id"]]
    return {
        "groups": groups,
        "unresolved": unresolved,
        "all_resolved": not unresolved,
    }


def assert_generation_allowed(conn, batch_id: str) -> None:
    batch = db.get_batch(conn, batch_id)
    if batch is None:
        raise SelectionError("批次不存在")
    if not batch["valid"]:
        errors = json.loads(batch["errors_json"])
        raise SelectionError(
            "批次存在阻断性错误，无法生成：" + "；".join(errors)
        )
    files = db.get_files(conn, batch_id)
    if sum(1 for f in files if f["valid"]) < config.MIN_FILES:
        raise SelectionError("有效 VCF 文件不足 2 份，无法生成")
    summary = selection_summary(conn, batch_id)
    if not summary["all_resolved"]:
        names = "、".join(summary["unresolved"])
        raise SelectionError(f"以下重复样本尚未选择唯一取用原件：{names}")
    if db.has_validating_replacement(conn, batch_id):
        raise SelectionError("换件校验进行中，请等待换件完成后再生成")
    active = db.get_active_generation(conn, batch_id)
    if active is not None:
        raise SelectionError("已有一个生成任务在进行中，请等待其结束")


def build_sources(conn, batch_id: str) -> list[SampleSource]:
    groups = db.get_dup_groups(conn, batch_id)
    sources: list[SampleSource] = []
    for group in groups:
        chosen = db.get_file(conn, group["chosen_file_id"])
        if chosen is None:
            raise SelectionError(
                f"样本 {group['sample_name']} 未选择取用原件"
            )
        sources.append(
            SampleSource(
                sample=group["sample_name"],
                file_path=chosen["stored_path"],
                file_name=chosen["original_name"],
            )
        )
    return sources


# --------------------------------------------------------------------------- #
# Background worker
# --------------------------------------------------------------------------- #
def _run_generation(gen_id: str, batch_id: str, version: int) -> None:
    with db.get_db() as conn:
        db.mark_generation(conn, gen_id, config.STATUS_RUNNING, started=True)
        try:
            sources = build_sources(conn, batch_id)
        except SelectionError as exc:
            db.mark_generation(
                conn, gen_id, config.STATUS_FAILED, error=str(exc), finished=True
            )
            return

    output_path = str(config.OUTPUT_DIR / f"cohort_{batch_id}_{gen_id}.vcf")

    def is_aborted() -> bool:
        with db.get_db() as check_conn:
            batch = db.get_batch(check_conn, batch_id)
            gen = db.get_generation(check_conn, gen_id)
            if batch is None or gen is None:
                return True
            if gen["status"] != config.STATUS_RUNNING:
                return True
            return batch["selection_version"] != version

    try:
        stats = merge_vcfs(sources, output_path, is_aborted=is_aborted)
    except MergeAborted:
        _safe_remove(output_path)
        with db.get_db() as conn:
            current = db.get_generation(conn, gen_id)
            if current is not None and current["status"] == config.STATUS_RUNNING:
                db.mark_generation(
                    conn,
                    gen_id,
                    config.STATUS_SUPERSEDED,
                    error="取用源选择已变更，该在途生成立即失效",
                    finished=True,
                )
        return
    except Exception as exc:  # noqa: BLE001 - failure reason must be surfaced
        _safe_remove(output_path)
        with db.get_db() as conn:
            db.mark_generation(
                conn,
                gen_id,
                config.STATUS_FAILED,
                error=f"合并失败：{exc}",
                finished=True,
            )
        return

    # Re-check before publishing: never publish a stale result.
    if is_aborted():
        _safe_remove(output_path)
        with db.get_db() as conn:
            current = db.get_generation(conn, gen_id)
            if current is not None and current["status"] == config.STATUS_RUNNING:
                db.mark_generation(
                    conn,
                    gen_id,
                    config.STATUS_SUPERSEDED,
                    error="取用源选择已变更，结果不予发布",
                    finished=True,
                )
        return

    logger.info(
        "generation %s merged %s sites for %s samples -> %s",
        gen_id, stats["sites"], stats["samples"], output_path,
    )
    with db.get_db() as conn:
        db.mark_generation(
            conn,
            gen_id,
            config.STATUS_DONE,
            output_path=output_path,
            error="",
            finished=True,
        )


def start_generation(batch_id: str) -> str:
    """Validate prerequisites synchronously, then spawn the worker thread."""
    with db.get_db() as conn:
        assert_generation_allowed(conn, batch_id)
        version = db.get_batch(conn, batch_id)["selection_version"]
        gen_id = db.create_generation(conn, batch_id, version)

    thread = threading.Thread(
        target=_run_generation,
        args=(gen_id, batch_id, version),
        name=f"merge-{gen_id[:8]}",
        daemon=True,
    )
    thread.start()
    return gen_id


def change_selection(batch_id: str, group_id: str, file_id: str | None) -> None:
    """Change a sample's source; invalidates all generations immediately.

    Allowed while a replacement validates: the swap recomputes groups inside
    its commit transaction and honours the latest selection.
    """
    with db.get_db() as conn:
        with db.transaction(conn):
            group = db.get_dup_group(conn, group_id)
            if group is None or group["batch_id"] != batch_id:
                raise SelectionError("重复样本分组不存在")
            updated = db.set_dup_group_choice(conn, group_id, file_id)
            if updated is None:
                raise SelectionError("选择无效")
            db._bump_selection_version_locked(conn, batch_id)

    # Remove output files of generations that were just superseded so no
    # half-product or stale result can ever be downloaded.
    _cleanup_invalid_outputs(batch_id)


def _cleanup_invalid_outputs(batch_id: str) -> None:
    with db.get_db() as conn:
        rows = conn.execute(
            "SELECT output_path FROM generations WHERE batch_id=? "
            "AND status=? AND output_path IS NOT NULL",
            (batch_id, config.STATUS_SUPERSEDED),
        ).fetchall()
    for row in rows:
        _safe_remove(row["output_path"])


def _safe_remove(path: str | None) -> None:
    if not path:
        return
    try:
        os.unlink(path)
    except FileNotFoundError:
        pass


# --------------------------------------------------------------------------- #
# Original-file replacement
# --------------------------------------------------------------------------- #
def _other_files_valid(conn, batch_id: str, target_file_id: str) -> str | None:
    """Every non-target original must already be valid, else replacement is
    refused (a replacement corrects one file, it does not repair a batch)."""
    broken = [
        row["original_name"]
        for row in db.get_files(conn, batch_id)
        if row["id"] != target_file_id and not row["valid"]
    ]
    if broken:
        return (
            "其余原件尚未通过单文件校验，不能换件："
            + "、".join(broken)
        )
    return None


def _stage_replacement_file(batch_id: str, filename: str, fileobj) -> str:
    safe_stem, suffix = _safe_stem(filename)
    name = f"{batch_id}_replace_{db.new_id()}_{safe_stem}{suffix}"
    config.TMP_DIR.mkdir(parents=True, exist_ok=True)
    path = config.TMP_DIR / name
    with open(path, "wb") as handle:
        shutil.copyfileobj(fileobj, handle)
    if path.stat().st_size > config.MAX_UPLOAD_BYTES:
        path.unlink(missing_ok=True)
        raise ReplacementError(f"文件 {filename} 超过大小限制")
    return str(path)


def submit_replacement(
    batch_id: str, target_file_id: str, filename: str, fileobj
) -> str:
    """Validate prerequisites synchronously, persist the staged file and
    return a replacement id whose validation runs in the background."""
    filename = (filename or "").strip()
    if not any(filename.lower().endswith(ext) for ext in config.ALLOWED_EXTENSIONS):
        raise ReplacementError(
            f"文件 {filename or '(未命名)'} 扩展名不支持，仅允许 .vcf / .vcf.gz"
        )
    if not target_file_id:
        raise ReplacementError("必须指定要替换的原件")

    with db.get_db() as conn:
        batch = db.get_batch(conn, batch_id)
        if batch is None:
            raise ReplacementError("批次不存在")
        target = db.get_file(conn, target_file_id)
        if target is None or target["batch_id"] != batch_id:
            raise ReplacementError("要替换的原件不属于该批次")
        reason = _other_files_valid(conn, batch_id, target_file_id)
        if reason:
            raise ReplacementError(reason)
        base_seq = batch["replace_seq"]

    staged_path = _stage_replacement_file(batch_id, filename, fileobj)
    try:
        with db.get_db() as conn:
            replacement_id = db.create_replacement(
                conn, batch_id, target_file_id, filename, staged_path, base_seq
            )
    except BaseException:
        _safe_remove(staged_path)
        raise

    thread = threading.Thread(
        target=_run_replacement,
        args=(replacement_id,),
        name=f"replace-{replacement_id[:8]}",
        daemon=True,
    )
    thread.start()
    return replacement_id


def _run_replacement(replacement_id: str) -> None:
    with db.get_db() as conn:
        row = db.get_replacement(conn, replacement_id)
        if row is None or row["status"] != config.REPLACE_VALIDATING:
            return

    delay = config.REPLACE_VALIDATE_DELAY_SECONDS
    if delay > 0:
        time.sleep(delay)

    try:
        new_validation = validate_vcf(row["temp_path"])
    except Exception as exc:  # noqa: BLE001 - reason must be surfaced
        _finish_replacement_failed(
            replacement_id, f"新原件校验时发生内部错误：{exc}"
        )
        return

    if not new_validation.valid:
        _finish_replacement_failed(
            replacement_id,
            f"新原件未通过单文件校验：{new_validation.error}",
        )
        return

    target_file_id = row["target_file_id"]
    batch_id = row["batch_id"]
    with db.get_db() as conn:
        batch = db.get_batch(conn, batch_id)
        target = db.get_file(conn, target_file_id)
        if batch is None or target is None:
            _finish_replacement_failed(replacement_id, "批次或原件已不存在")
            return
        require_valid = _other_files_valid(conn, batch_id, target_file_id)
        if require_valid:
            _finish_replacement_failed(replacement_id, require_valid)
            return

        # Whole-batch view with the new original substituted: every file is
        # individually valid here, so the cross-file checks are the only
        # remaining failure mode (reference mismatch, REF conflict, contig
        # order) and they must block the swap, not be committed.
        rows = []
        for other in db.get_files(conn, batch_id):
            if other["id"] == target_file_id:
                rows.append(
                    {
                        "id": target_file_id,
                        "original_name": row["original_name"],
                        "stored_path": row["temp_path"],
                        "valid": 1,
                        "samples_json": json.dumps(new_validation.samples),
                    }
                )
            else:
                rows.append(dict(other))
        reference_version, batch_errors = _cross_check_rows(rows)
        if batch_errors:
            _finish_replacement_failed(
                replacement_id,
                "换件后整批校验失败：" + "；".join(batch_errors),
            )
            return

        validation = dict(
            valid=True,
            error=None,
            reference_version=new_validation.reference_version,
            sample_count=len(new_validation.samples),
            record_count=new_validation.record_count,
            samples=new_validation.samples,
            contigs=new_validation.contigs,
        )
        stored_name = _replacement_stored_name(batch_id, row["original_name"])

    staged = Path(row["temp_path"])
    final_path = config.UPLOAD_DIR / stored_name
    try:
        os.replace(staged, final_path)
    except OSError as exc:
        _finish_replacement_failed(replacement_id, f"换件文件落盘失败：{exc}")
        return

    with db.get_db() as conn:
        committed = db.commit_replacement(
            conn,
            batch_id,
            replacement_id,
            target_file_id,
            row["original_name"],
            str(final_path),
            row["base_seq"],
            validation,
            reference_version,
            batch_errors,
        )
    if not committed:
        # A later replacement won; discard ours and make the UI explain.
        _safe_remove(str(final_path))
        with db.get_db() as conn:
            db.mark_replacement(
                conn,
                replacement_id,
                config.REPLACE_REJECTED,
                error=(
                    "换件校验期间已有另一次换件先提交，本次换件已被拒绝，"
                    "请刷新页面后重试"
                ),
            )
        _cleanup_invalid_outputs(batch_id)
        return

    _cleanup_invalid_outputs(batch_id)
    logger.info(
        "replacement %s committed for batch %s (file %s)",
        replacement_id, batch_id, target_file_id,
    )


def _replacement_stored_name(batch_id: str, new_name: str) -> str:
    safe_stem, suffix = _safe_stem(new_name)
    return f"{batch_id}_repl_{db.new_id()}_{safe_stem}{suffix}"


def _finish_replacement_failed(replacement_id: str, error: str) -> None:
    with db.get_db() as conn:
        row = db.get_replacement(conn, replacement_id)
        if row is None or row["status"] != config.REPLACE_VALIDATING:
            return
        db.mark_replacement(
            conn, replacement_id, config.REPLACE_FAILED, error=error
        )
        staged = row["temp_path"]
    # The batch keeps using the old original untouched.
    _safe_remove(staged)
