"""SQLite persistence layer.

Schema overview
---------------
batches            : one upload transaction (2-3 VCF files)
vcf_files          : uploaded files and their per-file validation outcome
dup_groups         : groups of files sharing a sample; exactly one must be
                     chosen as the unique source of that sample
generations        : merge runs, tagged with the ``selection_version`` they
                     were started with; any later selection change invalidates
                     them immediately
replacements       : attempts to replace one original with a corrected VCF
                     reissued by the lab. Each attempt is validated in the
                     background and committed atomically; a commit bumps
                     ``replace_seq`` so a concurrent attempt that started
                     earlier is rejected instead of overwriting it.

The schema deliberately keeps no per-locus data: VCF contents stay on disk
and the database only tracks workflow state.
"""
from __future__ import annotations

import json
import sqlite3
import threading
import uuid
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterator

from . import config

_local = threading.local()
_init_lock = threading.Lock()


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def new_id() -> str:
    return uuid.uuid4().hex


def _connect(db_path: Path | None = None) -> sqlite3.Connection:
    conn = sqlite3.connect(
        str(db_path or config.DB_PATH),
        timeout=30,
        check_same_thread=False,
        isolation_level=None,  # autocommit; explicit BEGIN for transactions
    )
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA foreign_keys=ON")
    conn.execute("PRAGMA busy_timeout=30000")
    conn.execute("PRAGMA synchronous=NORMAL")
    return conn


@contextmanager
def get_db() -> Iterator[sqlite3.Connection]:
    """Thread-local connection used by request handlers and workers."""
    conn = getattr(_local, "conn", None)
    if conn is None:
        conn = _connect()
        _local.conn = conn
    try:
        yield conn
    except sqlite3.Error:
        conn.rollback()
        raise


SCHEMA = """
CREATE TABLE IF NOT EXISTS batches (
    id                 TEXT PRIMARY KEY,
    created_at         TEXT NOT NULL,
    reference_version  TEXT,
    valid              INTEGER NOT NULL DEFAULT 1,
    -- Batch-level blocking reasons (reference mismatch, REF conflicts ...)
    errors_json        TEXT NOT NULL DEFAULT '[]',
    selection_version  INTEGER NOT NULL DEFAULT 1,
    -- Incremented on every committed original replacement. A replacement
    -- attempt records the value it observed when submitted; committing is
    -- only allowed while the value is unchanged (rejects stale attempts).
    replace_seq        INTEGER NOT NULL DEFAULT 0
);

CREATE TABLE IF NOT EXISTS vcf_files (
    id                 TEXT PRIMARY KEY,
    batch_id           TEXT NOT NULL REFERENCES batches(id) ON DELETE CASCADE,
    original_name      TEXT NOT NULL,
    stored_path        TEXT NOT NULL,
    uploaded_at        TEXT NOT NULL,
    valid              INTEGER NOT NULL DEFAULT 1,
    error              TEXT,
    reference_version  TEXT,
    sample_count       INTEGER NOT NULL DEFAULT 0,
    record_count       INTEGER NOT NULL DEFAULT 0,
    samples_json       TEXT NOT NULL DEFAULT '[]',
    contigs_json       TEXT NOT NULL DEFAULT '[]'
);

CREATE TABLE IF NOT EXISTS dup_groups (
    id                 TEXT PRIMARY KEY,
    batch_id           TEXT NOT NULL REFERENCES batches(id) ON DELETE CASCADE,
    sample_name        TEXT NOT NULL,
    -- Files that carry this sample (vcf_files id list, JSON array)
    file_ids_json      TEXT NOT NULL,
    -- Chosen source file id; NULL means no selection yet (blocks generation)
    chosen_file_id     TEXT
);

CREATE TABLE IF NOT EXISTS generations (
    id                 TEXT PRIMARY KEY,
    batch_id           TEXT NOT NULL REFERENCES batches(id) ON DELETE CASCADE,
    status             TEXT NOT NULL,
    selection_version  INTEGER NOT NULL,
    output_path        TEXT,
    error              TEXT,
    created_at         TEXT NOT NULL,
    started_at         TEXT,
    finished_at        TEXT
);

CREATE INDEX IF NOT EXISTS idx_files_batch ON vcf_files(batch_id);
CREATE INDEX IF NOT EXISTS idx_dup_batch ON dup_groups(batch_id);
CREATE INDEX IF NOT EXISTS idx_gen_batch ON generations(batch_id);

CREATE TABLE IF NOT EXISTS replacements (
    id                 TEXT PRIMARY KEY,
    batch_id           TEXT NOT NULL REFERENCES batches(id) ON DELETE CASCADE,
    target_file_id     TEXT NOT NULL,
    original_name      TEXT NOT NULL,
    temp_path          TEXT NOT NULL,
    status             TEXT NOT NULL,
    error              TEXT,
    -- Replace seq observed at submission time.
    base_seq           INTEGER NOT NULL,
    created_at         TEXT NOT NULL,
    finished_at        TEXT
);

CREATE INDEX IF NOT EXISTS idx_repl_batch ON replacements(batch_id);
"""


def init_db() -> None:
    with _init_lock:
        config.ensure_dirs()
        conn = _connect()
        try:
            conn.executescript(SCHEMA)
            _migrate(conn)
            # A replacement that was still validating when the process died can
            # never commit; mark it failed so it is never mistaken for success.
            conn.execute(
                "UPDATE replacements SET status=?, finished_at=?, "
                "error=COALESCE(error, ?) WHERE status=?",
                (
                    config.REPLACE_FAILED,
                    _now(),
                    "服务在换件校验完成前重启，该次换件判定失败",
                    config.REPLACE_VALIDATING,
                ),
            )
            conn.execute(
                "UPDATE generations SET status=?, finished_at=?, "
                "error=COALESCE(error, ?) WHERE status IN (?, ?)",
                (
                    config.STATUS_FAILED,
                    _now(),
                    "服务在合并完成前重启，该在途生成判定失败",
                    config.STATUS_PENDING,
                    config.STATUS_RUNNING,
                ),
            )
        finally:
            conn.close()


def _migrate(conn: sqlite3.Connection) -> None:
    """Additive migrations for databases created before replace_seq existed."""
    columns = {
        row["name"]
        for row in conn.execute("PRAGMA table_info(batches)").fetchall()
    }
    if "replace_seq" not in columns:
        conn.execute(
            "ALTER TABLE batches ADD COLUMN replace_seq INTEGER NOT NULL DEFAULT 0"
        )


# --------------------------------------------------------------------------- #
# Reads
# --------------------------------------------------------------------------- #
def get_batch(conn: sqlite3.Connection, batch_id: str) -> sqlite3.Row | None:
    return conn.execute(
        "SELECT * FROM batches WHERE id=?", (batch_id,)
    ).fetchone()


def get_files(conn: sqlite3.Connection, batch_id: str) -> list[sqlite3.Row]:
    return conn.execute(
        "SELECT * FROM vcf_files WHERE batch_id=? ORDER BY uploaded_at, id",
        (batch_id,),
    ).fetchall()


def get_file(conn: sqlite3.Connection, file_id: str) -> sqlite3.Row | None:
    return conn.execute(
        "SELECT * FROM vcf_files WHERE id=?", (file_id,)
    ).fetchone()


def get_dup_groups(conn: sqlite3.Connection, batch_id: str) -> list[sqlite3.Row]:
    return conn.execute(
        "SELECT * FROM dup_groups WHERE batch_id=? ORDER BY sample_name",
        (batch_id,),
    ).fetchall()


def get_dup_group(conn: sqlite3.Connection, group_id: str) -> sqlite3.Row | None:
    return conn.execute(
        "SELECT * FROM dup_groups WHERE id=?", (group_id,)
    ).fetchone()


def get_active_generation(
    conn: sqlite3.Connection, batch_id: str
) -> sqlite3.Row | None:
    return conn.execute(
        "SELECT * FROM generations WHERE batch_id=? AND status IN (?, ?) "
        "ORDER BY created_at DESC LIMIT 1",
        (batch_id, config.STATUS_PENDING, config.STATUS_RUNNING),
    ).fetchone()


def get_latest_generation(
    conn: sqlite3.Connection, batch_id: str
) -> sqlite3.Row | None:
    # rowid reflects insertion order; created_at has only second-level
    # resolution and several generations can share one timestamp.
    return conn.execute(
        "SELECT * FROM generations WHERE batch_id=? ORDER BY rowid DESC "
        "LIMIT 1",
        (batch_id,),
    ).fetchone()


def get_generation(conn: sqlite3.Connection, gen_id: str) -> sqlite3.Row | None:
    return conn.execute(
        "SELECT * FROM generations WHERE id=?", (gen_id,)
    ).fetchone()


def get_running_generations(conn: sqlite3.Connection) -> list[sqlite3.Row]:
    return conn.execute(
        "SELECT * FROM generations WHERE status IN (?, ?)",
        (config.STATUS_PENDING, config.STATUS_RUNNING),
    ).fetchall()


# --------------------------------------------------------------------------- #
# Writes
# --------------------------------------------------------------------------- #
def create_batch(conn: sqlite3.Connection) -> str:
    batch_id = new_id()
    conn.execute(
        "INSERT INTO batches(id, created_at) VALUES(?, ?)",
        (batch_id, _now()),
    )
    return batch_id


def create_file(
    conn: sqlite3.Connection,
    batch_id: str,
    original_name: str,
    stored_path: str,
) -> str:
    file_id = new_id()
    conn.execute(
        "INSERT INTO vcf_files(id, batch_id, original_name, stored_path, "
        "uploaded_at) VALUES(?, ?, ?, ?, ?)",
        (file_id, batch_id, original_name, stored_path, _now()),
    )
    return file_id


def update_file_validation(
    conn: sqlite3.Connection,
    file_id: str,
    *,
    valid: bool,
    error: str | None,
    reference_version: str | None,
    sample_count: int,
    record_count: int,
    samples: list[str],
    contigs: list[str],
) -> None:
    conn.execute(
        "UPDATE vcf_files SET valid=?, error=?, reference_version=?, "
        "sample_count=?, record_count=?, samples_json=?, contigs_json=? "
        "WHERE id=?",
        (
            1 if valid else 0,
            error,
            reference_version,
            sample_count,
            record_count,
            json.dumps(samples),
            json.dumps(contigs),
            file_id,
        ),
    )


def set_batch_reference(conn: sqlite3.Connection, batch_id: str, ref: str | None) -> None:
    conn.execute(
        "UPDATE batches SET reference_version=? WHERE id=?", (ref, batch_id)
    )


def reset_batch_errors(conn: sqlite3.Connection, batch_id: str) -> None:
    """Clear batch-level errors; used before recomputing them after a swap."""
    conn.execute(
        "UPDATE batches SET errors_json='[]', valid=1 WHERE id=?", (batch_id,)
    )


def add_batch_errors(conn: sqlite3.Connection, batch_id: str, errors: list[str]) -> None:
    row = get_batch(conn, batch_id)
    existing: list[str] = json.loads(row["errors_json"]) if row else []
    existing.extend(errors)
    conn.execute(
        "UPDATE batches SET errors_json=?, valid=? WHERE id=?",
        (json.dumps(existing), 0 if existing else 1, batch_id),
    )


def create_dup_group(
    conn: sqlite3.Connection,
    batch_id: str,
    sample_name: str,
    file_ids: list[str],
) -> str:
    group_id = new_id()
    conn.execute(
        "INSERT INTO dup_groups(id, batch_id, sample_name, file_ids_json, "
        "chosen_file_id) VALUES(?, ?, ?, ?, NULL)",
        (group_id, batch_id, sample_name, json.dumps(file_ids)),
    )
    return group_id


@contextmanager
def transaction(conn: sqlite3.Connection) -> Iterator[None]:
    """Exclusive write transaction on an autocommit connection."""
    conn.execute("BEGIN IMMEDIATE")
    try:
        yield
    except BaseException:
        conn.execute("ROLLBACK")
        raise
    else:
        conn.execute("COMMIT")


def bump_selection_version(conn: sqlite3.Connection, batch_id: str) -> int:
    """Invalidate every result/in-flight run tied to older selections."""
    with transaction(conn):
        return _bump_selection_version_locked(conn, batch_id)


def _bump_selection_version_locked(
    conn: sqlite3.Connection, batch_id: str, reason: str | None = None
) -> int:
    """Version bump for callers already inside :func:`transaction`."""
    row = get_batch(conn, batch_id)
    new_version = row["selection_version"] + 1
    conn.execute(
        "UPDATE batches SET selection_version=? WHERE id=?",
        (new_version, batch_id),
    )
    _invalidate_generations(
        conn,
        batch_id,
        reason or "取用源选择已变更，旧结果与在途生成立即失效",
    )
    return new_version


def _invalidate_generations(
    conn: sqlite3.Connection, batch_id: str, reason: str
) -> None:
    """Mark every completed or in-flight generation as superseded."""
    conn.execute(
        "UPDATE generations SET status=?, finished_at=?, error=? "
        "WHERE batch_id=? AND status IN (?, ?, ?)",
        (
            config.STATUS_SUPERSEDED,
            _now(),
            reason,
            batch_id,
            config.STATUS_PENDING,
            config.STATUS_RUNNING,
            config.STATUS_DONE,
        ),
    )


def set_dup_group_choice(
    conn: sqlite3.Connection, group_id: str, file_id: str | None
) -> sqlite3.Row | None:
    """Persist a chosen source; returns the updated row or None on bad input."""
    group = get_dup_group(conn, group_id)
    if group is None:
        return None
    if file_id is not None:
        allowed = json.loads(group["file_ids_json"])
        if file_id not in allowed:
            raise ValueError("chosen file is not an original source of this sample")
    conn.execute(
        "UPDATE dup_groups SET chosen_file_id=? WHERE id=?",
        (file_id, group_id),
    )
    return get_dup_group(conn, group_id)


def create_generation(
    conn: sqlite3.Connection, batch_id: str, selection_version: int
) -> str:
    gen_id = new_id()
    conn.execute(
        "INSERT INTO generations(id, batch_id, status, selection_version, "
        "created_at) VALUES(?, ?, ?, ?, ?)",
        (gen_id, batch_id, config.STATUS_PENDING, selection_version, _now()),
    )
    return gen_id


def mark_generation(
    conn: sqlite3.Connection,
    gen_id: str,
    status: str,
    *,
    output_path: str | None = None,
    error: str | None = None,
    started: bool = False,
    finished: bool = False,
) -> None:
    sets = ["status=?"]
    params: list[Any] = [status]
    if output_path is not None:
        sets.append("output_path=?")
        params.append(output_path)
    if error is not None:
        sets.append("error=?")
        params.append(error)
    if started:
        sets.append("started_at=?")
        params.append(_now())
    if finished:
        sets.append("finished_at=?")
        params.append(_now())
    params.append(gen_id)
    conn.execute(f"UPDATE generations SET {', '.join(sets)} WHERE id=?", params)


# --------------------------------------------------------------------------- #
# Replacement attempts
# --------------------------------------------------------------------------- #
def create_replacement(
    conn: sqlite3.Connection,
    batch_id: str,
    target_file_id: str,
    original_name: str,
    temp_path: str,
    base_seq: int,
) -> str:
    replacement_id = new_id()
    conn.execute(
        "INSERT INTO replacements(id, batch_id, target_file_id, original_name, "
        "temp_path, status, base_seq, created_at) "
        "VALUES(?, ?, ?, ?, ?, ?, ?, ?)",
        (
            replacement_id,
            batch_id,
            target_file_id,
            original_name,
            temp_path,
            config.REPLACE_VALIDATING,
            base_seq,
            _now(),
        ),
    )
    return replacement_id


def get_replacements(conn: sqlite3.Connection, batch_id: str) -> list[sqlite3.Row]:
    return conn.execute(
        "SELECT * FROM replacements WHERE batch_id=? ORDER BY rowid DESC",
        (batch_id,),
    ).fetchall()


def get_replacement(conn: sqlite3.Connection, replacement_id: str) -> sqlite3.Row | None:
    return conn.execute(
        "SELECT * FROM replacements WHERE id=?", (replacement_id,)
    ).fetchone()


def has_validating_replacement(conn: sqlite3.Connection, batch_id: str) -> bool:
    row = conn.execute(
        "SELECT 1 FROM replacements WHERE batch_id=? AND status=? LIMIT 1",
        (batch_id, config.REPLACE_VALIDATING),
    ).fetchone()
    return row is not None


def mark_replacement(
    conn: sqlite3.Connection,
    replacement_id: str,
    status: str,
    *,
    error: str | None = None,
    finished: bool = True,
) -> None:
    sets = ["status=?"]
    params: list[Any] = [status]
    if error is not None:
        sets.append("error=?")
        params.append(error)
    if finished:
        sets.append("finished_at=?")
        params.append(_now())
    params.append(replacement_id)
    conn.execute(
        f"UPDATE replacements SET {', '.join(sets)} WHERE id=?", params
    )


def update_file_identity(
    conn: sqlite3.Connection, file_id: str, original_name: str, stored_path: str
) -> None:
    conn.execute(
        "UPDATE vcf_files SET original_name=?, stored_path=? WHERE id=?",
        (original_name, stored_path, file_id),
    )


def build_group_plan(
    conn: sqlite3.Connection,
    batch_id: str,
    *,
    old_groups: list[sqlite3.Row] | None = None,
    target_file_id: str | None = None,
) -> list[tuple[str, list[str], str | None]]:
    """Derive the sample -> source groups from the current valid files.

    Rules, evaluated after the target original (if any) was already swapped:

    * one source -> auto-selected;
    * multiple sources, was also multi-source before the swap and its chosen
      original is not the replaced file -> keep that selection;
    * otherwise unresolved (needs confirmation).

    ``old_groups`` is the pre-swap group snapshot; pass the groups read inside
    the same transaction to avoid racing a concurrent selection change.
    """
    prior: dict[str, tuple[int, str | None]] = {}
    for group in old_groups or []:
        ids = json.loads(group["file_ids_json"])
        # A single-source group had no operator choice; keep it as such so a
        # swap that introduces a second source does not inherit a "selection".
        prior[group["sample_name"]] = (
            len(ids),
            group["chosen_file_id"] if len(ids) > 1 else None,
        )

    sample_to_files: dict[str, list[str]] = {}
    for row in get_files(conn, batch_id):
        if not row["valid"]:
            continue
        for sample in json.loads(row["samples_json"]):
            sample_to_files.setdefault(sample, []).append(row["id"])

    plan: list[tuple[str, list[str], str | None]] = []
    for sample, file_ids in sorted(sample_to_files.items()):
        if len(file_ids) == 1:
            plan.append((sample, file_ids, file_ids[0]))
            continue
        prior_count, prior_choice = prior.get(sample, (0, None))
        if (
            prior_count > 1
            and prior_choice is not None
            and prior_choice in file_ids
            and prior_choice != target_file_id
        ):
            plan.append((sample, file_ids, prior_choice))
        else:
            plan.append((sample, file_ids, None))
    return plan


def commit_replacement(
    conn: sqlite3.Connection,
    batch_id: str,
    replacement_id: str,
    target_file_id: str,
    original_name: str,
    stored_path: str,
    base_seq: int,
    validation: dict,
    reference_version: str | None,
    batch_errors: list[str],
) -> bool:
    """Swap the original file and rebuild derived state in one transaction.

    Returns False (and changes nothing) when another replacement committed
    after this attempt was submitted, i.e. the observed ``base_seq`` is stale.
    Among concurrent attempts whichever commits first wins; every other one
    observes a bumped ``replace_seq`` and is rejected.
    """
    with transaction(conn):
        row = get_batch(conn, batch_id)
        target = get_file(conn, target_file_id)
        if (
            row is None
            or row["replace_seq"] != base_seq
            or target is None
            or target["batch_id"] != batch_id
        ):
            return False
        # Snapshot groups in the same transaction as the swap, so a selection
        # submitted while we were validating is reflected in the new plan.
        old_groups = get_dup_groups(conn, batch_id)
        update_file_identity(conn, target_file_id, original_name, stored_path)
        update_file_validation(conn, target_file_id, **validation)
        conn.execute(
            "UPDATE batches SET replace_seq=? WHERE id=?",
            (base_seq + 1, batch_id),
        )
        set_batch_reference(conn, batch_id, reference_version)
        reset_batch_errors(conn, batch_id)
        if batch_errors:
            add_batch_errors(conn, batch_id, batch_errors)
        conn.execute("DELETE FROM dup_groups WHERE batch_id=?", (batch_id,))
        for sample, file_ids, choice in build_group_plan(
            conn, batch_id, old_groups=old_groups, target_file_id=target_file_id
        ):
            group_id = create_dup_group(conn, batch_id, sample, file_ids)
            if choice is not None:
                set_dup_group_choice(conn, group_id, choice)
        # A swap changes what a cohort file means: invalidate old outputs and
        # in-flight merges and move the selection version forward so any
        # worker still reading the previous version aborts.
        _bump_selection_version_locked(
            conn,
            batch_id,
            "原件已换件，旧合并结果与在途生成失效，需重新生成",
        )
        mark_replacement(
            conn, replacement_id, config.REPLACE_COMMITTED, finished=True
        )
        return True
