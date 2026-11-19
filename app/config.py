"""Application configuration and filesystem layout.

All runtime state lives under ``DATA_DIR`` so the whole directory can be
mounted as a Docker volume and survive container restarts.
"""
from __future__ import annotations

import os
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = Path(os.environ.get("VCF_DATA_DIR", BASE_DIR / "data"))

UPLOAD_DIR = DATA_DIR / "uploads"
OUTPUT_DIR = DATA_DIR / "outputs"
TMP_DIR = DATA_DIR / "tmp"
DB_PATH = DATA_DIR / "vcf_merge.db"

# Business rules from the specification.
MIN_FILES = 2
MAX_FILES = 3
MAX_UPLOAD_BYTES = int(os.environ.get("VCF_MAX_UPLOAD_BYTES", 200 * 1024 * 1024))
ALLOWED_EXTENSIONS = {".vcf", ".vcf.gz"}
ALLOWED_FORMAT_FIELDS = ("GT", "AD", "PL")
REQUIRED_FORMAT_FIELDS = ("GT",)
ALLELE_RE = __import__("re").compile(r"^[ACGT]+$")
DIPLOID_PLOIDY = 2

# Generation status values.
STATUS_PENDING = "pending"
STATUS_RUNNING = "running"
STATUS_DONE = "done"
STATUS_FAILED = "failed"
STATUS_SUPERSEDED = "superseded"
STATUS_ABORTED = "aborted"

TERMINAL_STATUSES = {STATUS_DONE, STATUS_FAILED, STATUS_SUPERSEDED, STATUS_ABORTED}

# Original-file replacement (lab reissues a corrected VCF) status values.
REPLACE_VALIDATING = "validating"
REPLACE_COMMITTED = "committed"
REPLACE_FAILED = "failed"
REPLACE_REJECTED = "rejected"

REPLACE_TERMINAL_STATUSES = {
    REPLACE_COMMITTED,
    REPLACE_FAILED,
    REPLACE_REJECTED,
}

# Optional artificial delay before a replacement validation runs. It makes the
# "selection may change while validating" and "a second replacement commits
# first" behaviours demonstrable from the web UI in Docker Compose. Default 0
# keeps tests and normal use fast.
REPLACE_VALIDATE_DELAY_SECONDS = float(
    os.environ.get("VCF_REPLACE_VALIDATE_DELAY_SECONDS", "0") or 0
)


def ensure_dirs() -> None:
    for directory in (UPLOAD_DIR, OUTPUT_DIR, TMP_DIR):
        directory.mkdir(parents=True, exist_ok=True)
