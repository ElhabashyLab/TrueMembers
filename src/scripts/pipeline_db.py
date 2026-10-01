import json
import sqlite3
from pathlib import Path
from typing import Optional


def _domains_to_json(domains: list[dict[str, tuple[int, int]]]) -> str:
    if not domains:
        return "[]"
    return json.dumps([{k: list(v) for k, v in d.items()} for d in domains])


def init_db(db_path: Path) -> None:
    db_path.parent.mkdir(parents=True, exist_ok=True)
    with sqlite3.connect(str(db_path)) as conn:
        conn.executescript(_SCHEMA)


def insert_member(db_path: Path, hit) -> None:
    with sqlite3.connect(str(db_path)) as conn:
        conn.execute(
            "INSERT OR REPLACE INTO members (unique_id, uniprot_id, region_start, region_end, "
            "sequence_path, structure_path, detection_method, qtmscore, is_member, domains) "
            "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            (
                hit.unique_id,
                hit.uniprot_id,
                hit.region_start,
                hit.region_end,
                str(hit.sequence_path),
                str(hit.structure_path),
                hit.detection_method.value,
                hit.qtmscore,
                int(hit.is_member),
                _domains_to_json(hit.domains),
            ),
        )


def insert_all_member(db_path: Path, hit) -> None:
    with sqlite3.connect(str(db_path)) as conn:
        conn.execute(
            "INSERT OR REPLACE INTO all_members (unique_id, uniprot_id, region_start, region_end, "
            "sequence_path, structure_path, detection_method, qtmscore, is_member, domains) "
            "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            (
                hit.unique_id,
                hit.uniprot_id,
                hit.region_start,
                hit.region_end,
                str(hit.sequence_path),
                str(hit.structure_path),
                hit.detection_method.value,
                hit.qtmscore,
                int(hit.is_member),
                _domains_to_json(hit.domains),
            ),
        )


def pick_unused_seeds(db_path: Path, limit: int) -> list[dict]:
    with sqlite3.connect(str(db_path)) as conn:
        conn.row_factory = sqlite3.Row
        rows = conn.execute(
            "SELECT * FROM members WHERE is_member = 1 "
            "AND unique_id NOT IN (SELECT unique_id FROM used_seeds) "
            "ORDER BY qtmscore DESC LIMIT ?",
            (limit,),
        ).fetchall()
    return [dict(r) for r in rows]


def insert_used_seed(db_path: Path, unique_id: str, uniprot_id: str,
                     iteration: int, new_members: Optional[int] = None) -> None:
    with sqlite3.connect(str(db_path)) as conn:
        conn.execute(
            "INSERT INTO used_seeds (unique_id, uniprot_id, iteration, new_members) "
            "VALUES (?, ?, ?, ?)",
            (unique_id, uniprot_id, iteration, new_members),
        )


def insert_sequence_identity(db_path: Path, seed_id: str, hit_id: str,
                             pident: float, detected: bool) -> None:
    with sqlite3.connect(str(db_path)) as conn:
        conn.execute(
            "INSERT OR REPLACE INTO sequence_identity (seed_id, hit_id, pident, detected) "
            "VALUES (?, ?, ?, ?)",
            (seed_id, hit_id, pident, int(detected)),
        )


_SCHEMA = """
CREATE TABLE IF NOT EXISTS members (
    unique_id TEXT PRIMARY KEY,
    uniprot_id TEXT NOT NULL,
    region_start INTEGER NOT NULL,
    region_end INTEGER NOT NULL,
    sequence_path TEXT NOT NULL,
    structure_path TEXT NOT NULL,
    detection_method TEXT NOT NULL,
    qtmscore REAL NOT NULL,
    is_member INTEGER NOT NULL DEFAULT 0,
    domains TEXT NOT NULL DEFAULT '[]'
);

CREATE TABLE IF NOT EXISTS all_members (
    unique_id TEXT PRIMARY KEY,
    uniprot_id TEXT NOT NULL,
    region_start INTEGER NOT NULL,
    region_end INTEGER NOT NULL,
    sequence_path TEXT NOT NULL,
    structure_path TEXT NOT NULL,
    detection_method TEXT NOT NULL,
    qtmscore REAL NOT NULL,
    is_member INTEGER NOT NULL DEFAULT 0,
    domains TEXT NOT NULL DEFAULT '[]'
);

CREATE TABLE IF NOT EXISTS used_seeds (
    unique_id TEXT NOT NULL,
    uniprot_id TEXT NOT NULL,
    iteration INTEGER NOT NULL,
    new_members INTEGER,
    PRIMARY KEY (unique_id, iteration)
);

CREATE TABLE IF NOT EXISTS sequence_identity (
    seed_id TEXT NOT NULL,
    hit_id TEXT NOT NULL,
    pident REAL NOT NULL,
    detected INTEGER NOT NULL DEFAULT 0,
    PRIMARY KEY (seed_id, hit_id)
);
"""