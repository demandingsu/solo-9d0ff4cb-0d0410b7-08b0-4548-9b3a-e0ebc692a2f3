"""SQLite 访问:每次操作独立连接,提交即持久化;WAL 模式支持读写并发。"""
import sqlite3
from contextlib import contextmanager

SCHEMA = """
CREATE TABLE IF NOT EXISTS runs (
  id             INTEGER PRIMARY KEY AUTOINCREMENT,
  run_number     TEXT NOT NULL UNIQUE,
  status         TEXT NOT NULL CHECK (status IN
                   ('processing', 'published', 'archiving', 'archived', 'failed')),
  error          TEXT,
  staging_dir    TEXT NOT NULL,
  r1_name        TEXT NOT NULL DEFAULT '',
  r2_name        TEXT NOT NULL DEFAULT '',
  i1_name        TEXT NOT NULL DEFAULT '',
  samples_name   TEXT NOT NULL DEFAULT '',
  total_pairs    INTEGER NOT NULL DEFAULT 0,
  created_at     TEXT NOT NULL,
  updated_at     TEXT NOT NULL,
  published_at   TEXT,
  archived_at    TEXT,
  archive_error  TEXT
);
CREATE TABLE IF NOT EXISTS sample_counts (
  run_id     INTEGER NOT NULL REFERENCES runs(id) ON DELETE CASCADE,
  sample_id  TEXT NOT NULL,                -- 'undetermined' 表示未分配组
  index_seq  TEXT NOT NULL DEFAULT '',
  pairs      INTEGER NOT NULL,
  corrected  INTEGER NOT NULL DEFAULT 0,   -- 其中经单碱基纠错的对数
  PRIMARY KEY (run_id, sample_id)
);
"""

# 旧库(只有 processing/published/failed 三态)迁移到带归档字段的新表
_MIGRATION = """
ALTER TABLE runs RENAME TO runs_old;
CREATE TABLE runs (
  id             INTEGER PRIMARY KEY AUTOINCREMENT,
  run_number     TEXT NOT NULL UNIQUE,
  status         TEXT NOT NULL CHECK (status IN
                   ('processing', 'published', 'archiving', 'archived', 'failed')),
  error          TEXT,
  staging_dir    TEXT NOT NULL,
  r1_name        TEXT NOT NULL DEFAULT '',
  r2_name        TEXT NOT NULL DEFAULT '',
  i1_name        TEXT NOT NULL DEFAULT '',
  samples_name   TEXT NOT NULL DEFAULT '',
  total_pairs    INTEGER NOT NULL DEFAULT 0,
  created_at     TEXT NOT NULL,
  updated_at     TEXT NOT NULL,
  published_at   TEXT,
  archived_at    TEXT,
  archive_error  TEXT
);
INSERT INTO runs (id, run_number, status, error, staging_dir, r1_name, r2_name,
                  i1_name, samples_name, total_pairs, created_at, updated_at,
                  published_at, archived_at, archive_error)
SELECT id, run_number, status, error, staging_dir, r1_name, r2_name, i1_name,
       samples_name, total_pairs, created_at, updated_at, published_at,
       NULL, NULL
FROM runs_old;
DROP TABLE runs_old;
"""


@contextmanager
def db_conn(path):
    """成功提交、异常回滚、用完关闭的连接上下文。"""
    conn = sqlite3.connect(str(path), timeout=30)
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA foreign_keys = ON")
    try:
        yield conn
        conn.commit()
    except Exception:
        conn.rollback()
        raise
    finally:
        conn.close()


def init_db(path):
    # 迁移需要临时关掉外键,单独使用一条裸连接
    conn = sqlite3.connect(str(path), timeout=30)
    conn.row_factory = sqlite3.Row
    try:
        conn.execute("PRAGMA journal_mode = WAL")
        cols = {r["name"] for r in conn.execute("PRAGMA table_info(runs)")}
        if cols and not {"archived_at", "archive_error"} <= cols:
            conn.execute("PRAGMA foreign_keys = OFF")
            conn.executescript(_MIGRATION)
            conn.execute("PRAGMA foreign_keys = ON")
        else:
            conn.executescript(SCHEMA)
        conn.commit()
    finally:
        conn.close()
