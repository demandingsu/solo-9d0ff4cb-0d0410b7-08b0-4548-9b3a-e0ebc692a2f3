"""归档功能自测:核对、打包、校验、幂等、失败回退、中断恢复、并发下载。"""
import sys
import time
import zipfile
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from app import archive as arc  # noqa: E402
from app.config import archives_dir, db_path, published_dir, staging_dir  # noqa: E402
from app.db import db_conn  # noqa: E402

SAMPLES_CSV = "sample_id,index\nS1,AAAACCCC\nS2,AAAACCGG\nS3,GGGGTTTT\n"

READS = [
    ("r001", "AAAACCCC"),
    ("r002", "AAAACCGG"),
    ("r003", "AAAACCCA"),
    ("r004", "AAAACCGC"),
    ("r005", "NNNNNNNN"),
    ("r006", "GGGGTTTA"),
    ("r007", "AAAACCCC"),
]

SEQ1, SEQ2 = "ACGT" * 10, "TGCA" * 10
QUAL = "I" * 40


def fq(records):
    return "".join(f"@{h}\n{s}\n+\n{q}\n" for h, s, q in records)


def make_inputs(reads=READS):
    r1, r2, i1 = [], [], []
    for n, (rid, idx) in enumerate(reads):
        head = f"INST:1:FC:1:1101:{1000 + n}:{2000 + n}"
        r1.append((f"{head} 1:N:0:0", SEQ1, QUAL))
        r2.append((f"{head} 2:N:0:0", SEQ2, QUAL))
        i1.append((f"{head} 1:N:0:0", idx, "I" * len(idx)))
    return fq(r1).encode(), fq(r2).encode(), fq(i1).encode(), SAMPLES_CSV.encode()


@pytest.fixture()
def client(tmp_path, monkeypatch):
    monkeypatch.setenv("DATA_DIR", str(tmp_path / "data"))
    from fastapi.testclient import TestClient
    from app.main import app

    with TestClient(app) as c:
        c.data_dir = tmp_path / "data"
        yield c


def submit_good(client, run_number="RUN1"):
    r1, r2, i1, csv = make_inputs()
    return client.post(
        "/runs",
        data={"run_number": run_number},
        files={
            "r1": ("r1.fastq", r1),
            "r2": ("r2.fastq", r2),
            "i1": ("i1.fastq", i1),
            "samples": ("samples.csv", csv),
        },
        follow_redirects=False,
    )


def publish(client, run_number="RUN1"):
    resp = submit_good(client, run_number)
    assert resp.status_code == 303
    rid = int(resp.headers["location"].rsplit("/", 1)[1])
    deadline = time.time() + 15
    while time.time() < deadline:
        st = client.get(f"/api/runs/{rid}").json()
        if st["status"] != "processing":
            break
        time.sleep(0.02)
    assert st["status"] == "published", st
    return rid, st


def post_archive(client, rid):
    return client.post(f"/runs/{rid}/archive", follow_redirects=False)


def wait_archived(client, rid, timeout=15):
    deadline = time.time() + timeout
    last = None
    while time.time() < deadline:
        last = client.get(f"/api/runs/{rid}").json()
        if last["status"] not in ("processing", "archiving"):
            return last
        time.sleep(0.02)
    raise TimeoutError(f"归档超时: {last}")


def expected_zip_entries(run_number):
    prefix = f"{run_number}/"
    names = [prefix + "samples.csv", prefix + "counts.csv", prefix + "demux.log"]
    for sid in ("S1", "S2", "S3", "undetermined"):
        names += [prefix + f"{sid}_R1.fastq", prefix + f"{sid}_R2.fastq"]
    return sorted(names)


# ------------------------------------------------------------ 正常归档

def test_archive_happy_path(client):
    rid, before = publish(client, "RUN1")
    run_dir = client.data_dir / "published" / "RUN1"
    loose_before = {p.name: p.read_bytes() for p in run_dir.iterdir()}

    resp = post_archive(client, rid)
    assert resp.status_code == 303 and resp.headers["location"] == f"/runs/{rid}"
    st = wait_archived(client, rid)

    assert st["status"] == "archived"
    assert st["total_pairs"] == 7  # 计数原样保留
    assert {c["sample_id"]: c["pairs"] for c in st["counts"]} == {
        "S1": 3, "S2": 1, "S3": 1, "undetermined": 2,
    }
    assert st["archived_at"]

    zip_path = client.data_dir / "archives" / "RUN1.zip"
    assert zip_path.is_file()
    # 散装双端文件已清理,卷空间腾出
    assert not run_dir.exists() or not any(run_dir.iterdir())
    # 暂存目录清空
    assert not any((client.data_dir / "staging").iterdir())

    with zipfile.ZipFile(zip_path) as zf:
        assert zf.testzip() is None
        assert sorted(zf.namelist()) == expected_zip_entries("RUN1")
        # 双端文件内容与归档前完全一致
        for name, data in loose_before.items():
            if name.endswith(".fastq"):
                assert zf.read(f"RUN1/{name}") == data
        counts = zf.read("RUN1/counts.csv").decode()
        assert "sample_id,index_seq,pairs,corrected" in counts
        assert "S1,AAAACCCC,3,1" in counts
        assert "undetermined,,2,0" in counts
        assert "TOTAL,,7," in counts
        samples = zf.read("RUN1/samples.csv").decode()
        assert "S1,AAAACCCC" in samples and "undetermined" not in samples
        assert "总读段对: 7" in zf.read("RUN1/demux.log").decode()

    # 归档后散装下载停止,归档 ZIP 可下载
    dl = client.get(f"/runs/{rid}/download/S1_R1.fastq")
    assert dl.status_code == 410
    z = client.get(f"/runs/{rid}/archive")
    assert z.status_code == 200 and z.headers["content-type"] == "application/zip"
    assert z.content == zip_path.read_bytes()

    page = client.get(f"/runs/{rid}")
    assert page.status_code == 200 and "已归档" in page.text
    assert "RUN1.zip" in page.text and "3" in page.text  # 页面仍显示原计数


def test_archive_idempotent_and_guards(client):
    rid, _ = publish(client, "RUN6")
    assert post_archive(client, rid).status_code == 303
    # 归档中的重复提交被拒绝,不会启动第二份
    mid = post_archive(client, rid)
    assert mid.status_code == 409 and "归档" in mid.json()["detail"]
    assert wait_archived(client, rid)["status"] == "archived"
    # 归档完成后重复提交:只产生一份 ZIP
    again = post_archive(client, rid)
    assert again.status_code == 303
    zips = list((client.data_dir / "archives").iterdir())
    assert len(zips) == 1 and zips[0].name == "RUN6.zip"
    # 未发布的运行不能归档:新造一个失败运行
    r1, r2, i1, csv = make_inputs()
    bad_r2 = r2[: r2.rfind(b"\n", 0, r2.rfind(b"\n"))]
    resp = client.post(
        "/runs", data={"run_number": "RUN7"},
        files={"r1": ("r1.fastq", r1), "r2": ("r2.fastq", bad_r2),
               "i1": ("i1.fastq", i1), "samples": ("samples.csv", csv)},
        follow_redirects=False,
    )
    fid = int(resp.headers["location"].rsplit("/", 1)[1])
    deadline = time.time() + 15
    while client.get(f"/api/runs/{fid}").json()["status"] == "processing":
        time.sleep(0.02)
    rej = post_archive(client, fid)
    assert rej.status_code == 409
    assert client.get(f"/runs/{fid}/archive").status_code == 404


def test_loose_download_still_works_during_archiving(client, monkeypatch):
    rid, _ = publish(client, "RUN8")
    run_dir = client.data_dir / "published" / "RUN8"
    original = (run_dir / "S1_R1.fastq").read_bytes()
    entered = __import__("threading").Event()
    proceed = __import__("threading").Event()

    from app.main import archive_run as real_archive_run

    def slow_archive(run_id):
        # 在状态已切到 archiving、文件尚未清理时放行检查
        entered.set()
        proceed.wait(timeout=10)
        real_archive_run(run_id)

    monkeypatch.setattr("app.main.archive_run", slow_archive)
    assert post_archive(client, rid).status_code == 303
    assert entered.wait(timeout=5)
    try:
        st = client.get(f"/api/runs/{rid}").json()
        assert st["status"] == "archiving"
        # 归档进行中的新请求:拿到的是完整旧文件,不是半成品
        dl = client.get(f"/runs/{rid}/download/S1_R1.fastq")
        assert dl.status_code == 200 and dl.content == original
        # 此时归档 ZIP 尚未对外提供
        assert client.get(f"/runs/{rid}/archive").status_code == 409
    finally:
        proceed.set()
    assert wait_archived(client, rid)["status"] == "archived"


# ------------------------------------------------------------ 失败场景

def test_verification_failure_keeps_run_intact(client, monkeypatch):
    rid, before = publish(client, "RUN9")
    run_dir = client.data_dir / "published" / "RUN9"
    loose_snapshot = {p.name: p.read_bytes() for p in run_dir.iterdir()}
    target = run_dir / "S1_R1.fastq"
    target.write_bytes(target.read_bytes() + b"@broken\nAC\n+\nII\n")  # 多出 1 对

    resp = post_archive(client, rid)
    assert resp.status_code == 303
    st = wait_archived(client, rid)
    assert st["status"] == "published"  # 安全回退,原状态恢复
    assert "归档失败" in st["archive_error"]
    # 散装文件一个没少(被注入的额外记录也原样保留),没有归档 ZIP
    assert {p.name for p in run_dir.iterdir()} == set(loose_snapshot)
    assert not (client.data_dir / "archives" / "RUN9.zip").exists()
    assert not any((client.data_dir / "staging").iterdir())
    # 普通下载照常工作
    assert client.get(f"/runs/{rid}/download/S1_R2.fastq").status_code == 200
    # 页面提示上次归档失败,且可重新归档(恢复文件后成功)
    target.write_bytes(loose_snapshot["S1_R1.fastq"])
    assert post_archive(client, rid).status_code == 303
    assert wait_archived(client, rid)["status"] == "archived"


def test_count_mismatch_fails_archive(client, tmp_path):
    rid, _ = publish(client, "RUN10")
    run_dir = client.data_dir / "published" / "RUN10"
    # 直接改库:入库对数与文件实际对数不符
    with db_conn(db_path()) as conn:
        conn.execute(
            "UPDATE sample_counts SET pairs=pairs+1 WHERE run_id=? AND sample_id='S1'",
            (rid,),
        )
    post_archive(client, rid)
    st = wait_archived(client, rid)
    assert st["status"] == "published"
    assert "入库对数" in st["archive_error"]
    assert (run_dir / "S1_R1.fastq").is_file()
    assert not (client.data_dir / "archives" / "RUN10.zip").exists()


def test_zip_build_failure_keeps_run_intact(client, monkeypatch):
    rid, _ = publish(client, "RUN11")
    run_dir = client.data_dir / "published" / "RUN11"
    snapshot = {p.name: p.read_bytes() for p in run_dir.iterdir()}

    def boom(*args, **kwargs):
        raise arc.ArchiveError("模拟磁盘写失败")

    monkeypatch.setattr(arc, "build_archive_zip", boom)
    post_archive(client, rid)
    st = wait_archived(client, rid)
    assert st["status"] == "published" and "模拟磁盘写失败" in st["archive_error"]
    assert {p.name: p.read_bytes() for p in run_dir.iterdir()} == snapshot
    assert not (client.data_dir / "archives" / "RUN11.zip").exists()


def test_corrupt_zip_rejected(client, monkeypatch):
    rid, _ = publish(client, "RUN12")
    real_verify = arc.verify_zip

    def corrupt_then_verify(zip_path, entries):
        with open(zip_path, "r+b") as fh:
            fh.seek(1024)
            fh.write(b"\xff" * 64)  # 破坏压缩数据,CRC 应失败
        real_verify(zip_path, entries)

    monkeypatch.setattr(arc, "verify_zip", corrupt_then_verify)
    post_archive(client, rid)
    st = wait_archived(client, rid)
    assert st["status"] == "published"
    assert "归档失败" in st["archive_error"]
    assert (client.data_dir / "published" / "RUN12" / "S1_R1.fastq").is_file()
    assert not (client.data_dir / "archives" / "RUN12.zip").exists()


# ------------------------------------------------------------ 中断恢复

def test_recovery_mid_archiving_with_loose_files(client):
    """归档线程在改名前被杀:状态 archiving、暂存半成品、散装文件完整。"""
    from app.main import recover_interrupted

    rid, _ = publish(client, "RUN13")
    run_dir = client.data_dir / "published" / "RUN13"
    snapshot = {p.name: p.read_bytes() for p in run_dir.iterdir()}
    staging = client.data_dir / "staging" / "archive-deadbeef"
    staging.mkdir(parents=True)
    (staging / "RUN13.zip").write_bytes(b"not a real zip")
    with db_conn(db_path()) as conn:
        conn.execute(
            "UPDATE runs SET status='archiving', staging_dir=? WHERE id=?",
            (str(staging), rid),
        )

    recover_interrupted()

    st = client.get(f"/api/runs/{rid}").json()
    assert st["status"] == "published" and st["archive_error"]
    assert {p.name: p.read_bytes() for p in run_dir.iterdir()} == snapshot
    assert not staging.exists()
    assert not (client.data_dir / "archives" / "RUN13.zip").exists()
    # 恢复后可以重新归档成功
    assert post_archive(client, rid).status_code == 303
    assert wait_archived(client, rid)["status"] == "archived"


def test_recovery_completes_archive_when_zip_ready(client):
    """改名后、清散装前被杀:完整 ZIP 已在归档区,重启后完成归档不丢数据。"""
    from app.main import recover_interrupted

    rid, _ = publish(client, "RUN14")
    run_dir = client.data_dir / "published" / "RUN14"
    loose = {p.name: p.read_bytes() for p in run_dir.iterdir()}
    with db_conn(db_path()) as conn:
        rows = conn.execute(
            "SELECT * FROM sample_counts WHERE run_id=? ORDER BY sample_id='undetermined', sample_id",
            (rid,),
        ).fetchall()
        total = conn.execute("SELECT total_pairs FROM runs WHERE id=?", (rid,)).fetchone()[0]
    tmp_zip = client.data_dir / "staging" / "archive-cafe" / "RUN14.zip"
    tmp_zip.parent.mkdir(parents=True)
    entries = arc.build_archive_zip(tmp_zip, "RUN14", run_dir, rows, total)
    arc.verify_zip(tmp_zip, entries)
    final_zip = client.data_dir / "archives" / "RUN14.zip"
    final_zip.parent.mkdir(parents=True, exist_ok=True)
    import os
    os.rename(tmp_zip, final_zip)
    with db_conn(db_path()) as conn:
        conn.execute(
            "UPDATE runs SET status='archiving', staging_dir=? WHERE id=?",
            (str(tmp_zip.parent), rid),
        )

    recover_interrupted()

    st = client.get(f"/api/runs/{rid}").json()
    assert st["status"] == "archived"
    assert final_zip.is_file()
    assert not run_dir.exists() or not any(run_dir.iterdir())
    with zipfile.ZipFile(final_zip) as zf:
        for name, data in loose.items():
            if name.endswith(".fastq"):
                assert zf.read(f"RUN14/{name}") == data


def test_recovery_rejects_corrupt_zip_even_if_central_dir_ok(client):
    """中央目录可读但压缩数据损坏的 ZIP 不算完整副本:散装在则回退,不删唯一副本。"""
    from app.main import recover_interrupted

    rid, _ = publish(client, "RUN14B")
    run_dir = client.data_dir / "published" / "RUN14B"
    snapshot = {p.name: p.read_bytes() for p in run_dir.iterdir()}
    with db_conn(db_path()) as conn:
        rows = conn.execute(
            "SELECT * FROM sample_counts WHERE run_id=? ORDER BY sample_id='undetermined', sample_id",
            (rid,),
        ).fetchall()
        total = conn.execute("SELECT total_pairs FROM runs WHERE id=?", (rid,)).fetchone()[0]
    tmp_zip = client.data_dir / "staging" / "archive-bad" / "RUN14B.zip"
    tmp_zip.parent.mkdir(parents=True)
    entries = arc.build_archive_zip(tmp_zip, "RUN14B", run_dir, rows, total)
    final_zip = client.data_dir / "archives" / "RUN14B.zip"
    final_zip.parent.mkdir(parents=True, exist_ok=True)
    import os
    os.rename(tmp_zip, final_zip)
    # 破坏一个 FASTQ 条目的压缩字节;中央目录仍在,但 CRC 必然失败
    with zipfile.ZipFile(final_zip, "r") as zf:
        info = zf.getinfo("RUN14B/S1_R1.fastq")
        offset = info.header_offset + 30
    with open(final_zip, "r+b") as fh:
        fh.seek(offset)
        fh.write(b"\xff" * 32)
    with zipfile.ZipFile(final_zip) as zf:  # 中央目录确实仍可打开
        assert "RUN14B/S1_R1.fastq" in zf.namelist()
    with db_conn(db_path()) as conn:
        conn.execute(
            "UPDATE runs SET status='archiving', staging_dir=? WHERE id=?",
            (str(tmp_zip.parent), rid),
        )

    recover_interrupted()

    st = client.get(f"/api/runs/{rid}").json()
    assert st["status"] == "published"  # 拒绝承认损坏 ZIP,安全回退
    assert {p.name: p.read_bytes() for p in run_dir.iterdir()} == snapshot
    assert not final_zip.exists()
    assert client.get(f"/runs/{rid}/download/S1_R1.fastq").status_code == 200


def test_recovery_archived_row_with_missing_zip_and_loose(client):
    """archived 行的 ZIP 丢失但散装完整:回退 published,绝不删唯一副本。"""
    from app.main import recover_interrupted

    rid, _ = publish(client, "RUN15")
    run_dir = client.data_dir / "published" / "RUN15"
    with db_conn(db_path()) as conn:
        conn.execute("UPDATE runs SET status='archived', archived_at=? WHERE id=?",
                     ("t", rid))
    recover_interrupted()
    st = client.get(f"/api/runs/{rid}").json()
    assert st["status"] == "published"
    assert (run_dir / "S1_R1.fastq").is_file()


def test_recovery_orphan_dirs_cleaned(client):
    rid, _ = publish(client, "RUN16")
    post_archive(client, rid)
    assert wait_archived(client, rid)["status"] == "archived"
    # 无主暂存目录与无主归档 ZIP
    stray_stage = client.data_dir / "staging" / "archive-stray"
    stray_stage.mkdir(parents=True)
    (stray_stage / "x").write_text("half")
    stray_zip = client.data_dir / "archives" / "GHOST.zip"
    stray_zip.write_bytes(b"PK\x03\x04junk")

    from app.main import recover_interrupted
    recover_interrupted()
    assert not stray_stage.exists() and not stray_zip.exists()
    # 正式归档仍在
    assert (client.data_dir / "archives" / "RUN16.zip").is_file()


# ------------------------------------------------------------ 流式与单元

def test_large_files_streamed_not_loaded_whole(client, monkeypatch):
    """大文件归档:分块搬运,峰值缓冲远小于文件总大小,且内容字节一致。"""
    rid, _ = publish(client, "BIG1")
    run_dir = client.data_dir / "published" / "BIG1"
    big = (b"@read%d\n" + b"ACGT" * 250 + b"\n+\n" + b"I" * 1000 + b"\n")
    payload = b"".join(big % n for n in range(400))  # 每组约 0.5 MiB
    for sid in ("S1", "S2", "S3", "undetermined"):
        (run_dir / f"{sid}_R1.fastq").write_bytes(payload)
        (run_dir / f"{sid}_R2.fastq").write_bytes(payload)
    # 散装文件被换成大文件后,同步入库计数与总对数,使核对仍闭合
    with db_conn(db_path()) as conn:
        conn.execute(
            "UPDATE sample_counts SET pairs=400, corrected=0 WHERE run_id=?", (rid,)
        )
        conn.execute("UPDATE runs SET total_pairs=1600 WHERE id=?", (rid,))
    total_size = sum(p.stat().st_size for p in run_dir.glob("*.fastq"))
    assert total_size > 4 * 1024 * 1024

    # 用带计数的 open 包装归档模块:每次 read 最多取走 1 MiB,证明没有整份读入
    import builtins
    import app.archive as archive_mod
    real_open = builtins.open
    max_read = 0

    class TrackingFile:
        def __init__(self, fh):
            self._fh = fh

        def read(self, size=-1):
            nonlocal max_read
            if isinstance(size, int) and size > 0:
                max_read = max(max_read, size)
            return self._fh.read(size)

        def write(self, data):
            return self._fh.write(data)

        def close(self):
            return self._fh.close()

        def __enter__(self):
            self._fh.__enter__()
            return self

        def __exit__(self, *exc):
            return self._fh.__exit__(*exc)

    def tracking_open(path, *args, **kwargs):
        fh = real_open(path, *args, **kwargs)
        return TrackingFile(fh) if "rb" in (args[0] if args else kwargs.get("mode", "")) else fh

    monkeypatch.setattr(archive_mod, "open", tracking_open, raising=False)
    post_archive(client, rid)
    st = wait_archived(client, rid)
    assert st["status"] == "archived"
    assert 0 < max_read <= arc.CHUNK

    with zipfile.ZipFile(client.data_dir / "archives" / "BIG1.zip") as zf:
        assert zf.testzip() is None
        for sid in ("S1", "S2", "S3", "undetermined"):
            assert zf.read(f"BIG1/{sid}_R1.fastq") == payload
            assert zf.read(f"BIG1/{sid}_R2.fastq") == payload


def test_verify_group_detects_pair_mismatch(tmp_path):
    d = tmp_path / "RUN"
    d.mkdir()
    (d / "S1_R1.fastq").write_text("@a\nACGT\n+\nIIII\n@b\nACGT\n+\nIIII\n")
    (d / "S1_R2.fastq").write_text("@a\nACGT\n+\nIIII\n@x\nACGT\n+\nIIII\n")
    with pytest.raises(arc.ArchiveError, match="读段号不一致"):
        arc.verify_group(d, "S1", 2)
    (d / "S1_R2.fastq").write_text("@a\nACGT\n+\nIIII\n")  # 记录数不同
    with pytest.raises(arc.ArchiveError, match="记录数不一致"):
        arc.verify_group(d, "S1", 2)
    # 两份文件都是 1 对,但入库计数声称 9 对
    (d / "S1_R1.fastq").write_text("@a\nACGT\n+\nIIII\n")
    with pytest.raises(arc.ArchiveError, match="入库对数"):
        arc.verify_group(d, "S1", 9)


def test_verify_zip_rejects_manipulated_entry(tmp_path):
    zip_path = tmp_path / "r.zip"
    with zipfile.ZipFile(zip_path, "w") as zf:
        zf.writestr("RUN/a.txt", b"hello")
    with pytest.raises(arc.ArchiveError, match="条目清单不一致"):
        arc.verify_zip(zip_path, ["RUN/other.txt"])
    with open(zip_path, "r+b") as fh:
        fh.seek(20)
        fh.write(b"\x00" * 8)
    with pytest.raises(arc.ArchiveError):
        arc.verify_zip(zip_path, ["RUN/a.txt"])
