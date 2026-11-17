"""泳道拆样台:上传 R1/R2/I1 与样本索引 CSV,按样本拆出成对读段。

发布模型:
- 处理中:上传件与半成品输出都在 staging/<token>/,不对外可见。
- 成功:输出目录整体 os.rename 到 published/<运行号>/,同事务写入计数,
  之后状态置为 published。任何一步失败都不会留下半套已发布文件。
- 失败/中断:状态置 failed,暂存与可能的半成品目录被清理;同一运行号可重试。
- 已发布:同一运行号拒绝改写。

归档模型(长期不用的已发布泳道压缩归档,腾出 FASTQ 卷空间):
- 仅 published 运行可归档;先把状态乐观置为 archiving,归档中普通双端下载
  继续有效,重复提交不会启动第二份归档。
- 归档线程在 staging/archive-<token>/ 里逐组核对 R1/R2 配对与入库对数,
  流式打包 ZIP 并重新校验 CRC;通过后原子改名到 archives/<运行号>.zip,
  状态置 archived,最后才删除 published 下的散装 FASTQ。
- 打包或校验失败:状态回到 published、记 archive_error,原始文件原样保留。
- 散装文件只在"归档 ZIP 已是唯一完整副本"之后才删除;已开始的下载持有
  已打开的文件描述符,Linux 下 unlink 不影响其读完整文件。
- 启动恢复处理归档与清理之间的残留:有完整 ZIP 就完成归档,散装文件完整
  就回退到 published,二者皆无才标 failed,绝不删掉唯一完整副本。
"""
import os
import re
import shutil
import sqlite3
import threading
import uuid
import zipfile
from contextlib import asynccontextmanager
from datetime import datetime, timezone
from pathlib import Path

from fastapi import FastAPI, File, Form, HTTPException, Request, UploadFile
from fastapi.responses import FileResponse, HTMLResponse, JSONResponse, RedirectResponse
from fastapi.templating import Jinja2Templates

from . import archive as arc
from .config import archives_dir, db_path, published_dir, staging_dir
from .db import db_conn, init_db
from .demux import demultiplex
from .samples import load_samples

RUN_NUMBER_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,63}$")
STATUS_LABELS = {
    "processing": "处理中",
    "published": "已发布",
    "archiving": "归档中",
    "archived": "已归档",
    "failed": "失败",
}
UNDETERMINED = "undetermined"

BASE_DIR = Path(__file__).resolve().parent
templates = Jinja2Templates(directory=str(BASE_DIR / "templates"))


def utcnow() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _remove(path) -> None:
    p = Path(path)
    if p.is_dir() and not p.is_symlink():
        shutil.rmtree(p, ignore_errors=True)
    else:
        try:
            p.unlink()
        except (FileNotFoundError, NotADirectoryError):
            pass


def _published_intact(run_number: str) -> bool:
    """散装发布目录是否仍完整(目录存在且含 demux.log)。"""
    d = published_dir() / run_number
    return d.is_dir() and (d / "demux.log").is_file()


def _archive_zip_fully_ok(run_number: str, zip_path: Path) -> bool:
    """启动恢复用的完整校验:中央目录 + 全量 CRC32 + 条目清单与入库计数闭合。

    只有通过这一层,ZIP 才被承认为"完整归档副本";否则按散装副本决定回退。
    """
    try:
        with zipfile.ZipFile(zip_path, "r") as zf:
            if zf.testzip() is not None:
                return False
            names = set(zf.namelist())
        with db_conn(db_path()) as conn:
            rows = conn.execute(
                "SELECT sc.sample_id FROM sample_counts sc "
                "JOIN runs r ON r.id = sc.run_id WHERE r.run_number=?",
                (run_number,),
            ).fetchall()
        expected = {
            f"{run_number}/samples.csv", f"{run_number}/counts.csv",
            f"{run_number}/demux.log",
        }
        for r in rows:
            expected.add(f"{run_number}/{r['sample_id']}_R1.fastq")
            expected.add(f"{run_number}/{r['sample_id']}_R2.fastq")
        return names == expected
    except (zipfile.BadZipFile, OSError):
        return False


def _finish_archived(conn, row) -> None:
    """完整归档 ZIP 已就位:置 archived,清散装文件与暂存目录。"""
    conn.execute(
        "UPDATE runs SET status='archived', archived_at=COALESCE(archived_at, ?), "
        "archive_error=NULL, error=NULL, updated_at=? WHERE id=?",
        (utcnow(), utcnow(), row["id"]),
    )
    _remove(published_dir() / row["run_number"])
    _remove(row["staging_dir"])


def _revert_to_published(conn, row, reason: str) -> None:
    """散装文件仍是完整副本:回退到 published,失败原因留档供页面展示。"""
    conn.execute(
        "UPDATE runs SET status='published', archive_error=?, updated_at=? WHERE id=?",
        (reason, utcnow(), row["id"]),
    )
    _remove(row["staging_dir"])
    _remove(archives_dir() / arc.archive_zip_name(row["run_number"]))


def recover_interrupted() -> None:
    """启动恢复:凡未走到稳定状态的运行一律处理掉半成品,
    保证对外只提供完整的散装文件或完整的归档 ZIP,绝不删掉唯一完整副本。"""
    with db_conn(db_path()) as conn:
        stale = conn.execute(
            "SELECT * FROM runs WHERE status = 'processing'"
        ).fetchall()
        for row in stale:
            conn.execute(
                "UPDATE runs SET status='failed', error=?, updated_at=? WHERE id=?",
                ("进程中断,未完成即停止;可用同一运行号重新提交", utcnow(), row["id"]),
            )
        archiving = conn.execute(
            "SELECT * FROM runs WHERE status = 'archiving'"
        ).fetchall()
        archived_rows = conn.execute(
            "SELECT * FROM runs WHERE status = 'archived'"
        ).fetchall()

        # 归档中断:按"现存完整副本"决定完成归档还是安全回退
        for row in archiving:
            zip_path = archives_dir() / arc.archive_zip_name(row["run_number"])
            loose_ok = _published_intact(row["run_number"])
            if _archive_zip_fully_ok(row["run_number"], zip_path):
                _finish_archived(conn, row)
            elif loose_ok:
                _revert_to_published(conn, row, "归档中断,已回退为已发布状态;可重新归档")
            else:
                conn.execute(
                    "UPDATE runs SET status='failed', error=?, "
                    "archive_error=?, updated_at=? WHERE id=?",
                    ("归档中断且归档与散装文件均不完整", "归档中断,无完整副本",
                     utcnow(), row["id"]),
                )
                _remove(row["staging_dir"])
                _remove(zip_path)

        # 已归档行:ZIP 必须完整;散装残留只在 ZIP 确认完整后才清
        for row in archived_rows:
            zip_path = archives_dir() / arc.archive_zip_name(row["run_number"])
            if _archive_zip_fully_ok(row["run_number"], zip_path):
                _finish_archived(conn, row)
            elif _published_intact(row["run_number"]):
                _revert_to_published(conn, row, "归档 ZIP 缺失或损坏,已回退为已发布状态;可重新归档")
            else:
                conn.execute(
                    "UPDATE runs SET status='failed', error=?, updated_at=? WHERE id=?",
                    ("归档 ZIP 缺失或损坏,且散装文件已不存在", utcnow(), row["id"]),
                )
                _remove(zip_path)

        published_nums = {
            r["run_number"]
            for r in conn.execute("SELECT run_number FROM runs WHERE status='published'")
        }
        known_nums = {r["run_number"] for r in conn.execute("SELECT run_number FROM runs")}

    for row in stale:
        _remove(row["staging_dir"])
        _remove(published_dir() / row["run_number"])
    if staging_dir().exists():
        for child in staging_dir().iterdir():  # 无主的暂存目录一并清理
            _remove(child)
    if published_dir().exists():
        for child in published_dir().iterdir():
            if child.name not in known_nums:
                _remove(child)
    if archives_dir().exists():  # 无数据库行引用的归档 ZIP 不对外提供
        for child in archives_dir().iterdir():
            if child.name not in {arc.archive_zip_name(n) for n in known_nums}:
                _remove(child)
    with db_conn(db_path()) as conn:  # 已发布但文件丢失:不允许提供残缺下载
        for num in published_nums:
            if not (published_dir() / num).is_dir():
                conn.execute(
                    "UPDATE runs SET status='failed', error=?, updated_at=? "
                    "WHERE run_number=? AND status='published'",
                    ("已发布文件缺失,运行标记为失败", utcnow(), num),
                )


def setup() -> None:
    staging_dir().mkdir(parents=True, exist_ok=True)
    published_dir().mkdir(parents=True, exist_ok=True)
    archives_dir().mkdir(parents=True, exist_ok=True)
    init_db(db_path())
    recover_interrupted()


@asynccontextmanager
async def lifespan(_app: FastAPI):
    setup()
    yield


app = FastAPI(title="泳道拆样台", lifespan=lifespan)


# ---------------------------------------------------------------- 后台处理

def write_log(path: Path, run, sample_list, stats) -> None:
    lines = [
        f"运行号: {run['run_number']}",
        f"提交时间: {run['created_at']}",
        f"完成时间: {utcnow()}",
        f"输入文件: R1={run['r1_name']} R2={run['r2_name']} "
        f"I1={run['i1_name']} 样本表={run['samples_name']}",
        "样本索引:",
    ]
    lines += [f"  {sid}\t{idx}" for sid, idx in sample_list]
    lines.append(f"总读段对: {stats['total_pairs']}")
    lines.append("分组对数(对数\t其中单碱基纠错):")
    for sid in stats["order"]:
        pairs, corrected = stats["counts"][sid]
        lines.append(f"  {sid}\t{pairs}\t{corrected}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def process_run(run_id: int) -> None:
    with db_conn(db_path()) as conn:
        run = conn.execute("SELECT * FROM runs WHERE id=?", (run_id,)).fetchone()
    if run is None:
        return
    run_number = run["run_number"]
    staging = Path(run["staging_dir"])
    uploads = staging / "uploads"
    out = staging / "out"
    final = published_dir() / run_number
    try:
        sample_list = load_samples(uploads / "samples.csv")
        out.mkdir(parents=True)
        stats = demultiplex(
            str(uploads / "r1.fastq"), str(uploads / "r2.fastq"),
            str(uploads / "i1.fastq"), sample_list, str(out),
        )
        write_log(out / "demux.log", run, sample_list, stats)
        if final.exists():
            raise RuntimeError(f"发布目录已存在: {final}")
        os.rename(out, final)  # 同一数据卷内原子移动:先完整落盘,再发布
        with db_conn(db_path()) as conn:
            updated = conn.execute(
                "UPDATE runs SET status='published', total_pairs=?, published_at=?, "
                "updated_at=? WHERE id=? AND status='processing'",
                (stats["total_pairs"], utcnow(), utcnow(), run_id),
            ).rowcount
            if updated != 1:
                raise RuntimeError("运行状态已变化,放弃发布")
            index_of = dict(sample_list)
            for sid in stats["order"]:
                pairs, corrected = stats["counts"][sid]
                conn.execute(
                    "INSERT INTO sample_counts(run_id, sample_id, index_seq, pairs, "
                    "corrected) VALUES (?,?,?,?,?)",
                    (run_id, sid, index_of.get(sid, ""), pairs, corrected),
                )
        _remove(staging)
    except Exception as exc:  # 任何失败:标失败 + 清半成品,不发布
        with db_conn(db_path()) as conn:
            conn.execute(
                "UPDATE runs SET status='failed', error=?, updated_at=? WHERE id=?",
                (str(exc), utcnow(), run_id),
            )
        _remove(staging)
        _remove(final)


# ---------------------------------------------------------------- 归档处理

def archive_run(run_id: int) -> None:
    """后台归档:核对 -> 暂存打包 -> CRC 校验 -> 原子发布 ZIP -> 清散装文件。

    任何一步失败都把运行安全回退到 published 并保留原始文件,不改计数。
    """
    with db_conn(db_path()) as conn:
        run = conn.execute("SELECT * FROM runs WHERE id=?", (run_id,)).fetchone()
    if run is None or run["status"] != "archiving":
        return
    run_number = run["run_number"]
    run_dir = published_dir() / run_number
    staging = Path(run["staging_dir"])
    final_zip = archives_dir() / arc.archive_zip_name(run_number)
    try:
        with db_conn(db_path()) as conn:
            count_rows = conn.execute(
                "SELECT * FROM sample_counts WHERE run_id=? "
                "ORDER BY sample_id = 'undetermined', sample_id",
                (run_id,),
            ).fetchall()
        # 1) 逐组核对 R1/R2 配对与入库对数(流式,不整份读入)
        arc.verify_run(run_dir, count_rows, run["total_pairs"])

        # 2) 暂存目录里流式打包,对外不可见
        staging.mkdir(parents=True, exist_ok=True)
        tmp_zip = staging / arc.archive_zip_name(run_number)
        entries = arc.build_archive_zip(
            tmp_zip, run_number, run_dir, count_rows, run["total_pairs"]
        )

        # 3) 重新流式校验整个 ZIP(CRC32 + 条目清单)
        arc.verify_zip(tmp_zip, entries)

        # 4) 原子改名到归档区;目标已存在说明是中断残留,且散装文件还在,
        #    先移走再改名,避免覆盖唯一副本
        if final_zip.exists():
            _remove(final_zip)
        os.rename(tmp_zip, final_zip)

        # 5) 状态先落 archived,然后才删散装文件:此刻 ZIP 已是完整副本
        with db_conn(db_path()) as conn:
            updated = conn.execute(
                "UPDATE runs SET status='archived', archived_at=?, "
                "archive_error=NULL, error=NULL, updated_at=? "
                "WHERE id=? AND status='archiving'",
                (utcnow(), utcnow(), run_id),
            ).rowcount
            if updated != 1:
                raise RuntimeError("运行状态已变化,放弃归档清理")
        for row in count_rows:  # 只删双端 FASTQ 与日志;计数与索引已在 ZIP 和库中
            for name in arc.group_files(run_number, row["sample_id"]):
                try:
                    (run_dir / name).unlink()
                except FileNotFoundError:
                    pass
        try:
            (run_dir / "demux.log").unlink()
        except FileNotFoundError:
            pass
        try:
            run_dir.rmdir()  # 目录空了才删;有残留(非本运行文件)则保留
        except OSError:
            pass
        _remove(staging)
    except Exception as exc:  # 失败:回退 published,原运行状态与文件不变
        with db_conn(db_path()) as conn:
            conn.execute(
                "UPDATE runs SET status='published', archive_error=?, updated_at=? "
                "WHERE id=? AND status='archiving'",
                (f"归档失败,散装文件保留: {exc}", utcnow(), run_id),
            )
        _remove(staging)
        if run_dir.is_dir():  # 散装完整副本还在,改名后的 ZIP 才允许清理
            _remove(final_zip)


# ---------------------------------------------------------------- 页面与接口

@app.exception_handler(HTTPException)
async def http_exception_handler(request: Request, exc: HTTPException):
    if "text/html" in request.headers.get("accept", ""):
        return templates.TemplateResponse(
            request, "error.html",
            {"status": exc.status_code, "detail": exc.detail},
            status_code=exc.status_code,
        )
    return JSONResponse({"detail": exc.detail}, status_code=exc.status_code)


@app.get("/", response_class=HTMLResponse)
def index(request: Request):
    with db_conn(db_path()) as conn:
        runs = conn.execute("SELECT * FROM runs ORDER BY id DESC").fetchall()
    return templates.TemplateResponse(
        request, "index.html", {"runs": runs, "labels": STATUS_LABELS}
    )


async def save_upload(upload: UploadFile, dest: Path) -> None:
    with open(dest, "wb") as fh:
        while chunk := await upload.read(1 << 20):
            fh.write(chunk)


@app.post("/runs")
async def create_run(
    run_number: str = Form(...),
    r1: UploadFile = File(...),
    r2: UploadFile = File(...),
    i1: UploadFile = File(...),
    samples: UploadFile = File(...),
):
    run_number = run_number.strip()
    if not RUN_NUMBER_RE.fullmatch(run_number):
        raise HTTPException(400, "运行号格式不正确:字母或数字开头,可含 - _ .,最长 64 字符")
    for label, up in (("R1", r1), ("R2", r2), ("I1", i1), ("样本表", samples)):
        if not up.filename:
            raise HTTPException(400, f"请上传{label}文件")

    token = uuid.uuid4().hex
    staging = staging_dir() / token
    uploads = staging / "uploads"
    names = {k: Path(u.filename or "").name for k, u in
             (("r1", r1), ("r2", r2), ("i1", i1), ("samples", samples))}
    try:
        with db_conn(db_path()) as conn:
            row = conn.execute(
                "SELECT id, status FROM runs WHERE run_number=?", (run_number,)
            ).fetchone()
            if row and row["status"] == "published":
                raise HTTPException(409, f"运行 {run_number} 已发布,不可改写")
            if row and row["status"] == "processing":
                raise HTTPException(409, f"运行 {run_number} 正在处理中,请稍后查看")
            if row:  # failed:允许同号重试
                conn.execute("DELETE FROM runs WHERE id=?", (row["id"],))
            cur = conn.execute(
                "INSERT INTO runs(run_number, status, staging_dir, r1_name, r2_name, "
                "i1_name, samples_name, created_at, updated_at) "
                "VALUES (?,?,?,?,?,?,?,?,?)",
                (run_number, "processing", str(staging), names["r1"], names["r2"],
                 names["i1"], names["samples"], utcnow(), utcnow()),
            )
            run_id = cur.lastrowid
    except sqlite3.IntegrityError:
        raise HTTPException(409, f"运行 {run_number} 已存在,请刷新后查看")

    try:
        uploads.mkdir(parents=True)
        await save_upload(r1, uploads / "r1.fastq")
        await save_upload(r2, uploads / "r2.fastq")
        await save_upload(i1, uploads / "i1.fastq")
        await save_upload(samples, uploads / "samples.csv")
        load_samples(uploads / "samples.csv")  # 先校验样本表,尽早报错
    except HTTPException:
        raise
    except Exception as exc:
        with db_conn(db_path()) as conn:
            conn.execute(
                "UPDATE runs SET status='failed', error=?, updated_at=? WHERE id=?",
                (f"上传或样本表错误: {exc}", utcnow(), run_id),
            )
        _remove(staging)
        raise HTTPException(400, f"上传或样本表错误: {exc}")

    threading.Thread(target=process_run, args=(run_id,), daemon=True).start()
    return RedirectResponse(f"/runs/{run_id}", status_code=303)


def get_run_or_404(run_id: int):
    with db_conn(db_path()) as conn:
        run = conn.execute("SELECT * FROM runs WHERE id=?", (run_id,)).fetchone()
    if run is None:
        raise HTTPException(404, "运行不存在")
    return run


def get_counts(run_id: int):
    with db_conn(db_path()) as conn:
        return conn.execute(
            "SELECT * FROM sample_counts WHERE run_id=? "
            "ORDER BY sample_id = 'undetermined', sample_id",
            (run_id,),
        ).fetchall()


@app.get("/runs/{run_id}", response_class=HTMLResponse)
def run_detail(request: Request, run_id: int):
    run = get_run_or_404(run_id)
    groups = [
        {
            "sample_id": c["sample_id"],
            "index_seq": c["index_seq"],
            "pairs": c["pairs"],
            "corrected": c["corrected"],
            "r1": f"{c['sample_id']}_R1.fastq",
            "r2": f"{c['sample_id']}_R2.fastq",
        }
        for c in get_counts(run_id)
    ]
    archive_name = arc.archive_zip_name(run["run_number"])
    return templates.TemplateResponse(
        request, "run.html",
        {"run": run, "groups": groups, "labels": STATUS_LABELS,
         "archive_name": archive_name},
    )


@app.get("/api/runs/{run_id}")
def run_status(run_id: int):
    run = get_run_or_404(run_id)
    counts = get_counts(run_id)
    return {
        "run_number": run["run_number"],
        "status": run["status"],
        "error": run["error"],
        "archive_error": run["archive_error"],
        "total_pairs": run["total_pairs"],
        "archived_at": run["archived_at"],
        "counts": [
            {"sample_id": c["sample_id"], "index_seq": c["index_seq"],
             "pairs": c["pairs"], "corrected": c["corrected"]}
            for c in counts
        ],
    }


@app.post("/runs/{run_id}/archive")
def request_archive(run_id: int):
    """提交归档:仅 published 可启动;archiving/archived 重复提交幂等。"""
    with db_conn(db_path()) as conn:
        run = conn.execute("SELECT * FROM runs WHERE id=?", (run_id,)).fetchone()
        if run is None:
            raise HTTPException(404, "运行不存在")
        if run["status"] == "archived":  # 重复提交:一份归档,直接回运行页
            return RedirectResponse(f"/runs/{run_id}", status_code=303)
        if run["status"] == "archiving":
            raise HTTPException(409, "该运行正在归档中,请稍后查看")
        if run["status"] != "published":
            raise HTTPException(409, "只有已发布的运行可以归档")
        staging = staging_dir() / f"archive-{uuid.uuid4().hex}"
        updated = conn.execute(
            "UPDATE runs SET status='archiving', archive_error=NULL, "
            "staging_dir=?, updated_at=? WHERE id=? AND status='published'",
            (str(staging), utcnow(), run_id),
        ).rowcount
        if updated != 1:  # 并发提交:只允许一个归档任务
            raise HTTPException(409, "该运行正在归档中,请稍后查看")
    threading.Thread(target=archive_run, args=(run_id,), daemon=True).start()
    return RedirectResponse(f"/runs/{run_id}", status_code=303)


@app.get("/runs/{run_id}/archive")
def download_archive(run_id: int):
    """归档 ZIP 下载:archived 提供完整 ZIP;archiving 期间暂不提供半成品。"""
    run = get_run_or_404(run_id)
    if run["status"] == "archiving":
        raise HTTPException(409, "归档正在进行中,完成后才能下载归档 ZIP")
    if run["status"] != "archived":
        raise HTTPException(404, "该运行尚未归档,没有归档 ZIP")
    path = archives_dir() / arc.archive_zip_name(run["run_number"])
    if not path.is_file():  # 落库与文件不一致时拒绝残缺下载
        raise HTTPException(404, "归档文件不存在")
    return FileResponse(
        path, filename=arc.archive_zip_name(run["run_number"]),
        media_type="application/zip",
    )


@app.get("/runs/{run_id}/download/{filename}")
def download(run_id: int, filename: str):
    run = get_run_or_404(run_id)
    if run["status"] == "archived":
        raise HTTPException(410, "该运行已归档,散装双端文件已停止下载,请下载归档 ZIP")
    if run["status"] == "archiving":
        # 归档进行中:旧下载继续有效,新请求拿到的仍是完整散装文件
        pass
    elif run["status"] != "published":
        raise HTTPException(409, "运行未发布,没有完整文件可下载")
    valid = {"demux.log"}
    for c in get_counts(run_id):
        valid.add(f"{c['sample_id']}_R1.fastq")
        valid.add(f"{c['sample_id']}_R2.fastq")
    if filename not in valid:
        raise HTTPException(404, "文件不存在")
    path = published_dir() / run["run_number"] / filename
    if not path.is_file():
        raise HTTPException(404, "文件不存在")
    return FileResponse(path, filename=filename, media_type="application/octet-stream")
