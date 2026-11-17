"""归档核心:逐组核对 R1/R2 配对与入库对数,流式打包并校验 ZIP。

可靠性约束:
- 打包与校验都在调用方提供的暂存目录里进行,对外不可见;成功后由调用方
  原子改名到归档区,任何异常都不改变已发布的原始文件。
- 全程流式:FASTQ 逐记录迭代、ZIP 按块写入,大文件不整份读入内存。
- 校验三层:
  1. 每组 R1/R2 文件齐全;
  2. R1/R2 逐对读段号一致、结构完整,且记录数与数据库入库对数一致;
  3. 打包后重新流式读取 ZIP 全部条目做 CRC32 校验,并核对条目清单。
"""
from __future__ import annotations

import csv
import io
import zipfile

from .demux import FastqFormatError, iter_fastq, read_id

CHUNK = 1 << 20  # 1 MiB:ZIP 内大文件按块搬运的缓冲大小


class ArchiveError(Exception):
    """归档前置核对、打包或 ZIP 校验失败。"""


def group_files(run_number: str, sample_id: str) -> tuple[str, str]:
    return f"{sample_id}_R1.fastq", f"{sample_id}_R2.fastq"


def verify_group(run_dir, sample_id: str, expected_pairs: int) -> None:
    """核对单个分组:R1/R2 齐全、记录数等于入库对数、逐对读段号一致。"""
    r1_name, r2_name = group_files(run_dir, sample_id)
    r1_path, r2_path = run_dir / r1_name, run_dir / r2_name
    if not r1_path.is_file():
        raise ArchiveError(f"分组 {sample_id}:缺少 R1 文件 {r1_name}")
    if not r2_path.is_file():
        raise ArchiveError(f"分组 {sample_id}:缺少 R2 文件 {r2_name}")
    pairs = 0
    try:
        for rec1, rec2 in zip(
            iter_fastq(str(r1_path), f"{sample_id}/R1"),
            iter_fastq(str(r2_path), f"{sample_id}/R2"),
            strict=True,
        ):
            pairs += 1
            id1 = read_id(rec1[0])
            id2 = read_id(rec2[0])
            if not id1 or id1 != id2:
                raise ArchiveError(
                    f"分组 {sample_id}:第 {pairs} 对读段号不一致 R1={rec1[0]!r} R2={rec2[0]!r}"
                )
    except FastqFormatError as exc:
        raise ArchiveError(f"分组 {sample_id}:FASTQ 结构损坏:{exc}") from exc
    except ValueError as exc:  # zip(strict=True) 记录数不一致
        raise ArchiveError(f"分组 {sample_id}:R1/R2 记录数不一致") from exc
    if pairs != expected_pairs:
        raise ArchiveError(
            f"分组 {sample_id}:文件对数 {pairs} 与入库对数 {expected_pairs} 不一致"
        )


def samples_csv_text(count_rows) -> str:
    """样本索引:不含未分配组;首行表头,与上传样本表同格式。"""
    buf = io.StringIO()
    writer = csv.writer(buf, lineterminator="\n")
    writer.writerow(["sample_id", "index"])
    for row in count_rows:
        if row["sample_id"] != "undetermined":
            writer.writerow([row["sample_id"], row["index_seq"]])
    return buf.getvalue()


def counts_csv_text(count_rows, total_pairs: int) -> str:
    """计数清单:每组对数与纠错对数,外加总对数,供日后追查。"""
    buf = io.StringIO()
    writer = csv.writer(buf, lineterminator="\n")
    writer.writerow(["sample_id", "index_seq", "pairs", "corrected"])
    for row in count_rows:
        writer.writerow([row["sample_id"], row["index_seq"],
                         row["pairs"], row["corrected"]])
    writer.writerow(["TOTAL", "", total_pairs, ""])
    return buf.getvalue()


def build_archive_zip(
    zip_path, run_number: str, run_dir, count_rows, total_pairs: int
) -> list[str]:
    """把全部双端文件、样本索引、计数与日志流式装入 zip_path。

    返回 ZIP 内条目名清单(已按名排序)。调用前应已通过 verify_run 核对。
    """
    prefix = f"{run_number}/"
    entries: list[str] = []
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED, allowZip64=True) as zf:
        # 小的元数据直接写入字符串;大 FASTQ 按块搬运,不整份入内存
        for name, text in (
            ("samples.csv", samples_csv_text(count_rows)),
            ("counts.csv", counts_csv_text(count_rows, total_pairs)),
        ):
            arc = prefix + name
            zf.writestr(arc, text.encode("utf-8"))
            entries.append(arc)

        fastq_names: list[str] = []
        for row in count_rows:
            for name in group_files(run_number, row["sample_id"]):
                src = run_dir / name
                if not src.is_file():  # 兜底:verify_run 已查过,这里防并发缺失
                    raise ArchiveError(f"打包时文件消失:{name}")
                fastq_names.append(name)
        for name in fastq_names:
            arc = prefix + name
            entries.append(arc)
            with open(run_dir / name, "rb") as src, zf.open(arc, "w", force_zip64=True) as dst:
                while True:
                    chunk = src.read(CHUNK)
                    if not chunk:
                        break
                    dst.write(chunk)

        log_name = "demux.log"
        log_path = run_dir / log_name
        if not log_path.is_file():
            raise ArchiveError("缺少处理日志 demux.log")
        arc = prefix + log_name
        entries.append(arc)
        with open(log_path, "rb") as src, zf.open(arc, "w", force_zip64=True) as dst:
            while True:
                chunk = src.read(CHUNK)
                if not chunk:
                    break
                dst.write(chunk)
    return sorted(entries)


def verify_zip(zip_path, expected_entries: list[str]) -> None:
    """重新流式读取整个 ZIP:中央目录齐全、每条目 CRC32 通过。"""
    try:
        with zipfile.ZipFile(zip_path, "r") as zf:
            bad = zf.testzip()  # 逐条目流式校验 CRC,返回第一个坏条目或 None
            if bad is not None:
                raise ArchiveError(f"归档 ZIP 校验失败,损坏条目: {bad}")
            names = sorted(zf.namelist())
    except zipfile.BadZipFile as exc:
        raise ArchiveError(f"归档 ZIP 已损坏: {exc}") from exc
    if names != sorted(expected_entries):
        raise ArchiveError(
            f"归档 ZIP 条目清单不一致: 期望 {sorted(expected_entries)}, 实际 {names}"
        )


def verify_run(run_dir, count_rows, total_pairs: int) -> int:
    """归档前完整核对:逐组 R1/R2 配对与入库对数、日志齐全、总对数闭合。"""
    if not run_dir.is_dir():
        raise ArchiveError(f"发布目录不存在: {run_dir}")
    summed = 0
    for row in count_rows:
        verify_group(run_dir, row["sample_id"], row["pairs"])
        summed += row["pairs"]
    if summed != total_pairs:
        raise ArchiveError(
            f"各组对数之和 {summed} 与运行总对数 {total_pairs} 不一致"
        )
    if not (run_dir / "demux.log").is_file():
        raise ArchiveError("缺少处理日志 demux.log")
    return summed


def archive_zip_name(run_number: str) -> str:
    return f"{run_number}.zip"
