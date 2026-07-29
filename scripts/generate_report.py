#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import json
import os
from typing import Dict, Tuple, List

def read_kv_file(path: str) -> Dict[str, str]:
    d: Dict[str, str] = {}
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return d
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            k, v = (line.split("\t", 1) + [""])[:2]
            d[k] = v
    return d

def safe_float(x: str, default: float = 0.0) -> float:
    try:
        return float(x)
    except Exception:
        return default

def safe_int(x: str, default: int = 0) -> int:
    try:
        return int(float(x))
    except Exception:
        return default

def load_fastp_summary(fastp_json: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    if not fastp_json or not os.path.exists(fastp_json):
        return out
    try:
        with open(fastp_json, "r", encoding="utf-8") as f:
            fp = json.load(f)
        s = fp.get("summary", {}) or {}
        b = s.get("before_filtering", {}) or {}
        a = s.get("after_filtering", {}) or {}
        out["before_reads"] = str(b.get("total_reads", "NA"))
        out["before_bases"] = str(b.get("total_bases", "NA"))
        out["after_reads"] = str(a.get("total_reads", "NA"))
        out["after_bases"] = str(a.get("total_bases", "NA"))
        out["q20_rate"] = str(a.get("q20_rate", "NA"))
        out["q30_rate"] = str(a.get("q30_rate", "NA"))
    except Exception:
        return {}
    return out

def cs_corrected_stats(cs_kv: Dict[str, str]) -> Tuple[int, int, int, int]:
    sub = safe_int(cs_kv.get("sub", "0"))
    ins = safe_int(cs_kv.get("ins", "0"))
    dele = safe_int(cs_kv.get("del", "0"))
    corrected = sub + ins + dele
    return corrected, sub, ins, dele

def pct(x: float) -> str:
    return f"{x*100:.2f}%"

def write_text(path: str, lines: List[str]):
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")

def build_warning_messages(cov: Dict[str, str]) -> Tuple[str, List[str]]:
    """
    Use coverage_qc.py outputs:
      avg_depth, total_pos, low_dp, low_pos, low_frac, zero_pos, zero_frac, th_low_frac, th_zero_frac, status
    """
    status = (cov.get("status", "OK") or "OK").strip()
    avg_depth = safe_float(cov.get("avg_depth", "0"))
    total = safe_int(cov.get("total_pos", "0"))
    low_dp = safe_int(cov.get("low_dp", "10"))
    low_pos = safe_int(cov.get("low_pos", "0"))
    low_frac = safe_float(cov.get("low_frac", "0"))
    zero_pos = safe_int(cov.get("zero_pos", "0"))
    zero_frac = safe_float(cov.get("zero_frac", "0"))
    th_low = safe_float(cov.get("th_low_frac", "0.30"))
    th_zero = safe_float(cov.get("th_zero_frac", "0.05"))

    msgs = []
    msgs.append(f"平均测序深度（Average depth）: {avg_depth:.2f}")
    if total > 0:
        msgs.append(f"深度小于 {low_dp} 的碱基占比: {pct(low_frac)} ({low_pos}/{total})")
        msgs.append(f"深度等于 0 的碱基占比: {pct(zero_frac)} ({zero_pos}/{total})")
    else:
        msgs.append("覆盖度统计失败：total_pos=0（depth 文件为空或比对异常）。")

    # Thresholds in human readable form
    msgs.append(f"预警阈值: 当“深度<{low_dp}的碱基占比 ≥ {pct(th_low)}”或“深度=0的碱基占比 ≥ {pct(th_zero)}”触发预警")

    # Suggestions integrated as bullets (no separate section)
    # Keep them short and actionable
    if status == "OK":
        msgs.append("建议: 覆盖度整体良好；若仍出现序列差异，优先检查 indel 富集区域的 pileup 支持度。")
    else:
        # Distinguish likely causes
        if zero_frac >= th_zero:
            msgs.append("建议: 存在明显零覆盖区，提示可能仅扩增到部分基因组；全长共识在缺失区域不可靠，建议标注为 partial genome。")
        if low_frac >= th_low:
            msgs.append("建议: 大范围低覆盖会使共识更易受噪音/indel 影响；建议提高目标区覆盖或对低覆盖区域做更严格变异过滤后再生成共识。")

    return status, msgs

def build_mask_message(mask: Dict[str, str]) -> str:
    if not mask:
        return ""
    original_len = safe_int(mask.get("original_length", "0"))
    final_len = safe_int(mask.get("final_length", "0"))
    trim_left = safe_int(mask.get("trim_left", "0"))
    trim_right = safe_int(mask.get("trim_right", "0"))
    internal_n = safe_int(mask.get("internal_zero_masked", "0"))
    return (
        "零覆盖处理（Zero-coverage handling）: "
        f"terminal_trim left={trim_left} bp, right={trim_right} bp, "
        f"internal_N={internal_n} bp; length {original_len}->{final_len} bp"
    )

def build_ref_fill_message(ref_fill: Dict[str, str], ref_filled_fasta: str = "") -> str:
    if not ref_fill:
        return ""
    original_len = safe_int(ref_fill.get("original_length", "0"))
    ref_len = safe_int(ref_fill.get("best_ref_length", "0"))
    output_len = safe_int(ref_fill.get("output_length", "0"))
    filled = safe_int(ref_fill.get("filled_zero_with_ref", "0"))
    unfilled = safe_int(ref_fill.get("zero_without_ref", "0"))
    nonzero_kept = safe_int(ref_fill.get("nonzero_kept", "0"))
    message = (
        "参考填补版本 / Best-ref-filled zero-coverage version: "
        f"文件 / FASTA={ref_filled_fasta or 'NA'}; "
        f"用初始 best ref 同坐标碱基填补 depth=0 位点 / depth=0 positions filled from initial best ref={filled}; "
        f"超出 best ref 长度未填补 / zero-depth positions beyond best ref={unfilled}; "
        f"非零覆盖位点保持 consensus / nonzero-depth consensus bases kept={nonzero_kept}; "
        f"长度 / length {original_len}->{output_len} bp; best_ref_length={ref_len} bp"
    )
    return message

def build_ref_fill_cn_lines(ref_fill: Dict[str, str], ref_filled_fasta: str = "") -> List[str]:
    if not ref_fill:
        return ["- 未生成参考填补版本。"]
    original_len = safe_int(ref_fill.get("original_length", "0"))
    ref_len = safe_int(ref_fill.get("best_ref_length", "0"))
    output_len = safe_int(ref_fill.get("output_length", "0"))
    filled = safe_int(ref_fill.get("filled_zero_with_ref", "0"))
    unfilled = safe_int(ref_fill.get("zero_without_ref", "0"))
    nonzero_kept = safe_int(ref_fill.get("nonzero_kept", "0"))
    return [
        f"- FASTA 文件: {ref_filled_fasta or 'NA'}",
        f"- 用初始 best ref 同坐标碱基填补 depth=0 位点: {filled}",
        f"- 超出 best ref 长度而未填补的 depth=0 位点: {unfilled}",
        f"- 保持 consensus 的非零覆盖位点: {nonzero_kept}",
        f"- 序列长度: {original_len}->{output_len} bp；best ref 长度: {ref_len} bp",
        "- 规则: 仅替换 depth=0 的位点；depth>0 的 consensus 位点不因该版本改变。",
    ]

def build_ref_fill_en_lines(ref_fill: Dict[str, str], ref_filled_fasta: str = "") -> List[str]:
    if not ref_fill:
        return ["- No best-reference-filled version was generated."]
    original_len = safe_int(ref_fill.get("original_length", "0"))
    ref_len = safe_int(ref_fill.get("best_ref_length", "0"))
    output_len = safe_int(ref_fill.get("output_length", "0"))
    filled = safe_int(ref_fill.get("filled_zero_with_ref", "0"))
    unfilled = safe_int(ref_fill.get("zero_without_ref", "0"))
    nonzero_kept = safe_int(ref_fill.get("nonzero_kept", "0"))
    return [
        f"- FASTA: {ref_filled_fasta or 'NA'}",
        f"- Depth=0 positions filled from the initial best reference at the same coordinates: {filled}",
        f"- Depth=0 positions not filled because they were beyond the best-reference length: {unfilled}",
        f"- Nonzero-depth consensus positions kept unchanged: {nonzero_kept}",
        f"- Sequence length: {original_len}->{output_len} bp; best-reference length: {ref_len} bp",
        "- Rule: only depth=0 positions are replaced; depth>0 consensus positions are not changed in this version.",
    ]

def build_canu_extension_message(extension: Dict[str, str], extended_fasta: str = "") -> str:
    if not extension:
        return ""
    status = extension.get("status", "unknown")
    left_status = extension.get("left_status", "no_candidate")
    right_status = extension.get("right_status", "no_candidate")
    left_retained = safe_int(extension.get("left_retained_bp", "0"))
    right_retained = safe_int(extension.get("right_retained_bp", "0"))
    message = (
        "Canu 末端扩展 / Canu-supported terminal extension: "
        f"status={status}; left={left_status} retained={left_retained} bp; "
        f"right={right_status} retained={right_retained} bp"
    )
    if status == "generated" and extended_fasta:
        message += f"; FASTA={extended_fasta}"
    return message

def build_canu_extension_end_line(extension: Dict[str, str], end: str, label: str) -> str:
    status = extension.get(f"{end}_status", "no_candidate")
    retained = safe_int(extension.get(f"{end}_retained_bp", "0"))
    contig = extension.get(f"{end}_contig", "")
    if not contig:
        return f"- {label}: {status}; retained={retained} bp"
    anchor = safe_int(extension.get(f"{end}_anchor_bp", "0"))
    identity = extension.get(f"{end}_identity", "NA")
    mapq = safe_int(extension.get(f"{end}_mapq", "0"))
    candidate = safe_int(extension.get(f"{end}_candidate_overhang_bp", "0"))
    return (
        f"- {label}: {status}; contig={contig}; anchor={anchor} bp; identity={identity}; "
        f"MAPQ={mapq}; candidate={candidate} bp; retained={retained} bp"
    )

def build_canu_extension_cn_lines(extension: Dict[str, str], extended_fasta: str = "") -> List[str]:
    status = extension.get("status", "unknown")
    lines = [f"- 状态: {status}"]
    if status == "generated" and extended_fasta:
        lines.append(f"- 独立 FASTA 文件: {extended_fasta}")
    else:
        lines.append("- 未生成独立 FASTA 文件。")
    lines.append(build_canu_extension_end_line(extension, "left", "左端"))
    lines.append(build_canu_extension_end_line(extension, "right", "右端"))
    lines.append(
        f"- 验证规则: 仅保留与原共识边界连续相接且 depth>={safe_int(extension.get('min_validation_depth', '3'))} 的延伸碱基。"
    )
    return lines

def build_canu_extension_en_lines(extension: Dict[str, str], extended_fasta: str = "") -> List[str]:
    status = extension.get("status", "unknown")
    lines = [f"- Status: {status}"]
    if status == "generated" and extended_fasta:
        lines.append(f"- Separate FASTA: {extended_fasta}")
    else:
        lines.append("- No separate FASTA was generated.")
    lines.append(build_canu_extension_end_line(extension, "left", "Left"))
    lines.append(build_canu_extension_end_line(extension, "right", "Right"))
    lines.append(
        f"- Validation rule: retain only extension bases continuously adjacent to the primary consensus with depth>={safe_int(extension.get('min_validation_depth', '3'))}."
    )
    return lines

def build_degenerate_message(degen: Dict[str, str]) -> str:
    if not degen:
        return ""
    total = safe_int(degen.get("degenerate_total", "0"))
    if total == 0:
        return (
            "简并碱基处理 / IUPAC degenerate-base handling: "
            "未检测到 R/Y/S/W/K/M/B/D/H/V 简并碱基，未进行替换 / "
            "no R/Y/S/W/K/M/B/D/H/V degenerate bases detected; no replacement was needed"
        )
    resolved_acgt = safe_int(degen.get("resolved_to_acgt", "0"))
    resolved_iupac = safe_int(degen.get("resolved_to_iupac", "0"))
    masked = safe_int(degen.get("masked_to_n", "0"))
    no_pileup = safe_int(degen.get("no_pileup", "0"))
    incompatible_top = safe_int(degen.get("incompatible_top", "0"))
    tied_top = safe_int(degen.get("tied_top", "0"))
    unrepresentable_tie = safe_int(degen.get("unrepresentable_tie", "0"))
    return (
        "简并碱基处理（IUPAC degenerate-base handling）: "
        f"total={total}, resolved_ACGT={resolved_acgt}, resolved_IUPAC={resolved_iupac}, masked_N={masked}; "
        "rule=use observed sequencing bases; only zero A/C/G/T observations become N; "
        f"notes no_pileup={no_pileup}, incompatible_top={incompatible_top}, "
        f"tied_top={tied_top}, unrepresentable_tie={unrepresentable_tie}"
    )

def build_indel_message(indel: Dict[str, str], indel_tsv: str = "") -> Tuple[bool, str]:
    if not indel:
        return False, ""
    flagged = safe_int(indel.get("flagged_sites", "0"))
    low_depth = safe_int(indel.get("low_depth", "0"))
    del_dominant = safe_int(indel.get("del_dominant", "0"))
    mixed_base = safe_int(indel.get("mixed_base", "0"))
    ins_supported = safe_int(indel.get("ins_supported", "0"))
    del_event = safe_int(indel.get("del_event_supported", "0"))
    message = (
        "Indel QC（indel/mixed-base warning）: "
        f"flagged_sites={flagged}, del_dominant={del_dominant}, "
        f"del_event_supported={del_event}, ins_supported={ins_supported}, "
        f"mixed_base={mixed_base}, low_depth={low_depth}; "
        "rule=report only, consensus sequence is not changed by this QC step"
    )
    if indel_tsv:
        message += f"; detail_tsv={indel_tsv}"
    if flagged > 0:
        message += "; suggestion=review flagged indel-complex or mixed-base loci before downstream use"
    else:
        message += "; no flagged indel-complex sites"
    return flagged > 0, message

def main():
    ap = argparse.ArgumentParser(description="Generate bilingual PRRSV pipeline report and warning.txt")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--input-fastq", required=True)
    ap.add_argument("--final-fasta", required=True)
    ap.add_argument("--best-ref-id", required=True)
    ap.add_argument("--fastp-html", required=True)
    ap.add_argument("--fastp-json", default="")
    ap.add_argument("--coverage-kv", required=True)
    ap.add_argument("--cs-round1-kv", required=True)
    ap.add_argument("--cs-round2-kv", required=True)
    ap.add_argument("--mask-stats-kv", default="")
    ap.add_argument("--degenerate-stats-kv", default="")
    ap.add_argument("--indel-qc-kv", default="")
    ap.add_argument("--indel-qc-tsv", default="")
    ap.add_argument("--ref-fill-stats-kv", default="")
    ap.add_argument("--ref-filled-fasta", default="")
    ap.add_argument("--canu-extension-stats-kv", default="")
    ap.add_argument("--canu-extended-fasta", default="")
    ap.add_argument("--out-cn", required=True)
    ap.add_argument("--out-en", required=True)
    ap.add_argument("--out-warning", required=True)
    args = ap.parse_args()

    cov = read_kv_file(args.coverage_kv)
    cs1 = read_kv_file(args.cs_round1_kv)
    cs2 = read_kv_file(args.cs_round2_kv)
    mask = read_kv_file(args.mask_stats_kv)
    degen = read_kv_file(args.degenerate_stats_kv)
    indel = read_kv_file(args.indel_qc_kv)
    ref_fill = read_kv_file(args.ref_fill_stats_kv)
    canu_extension = read_kv_file(args.canu_extension_stats_kv)
    fp = load_fastp_summary(args.fastp_json)

    status, warn_msgs = build_warning_messages(cov)
    mask_msg = build_mask_message(mask)
    if mask_msg:
        warn_msgs.append(mask_msg)
    ref_fill_msg = build_ref_fill_message(ref_fill, args.ref_filled_fasta)
    canu_extension_msg = build_canu_extension_message(canu_extension, args.canu_extended_fasta)
    degen_msg = build_degenerate_message(degen)
    if degen_msg:
        warn_msgs.append(degen_msg)
    _, indel_msg = build_indel_message(indel, args.indel_qc_tsv)

    corr1, sub1, ins1, del1 = cs_corrected_stats(cs1)
    corr2, sub2, ins2, del2 = cs_corrected_stats(cs2)

    # warning.txt (keep both concise + readable)
    warn_lines = [f"Sample\t{args.sample}", f"Status\t{status}"]
    for m in warn_msgs:
        warn_lines.append(f"Message\t{m}")
    if ref_fill_msg:
        warn_lines.append(f"Message\t{ref_fill_msg}")
    if canu_extension_msg:
        warn_lines.append(f"Message\t{canu_extension_msg}")
    write_text(args.out_warning, warn_lines)

    # CN report
    cn = []
    cn.append("PRRSV 自动组装与共识流程报告（中文）")
    cn.append("=" * 72)
    cn.append(f"样本ID: {args.sample}")
    cn.append(f"输入FASTQ: {args.input_fastq}")
    cn.append(f"最终共识序列: {args.final_fasta}")
    cn.append(f"参考序列ID: {args.best_ref_id}")
    cn.append("")
    cn.append("一、fastp 报告")
    cn.append(f"- HTML 报告: {args.fastp_html}")
    if fp:
        cn.append(f"- 过滤前 reads/bases: {fp.get('before_reads','NA')} / {fp.get('before_bases','NA')}")
        cn.append(f"- 过滤后 reads/bases: {fp.get('after_reads','NA')} / {fp.get('after_bases','NA')}")
        cn.append(f"- Q20/Q30 比例(过滤后): {fp.get('q20_rate','NA')} / {fp.get('q30_rate','NA')}")
    else:
        cn.append("- fastp.json 未提供或解析失败（不影响主流程）。")
    cn.append("")
    cn.append("二、覆盖度分布与临床预警")
    cn.append(f"- 预警状态: {status}")
    for m in warn_msgs:
        cn.append(f"- {m}")
    if indel_msg:
        cn.append(f"- Indel QC 复核信息（不影响预警状态）: {indel_msg}")
    cn.append("")
    cn_section_labels = {3: "三", 4: "四", 5: "五"}
    cn_section = 3
    if ref_fill:
        cn.append(f"{cn_section_labels[cn_section]}、参考填补版本")
        for m in build_ref_fill_cn_lines(ref_fill, args.ref_filled_fasta):
            cn.append(m)
        cn.append("")
        cn_section += 1
    if canu_extension:
        cn.append(f"{cn_section_labels[cn_section]}、Canu 末端扩展（独立版本）")
        for m in build_canu_extension_cn_lines(canu_extension, args.canu_extended_fasta):
            cn.append(m)
        cn.append("")
        cn_section += 1
    cn.append(f"{cn_section_labels[cn_section]}、bcftools 共识矫正统计（近似，基于 minimap2 cs 标签）")
    cn.append(f"- Round1（ref -> round1）: corrected={corr1} (sub={sub1}, ins_bases={ins1}, del_bases={del1})")
    cn.append(f"- Round2（round1 -> round2）: corrected={corr2} (sub={sub2}, ins_bases={ins2}, del_bases={del2})")
    write_text(args.out_cn, cn)

    # EN report (match the same “non-confusing” style)
    en = []
    en.append("PRRSV Auto-assembly & Consensus Report (English)")
    en.append("=" * 72)
    en.append(f"Sample: {args.sample}")
    en.append(f"Input FASTQ: {args.input_fastq}")
    en.append(f"Final consensus FASTA: {args.final_fasta}")
    en.append(f"Selected reference ID: {args.best_ref_id}")
    en.append("")
    en.append("1) fastp report")
    en.append(f"- HTML report: {args.fastp_html}")
    if fp:
        en.append(f"- Before filtering reads/bases: {fp.get('before_reads','NA')} / {fp.get('before_bases','NA')}")
        en.append(f"- After filtering reads/bases: {fp.get('after_reads','NA')} / {fp.get('after_bases','NA')}")
        en.append(f"- Q20/Q30 rate (after): {fp.get('q20_rate','NA')} / {fp.get('q30_rate','NA')}")
    else:
        en.append("- fastp.json not provided or failed to parse (pipeline still valid).")
    en.append("")
    en.append("2) Coverage distribution & clinical warning")
    en.append(f"- Status: {status}")
    # Provide English equivalents for the two key ratios while keeping the CN text understandable if mixed
    # We reuse the same messages (they include English average depth already).
    for m in warn_msgs:
        en.append(f"- {m}")
    if indel_msg:
        en.append(f"- Indel QC review information (does not affect warning status): {indel_msg}")
    en.append("")
    en_section = 3
    if ref_fill:
        en.append(f"{en_section}) Best-reference-filled zero-coverage version")
        for m in build_ref_fill_en_lines(ref_fill, args.ref_filled_fasta):
            en.append(m)
        en.append("")
        en_section += 1
    if canu_extension:
        en.append(f"{en_section}) Canu-supported terminal extension (separate version)")
        for m in build_canu_extension_en_lines(canu_extension, args.canu_extended_fasta):
            en.append(m)
        en.append("")
        en_section += 1
    en.append(f"{en_section}) Estimated corrected bases after bcftools consensus (approx. via minimap2 cs-tag)")
    en.append(f"- Round1 (ref -> round1): corrected={corr1} (sub={sub1}, ins_bases={ins1}, del_bases={del1})")
    en.append(f"- Round2 (round1 -> round2): corrected={corr2} (sub={sub2}, ins_bases={ins2}, del_bases={del2})")
    write_text(args.out_en, en)

if __name__ == "__main__":
    main()
