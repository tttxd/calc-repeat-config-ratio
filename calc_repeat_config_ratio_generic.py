#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function

import argparse
import csv
import re
import subprocess
from collections import defaultdict

REF_CONSUME = set(["M", "D", "N", "=", "X"])


def normalize_ref_name(name):
    """
    统一参考序列名称：
    1. 只保留空格前第一段
    2. 去掉末尾 .fa / .fasta / .fna
    """
    name = name.strip().split()[0]
    name = re.sub(r"\.(fa|fasta|fna)$", "", name, flags=re.IGNORECASE)
    return name


def parse_ref_list(text):
    """
    将逗号分隔的路径名转成列表，并规范化名称
    """
    refs = []
    for x in text.split(","):
        x = x.strip()
        if x:
            refs.append(normalize_ref_name(x))
    return refs


def cigar_ref_len(cigar):
    """
    根据 CIGAR 计算 alignment 在 reference 上覆盖的长度
    只统计消耗 reference 的操作：M, D, N, =, X
    """
    if cigar == "*" or (not cigar):
        return 0
    total = 0
    for length, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        if op in REF_CONSUME:
            total += int(length)
    return total


def main():
    parser = argparse.ArgumentParser(
        description="Generic script for estimating proportions of two repeat-mediated configurations from long-read BAM using samtools only."
    )
    parser.add_argument("-b", "--bam", required=True, help="Input BAM file")
    parser.add_argument("--config1", required=True,
                        help="Comma-separated path names for configuration 1, e.g. pathA,pathB")
    parser.add_argument("--config2", required=True,
                        help="Comma-separated path names for configuration 2, e.g. pathC,pathD")
    parser.add_argument("--repeat-len", type=int, required=True, help="Repeat length")
    parser.add_argument("--left-flank", type=int, required=True, help="Left flank length")
    parser.add_argument("--right-flank", type=int, required=True, help="Right flank length")
    parser.add_argument("--anchors", default="20,50,100",
                        help="Comma-separated anchor thresholds, default: 20,50,100")
    parser.add_argument("--mapq", type=int, default=30, help="Minimum MAPQ, default: 30")
    parser.add_argument("--samtools", default="samtools", help="Path to samtools, default: samtools")
    parser.add_argument("--label", default="repeat", help="Label for this repeat, default: repeat")
    parser.add_argument("--out", default="repeat_ratio_multi_anchor.tsv", help="Output TSV file")
    args = parser.parse_args()

    config1_refs = parse_ref_list(args.config1)
    config2_refs = parse_ref_list(args.config2)

    if len(config1_refs) == 0:
        raise ValueError("`--config1` is empty")
    if len(config2_refs) == 0:
        raise ValueError("`--config2` is empty")

    overlap = set(config1_refs) & set(config2_refs)
    if overlap:
        raise ValueError("The same reference(s) appear in both config1 and config2: %s" %
                         ",".join(sorted(overlap)))

    all_target_refs = set(config1_refs + config2_refs)

    anchors = sorted(set([int(x) for x in args.anchors.split(",") if x.strip()]))
    if len(anchors) == 0:
        raise ValueError("No valid anchors provided")

    # 1-based repeat boundaries on the path reference
    repeat_start_1 = args.left_flank + 1
    repeat_end_1 = args.left_flank + args.repeat_len

    # 每个 anchor 下：read -> 支持了哪些路径
    read_to_refs_by_anchor = {}
    for a in anchors:
        read_to_refs_by_anchor[a] = defaultdict(set)

    cmd = [
        args.samtools, "view",
        "-q", str(args.mapq),
        "-F", "2308",   # 4(unmapped) + 256(secondary) + 2048(supplementary)
        args.bam
    ]

    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    total_lines = 0
    used_lines = 0

    for raw in p.stdout:
        if not raw:
            continue
        if not isinstance(raw, str):
            raw = raw.decode("utf-8", "replace")
        line = raw.rstrip("\n")
        if not line:
            continue

        total_lines += 1
        fields = line.split("\t")
        if len(fields) < 6:
            continue

        qname = fields[0]
        rname = normalize_ref_name(fields[2])
        pos1 = int(fields[3])   # 1-based inclusive
        cigar = fields[5]

        if rname not in all_target_refs:
            continue

        ref_len = cigar_ref_len(cigar)
        if ref_len <= 0:
            continue

        end1 = pos1 + ref_len - 1   # 1-based inclusive

        # 必须覆盖整个 repeat 区间
        if pos1 > repeat_start_1:
            continue
        if end1 < repeat_end_1:
            continue

        # 两侧 flank 实际支持长度
        left_supported = args.left_flank - (pos1 - 1)
        right_supported = end1 - repeat_end_1
        max_supported_anchor = min(left_supported, right_supported)

        if max_supported_anchor < min(anchors):
            continue

        used_lines += 1

        for a in anchors:
            if max_supported_anchor >= a:
                read_to_refs_by_anchor[a][qname].add(rname)

    stderr_txt = p.stderr.read()
    if stderr_txt and not isinstance(stderr_txt, str):
        stderr_txt = stderr_txt.decode("utf-8", "replace")
    ret = p.wait()
    if ret != 0:
        raise RuntimeError("samtools view failed:\n%s" % stderr_txt)

    # 动态表头
    header = [
        "label",
        "anchor",
        "mapq_cutoff",
        "repeat_len",
        "left_flank",
        "right_flank",
    ]

    ordered_refs = config1_refs + config2_refs
    header.extend(ordered_refs)

    header.extend([
        "config1_support",
        "config2_support",
        "config1_proportion",
        "config2_proportion",
        "config1_to_config2_ratio",
        "ambiguous_reads_removed",
        "unique_supporting_reads_total",
    ])

    with open(args.out, "w") as fout:
        writer = csv.writer(fout, delimiter="\t")
        writer.writerow(header)

        for a in anchors:
            per_read = read_to_refs_by_anchor[a]
            unique_counts = defaultdict(int)
            ambiguous_reads = 0

            for qname, refs in per_read.items():
                if len(refs) == 1:
                    ref = list(refs)[0]
                    unique_counts[ref] += 1
                else:
                    ambiguous_reads += 1

            config1_support = sum([unique_counts[x] for x in config1_refs])
            config2_support = sum([unique_counts[x] for x in config2_refs])
            total_support = config1_support + config2_support

            if total_support > 0:
                p1 = float(config1_support) / total_support
                p2 = float(config2_support) / total_support
                ratio = "Inf" if config2_support == 0 else ("%.10f" % (float(config1_support) / config2_support))
            else:
                p1 = 0.0
                p2 = 0.0
                ratio = "NA"

            row = [
                args.label,
                a,
                args.mapq,
                args.repeat_len,
                args.left_flank,
                args.right_flank,
            ]

            for ref in ordered_refs:
                row.append(unique_counts[ref])

            row.extend([
                config1_support,
                config2_support,
                "%.10f" % p1,
                "%.10f" % p2,
                ratio,
                ambiguous_reads,
                total_support,
            ])

            writer.writerow(row)

    print("Done.")
    print("Label:", args.label)
    print("Total SAM records scanned:", total_lines)
    print("Primary records passing repeat-span filter:", used_lines)
    print("Output:", args.out)

    print("\nSummary:")
    for a in anchors:
        per_read = read_to_refs_by_anchor[a]
        unique_counts = defaultdict(int)
        ambiguous_reads = 0

        for qname, refs in per_read.items():
            if len(refs) == 1:
                ref = list(refs)[0]
                unique_counts[ref] += 1
            else:
                ambiguous_reads += 1

        config1_support = sum([unique_counts[x] for x in config1_refs])
        config2_support = sum([unique_counts[x] for x in config2_refs])
        total_support = config1_support + config2_support

        if total_support > 0:
            p1 = float(config1_support) / total_support
            p2 = float(config2_support) / total_support
        else:
            p1 = 0.0
            p2 = 0.0

        print("anchor=%s: config1=%s, config2=%s, P(config1)=%.6f, P(config2)=%.6f, ambiguous_removed=%s"
              % (a, config1_support, config2_support, p1, p2, ambiguous_reads))


if __name__ == "__main__":
    main()
