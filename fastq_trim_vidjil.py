# !/usr/bin/python
# !/bin/sh
# !/usr/tcl/bin
# -*- coding: utf-8 -*-
import argparse
import logging as log
import os
import subprocess


def main(args):
    sample_list = ["X32029"]
    primers = {}

    for sample in sample_list:
        fastq1 = f"{sample}_filtered.R1.fastq.gz"
        fastq2 = f"{sample}_filtered.R2.fastq.gz"
        with open(args.primers) as file:
            for line in file:
                line = line.split()
                adapter_name = line[0]
                adapter_fwd = line[1]
                adapter_rev = line[2].split(",")[0]
                adapter_rev2 = line[2].split(",")[1]
                subprocess.run(
                    [
                        "cutadapt",
                        "-j 7",
                        "-g",
                        f"{adapter_fwd};min_overlap={len(adapter_fwd)}",
                        "-g",
                        f"{adapter_fwd};min_overlap={len(adapter_fwd)}",
                        "-G",
                        f"{adapter_rev};min_overlap={len(adapter_rev)}",
                        "-G",
                        f"{adapter_rev2};min_overlap={len(adapter_rev2)}",
                        "--pair-adapters",
                        "-e",
                        "0",
                        "--untrimmed-output",
                        f"{adapter_name}_untrimmed.R1.fastq.gz",
                        "--untrimmed-paired-output",
                        f"{adapter_name}_untrimmed.R2.fastq.gz",
                        "-o",
                        adapter_name + ".R1.fastq.gz",
                        "-p",
                        adapter_name + ".R2.fastq.gz",
                        fastq1,
                        fastq2,
                    ]
                )


def parse_args():
    parser = argparse.ArgumentParser(
        prog="fastq_trim",
        description="Hello",
    )
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        default="*.fastq.gz",
        nargs="+",
        help="Input fastq, formatted in Sample_SX_RX_001.f(ast)q.gz or Sample.RX.f(ast)q.gz, each R1 need its R2",
    )
    parser.add_argument(
        "-p",
        "--primers",
        type=str,
        required=True,
        help="Tabulated file containing 5'->3' primers : name  primer  reverse",
    )
    parser.add_argument("-o", "--output", type=str, default=".", help="Output path")

    args = parser.parse_args()
    main(args)


if __name__ == "__main__":
    parse_args()
