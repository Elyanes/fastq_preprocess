# !/usr/bin/python
# !/bin/sh
# !/usr/tcl/bin
# -*- coding: utf-8 -*-
import argparse
import os
import subprocess
import glob
import shutil


def get_samples(args):
    R1_files = []
    R2_files = []
    sample_names = []
    for i in args.input:
        sample_name = os.path.basename(i).split(".")[0]
        if sample_name not in sample_names:
            sample_names.append(sample_name)
        if "R1" in i:
            R1_files.append(i)
        if "R2" in i:
            R2_files.append(i)

    return (R1_files, R2_files, sample_names)


def filter_fastq(args, sample_list):
    sample_names = []
    filtered_R1 = []
    filtered_R2 = []
    loop = 0
    for i in sample_list[2]:
        if not os.path.isdir(os.path.join(args.output, i)):
            os.mkdir(os.path.join(args.output, i))
        subprocess.run(
            [
                "fastp",
                "-q",
                str(args.quality_filter),
                "-A",
                "-i",
                sample_list[0][loop],
                "-I",
                sample_list[1][loop],
                "-o",
                os.path.join(args.output, i, i + "_filtered.R1.fastq.gz"),
                "-O",
                os.path.join(args.output, i, i + "_filtered.R2.fastq.gz"),
                "-h",
                os.path.join(args.output, i, i + "_fastp_log.html"),
                "-j",
                os.path.join(args.output, i, i + "_fastp_log.json"),
            ]
        )
        loop += 1
        sample_names.append(i)
        filtered_R1.append(os.path.join(args.output, i, i + "_filtered.R1.fastq.gz"))
        filtered_R2.append(os.path.join(args.output, i, i + "_filtered.R2.fastq.gz"))

    return (filtered_R1, filtered_R2, sample_names)


def concat_results(args, adapter_names, sample_names):
    for i in sample_names:
        concatenate_result_R1 = os.path.join(
            args.output, i, f"{i}_processed.R1.fastq.gz"
        )
        concatenate_result_R2 = os.path.join(
            args.output, i, f"{i}_processed.R2.fastq.gz"
        )
        R1_files = []
        R2_files = []
        for j in adapter_names:
            R1_files.append(os.path.join(args.output, i, f"{i}_{j}.R1.fastq.gz"))
            R2_files.append(os.path.join(args.output, i, f"{i}_{j}.R2.fastq.gz"))
        with open(concatenate_result_R1, "wb") as write_file:
            for R1_file in R1_files:
                with open(R1_file, "rb") as read_file:
                    shutil.copyfileobj(read_file, write_file)
                os.remove(R1_file)
        with open(concatenate_result_R2, "wb") as write_file:
            for R2_file in R2_files:
                with open(R2_file, "rb") as read_file:
                    shutil.copyfileobj(read_file, write_file)
                os.remove(R2_file)


def main(args):
    sample_list = get_samples(args)
    filtered_fastq = filter_fastq(args, sample_list)
    loop = 0
    adapter_names = []

    with open(args.primers) as file:
        for line in file:
            cmd_adapter_fwd = ""
            cmd_adapter_rev = ""
            line = line.split()
            adapter_name = line[0]
            adapter_names.append(adapter_name)

            adapters_fwd = line[1].split(",")
            for adapter_fwd in adapters_fwd:
                cmd_adapter_fwd = (
                    cmd_adapter_fwd
                    + " -g, "
                    + adapter_fwd
                    + ";min_overlap="
                    + str(len(adapter_fwd))
                    + ","
                )
            adapters_rev = line[2].split(",")
            for adapter_rev in adapters_rev:
                cmd_adapter_rev = (
                    cmd_adapter_rev
                    + " -G, "
                    + adapter_rev
                    + ";min_overlap="
                    + str(len(adapter_rev))
                    + ","
                )

    for sample in filtered_fastq[2]:
        R1_fastq = filtered_fastq[0][loop]
        R2_fastq = filtered_fastq[1][loop]
        for adapter_name in adapter_names:
            subprocess.run(
                [
                    "cutadapt",
                    "-j 7",
                    "-g",
                    f"{adapter_fwd};min_overlap={len(adapter_fwd)}",
                    "-G",
                    f"{adapter_rev};min_overlap={len(adapter_rev)}",
                    "--pair-adapters",
                    "-e",
                    "0",
                    "--untrimmed-output",
                    f"{args.output}/{sample}/{sample}_{adapter_name}_untrimmed.R1.fastq.gz",
                    "--untrimmed-paired-output",
                    f"{args.output}/{sample}/{sample}_{adapter_name}_untrimmed.R2.fastq.gz",
                    "-o",
                    f"{args.output}/{sample}/{sample}_{adapter_name}.R1.fastq.gz",
                    "-p",
                    f"{args.output}/{sample}/{sample}_{adapter_name}.R2.fastq.gz",
                    f"{R1_fastq}",
                    f"{R2_fastq}",
                ]
            )
            if args.keep_untrimmed is False:
                untrimmed_files = glob.glob(
                    f"{args.output}/{sample}/{sample}_{adapter_name}_untrimmed.*"
                )
                for file in untrimmed_files:
                    os.remove(file)
        loop += 1

    concat_results(args, adapter_names, filtered_fastq[2])


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
        help="Tabulated file containing 5'->3' primers in this order : name  primer  reverse",
    )
    parser.add_argument(
        "-qf",
        "--quality_filter",
        type=int,
        default=0,
        help="Specify a quality value to filter your fastq with fastp, default is 0",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default=".",
        help="Output path, by default is the launch folder",
    )
    parser.add_argument(
        "-ku",
        "--keep_untrimmed",
        action="store_true",
        help="If used, keeping untrimmed reads in specific files",
    )

    args = parser.parse_args()
    main(args)


if __name__ == "__main__":
    parse_args()
