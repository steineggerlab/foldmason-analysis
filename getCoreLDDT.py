#!/usr/bin/env python3

"""
python3 getCoreLDDT.py msa2lddt.html

Pulls out JSON section, finds gapless columns in aligned sequences,
then computes "core" LDDT from per-column scores of gapless columns
"""

import re
import json
import argparse


def find_gapless_columns(sequences):
    length = len(sequences[0])
    columns = []
    for i in range(length):
        if any(sequence[i] == '-' for sequence in sequences):
            continue
        columns.append(i)
    return columns


def main(html):
    with open(html) as fp:
        content = fp.read()
    match = re.search("({\"entries.*}})</script", content).group(1).strip().replace("\0", "")
    try:
        data = json.loads(match)
    except:
        print(match)
        raise
    count = len(data["entries"])
    columns = find_gapless_columns([entry["aa"] for entry in data["entries"]])
    if len(columns) > 0:
        lddt = sum(data["scores"][idx] for idx in columns) / len(columns)
    else:
        lddt = 0
    *_, family, tool = html.split("/")
    tool = tool.replace(".html", "")
    print(f"{family}\t{tool}\tlddt_core\t{lddt}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("html")
    args = parser.parse_args()
    main(args.html)
