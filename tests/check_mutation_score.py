"""Fail if the mutation score drops below a threshold (CI ratchet).

Usage (after `mutmut run` and `mutmut export-cicd-stats`):
    python tests/check_mutation_score.py --min 68

Score = (killed + timeout) / (total - skipped). Writes a Markdown summary to
$GITHUB_STEP_SUMMARY when running in GitHub Actions.
"""
import argparse
import collections
import json
import os
import subprocess
import sys


def surviving_by_function():
    out = subprocess.run([sys.executable, "-m", "mutmut", "results"],
                         capture_output=True, text=True).stdout
    counts = collections.Counter()
    for line in out.splitlines():
        name, _, status = line.strip().rpartition(": ")
        if status == "survived":
            counts[name.split("__mutmut_")[0].replace(".x_", ".")] += 1
    return counts


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--min", type=float, required=True, help="minimal score in percent")
    p.add_argument("--stats", default="mutants/mutmut-cicd-stats.json")
    args = p.parse_args()

    s = json.load(open(args.stats))
    tested = s["total"] - s["skipped"]
    killed = s["killed"] + s["timeout"]
    score = 100.0 * killed / tested if tested else 0.0
    ok = score >= args.min

    lines = [
        "## Mutation testing",
        "",
        f"**Score: {score:.1f}%** (minimum {args.min:.1f}%) {'✅' if ok else '❌'}",
        "",
        "| killed | timeout | survived | no tests | suspicious | total |",
        "|---:|---:|---:|---:|---:|---:|",
        f"| {s['killed']} | {s['timeout']} | {s['survived']} | {s['no_tests']} "
        f"| {s['suspicious']} | {s['total']} |",
    ]
    survivors = surviving_by_function()
    if survivors:
        lines += ["", "<details><summary>Surviving mutants per function</summary>", "",
                  "| function | survived |", "|---|---:|"]
        lines += [f"| `{f}` | {n} |" for f, n in survivors.most_common()]
        lines += ["", "Inspect one with `mutmut show <name>` or browse with `mutmut browse`.",
                  "</details>"]
    report = "\n".join(lines) + "\n"
    print(report)
    if os.getenv("GITHUB_STEP_SUMMARY"):
        with open(os.environ["GITHUB_STEP_SUMMARY"], "a") as f:
            f.write(report)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
