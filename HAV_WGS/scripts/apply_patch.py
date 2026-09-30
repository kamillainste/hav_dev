"""Apply a unified diff's -/+ line pairs to one file as exact replacements.

Usage: python apply_patch.py <diff> <target_file>

Used by setup.sh instead of `patch`, which is not always installed. Every
removed line must occur exactly once in the target, otherwise nothing is
written and the script exits 1.
"""
import sys

diff, target = sys.argv[1], sys.argv[2]
lines = open(diff).read().splitlines()
old = [l[1:] for l in lines if l.startswith("-") and not l.startswith("---")]
new = [l[1:] for l in lines if l.startswith("+") and not l.startswith("+++")]
if not old or len(old) != len(new):
    sys.exit(f"malformed diff: {diff}")

text = open(target).read()
for o, n in zip(old, new):
    if text.count(o + "\n") != 1:
        sys.exit(f"line not found exactly once in {target}: {o.strip()}")
    text = text.replace(o + "\n", n + "\n")
open(target, "w").write(text)
print(f"applied {len(old)} line change(s) to {target}")
