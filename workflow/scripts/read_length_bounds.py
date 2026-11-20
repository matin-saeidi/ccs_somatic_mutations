#!/usr/bin/env python3
"""
Read-length bounds of the read filter: the mean +/- 2 s.d. of the lengths of
the strand consensus reads of one library (Methods, Read filters;
Supplementary Fig. 2).

Input:  <sample>.np_len.txt, with a header and the columns
        zmw  strand  np  read_length
Output: one line "<min>\t<max>" to the output path (or stdout)

Usage:  python read_length_bounds.py <sample>.np_len.txt <sample>.read_length_bounds.txt
"""
import sys, statistics

src = sys.argv[1]
lengths = []
with open(src) as f:
    next(f)  # header
    for line in f:
        parts = line.split()
        if len(parts) >= 4:
            try:
                lengths.append(int(parts[3]))
            except ValueError:
                continue

if not lengths:
    sys.exit(f"No read lengths parsed from {src}")

mu = statistics.fmean(lengths)
sd = statistics.stdev(lengths)         # sample s.d.
lo = max(0, int(round(mu - 2 * sd)))
hi = int(round(mu + 2 * sd))

line = f"{lo}\t{hi}\n"
if len(sys.argv) > 2:
    with open(sys.argv[2], "w") as g:
        g.write(line)
else:
    sys.stdout.write(line)
