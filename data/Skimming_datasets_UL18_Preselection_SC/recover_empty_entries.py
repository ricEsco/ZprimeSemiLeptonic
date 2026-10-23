#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Re-check every '<!--EMPTY ... -->' entry in dataset XML files written by readaMCatNloEntries.py.

readaMCatNloEntries.py comments a file out as EMPTY when it fails to open it within a 10 s timeout,
which can happen for perfectly good files (cold ROOT start-up when 10 worker processes open files at the same time). 
This script opens each commented-out file directly (no timeout) and:
  - un-comments it if its AnalysisTree has > 0 entries
  - leaves it commented if it really is empty, has no AnalysisTree, or cannot be opened
  - rewrites the '<!-- < NumberEntries=... -->' summary line with the corrected total

usage: python recover_empty_entries.py MC_*.xml DATA_*.xml          (add --dry-run to only report)
"""
from __future__ import print_function
import sys, re
import ROOT
ROOT.gROOT.SetBatch(True)

EMPTY_RE = re.compile(r'^<!--EMPTY\s+(<In FileName="([^"]+)"[^>]*/>)\s*-->\s*$')
IN_RE    = re.compile(r'^<In FileName="([^"]+)"')

def count_entries(fname):
    """Number of AnalysisTree entries, or None if the file/tree cannot be read."""
    try:
        f = ROOT.TFile.Open(fname)     # ROOT 6.14 returns a zombie on failure, newer PyROOT raises
    except Exception:
        return None
    if not f or f.IsZombie():
        return None
    tree = f.Get("AnalysisTree")
    n = tree.GetEntriesFast() if tree else None
    f.Close()
    return n

def fix_xml(xml, dry_run=False):
    lines = open(xml).read().splitlines()
    out, total, recovered, still_empty, unreadable = [], 0, [], [], []
    for line in lines:
        m_empty = EMPTY_RE.match(line.strip())
        m_in    = IN_RE.match(line.strip())
        if m_empty:
            entry, fname = m_empty.group(1), m_empty.group(2)
            n = count_entries(fname)
            if n is None:
                unreadable.append(fname); out.append(line)
            elif n == 0:
                still_empty.append(fname); out.append(line)
            else:
                recovered.append((fname, n)); out.append(entry); total += n
        elif m_in:
            n = count_entries(m_in.group(1))
            if n is None:
                print("  WARNING: active entry cannot be read:", m_in.group(1))
            else:
                total += n
            out.append(line)
        elif line.startswith('<!-- < NumberEntries'):
            continue                                   # rewritten below
        else:
            out.append(line)
    out.append('<!-- < NumberEntries="%d" Method=fast /> -->' % total)

    print("%s: recovered %d, truly empty %d, unreadable %d, total entries now %d"
          % (xml, len(recovered), len(still_empty), len(unreadable), total))
    for fname, n in recovered:  print("   recovered  %10d  %s" % (n, fname))
    for fname in still_empty:   print("   empty (kept commented)  %s" % fname)
    for fname in unreadable:    print("   UNREADABLE (kept commented)  %s" % fname)
    if not dry_run:
        open(xml, "w").write("\n".join(out) + "\n")

if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if a != "--dry-run"]
    for xml in args:
        fix_xml(xml, dry_run=("--dry-run" in sys.argv))