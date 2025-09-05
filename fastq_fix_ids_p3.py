#!/usr/bin/env python3
"""
Normalize FASTQ read headers:
- If a header line contains spaces (CASAVA 1.8+ style like '@... 1:N:0:INDEX'),
  convert to '@<first_token>/<side>' (dropping the trailing fields by default).
- If a header already ends with '/1' or '/2', leave it as-is.
- Side is auto-detected from header (when available) or from filename (R1/R2 or _1/_2).

Usage:
  fastq_fix_ids.py list.txt [--gzip | --no-gzip] [--keep-rest]
                            [--outdir OUTDIR] [--default-side 1]
                            [--dry-run]

The list file should contain one FASTQ(.gz) path per line.
"""

from __future__ import annotations
import argparse
import gzip
from pathlib import Path
import re
import sys
from typing import IO, Tuple, Optional

R1_PAT = re.compile(r'(?i)(?:^|[_.-])R?1(?:[_.-]|$)')
R2_PAT = re.compile(r'(?i)(?:^|[_.-])R?2(?:[_.-]|$)')

def open_in(path: Path) -> IO[str]:
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open("rt")

def open_out(path: Path, use_gzip: bool) -> IO[str]:
    if use_gzip:
        # .gz extension is ensured later
        return gzip.open(path, "wt", compresslevel=6)
    return path.open("wt")

def guess_side_from_header(h: str) -> Optional[str]:
    """
    Header may look like:
      '@ID 1:N:0:INDEX'  -> returns '1'
      '@ID 2:N:0:INDEX'  -> returns '2'
      '@ID/1'            -> returns '1'
      '@ID/2'            -> returns '2'
    """
    if h.endswith("/1") or h.endswith("/2"):
        return h[-1]
    parts = h.split()
    if len(parts) >= 2 and parts[1] and parts[1][0] in ("1", "2"):
        return parts[1][0]
    return None

def guess_side_from_filename(p: Path) -> Optional[str]:
    name = p.name
    if R1_PAT.search(name) and not R2_PAT.search(name):
        return "1"
    if R2_PAT.search(name) and not R1_PAT.search(name):
        return "2"
    # ambiguous or not found
    return None

def normalize_header(raw_header: str, side: str, keep_rest: bool) -> Tuple[str, bool]:
    """
    Returns (new_header_line, changed?)
    - Only modifies the header (line 1 of each 4-line FASTQ record).
    - Preserves '@' at the start.
    """
    assert raw_header.startswith("@"), "FASTQ header must start with '@'"
    header = raw_header.rstrip("\n")

    # Already suffixed with /1 or /2? If so, leave unchanged.
    if header.endswith("/1") or header.endswith("/2"):
        return header + "\n", False

    # If it has a space, take the first token; drop or keep the rest.
    parts = header.split(maxsplit=1)
    first = parts[0]
    rest = parts[1] if len(parts) == 2 else ""

    if keep_rest and rest:
        # Append /side but keep the remainder after a space.
        # Note: this duplicates the read number information; that’s intentional when keep_rest=True.
        new_h = f"{first}/{side} {rest}"
    else:
        new_h = f"{first}/{side}"

    if new_h == header:
        return header + "\n", False
    return new_h + "\n", True

def out_path_for(inpath: Path, outdir: Optional[Path], use_gzip: bool) -> Path:
    """
    Produce '<basename>.checked.fastq(.gz)' next to input (or in outdir).
    """
    base = inpath.name
    # strip one or two suffixes if .gz
    stem = base[:-3] if base.endswith(".gz") else base
    # ensure we end with .fastq or .fq in name before appending suffix
    if stem.endswith(".fastq") or stem.endswith(".fq"):
        stem_nosuf = stem.rsplit(".", 1)[0]
    else:
        stem_nosuf = stem
    outname = f"{stem_nosuf}.checked.fastq"
    if use_gzip:
        outname += ".gz"
    return (outdir or inpath.parent) / outname

def process_one(inpath: Path, use_gzip: bool, keep_rest: bool,
                default_side: str, dry_run: bool=False) -> Tuple[Path, bool, int]:
    """
    Returns (outpath, changed_any, records_processed)
    """
    # Peek the first header to infer side if possible
    with open_in(inpath) as fh:
        first_header = fh.readline().rstrip("\n")
    if not first_header.startswith("@"):
        raise ValueError(f"{inpath}: not a FASTQ file (first line does not start with '@').")

    side = guess_side_from_header(first_header) or guess_side_from_filename(inpath) or default_side

    outpath = out_path_for(inpath, None, use_gzip)
    changed_any = False
    nrec = 0

    if dry_run:
        return outpath, True, 0

    with open_in(inpath) as fin, open_out(outpath, use_gzip) as fout:
        while True:
            h = fin.readline()
            if not h:
                break  # EOF
            s = fin.readline()
            plus = fin.readline()
            q = fin.readline()
            if not q:
                raise ValueError(f"{inpath}: truncated FASTQ record at record {nrec+1}.")

            new_h, changed = normalize_header(h, side, keep_rest)
            if changed:
                changed_any = True
            fout.write(new_h)
            fout.write(s)
            fout.write(plus)
            fout.write(q)
            nrec += 1

    return outpath, changed_any, nrec

def main():
    ap = argparse.ArgumentParser(description="Normalize FASTQ headers to include '/1' or '/2'.")
    ap.add_argument("listfile", help="Text file with one FASTQ(.gz) path per line")
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--gzip", dest="gzip_out", action="store_true", help="Write .gz output")
    g.add_argument("--no-gzip", dest="gzip_out", action="store_false", help="Write uncompressed output")
    ap.set_defaults(gzip_out=True)
    ap.add_argument("--keep-rest", action="store_true",
                    help="Keep header text after the first space (default: drop it, like the original script)")
    ap.add_argument("--outdir", type=Path, default=None, help="Write outputs here (default: alongside inputs)")
    ap.add_argument("--default-side", choices=["1","2"], default="1",
                    help="Fallback side if it cannot be inferred from header/filename (default: 1)")
    ap.add_argument("--dry-run", action="store_true", help="Just show what would be written, do not create files")
    args = ap.parse_args()

    listfile = Path(args.listfile)
    if not listfile.exists():
        sys.exit(f"List file not found: {listfile}")

    with listfile.open() as lf:
        paths = [Path(line.strip()) for line in lf if line.strip()]

    if not paths:
        sys.exit("No inputs found in the list.")

    for p in paths:
        if not p.exists():
            print(f"[skip] Missing: {p}", file=sys.stderr)
            continue
        # Decide compression per user flag
        outp = out_path_for(p, args.outdir, args.gzip_out)
        if args.dry_run:
            print(f"[dry-run] {p} -> {outp}")
            continue
        outp.parent.mkdir(parents=True, exist_ok=True)
        outp, changed, nrec = process_one(p, args.gzip_out, args.keep_rest, args.default_side)
        print(f"[ok] {p.name} -> {outp.name} ({nrec} reads){' [modified]' if changed else ' [unchanged]'}")

if __name__ == "__main__":
    main()

