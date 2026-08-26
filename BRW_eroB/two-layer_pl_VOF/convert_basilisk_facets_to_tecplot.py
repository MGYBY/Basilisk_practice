#!/usr/bin/env python3
"""Convert 2-D Basilisk output_facets files to Tecplot FELINESEG ASCII files."""

from pathlib import Path
import argparse

def read_segments(path):
    segments = []
    current = []
    with Path(path).open("r", encoding="utf-8", errors="replace") as f:
        for raw in f:
            line = raw.strip()
            if not line:
                if current:
                    segments.append(current)
                    current = []
                continue
            values = [float(v) for v in line.split()]
            if len(values) != 2:
                raise ValueError(f"Expected two columns, got {len(values)}: {line}")
            current.append(values)
    if current:
        segments.append(current)

    for i, seg in enumerate(segments, start=1):
        if len(seg) != 2:
            raise ValueError(
                f"Facet {i} contains {len(seg)} points; expected 2 for a 2-D output_facets file."
            )
    return segments

def convert(src, dst, zone_name):
    segments = read_segments(src)
    with Path(dst).open("w", encoding="utf-8") as f:
        f.write('TITLE = "Basilisk output_facets interface"\n')
        f.write('VARIABLES = "X" "Y"\n')
        f.write(
            f'ZONE T="{zone_name}", N={2*len(segments)}, E={len(segments)}, '
            'DATAPACKING=POINT, ZONETYPE=FELINESEG\n'
        )

        for segment in segments:
            for x, y in segment:
                f.write(f"{x:.16g} {y:.16g}\n")

        for i in range(len(segments)):
            f.write(f"{2*i+1} {2*i+2}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input")
    parser.add_argument("output")
    parser.add_argument("--zone-name", default="Basilisk_interface")
    args = parser.parse_args()
    convert(args.input, args.output, args.zone_name)
