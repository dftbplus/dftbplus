#!/usr/bin/env python3

with open("dftb_in.hsd.all") as inp:
    for ii, chunk in enumerate(inp.read().split("%")):
        with open(f"dftb_in.hsd.{ii:02d}", "w") as out:
            out.write("<<+ 'dftb_in.hsd.common'\n")
            out.write(chunk)
