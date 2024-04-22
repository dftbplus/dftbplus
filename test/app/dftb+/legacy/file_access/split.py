#!/usr/bin/env python3

with open("dftb_in.hsd.all") as inp:
    for ii, chunk in enumerate(inp.read().split("%")):
        # We use 2 digits suffixes, so only the first 100 entries will be correctly ordered with
        # alphabetical sorting (as used by ls command)
        if ii >= 100:
            raise ValueError("Obtained more than 100 fragments")
        with open(f"dftb_in.hsd.{ii:02d}", "w") as out:
            out.write("<<+ 'dftb_in.hsd.common'\n")
            out.write(chunk)
