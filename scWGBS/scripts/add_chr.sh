#!/bin/bash

for f in /crex/proj/uppstore2017165/nobackup/alva/scWGBS/K562_out_GB_hg38/CpG_coordinates_in_regions/*; do
   cp "$f" "$f~" &&    
   gzip -cd "$f~" | sed '2,$s/^/chr/' | gzip > "$f"
done

