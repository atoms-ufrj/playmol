#!/bin/bash
if [ "$#" -ne 1 ]; then
  echo "usage: pdb2xyz <pdb-file-name>"
  exit
fi
if [ ! -f "$1" ]; then
    echo "error: file '$1' was not found"
    exit
fi
grep '^ATOM\|^HETATM' $1 | wc -l
echo ""
grep '^ATOM\|^HETATM' $1 | cut -c13-16,31-38,39-46,47-54

