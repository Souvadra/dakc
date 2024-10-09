#!/bin/bash

for file in reference_files/*.fa; do
    filename=$(basename "$file")
    dirname="${filename%.*}"
    mkdir -p "fqfiles/$dirname"
    echo "art_illumina -ss HS25 -sam -i $file -l 150 -f 50 -o fqfiles/$dirname"
    art_illumina -ss HS25 -sam -i "$file" -l 150 -f 50 -o "fqfiles/$dirname"
done
