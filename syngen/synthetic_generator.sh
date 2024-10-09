#!/bin/bash

chmod +x fq2famaker.sh
chmod +x fq2txtmaker.sh

for file in *.fa; do
    filename=$(basename "$file")
    dirname="${filename%.*}"
    mkdir -p "fqfiles/$dirname"
    echo "art_illumina -ss HS25 -sam -i $file -l 150 -f 50 -o fqfiles/$dirname"
    art_illumina -ss HS25 -sam -i "$file" -l 150 -f 50 -o "fqfiles/$dirname"
    cd fqfiles
    mv synthetic*.fq $dirname
    mv synthetic*.sam $dirname
    mv synthetic*.aln $dirname
    cp ../fq2famaker.sh $dirname
    cp ../fq2txtmaker.sh $dirname
    cd $dirname
    ./fq2famaker.sh
    ./fq2txtmaker.sh
    rm *.sh
    cd ../../
done
