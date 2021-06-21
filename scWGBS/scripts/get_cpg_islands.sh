#!/bin/bash

#download and format cpg islands 

wget -qO- http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cpgIslandExt.txt.gz | gunzip -c > hg38.cpg.txt
awk '{print $2 " " $3 " " $4}' hg38.cpg.txt > hg38.cpg.positions.txt
sed -i 's/chr//g' hg38.cpg.positions.txt
sed  -i '1i chr start end' hg38.cpg.positions.txt

