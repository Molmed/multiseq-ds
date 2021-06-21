#!/bin/bash

#download and fornat cpg islands 

wget -qO- http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cpgIslandExt.txt.gz | gunzip -c > hg38.cpg.txt
awk '{print $2 " " $3 " " $4}' hg38.cpg.txt > hg38.cpg.positions.txt
sed -i 's/chr//g' hg38.cpg.positions.txt
sed  -i '1i chr start end' hg38.cpg.positions.txt

wget -qO- http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/cpgIslandExt.txt.gz | gunzip -c > hg19.cpg.txt
awk '{print $2 " " $3 " " $4}' hg19.cpg.txt > hg19.cpg.positions.txt
sed  -i '1i chr start end' hg19.cpg.positions.txt
