wget -qO- http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/knownGene.txt.gz | gunzip -c > knownGene.hg38.txt
awk '{print $2 " " $4 " " $5 }' knownGene.hg38.txt > hg38.gene.bodies.txt
sed  -i '1i chr start end' hg38.gene.bodies.txt

