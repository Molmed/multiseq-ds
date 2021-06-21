
#create TSS file from gene bodies file. TSS+-1000 bp

awk 'NR>1 {print $1 " " ($2 - 1000) " " ($2 + 1000)}' files/hg38.gene.bodies.txt > files/hg38.TSS.txt
sed  -i '1i chr start end' files/hg38.TSS.txt

