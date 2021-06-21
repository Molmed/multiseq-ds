
rm barcode_stats.csv

for FILE in /crex/proj/uppstore2017165/FU180/210428_A00605_0213_AH72YVDRXY/Sample_FU180-MULTI-seq-Barcodes/*
do 
  COUNT_1=$(gzip -cd $FILE | grep -o -i ^GGAGAAGA | wc -l)
  COUNT_2=$(gzip -cd $FILE | grep -o -i ^CCAACCGG | wc -l)
  DEV=4 
  TOT=$(gzip -cd $FILE | wc -l)
  TOT_SEQ=$((TOT/DEV))
  printf "$FILE GGAGAAGA $COUNT_1 $TOT_SEQ\n" >> barcode_stats.csv
  printf "$FILE CCAACCGG $COUNT_2 $TOT_SEQ\n" >> barcode_stats.csv   
done

sed  -i '1i file barcode count total' barcode_stats.csv
