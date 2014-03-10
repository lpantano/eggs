
set -e

SB_HOME=~/scripts/git_seqbuster
SB_DB="$SB_HOME/modules/miraligner/DB"
echo "download miRBase files"
wget -q ftp://mirbase.org/pub/mirbase/CURRENT/hairpin.fa.gz
wget -q ftp://mirbase.org/pub/mirbase/CURRENT/miRNA.str.gz

gunzip hairpin.fa.gz miRNA.str.gz

awk '{if ($0~/>/){print name"\t"seq;name=$1;seq=""}else{seq=seq$0}}' hairpin.fa | grep "^>hsa" | sed "s/U/T/g" | sed 's/\t/\n/' > hairpin.hsa.fa

python "$SB_HOME"/misc/miRNA.simulator.py -f hairpin.hsa.fa -m miRNA.str -n 10 -s hsa > sim.20.hsa.fa

mkdir -p miraligner
cd miraligner
 java -jar  "$SB_HOME"/modules/miraligner/miraligner.jar  -sub 1 -trim 3 -add  3 -s hsa -i ../sim.20.hsa.fa -db $SB_DB -o sim.20.hsa -pre
cd ..

mkdir -p bowtie2
cd bowtie2
 bowtie2-build ../hairpin.hsa.fa mirbowtie2
 bowtie2 -D 20 -R 3 -N 0 -L 10 -i S,1,0.50 -f -x mirbowtie2 -U ../sim.20.hsa.fa -S sim.20.hsa.sam
cd ..

mkdir -p novo
cd novo
 novoindex -k 14 -s 2 -m mirnovo ../hairpin.hsa.fa  
 novoalign -d mirnovo -F FA -f ../sim.20.hsa.fa -a -l 15 -t 30 -o SAM >sim.20.novo.sam
cd ..

mkdir -p gem
cd gem
 gem-indexer -i ../hairpin.hsa.fa  -o mirgem --complement 'yes' 
 gem-mapper -e 0.1 -m 0.1 -I mirgem.gem -i ../sim.20.hsa.fa -o sim.20.hsa 
 gem-2-sam   -I mirgem.gem -i sim.20.hsa.map -o sim.20.hsa.sam
cd ..

python check.align.py > stats

echo "now time to run stats.rmd"