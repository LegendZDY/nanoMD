
conda activate mines

cat m6A.bed|sort -k1,1 -k2,2n > m6a.sorted.bed

perl ~/soft/metaPlotR/annotate_bed_file.pl --bed ./m6a.sorted.bed --bed2 ~/Refernce/human_ensembl_genome/91/hg38_annot.sorted.bed > annot_m6a.sorted.bed

perl ~/soft/metaPlotR/rel_and_abs_dist_calc.pl --bed ./annot_m6a.sorted.bed --regions ~/Refernce/human_ensembl_genome/91/region_sizes.txt > m6a.dist.measures.txt

cat Treat.mod.bed | awk -F "\t" '$4 > 0.95 {file="Treat/" $7 ".bed"; print $0 > file}'

cat Control.mod.bed | awk -F "\t" '$4 > 0.95 {file="Control/" $7 ".bed"; print $0 > file}'

cat Treat.mod.bed | awk -F "\t" '$4 >= 0.9 {file="treat_" $7 ".bed"; print $0 > file}'

cat Control.mod.bed | awk -F "\t" '$4 >= 0.9 {file="control_" $7 ".bed"; print $0 > file}'

awk '!seen[$4]++ {first[$4]=$0} {count[$4]++} END {for (key in count) print first[key], count[key]}' treat_m6A.bed > treat_m6A.uniq.bed

nohup python ./split.py -i Treat.mod.bed -p treat -v 0.6 &
nohup python ./split.py -i WT.mod.bed -p wt -v 0.6 &

R CMD INSTALL legendBaseModel_0.0.6.tar.gz