
# conda activate mines

# cat m6A.bed|sort -k1,1 -k2,2n > m6a.sorted.bed

# perl ~/soft/metaPlotR/annotate_bed_file.pl --bed ./m6a.sorted.bed --bed2 ~/Refernce/human_ensembl_genome/91/hg38_annot.sorted.bed > annot_m6a.sorted.bed

# perl ~/soft/metaPlotR/rel_and_abs_dist_calc.pl --bed ./annot_m6a.sorted.bed --regions ~/Refernce/human_ensembl_genome/91/region_sizes.txt > m6a.dist.measures.txt

