######Complete######
 learn=("BU-2" "BU-3" "BU-4" "BU-5" "BU-6" "BU-7 " "BU-8" "BU-9" "BU-10" "BU-11")
for i in ${learn[@]}
 do
 x=$i
 prokka --outdir $x --prefix $x --locustag $x \
 --addgenes --addmrna --compliant --centre CDC \
 --genus Bacteroides --species "Bacteroides uniformis" --strain $x \
 --kingdom Bacteria --usegenus --rfam --rnammer --force ./$x.fna
 done
######Scoffold######
learn=("strains_names")
for i in ${learn[@]}
 do
 x=$i
 prokka --outdir $x --prefix $x --locustag $x \
 --addgenes --addmrna --compliant --centre CDC \
 --genus Bacteroides --species "Bacteroides uniformis" --strain $x \
 --kingdom Bacteria --usegenus --rfam --cpus 4 --rnammer --force ./$x.fna
 done
######Roary######
roary -h
conda install -y R
R
roary -e --mafft -p 8 -r  -i 95 -cd 100 -f out_result roary_gff/*.gff 
cd out_result
query_pan_genome -a intersection   -o core_genes_list.xls ../roary_gff/*.gff  
query_pan_genome -a union  -o union_list.xls ../roary_gff/*.gff 
perl -p -i -e 's/: /\t/g' union_list.xls core_genes_list.xls
