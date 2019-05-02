
cat  Rattus_norvegicus.Rnor_6.0.94.chr.gtf |sed "s/;/\t/" |cut -f 1,3,4,5,9|sed "s/gene_id \|\"//g" |grep "gene\|exon"|grep -v "^MT" |grep -v "#"|sed "s/^/chr/" >ensembl_gene_model.tab 

