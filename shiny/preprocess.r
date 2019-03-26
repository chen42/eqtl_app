library('dplyr')
chrstat<-read.table(file="./rn6_chr_length", header=F)
names(chrstat)<-c("chr","length")
chrstat$lengthadj<-c(0,cumsum(as.numeric(chrstat$length))[-21])
chrstat$chrmid<-chrstat$length/2+chrstat$lengthadj
head(chrstat)

genepos<-read.table(file="./rn6_gene_pos.tab", sep="\t", header=F)
names(genepos)<-c("gene","chr","start","end")
genepos<-merge(genepos,chrstat[,c("chr","lengthadj")], by="chr")
genepos$g_start<-genepos$start+genepos$lengthadj
genepos$g_end<-genepos$end+genepos$lengthadj
head(genepos)
genepos<-genepos[,c("chr","gene","g_start","g_end", "start", "end")]
names(genepos)<-paste("rn6",names(genepos),sep="_")
tail(genepos)

# for trans, p-values cut at -logp>4.9
filename="./eqtl_data.tab"
df0<-read.table(file=filename, header=F, sep="\t")
names(df0)<-c("region","gene","qtl_chr","qtl_bp","qtl_p")
df1<-merge(df0,chrstat,by.x="qtl_chr", by.y='chr')
df1$cumlength<-df1$qtl_bp+df1$lengthadj
df1<-df1[order(df1$gene,df1$cumlength),]
trans<-merge(df1,genepos, by.x="gene", by.y="rn6_gene")
trans$cistrans<-"inbetween"
transidx<- abs(trans$cumlength-trans$rn6_g_start)>1e+6 | abs(trans$cumlength-trans$rn6_g_end) > 1e+6
trans$cistrans[transidx]<-"trans"
#cisidx<- abs(trans$cumlength-trans$rn6_g_start)<1e+6 | abs(trans$cumlength-trans$rn6_g_end) < 1e+6
#trans$cistrans[cisidx]<-"cis"
trans$logp<- -log10(trans$qtl_p)
transeqtl<-subset(trans, logp>4.9 & cistrans=="trans") # threshold is 5.6 but keep a little more data points for the plot
dim(transeqtl)
# 348292 16

## for cis 
filename="./five_regions_cis_candidates.tab"
df0<-read.table(file=filename, header=F, sep="\t")
names(df0)<-c("region","gene","qtl_chr","qtl_bp","qtl_p")
head(df0)
head(chrstat)
df1<-merge(df0,chrstat,by.x="qtl_chr", by.y='chr')
df1$cumlength<-df1$qtl_bp+df1$lengthadj
df1<-df1[order(df1$gene,df1$cumlength),]
cisCandidates<-merge(df1,genepos, by.x="gene", by.y="rn6_gene")
cisCandidates$cistrans<-"cis"
cisCandidates$logp<- -1 * log10(cisCandidates$qtl_p)

# for cis-, correct p-values
eigenmt<-read.table(file="./fpkm_eigenmt_hits.txt", header=F,sep="\t")
names(eigenmt)<-c("region","chr_loc","gene","p","padj","tests")
head(eigenmt)
# eigenMT calculates the Meff (effective multi-test) for each gene, then do Bonferroni. Here I am recapture this tests and apply this to all SNPs within the cis- of each gene. The tests are slightly different (with one or two) between different brains regions. So I am taking the mean
gene_testcounts<-aggregate(tests~gene, FUN=function(x) round(mean(x),0), data=eigenmt)
ciseqtl0<-merge(cisCandidates,gene_testcounts, by.x="gene", by.y="gene")
ciseqtl0$qtl_p<-ciseqtl0$qtl_p*ciseqtl0$tests
#ciseqtl0<-ciseqtl0[,1:16]
#names(ciseqtl0)
#names(transeqtl)
min_p<-aggregate(qtl_p~gene+region, FUN=min,data=ciseqtl0)
min_p_sig<-subset(min_p, qtl_p<0.05)
names(min_p_sig)[3]<-"min_p"
head(min_p)
ciseqtl<-merge(ciseqtl0, min_p_sig, by.x=c("gene","region"), by.y=c("gene","region"))[,1:16]
#dim(ciseqtl0)
#dim(ciseqtl_with_min_p)
#head(ciseqtl_with_min_p)
#tail(ciseqtl_with_min_p,10)

all_data<-rbind(ciseqtl,transeqtl)
names(all_data)
dat<-all_data[,c("gene","qtl_chr","region","qtl_bp","rn6_chr","rn6_start","rn6_g_start", "cistrans", "logp", "cumlength")]

#fdr from another script (fdr.r)
fdr5pct<-c(0.00481,0.00541,0.00594,0.00549,0.00443)
names(fdr5pct)<-c("AC","IL","PL","OF","LH")




# gene symbl
temp<-read.table(file="./ensembl_id2symb.tab", header=F)
symb<-temp[,2]
names(symb)<-temp[,1]
head(symb)
symb["ENSRNOG00000052790"]

# gene model
gmodel<-read.table(file="./ensembl_gene_model.tab", header=F)
names(gmodel)<-c("gm_chr","gm_type","gm_start","gm_end","gene")

# gaps
gaps<-read.table(file="./gap.txt", header=F)[,c(2,3,4)]
names(gaps)<-c("gap_chr","gap_start","gap_end")
head(gaps)

# svs
svs<-read.table(file="./svs.tab", header=F)
names(svs)<-c("sv_set","sv_chr","sv_start","sv_end", "sv_score", "sv_type")
head(svs)

save(file="eqtl.Rdata", dat, chrstat, symb, gaps, svs, gmodel, fdr5pct)

hide<-function(){# define peaks?
	distance<- c(0,df1$cumlength) - c(df1$cumlength,2800000000)
	dist0<- distance[distance>0 ]
	dist0
	#hist(dist0, breaks=100, col="lightgreen")
	span=1e+6 #18241 eqtls
	span=5e+6 #16457 eqtls ## merge SNPs within the span into one eqtl
	moveDownOne<-c(0,df1$cumlength[c(1:(dim(df1)[1]-1))])
	idx<-abs(df1$cumlength-moveDownOne) > span
	df1$peakcnt<-cumsum(idx)
	head(df1)
	df1<-df1[order(df1$gene,df1$peakcnt,df1$p_value),]
	eqtl<-df1%>% group_by(peakcnt)%>%slice(which.min(p_value)) 
	names(eqtl)<-paste("eqtl",names(eqtl),sep="_")
	dim(eqtl)
	head(eqtl)
	dim(eqtl)
	#write.table(file="matrix_eqtl_pl.csv",as.data.frame(eqtl),sep="\t")
}
