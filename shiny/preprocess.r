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

filename="./eqtl_data.tab"

df0<-read.table(file=filename, header=F, sep="\t")
names(df0)<-c("region","gene","chr","bp","p")
head(df0)
length(unique(df0$gene))
head(df0)
dim(df0)
df1<-merge(df0,chrstat,by="chr")
df1$cumlength<-df1$bp+df1$lengthadj
df1<-df1[order(df1$gene,df1$cumlength),]
head(df1)

all_data<-merge(df1,genepos, by.x="gene", by.y="rn6_gene")
head(all_data)
dim(all_data)
all_data$cistrans<-"inbetween"
transidx<- abs(all_data$cumlength-all_data$rn6_g_start)>1e+6 | abs(all_data$cumlength-all_data$rn6_g_end) > 1e+6
cisidx<- abs(all_data$cumlength-all_data$rn6_g_start)<1e+6 | abs(all_data$cumlength-all_data$rn6_g_end) < 1e+6
all_data$cistrans[transidx]<-"trans"
all_data$cistrans[cisidx]<-"cis"
all_data$logp<- -log10(all_data$p)


temp<-read.table(file="./ensembl_id2symb.tab", header=F)
symb<-temp[,2]
names(symb)<-temp[,1]
head(symb)
symb["ENSRNOG00000052790"]


names(all_data)
dat<-all_data[,c("gene","chr","region","bp","rn6_start","rn6_g_start", "cistrans", "logp", "cumlength")]
save(file="eqtl.Rdata", dat, chrstat, symb)

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
