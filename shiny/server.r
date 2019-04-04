library("dplyr")
library("ggplot2")
library("plotly")

load(file="eqtl.Rdata")
bp_scale=1e+6
dat$qtl_bp<-dat$qtl_bp/bp_scale
dat$rn6_start<-dat$rn6_start/bp_scale
gaps$gap_start<-gaps$gap_start/bp_scale
gaps$gap_end<-gaps$gap_end/bp_scale
svs$sv_start<-svs$sv_start/bp_scale
svs$sv_end<-svs$sv_end/bp_scale
gmodel$gm_start<-gmodel$gm_start/1e+6
gmodel$gm_end<-gmodel$gm_end/1e+6

genomesize=sum(as.numeric(chrstat$length[1:21]))
regions<-c("all five brain regions", "accumbene core", "lateral habenula", "prelimbic cortex", "infralimbic cortex", "orbitofrontal cortex")
names(regions)<-c("All", "AC", "LH", "PL", "IL", "OF")

server <- function(input, output, session) {
	font_size=3 # labels in the figures
	dat_f<-reactive({	
		## hide some extra points
		hide<-abs(dat$rn6_start-input$loc/bp_scale)>input$win
		dat$rn6_start[hide]<-NA
		## filter by chr and location
		dat0<- dat %>% 
				filter( rn6_chr == input$chr & qtl_chr==input$chr & 
					abs(qtl_bp*bp_scale - input$loc) < input$win*1e+6   
					) %>% 
				droplevels() %>%  
				arrange(qtl_bp)
		# filter by brain region	
		if (input$region != "All")  {
			dat0<- dat0 %>% 
				filter(region==input$region) %>% 
				droplevels() %>%  
				arrange(qtl_bp)
		}
		if (input$cistrans != "Both")  {
			dat0<- dat0 %>% 
				filter(cistrans==input$cistrans) %>% 
				droplevels() %>%  
				arrange(qtl_bp)
		}
		data.frame(dat0)
	})
	## generate list of genes in the plot
	selections <- reactive({dat_f() %>% select(gene) %>%  distinct %>% arrange(gene)})  
	gaps_f<-reactive({
		bprange=list()
		bprange<-range(dat_f()$qtl_bp)
		gaps<- gaps %>% filter(gap_chr==input$chr  &  gap_end < bprange[2] & gap_start > bprange[1]) 
		data.frame(gaps)
	}) 
	svs_f<-reactive({
		bprange=list()
		bprange<-range(dat_f()$qtl_bp)
		svs<- svs %>% filter(sv_chr==input$chr  &  sv_end < bprange[2] & sv_start > bprange[1]) 
		data.frame(svs)
	}) 
	gene_f<-reactive({
		bprange=list()
		bprange<-range(dat_f()$qtl_bp)
		gmodel<- gmodel %>% filter(gm_chr==input$chr & gm_type=="gene"  &  gm_end < bprange[2] & gm_start > bprange[1] & gene %in% dat_f()$gene) 
		data.frame(gmodel)
	}) 
	exon_f<-reactive({
		bprange=list()
		bprange<-range(dat_f()$qtl_bp)
		gmodel<- gmodel %>% filter(gm_chr==input$chr & gm_type=="exon"  &  gm_end < bprange[2] & gm_start > bprange[1] & gene %in% dat_f()$gene) 
		data.frame(gmodel)
	}) 
	
#	type <- reactive({dat_f() %>% select(cistrans) %>%  distinct %>% arrange(cistrans)})
	dotsize <- reactive({
		dot_cnt<-dim(dat_f())[1]
		if (dot_cnt > 500){
			size=2.3
		} else if  (dot_cnt > 100) { 
			size=3
		} else {
			size=4	
		}
		size
	})  
	observe({
		updateSelectInput(session, "geneList", label="Ensembl ID", choices=selections())
	})
	## plot all SNPs in the selected region 
	output$first<-renderPlotly({
		if(input$cistrans=='trans') {
			shp=c(24)
		} else {
			shp=c(21,24)
		}
		p1<-ggplot( data=dat_f(), aes(x=qtl_bp, y=logp))+
			geom_point(aes(fill=gene,shape=cistrans), alpha=1, size=dotsize()) +
#			scale_color_discrete(show_guide=FALSE)+
			scale_shape_manual(values=shp)+
			scale_x_continuous(name="Mb") +
			scale_y_continuous(name="-log10 p-value") +
			#geom_point(data=dat_f(), aes(x=rn6_start, y=4.5, color=gene), shape=4, size=2) + 
			geom_hline(yintercept=c(1.3,2.3,5.6), color="grey", size=0.2, linetype="dotted")+
			annotate("text",x=min(dat_f()$qtl_bp-0.15), y=c(1.3,2.3,5.6)+0.4, label=c("cis-\nadj.p=0.05","cis-\nFDR=0.05","trans-\np=2.5e-06"), size=font_size)+
			theme(legend.position='none', 
				  panel.grid.major.x=element_blank(),
				  panel.grid.minor.x=element_blank(),
				  panel.grid.major.y=element_blank(),
				  panel.grid.minor.y=element_blank()
				  )
			y_sv=0.1
			y_gap=0.4
			y_gene=0.7
			if (dim(svs_f())[1]) {	
				p1<-p1+
				geom_segment(data=svs_f(), aes(x=sv_start, xend=sv_end, y=y_sv, yend=y_sv, color=sv_type, info1=sv_score), size=1.5)+
				annotate("text", x=min(dat_f()$qtl_bp)-0.2, y=y_sv,label= "SV", size=font_size, hjust=0)
			}
			if (dim(gaps_f())[1]) {	
				p1<-p1+
				geom_segment(data=gaps_f(), aes(x=gap_start, xend=gap_end, y=y_gap, yend=y_gap), size=2) +
				annotate("text", x=min(dat_f()$qtl_bp)-0.2, y=y_gap,label= "gap", size=font_size, hjust=0)
			}
			if (dim(gene_f())[1]) {	
				p1<-p1+
				geom_segment(data=gene_f(), aes(x=gm_start, xend=gm_end, y=y_gene, yend=y_gene, color=gene), size=1) +
				annotate("text", x=min(dat_f()$qtl_bp)-0.2, y=y_gene,label= "gene", size=font_size, hjust=0)

			}
			if (dim(exon_f())[1]) {	
				p1<-p1+geom_segment(data=exon_f(), aes(x=gm_start, xend=gm_end, y=y_gene, yend=y_gene, color=gene), size=2)
			}
		ggplotly(p1) %>% config(displayModeBar=T)
	})
	type=c("both cis- and trans-eQTL", "cis-eQTL", "trans-eQTL")
	names(type)<-c("Both","cis","trans")
	## Manhattan plot for one gene, all brain regions
	output$regionText<-renderText({paste("Focusing on ", input$chr, ", ", round(input$loc/bp_scale, 2) , " ± ", input$win, " Mb. Displaying ", type[input$cistrans], " in ", regions[input$region], ".", sep="" ) })
	output$legend<-renderText({paste("Only showing SNPs with -log10(P) > 4.9 (i.e., p < 1.25e-5). Colors: genes; Shape: cis- vs trans-", sep="") })
	dat_m<-eventReactive(input$submitButton, { 
		if (nchar(input$geneSymb) >1) {
			idx<-toupper(input$geneSymb) == toupper(symb)
			if (sum(idx)==1){
				geneID=names(symb[idx])
				dat0<-dat %>% 
				filter(gene==geneID) %>%
				droplevels() %>%  
				arrange(qtl_bp)
			}
		} else {
			dat0<- dat %>%  filter(gene == input$geneList) %>% 
			droplevels()  %>% 
			arrange(qtl_bp)
		}
		updateTextInput(session, "geneSymb",  value="")
		dat0
	})

	geneInfo<-eventReactive(input$submitButton, {
		if (nchar(input$geneSymb) >1) {
			idx<-toupper(input$geneSymb) == toupper(symb)
			if (sum(idx)==1){
				geneSymbol= input$geneSymb
				geneID=names(symb[idx])
			}
		} else {
			geneSymbol=as.character(symb[input$geneList])
			geneID=input$geneList
		}
		c(geneSymbol, geneID) 		
	})
	output$second<-renderPlotly({
		output$manhText<-renderText({paste("Abbr. Manhattan plot for ", geneInfo()[1], " (", geneInfo()[2], "), only showing SNPs with -log10(P) > 4.9 for trans-eQTL and adj.p<0.05 for cis-eQTL.", sep="") })  
		p<-ggplot(data=dat_m(), aes(x=cumlength, y=logp)) + 
			geom_point(aes(color=region, shape=cistrans), size=1, alpha=0.7) + 
			scale_x_continuous(name='Chr', label = gsub("chr", "", chrstat$chr), breaks= chrstat$chrmid, limits=c(0,genomesize)) +
			scale_y_continuous(expand = c(0, 0), name="-log10 p-value", limits=c(0.3,max(dat_m()$logp)+0.5)) + 
			geom_hline(yintercept=c(1.3,2.3,5.6), color="grey", size=0.2, linetype="dotted")+
			annotate("text",x=5e+7, y=c(1.3,2.3,5.6)+0.2, label=c("cis-\nadj.p=0.05","cis-\nFDR=0.05","trans-\np=2.5e-06"), size=font_size)+
			geom_vline(xintercept=chrstat$lengthadj, size=0.1, color="white")+
			geom_point(aes(x=rn6_g_start, y=0.5), shape=4, size=2, color="black") + 
			guides(color=FALSE)+
			theme(legend.position="right",
			 panel.grid.major.x = element_blank(),
			 panel.grid.major.y = element_blank(),
			 panel.grid.minor.y = element_blank(),
			 panel.grid.minor.x = element_blank()
			)
			ggplotly(p) %>% config(displayModeBar=F)
	})

}


