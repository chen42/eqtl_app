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

genomesize=sum( as.numeric(chrstat$length[1:21]))
regions<-c("all five brain regions", "accumbene core", "lateral habenula", "prelimbic cortex", "infralimbic cortex", "orbitofrontal cortex")
names(regions)<-c("All", "AC", "LH", "PL", "IL", "OF")

#head(dat)
server <- function(input, output, session) {
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
	
#	type <- reactive({dat_f() %>% select(cistrans) %>%  distinct %>% arrange(cistrans)})
	dotsize <- reactive({
		dot_cnt<-dim(dat_f())[1]
		if (dot_cnt > 500){
			size=0.6
		} else if  (dot_cnt > 100) { 
			size=1
		} else {
			size=2	
		}
		size
	})  
	observe({
		updateSelectInput(session, "geneList", label="Ensembl ID", choices=selections())
#		updateSelectInput(session, "cistrans", label="type", choices=list("Both", type()))
	})
	## plot all SNPs in the selected region 
	output$first<-renderPlotly({
		p1<-ggplot()+
#		dat_f(), aes(x=qtl_bp, y=logp, color=gene, shape=cistrans)) +
			geom_point(data=dat_f(), aes(x=qtl_bp, y=logp, color=gene, shape=cistrans), alpha=0.6, size=dotsize()) +
			scale_x_continuous(name="M bp") +
			scale_y_continuous(name="-Log10(P)") +
			geom_point(data=dat_f(), aes(x=rn6_start, y=4.5, color=gene), shape=4, size=2) + 
			geom_hline(yintercept=5.6, color="red", size=0.2)+
			theme(legend.position='none', 
				  panel.grid.major.x=element_blank(),
				  panel.grid.minor.x=element_blank()
				  )
			if (dim(svs_f())[1]) {	
				p1<-p1+geom_segment(data=svs_f(), aes(x=sv_start, xend=sv_end, y=4, yend=4, color=sv_type, info1=sv_score), size=2)
			}
			if (dim(gaps_f())[1]) {	
				p1<-p1+geom_segment(data=gaps_f(), aes(x=gap_start, xend=gap_end, y=4.3, yend=4.3), size=2)
			}
		ggplotly(p1) %>% config(displayModeBar=F)
	})
	type=c("both cis- and trans-eQTL", "cis-eQTL", "trans-eQTL")
	names(type)<-c("Both","cis","trans")
	## Manhattan plot for one gene, all brain regions
	output$regionText<-renderText({paste("Focusing on ", input$chr, ", ", round(input$loc/bp_scale, 2) , " ± ", input$win, " M bp. Displaying ", type[input$cistrans], " in ", regions[input$region], ".", sep="" ) })
	output$legend<-renderText({paste("Only showing SNPs with -log10(P) > 4.9 (i.e., p < 1.25e-5). Colors: genes; Shape: cis- vs trans-;  Hovering mouse over the points for more info. X: location of the gene, ", sep="") })
	dat_m<-reactive({ 
		dat0<- dat %>%  filter(gene == input$geneList) %>% 
		droplevels()  %>% 
		arrange(qtl_bp)
	})
	output$second<-renderPlotly({
		p<-ggplot(dat_m(), aes(x=cumlength, y=logp, color=region)) + 
			geom_point(size=1, alpha=0.6) + 
			scale_x_continuous(name='Chr', label = gsub("chr", "", chrstat$chr), breaks= chrstat$chrmid, limits=c(0,genomesize)) +
			scale_y_continuous(expand = c(0, 0), name="-Log10(P)", limits=c(4.3,max(dat_m()$logp)+0.5)) + 
			geom_hline(yintercept=5.6, color="red", size=0.2)+
			geom_vline(xintercept=chrstat$lengthadj, size=0.1, color="white")+
			geom_point(aes(x=rn6_g_start, y=4.5), shape=4, size=2, color="black") + 
			theme(legend.position="right",
			 panel.grid.major.x = element_blank(),
			 panel.grid.major.y = element_blank(),
			 panel.grid.minor.y = element_blank(),
			 panel.grid.minor.x = element_blank()
			)
			ggplotly(p) %>% config(displayModeBar=F)
	})
	output$manhText<-renderText({paste("Abbr. Manhattan plot for ", symb[input$geneList], " (", input$geneList, "), only showing SNPs with -log10(P) > 4.9 (i.e., p < 1.25e-5).", sep="") })  
	## search by gene ID and symb

}


