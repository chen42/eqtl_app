library("dplyr")
library("ggplot2")
library("plotly")

load(file="eqtl.Rdata")
bp_scale=1e+6
dat$bp<-dat$bp/bp_scale
dat$rn6_start<-dat$rn6_start/bp_scale
genomesize=sum(chrstat$length[1:21])
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
				filter( chr == input$chr & 
					abs(bp*bp_scale - input$loc) < input$win*1e+6   
					) %>% 
				droplevels() %>%  
				arrange(bp)
		# filter by brain region	
		if (input$region != "All")  {
			dat0<- dat0 %>% 
				filter(region==input$region) %>% 
				droplevels() %>%  
				arrange(bp)
		}
		if (input$cistrans != "Both")  {
			dat0<- dat0 %>% 
				filter(cistrans==input$cistrans) %>% 
				droplevels() %>%  
				arrange(bp)
		}
		dat0	
	})
	## generate list of genes in the plot
	selections <- dat_f %>% select(gene) %>%  distinct %>% arrange(gene)  
	type <- dat_f %>% select(cistrans) %>%  distinct %>% arrange(cistrans)  
	observe({
		updateSelectInput(session, "geneList", label="Ensembl ID", choices=selections())
#		updateSelectInput(session, "cistrans", label="type", choices=list("Both", type()))
		 
	})
	## plot all SNPs in the selected region 
	output$first<-renderPlotly({
		p1<-ggplot(dat_f(), aes(x=bp, y=logp, color=gene, shape=cistrans)) +
			geom_point(alpha=0.6, size=2)   +
			scale_x_continuous(name="M bp") +
			scale_y_continuous(name="-Log10(P)") +
			geom_point(aes(x=rn6_start, y=4.5), shape=4, size=2) + 
			geom_hline(yintercept=5.6, color="red", size=0.2)+
			theme(legend.position='none', 
				  panel.grid.major.x=element_blank(),
				  panel.grid.minor.x=element_blank()
				  )
		ggplotly(p1) %>% config(displayModeBar=F)
	})
	type=c("both cis- and trans-eQTL", "cis-eQTL", "trans-eQTL")
	names(type)<-c("Both","cis","trans")
	## Manhattan plot for one gene, all brain regions
	output$regionText<-renderText({paste("Focusing on ", input$chr, ", ", round(input$loc/bp_scale, 2) , " ± ", input$win, " M bp. Displaying ", type[input$cistrans], " in ", regions[input$region], ".", sep="" ) })
	output$legend<-renderText({paste("Only showing SNPs with -log10(P) > 4.9 (i.e., p < 1.25e-5). Colors: genes; Shape: cis- vs trans-;  Hovering mouse over the points for more info. X: location of the gene", sep="") })
	dat_m<-reactive({ 
		dat0<- dat %>%  filter(gene == input$geneList) %>% 
		droplevels()  %>% 
		arrange(bp)
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


