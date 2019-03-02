library("ggvis")
library("dplyr")
library("ggplot2")

load(file="eqtl.Rdata")
bp_scale=1e+6
dat$bp<-dat$bp/bp_scale
dat$rn6_start<-dat$rn6_start/bp_scale
genomesize=sum(chrstat$length)

#head(dat)
server <- function(input, output, session) {
	dat_f<-reactive({	
		hide<-abs(dat$rn6_start-input$loc/bp_scale)>input$win
		dat$rn6_start[hide]<-NA
		if (input$region != "All")  {
			dat0<- dat %>% 
				filter( chr == input$chr & 
					abs(bp*bp_scale - input$loc) < input$win*1e+6 & 
				 region==input$region) %>% 
				droplevels() %>%  
				arrange(bp)
		} else {
			dat0<- dat %>% 
				filter( chr == input$chr & 
					abs(bp*bp_scale - input$loc) < input$win*1e+6   
					) %>% 
				droplevels() %>%  
				arrange(bp)
		}	
	})
	## generate list of genes in the plot
	selections <- dat_f %>% select(gene) %>%  distinct  
	observe({
		updateSelectInput(session, "geneList", label="Ensembl ID", choices=selections())
	})

	output$region<-renderPlot({
		ggplot(dat_f(), aes(x=bp, y=logp, color=gene, shape=cistrans)) +
			geom_point(alpha=0.5, size=4)   +
			scale_x_continuous(name="M bp") +
			scale_y_continuous(name="-Log10(P)") +
			geom_point(aes(x=rn6_start, y=4.5), shape=4, size=4) + 
			geom_hline(yintercept=5.6, color="red", linetype="dashed")+
			theme_bw() + theme( legend.position="right",
			 panel.grid.major.x = element_blank(),
			 panel.grid.minor.x = element_blank()
			)
	})
	
	output$regionText<-renderText({paste("Showing all SNPs with p<1.25e-5 located on ", input$chr, ", ", round(input$loc/bp_scale, 2) , "±", input$win, " M bp. X: location of the gene", sep="" ) })
	## Manhattan plot
	dat_m<-reactive({ 
		dat0<- dat %>%  filter(gene == input$geneList) %>% 
		droplevels()  %>% 
		arrange(bp)
	})

	output$manh<-renderPlot({
		ggplot(dat_m(), aes(x=cumlength, y=logp, color=region)) + 
			geom_point(size=4, alpha=0.6) + 
			scale_x_continuous(name=NULL, label = gsub("chr", "", chrstat$chr), breaks= chrstat$chrmid, limits=c(0,genomesize)) +
			scale_y_continuous(expand = c(0, 0), name="-Log10(P)", limits=c(4.3,max(dat_m()$logp)+0.5)) + 
			geom_hline(yintercept=5.6, color="red", linetype="dashed")+
			geom_point(aes(x=rn6_g_start, y=4.5), shape=4, size=4, color="black") + 
			theme_bw() + theme( 
			 legend.position="right",
			 panel.grid.major.x = element_blank(),
			 panel.grid.minor.x = element_blank()
			)
	})
	output$manhText<-renderText({paste("Abbr. Manhattan plot for ", input$geneList, " (", symb[input$geneList], "), showing all SNPs with p<1.25e-5.", sep="") })  
	output$legend<-renderText({paste("AC: Accumbens Core, IL: Infralimbic Cortex, PL: Prelimbic Cortex, OF: Orbitofrontal Cortex, LH: Lat. Habenula", sep="") })

}


