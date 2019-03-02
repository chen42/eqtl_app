library("ggvis")
library("dplyr")
library("ggplot2")

load(file="eqtl.Rdata")
bp_scale=1e+6
dat$bp<-dat$bp/bp_scale
dat$rn6_start<-dat$rn6_start/bp_scale

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

	ggvis(dat_f, ~bp, ~logp, opacity :=0.5, fill =~gene ) %>% 
		#layer_points(shape=~cistrans) %>% 
		layer_points() %>% 
		layer_lines(y=5.6, stroke:="red") %>% 
		add_axis("x", title="M bp") %>% 
		add_axis("y", title="-Log10(p)") %>% 
		layer_points( x = ~rn6_start, y = 4.5, opacity :=.6, shape :="diamond") %>%
		bind_shiny("ggvis", "ggvis_ui")

	## Manhattan plot
	dat_m<-reactive({ 
		dat0<- dat %>%  filter(gene == input$geneList) %>% 
		droplevels()  %>% 
		arrange(bp)
	})
	# trouble shoting
	output$selected_var<-renderText({paste("Showing", names(dat)) })
#	output$selected_var<-renderText({paste(names(dat_m())) })
	output$manh<-renderPlot({
		#dat<-subset(dat, gene=='ENSRNOG00000018217')
		axisdf<- dat %>% group_by(chr) %>% summarize(center=(max(cumlength)+ min(cumlength))/2) 
		ggplot(dat, aes(x=cumlength, y=logp)) + 
			geom_point(aes(color=as.factor(chr), alpha=0.5, size=1.2, shape=as.factor(region))) + 
			scale_color_manual(values=rep(c("grey", "skyblue"), 22))+ 
			#geom_vline(xintercept=man_vline, color="grey80", linetype="dashed")+
			scale_x_continuous(name=NULL, label = axisdf$chr, breaks= axisdf$center, limits=c(1,2778702144)) +
			scale_y_continuous(expand = c(0, 0), name="-Log10(p)") + 
			theme_bw() + theme( 
			 legend.position="none",
			 panel.border = element_blank(),
			 panel.grid.major.x = element_blank(),
			 panel.grid.minor.x = element_blank()
			)
	})



}


