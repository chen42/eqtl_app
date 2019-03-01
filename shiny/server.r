library("ggvis")
library("dplyr")

load(file="eqtl.Rdata")
bp_scale=1e+6
dat$bp<-dat$bp/bp_scale
dat$rn6_start<-dat$rn6_start/bp_scale

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

	selections <- dat_f %>% select(gene) %>%  distinct  
	observe({
		updateSelectInput(session, "geneList", label="Ensembl ID", choices=selections())
	})


	output$selected_var<-renderText({paste("Showing", str(input$chr)) })
	ggvis(dat_f, ~bp, ~logp, opacity :=0.5, fill =~gene ) %>% 
			#layer_points(shape=~cistrans) %>% 
			layer_points() %>% 
			layer_lines(y=5.6, stroke:="red" ) %>% 
			add_axis("x", title="M bp") %>% 
			add_axis("y", title="-Log10(p)") %>% 
			layer_points( x = ~rn6_start, y = 4.5, opacity :=.6, shape :="diamond") %>%
			bind_shiny("ggvis", "ggvis_ui")
}

			#add_axis("x", title="M bp", properties=axis_props(grid=list(stroke="#00000010"))) %>% 
			#add_axis("y", title="-Log10(p)", properties=axis_props(grid=list(stroke="#00000010"))) %>% 
