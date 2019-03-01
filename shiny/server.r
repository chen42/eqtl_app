library("ggvis")
library("dplyr")


bp_scale=1e+6
dat<-read.table(file="eqtl_data.tab", sep="\t",header=F)
names(dat)<-c("region", "Gene_ID", "chr", "bp", "p")
dat$logp<- -log10(dat$p)
dat$bp<-dat$bp/bp_scale

server <- function(input, output, session) {
	dat_f<-reactive({	
			dat %>% 
				filter( chr == input$chr & 
					abs(bp*bp_scale - input$loc) < input$win*1e+6 &  
					region == input$region
					) %>% 
				droplevels()   
	})

	selections <- dat_f %>% select(Gene_ID) %>%  distinct  

	observe({
		updateSelectInput(session, "geneList", label="Ensembl ID", choices=selections())
	})
	mark <- reactive({data.frame(X=input$loc/bp_scale, Y=4.8)})
	#output$selected_var<-renderText({paste(str(selections)) })
	ggvis(dat_f, ~bp, ~logp, opacity :=0.3, fill =~Gene_ID  ) %>% 
			layer_points() %>% 
			layer_lines(y=5.6, stroke:="red" ) %>% 
			add_axis("x", title="M bp") %>% 
			add_axis("y", title="-Log10(p)") %>% 
			layer_points( data = mark, x = ~X, y = ~Y, fill :="black", opacity :=1, shape :="diamond") %>%
			bind_shiny("ggvis", "ggvis_ui")
}

			#add_axis("x", title="M bp", properties=axis_props(grid=list(stroke="#00000010"))) %>% 
			#add_axis("y", title="-Log10(p)", properties=axis_props(grid=list(stroke="#00000010"))) %>% 
