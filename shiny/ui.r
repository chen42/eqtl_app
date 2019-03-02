library("ggvis")
ui<-fluidPage(
    titlePanel("Brain eQTL of HS rats"),
	sidebarLayout(
		sidebarPanel(
			width=3,
			fluidRow( column(5, textInput("chr", "Chr", "chr20", width=120)),
					  column(7, numericInput("loc", "Location, bp", 43803758, width=120, step=1000 ))
			), 
			fluidRow(

					 column(5, numericInput("win", "Window, Mbp",2 , width=120)),
					 column(7, selectInput("region", "Brain Region", 
										   c("All" = "All",
											 "Accumbens core" = "AC",
											 "Infralimbic cortex" = "IL",
											 "Orbitofrontal cortex" = "OF",
											 "Prelimbic cortex" = "PL",
											 "Lat. Habenula"="LH"), 
										   tableOutput("output")
										   )
					 		)
			),
					 hr(), 
			fluidRow( 
					 selectInput("geneList", "Abbreviated Manhattan Plot",
								  c("All"),
								  selected=c("All"),

								 )
			)
		),
		mainPanel(
#			uiOutput("ggvis_ui"),
			ggvisOutput("ggvis"),
#			uiOutput("geneplot"),
			ggvisOutput("geneplot"),
			#textOutput("selected_var"),
			helpText("You can drag the lower right corner of the figure to adjust its size.")
		)
	)
)
