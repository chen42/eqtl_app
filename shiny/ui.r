library("ggvis")
ui<-fluidPage(
    titlePanel("Brain eQTL of HS rats"),
	sidebarLayout(
		sidebarPanel(
			width=3,
			fluidRow( column(5, textInput("chr", "Chr", "chr20", width=120)),
					  column(7, numericInput("loc", "Location, bp", 43803758, width=120 ))
			), 
			fluidRow(

					 column(5, numericInput("win", "Window, Mbp",2 , width=120)),
					 column(7, selectInput("region", "Brain Region", 
										   c("Accumbens core" = "AC",
											 "Infralimbic cortex" = "IL",
											 "Orbitofrontal cortex" = "OF",
											 "Prelimbic cortex" = "PL",
											 "Lat. Habenula"="LH"), 
										   tableOutput("output")
										   )
					 		)
			),
			fluidRow( 
					 selectInput("geneList", "Limit to Gene",
								  c("All")
								 )
			)
		),
		mainPanel(
#			textOutput("selected_var"),
			uiOutput("ggvis_ui"),
			ggvisOutput("ggvis"),
			helpText("You can drag the lower right corner of the figure to adjust its size.")
		)
	)
)
