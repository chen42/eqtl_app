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
			fluidRow( 
					 selectInput("geneList", "Focus on Gene",
								  c("All"),
								  selected=c("All"),

								 )
			)
		),
		mainPanel(
			uiOutput("ggvis_ui"),
			ggvisOutput("ggvis"),
			#textOutput("selected_var"),
			helpText("You can drag the lower right corner of the figure to adjust its size.")
		)
	)
)
