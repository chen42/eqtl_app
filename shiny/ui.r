library("plotly")

ui<-fluidPage(
    titlePanel("Brain eQTL of HS rats"),
	sidebarLayout(
		sidebarPanel(
			width=3,
			fluidRow( column(4, textInput("chr", "Chr", "chr1", width=130)),
					  column(5, numericInput("loc", "bp", 175030758, width=150, step=1000 )),
					  column(3, numericInput("win", "±, M bp",2 , width=120))
			), 
			fluidRow(
					 column(4, selectInput("cistrans","Type", c("Both"="Both", "cis-"="cis", "trans-"="trans"))),
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
					 selectInput("geneList", "Abbreviated Manhattan Plot",
								  c("All"),
								  selected=c("All"),
								width=250
								 )
			),
			hr(), 
			fluidRow(textInput("query", "Search by gene name", "Not working yet", width=250))
	 	),
		mainPanel(
			strong(textOutput("regionText")),
			textOutput('legend'),
			plotlyOutput('first'),
			hr(),
			strong(textOutput("manhText")),
			plotlyOutput('second')
		)
	)
)


