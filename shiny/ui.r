library("plotly")

ui<-fluidPage(
    titlePanel("Brain eQTL of HS rats"),
	sidebarLayout(
		sidebarPanel(
			width=3,
			fluidRow(column(3, offset=0, style='padding:5px;',  selectInput("chr", "Chr", 
											c("chr1"="chr1",
											  "chr2"="chr2", 
											  "chr3"="chr3", 
											  "chr4"="chr4", 
											  "chr5"="chr5", 
											  "chr6"="chr6", 
											  "chr7"="chr7", 
											  "chr8"="chr8", 
											  "chr9"="chr9", 
											  "chr10"="chr10", 
											  "chr11"="chr11", 
											  "chr12"="chr12", 
											  "chr13"="chr13", 
											  "chr14"="chr14", 
											  "chr15"="chr15", 
											  "chr16"="chr16", 
											  "chr17"="chr17", 
											  "chr18"="chr18", 
											  "chr19"="chr19", 
											  "chr20"="chr20" 
											  ),
											  selected=c("chr1"),
											  width=80)
							 ),
					column(5, offset=0, style='padding:5px;', numericInput("loc", "bp", 175030758, width=150, step=10000 )),
					column(3, offset=0, style='padding:5px;', numericInput("win", "± M bp",2 , width=60))
			), 
			fluidRow(
					column(4, offset=0, style='padding:5px;', selectInput("cistrans","Type", c("Both"="Both", "cis-"="cis", "trans-"="trans"))),
					column(6, offset=0, style='padding:5px;', selectInput("region", "Brain Region", 
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


