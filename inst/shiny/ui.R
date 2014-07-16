library(shiny)
shinyUI(navbarPage("Shiny Distance",
                   tabPanel("Features", titlePanel("Features"),
                            sidebarLayout(
                              sidebarPanel(
                                fluidRow(
                                  column(6,
                                      radioButtons("type",strong("Sample type"),choices=c("Line","Point"),selected="Line"), hr(),
                                      radioButtons("numobs",strong("Configuration"),choices=c("Single","Double"),selected="Single")),
                                  column(6, 
                                      radioButtons("obstype",strong("Observation type"),choices=c("Single","Cluster"),selected="Single"), hr(),
                                      conditionalPanel(
                                                 condition = "input.numobs == 'Double'",
                                                  radioButtons("config","",choices=c("Independent","Trial","Removal"),selected="Independent")
                                      ))
                                
                              ),hr(),
                                fluidRow(radioButtons("disttype",strong("Distance Measurement"),choices=c("Perpendicular","Radial&Angle"),selected="Perpendicular")
                                )
                               ),
                              
                              # Show nothing; only input on this page
                              mainPanel()
                            )),
                   tabPanel("Data",titlePanel("Data"),
                            sidebarLayout(
                              sidebarPanel(
                                fileInput("regionfile",strong("Regions")),
                                fileInput("samplefile",strong("Samples")),
                                fileInput("datafile",strong("Observations")),
                                fileInput("linkfile",strong("Obs-Sample-Region Links")),
                                checkboxInput("binned","Binned Distances",value=FALSE),
                                numericInput("width","Width",value=0,min=0),
                                conditionalPanel(
                                  condition = "input.binned",
                                  textInput("cutpoints", "Bin Cut Points",value="")),
                                conditionalPanel(
                                  condition = "input.binned",
                                  actionButton("dobins", "Create bins")),
                                hr()
                              ),
                              
                              # Show a plot of the generated distribution
                              mainPanel(dataTableOutput("datatable") )
                            )),
                   tabPanel("Detection Function",titlePanel("Detection Function"),
                            sidebarLayout(
                              sidebarPanel(
                                conditionalPanel(condition = "input.numobs == 'Double'",
                                                radioButtons("indeptype",strong("Independence Assumption"),choices=c("Full","Point"),selected="Point")
                                ),                                
                                conditionalPanel(condition = "input.indeptype == 'Point' | input.numobs=='Single'",
                                    hr(), strong("Distance Model (dsmodel)"),
                                    selectInput("keyfct","Key Function",choices=c("Uniform","Half-normal","Hazard rate","Gamma"),selected="Half-normal"),
                                    selectInput("adjfct","Adjustment Function",choices=c("None","Cosine","Simple Polynomial","Hermite polynomial"),selected="None"),
                                    conditionalPanel(
                                        condition = "input.adjfct != 'None'",
                                        textInput("order", "Order")),
                                        conditionalPanel(
                                            condition = "input.adjfct == 'None' & input.keyfct!='Uniform'",
                                            textInput("scale.formula", "Scale Formula",value="~1"))
                                        ),
                                conditionalPanel(condition = "input.numobs=='Double'",
                                                 hr(), strong("Mark-Recapture Model (mrmodel)"),
                                                 textInput("mr.formula", "MR Formula",value="~distance")),                              
                                actionButton("goButton", "Fit model") ,    
                                hr(),
                                conditionalPanel(
                                  condition = "!input.binned",
                                  sliderInput("bins",
                                              "Number of bins:",
                                              min = 1,
                                              max = 50,
                                              value = 4))
                              ),
                              
                              mainPanel(
                              tabsetPanel(type = "tabs", 
                                          tabPanel("Summary",verbatimTextOutput("fit")),
                                          tabPanel("Histogram Plots", plotOutput("distPlot")),  
                                          tabPanel("GOF", verbatimTextOutput("gofvalues")),
                                                    #htmlOutput("P"),
                                                    #hr(),
                                                    #dataTableOutput("gof"))
                                          tabPanel("Q-Q Plot", htmlOutput("qqstats"),plotOutput("qqplot"))
                              )
                            )
                            )),
                   tabPanel("Abundance",titlePanel("Abundance"),
                              verbatimTextOutput("abundance")),
                   tabPanel("Record",titlePanel("Record"),
                            sidebarLayout(
                              sidebarPanel(actionButton("stopButton", "Record and stop")  ),mainPanel())
                   )
                   
))
