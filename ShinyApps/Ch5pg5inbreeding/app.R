#==============================================================================#
# Falconer Chapter 5: Inbreeding Coefficients
#==============================================================================#

library(shiny)
library(pedigreemm)
library(ggplot2)

# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Application title
  titlePanel("Calculate Inbreeding Coefficients"),
  
  HTML("
    <p>
    This was developed by Austin Putz at Iowa State University. 
    Please contact me with problems (aputz@iastate.edu). 
    </p>
    <p>
    The use is to calculate inbreeding coefficients for a small pedigree. 
    Mostly for learning as R cannot handle huge pedigrees. 
    </p>
    "),
  
  tags$h2("Input"),
  
  # take file input
    fluidRow(column(6,
  wellPanel(
  fileInput('file1', 'Choose CSV File',
                accept=c('text/csv',
                         'text/comma-separated-values,text/plain',
                         '.csv'))))
  ), # END wellPanel()

  tags$h2("Output"),
  
  # Show the raw pedigree
  tags$h3("Raw Pedigree (0 is missing but replaced in R)"),
  fluidRow(column(6,
  dataTableOutput("output_data"))),
  
  actionButton("button", "Run Inbreeding Coefficients", icon("paper-plane"), 
    style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
  
  tags$h3("Dataset with inbreeding coefficients"),
  fluidRow(column(4,
  dataTableOutput("output_data2"))),
  
  tags$h3("Histogram of Inbreeding Coefficients"),
  fluidRow(column(8,
  plotOutput("coefs_plot")))

)

# Define server logic required to draw a histogram
server <- function(input, output) {
   
  myData <- reactive({
    
    # set the input file
    inFile <- input$file1
    
    # read the file in
    if (is.null(inFile)) return(NULL)
    data <- read.csv(inFile$datapath, header=TRUE)
    
    # set the 0's to NA for the pedigreemm package
    data[data==0] <- NA
    data
    
  })
  
  # output the raw pedigree
  output$output_data <- renderDataTable({
    myData()
  })
  
  observeEvent(input$button, {
  inbred <- reactive({

    new_data <- myData()
    
    new_ped  <- editPed(sire=new_data[,2], 
                         dam=new_data[,3], 
                         label=new_data[,1])
    new_ped2 <- pedigree(new_ped$sire, new_ped$dam, new_ped$label)
    ind <- new_ped2@label
    inbred_coefs <- inbreeding(new_ped2)
    
    data2 <- data.frame(Animal=ind, Inbreeding=inbred_coefs)
    data2 <- data2[order(data2$Animal), ]
    data2
    
  })
  output$output_data2 <- renderDataTable({
    inbred()
  })
  output$coefs_plot <- renderPlot({
    
    # plot a histogram of inbreeding coefficients
    inbred.data <- inbred()
    ggplot(inbred.data, aes(x=Inbreeding)) +
      geom_histogram(fill="dodgerblue3", color="grey50", bins=25) +
      ggtitle("Histogram of Inbreeding Coefficients (F)") +
      xlab("Inbreeding Coefficient (F)") +
      theme(text = element_text(size=16), 
            plot.title=element_text(hjust=0.5))
      
  })
  })
}

# Run the application 
shinyApp(ui = ui, server = server)

