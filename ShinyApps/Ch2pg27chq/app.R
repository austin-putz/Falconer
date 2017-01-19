#==============================================================================#
# Falconer (pg 27) Delta q under Selection
#==============================================================================#

library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("Delta q under Selection (Falconer pg 27)"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
         sliderInput("s",
                     "Selection pressure:",
                     min = 0,
                     max = 1,
                     value = 0.01)
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
         plotOutput("deltaqPlot")
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
   
   output$deltaqPlot <- renderPlot({
      
      # sequence for m
      m.seq <- seq(0,1,by=0.01)
      
      deltaq <- (input$s * m.seq^2*(1-m.seq)) / (1 - input$s * m.seq^2)
      
      # draw the histogram with the specified number of bins
      plot(x=m.seq, y=deltaq, type="l",
        main = "Change in q (delta q) over values of q",
        ylab="delta q", xlab="q", col="steelblue")
      
   })
}

# Run the application 
shinyApp(ui = ui, server = server)

