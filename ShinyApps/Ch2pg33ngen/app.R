#==============================================================================#
# Falconer (pg 32) Number of generations required for q to decrease
#==============================================================================#

library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("Falconer (pg 32): Number of generations required"),
   
   # Sidebar with inputs
   sidebarLayout(
      sidebarPanel(
         numericInput(inputId="q0",
                     label="Initial q (q0):",
                     value = 0.50,
                     min = 0,
                     max = 1, 
                     step=0.01),
        sliderInput(inputId="gen",
          label="Number of generations:",
          min=1, max=100, value=10,
          step=1)
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
         plotOutput("qtPlot")
      )
   )
)

server <- function(input, output) {
   
   output$qtPlot <- renderPlot({
      
      # set t and qt
      t <- seq(0, input$gen, by=1)
      qt <- input$q0 / (1 + t * input$q0)
      
      # plot t vs qt
      plot(x=t, y=qt, 
            main="q at generation t with selection against q", 
            xlab="Generation (t)", 
            ylab="q (at gen t)", 
            type="l",
            pch=16)
      
   })
   
}

# Run the application 
shinyApp(ui = ui, server = server)

