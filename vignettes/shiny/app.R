#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

source("../rogan-gladen-graph.R")

# Define UI for application that draws a histogram
ui <- fluidPage(
    # Application title
    titlePanel("Rogan Gladen Plot"),
    # Sidebar with a slider inputs
    sidebarLayout(
        fixedPanel(
            sliderInput("sens",
                        "Test Sensitivity:",
                        min = 0.5,
                        max = 1,
                        # step = 0.005,
                        value = 0.75),
            sliderInput("spec",
                    "Test Specificity",
                    min = 0.9,
                    max = 1,
                    # step = 0.001,
                    value = 0.99),
            top = 0,left = 0,right = 399,bottom = 800
        ),
        # Show a plot of the generated distribution
        fixedPanel(
           plotOutput("Plot"),
           top = 0,left = 400,right = 1400,bottom = 1000
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    output$Plot <- renderPlot({
        rogan_gladen_plot(spec=input$spec,sens=input$sens)},
        width=1000, height=1000
      )
}

# Run the application 
shinyApp(ui = ui, server = server)
