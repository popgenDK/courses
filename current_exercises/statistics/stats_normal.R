# Define UI for application
ui <- fluidPage(

  # Application title
  titlePanel("Normal Distribution Functions"),

  # Sidebar with inputs for the normal distribution parameters
  sidebarLayout(
    sidebarPanel(
      numericInput("mean", "Mean (μ):", value = 0, step = 0.1),
      numericInput("sd", "Standard Deviation (σ):", value = 1, min = 0.1, step = 0.1),
      numericInput("x", "Value (x) for dnorm and pnorm:", value = 0, step = 0.1),
      numericInput("q", "Quantile (q) for qnorm:", value = 0.5, min = 0, max = 1, step = 0.01)
    ),

    # Show a plot of the distribution and the results of the functions
    mainPanel(
      plotOutput("distPlot"),
      verbatimTextOutput("dnormResult"),
      verbatimTextOutput("pnormResult"),
      verbatimTextOutput("qnormResult")
    )
  )
)

# Define server logic
server <- function(input, output) {

  # Reactive expression to calculate the normal distribution
  normal <- reactive({
    mean <- input$mean
    sd <- input$sd
    x <- input$x
    q <- input$q

    list(
      dnorm = dnorm(x, mean, sd),
      pnorm = pnorm(x, mean, sd),
      qnorm = qnorm(q, mean, sd)
    )
  })

  # Output the results of the normal functions
  output$dnormResult <- renderText({
    paste("dnorm(", input$x, ",", input$mean, ",", input$sd, ") =", normal()$dnorm)
  })

  output$pnormResult <- renderText({
    paste("pnorm(", input$x, ",", input$mean, ",", input$sd, ") =", normal()$pnorm)
  })

  output$qnormResult <- renderText({
    paste("qnorm(", input$q, ",", input$mean, ",", input$sd, ") =", normal()$qnorm)
  })

  # Plot the normal distribution
  output$distPlot <- renderPlot({
    mean <- input$mean
    sd <- input$sd
    x <- seq(mean - 4*sd, mean + 4*sd, length.out = 100)
    chosen_x <- input$x

    # Create a plot
    plot(x, dnorm(x, mean, sd), type = "l", col = "blue", xlab = "x", ylab = "Density", main = "Normal Distribution")
    abline(v = chosen_x, col = "red", lwd = 2)
    points(chosen_x, dnorm(chosen_x, mean, sd), col = "red", pch = 19)
  })
}

