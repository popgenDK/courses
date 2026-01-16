
# Define UI for application
ui <- fluidPage(

  # Application title
  titlePanel("Binomial Distribution Functions"),

  # Sidebar with inputs for the binomial distribution parameters
  sidebarLayout(
    sidebarPanel(
      numericInput("size", "Size (n):", value = 10, min = 1),
      numericInput("prob", "Probability (p):", value = 0.5, min = 0, max = 1),
      numericInput("x", "Value (x) for dbinom and pbinom:", value = 5, min = 0),
      numericInput("q", "Quantile (q) for qbinom:", value = 0.5, min = 0, max = 1)
    ),

    # Show a plot of the distribution and the results of the functions
    mainPanel(
      plotOutput("distPlot"),
      verbatimTextOutput("dbinomResult"),
      verbatimTextOutput("pbinomResult"),
      verbatimTextOutput("qbinomResult")
    )
  )
)

# Define server logic
server <- function(input, output) {

  # Reactive expression to calculate the binomial distribution
  binomial <- reactive({
    size <- input$size
    prob <- input$prob
    x <- input$x
    q <- input$q

    list(
      dbinom = dbinom(x, size, prob),
      pbinom = pbinom(x, size, prob),
      qbinom = qbinom(q, size, prob)
    )
  })

  # Output the results of the binomial functions
  output$dbinomResult <- renderText({
    paste("dbinom(", input$x, ",", input$size, ",", input$prob, ") =", binomial()$dbinom)
  })

  output$pbinomResult <- renderText({
    paste("pbinom(", input$x, ",", input$size, ",", input$prob, ") =", binomial()$pbinom)
  })

  output$qbinomResult <- renderText({
    paste("qbinom(", input$q, ",", input$size, ",", input$prob, ") =", binomial()$qbinom)
  })

  # Plot the binomial distribution
  output$distPlot <- renderPlot({
    size <- input$size
    prob <- input$prob
    x <- 0:size
    chosen_x <- input$x

    # Create a bar plot
    barplot(dbinom(x, size, prob), names.arg = x, xlab = "x", ylab = "Probability", main = "Binomial Distribution",
            col = ifelse(x == chosen_x, "red", "blue"))
  })
}
