library(shiny)

# Needleman-Wunsch function
needleman_wunsch <- function(seq1, seq2, match = 1, mismatch = -1, gap = -1) {
  n <- nchar(seq1)
  m <- nchar(seq2)

  # Initialize the scoring matrix
  score_matrix <- matrix(0, nrow = n + 1, ncol = m + 1)
  score_matrix[, 1] <- (0:n) * gap
  score_matrix[1, ] <- (0:m) * gap

  # Fill the scoring matrix
  for (i in 2:(n + 1)) {
    for (j in 2:(m + 1)) {
      if (substr(seq1, i - 1, i - 1) == substr(seq2, j - 1, j - 1)) {
        match_mismatch <- match
      } else {
        match_mismatch <- mismatch
      }
      score_matrix[i, j] <- max(
        score_matrix[i - 1, j] + gap,
        score_matrix[i, j - 1] + gap,
        score_matrix[i - 1, j - 1] + match_mismatch
      )
    }
  }

  # Traceback to find the alignment
  align1 <- ""
  align2 <- ""
  i <- n + 1
  j <- m + 1

  while (i > 1 || j > 1) {
    if (i > 1 && j > 1 && score_matrix[i, j] == score_matrix[i - 1, j - 1] +
        ifelse(substr(seq1, i - 1, i - 1) == substr(seq2, j - 1, j - 1), match, mismatch)) {
      align1 <- paste0(substr(seq1, i - 1, i - 1), align1)
      align2 <- paste0(substr(seq2, j - 1, j - 1), align2)
      i <- i - 1
      j <- j - 1
    } else if (i > 1 && score_matrix[i, j] == score_matrix[i - 1, j] + gap) {
      align1 <- paste0(substr(seq1, i - 1, i - 1), align1)
      align2 <- paste0("-", align2)
      i <- i - 1
    } else {
      align1 <- paste0("-", align1)
      align2 <- paste0(substr(seq2, j - 1, j - 1), align2)
      j <- j - 1
    }
  }

  return(list(
    score = score_matrix[n + 1, m + 1],
    align1 = align1,
    align2 = align2,
    score_matrix = score_matrix
  ))
}

# UI
ui <- fluidPage(
  titlePanel("Needleman-Wunsch Algorithm Demo"),
  sidebarLayout(
    sidebarPanel(
      textInput("seq1", "Sequence 1:", value = "ATCG"),
      textInput("seq2", "Sequence 2:", value = "ATGC"),
      numericInput("match", "Match score:", value = 1),
      numericInput("mismatch", "Mismatch score:", value = -1),
      numericInput("gap", "Gap penalty:", value = -1),
      actionButton("align", "Align Sequences")
    ),
    mainPanel(
      h3("Alignment Result"),
      verbatimTextOutput("alignment"),
      h3("Scoring Matrix"),
      tableOutput("matrix")
    )
  )
)

# Server
server <- function(input, output) {
  observeEvent(input$align, {
    result <- needleman_wunsch(
      toupper(input$seq1),
      toupper(input$seq2),
      match = input$match,
      mismatch = input$mismatch,
      gap = input$gap
    )
    output$alignment <- renderPrint({
      cat("Score:", result$score, "\n")
      cat("Alignment:\n")
      cat(result$align1, "\n")
      cat(result$align2, "\n")
    })
    output$matrix <- renderTable({
      n <- nchar(input$seq1)
      m <- nchar(input$seq2)
      score_matrix <- result$score_matrix[-1, -1, drop = FALSE]
      colnames(score_matrix) <- strsplit(toupper(input$seq2), "")[[1]]
      rownames(score_matrix) <- strsplit(toupper(input$seq1), "")[[1]]
      score_matrix
    }, digits = 2)
  })
}

# Run the app
shinyApp(ui = ui, server = server)
