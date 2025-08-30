# app.R — Needleman–Wunsch (Global Alignment) visualizer for proteins using BLOSUM50
# ------------------------------------------------------------------------------
# Features
# - Input two protein sequences (single-letter codes)
# - Choose gap model: Linear or Affine (Gotoh)
# - Uses Biostrings::BLOSUM50 substitution matrix
# - Shows DP matrix (and M/Ix/Iy matrices for affine), optimal alignment, score
# - Highlights traceback path in the score matrix
# ------------------------------------------------------------------------------

if (!requireNamespace("shiny", quietly = TRUE)) install.packages("shiny")
if (!requireNamespace("DT", quietly = TRUE)) install.packages("DT")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("Biostrings", quietly = TRUE)) BiocManager::install("Biostrings", update = FALSE, ask = FALSE)

library(shiny)
library(DT)
library(Biostrings)

data("BLOSUM50", package = "Biostrings")
SUBMAT <- as.matrix(BLOSUM50)

AA_ALLOWED <- rownames(SUBMAT)  # expected: 20 standard residues

# --------------------------- Helpers ----------------------------------------
clean_seq <- function(x) {
  x <- toupper(gsub("[^A-Za-z]", "", x))
  chars <- strsplit(x, "", fixed = TRUE)[[1]]
  if (length(chars) == 0) return("")
  # Filter to allowed AAs only; warn via attribute
  ok <- chars %in% AA_ALLOWED
  attr(x, "dropped") <- paste0(chars[!ok], collapse = "")
  paste0(chars[ok], collapse = "")
}

score_pair <- function(a, b, submat) {
  if (!(a %in% rownames(submat)) || !(b %in% colnames(submat))) return(-10)
  submat[a, b]
}

# ------------------------ Needleman–Wunsch (Linear) -------------------------
NW_linear <- function(s1, s2, submat, gap) {
  # gap is positive (penalty magnitude). Internally add as negative.
  gap_pen <- -abs(gap)
  a <- strsplit(s1, "")[[1]]
  b <- strsplit(s2, "")[[1]]
  m <- length(a); n <- length(b)
  F <- matrix(-Inf, nrow = m + 1, ncol = n + 1)
  ptr <- matrix("", nrow = m + 1, ncol = n + 1)  # D, U, L, or "" for start

  F[1, 1] <- 0
  for (i in 2:(m + 1)) { F[i, 1] <- F[i - 1, 1] + gap_pen; ptr[i, 1] <- "U" }
  for (j in 2:(n + 1)) { F[1, j] <- F[1, j - 1] + gap_pen; ptr[1, j] <- "L" }

  for (i in 2:(m + 1)) {
    for (j in 2:(n + 1)) {
      diag <- F[i - 1, j - 1] + score_pair(a[i - 1], b[j - 1], submat)
      up   <- F[i - 1, j] + gap_pen  # gap in s2 (insert in s2)
      left <- F[i, j - 1] + gap_pen  # gap in s1
      best <- max(diag, up, left)
      F[i, j] <- best
      ptr[i, j] <- if (best == diag) "D" else if (best == up) "U" else "L"
    }
  }

  # Traceback
  i <- m + 1; j <- n + 1
  aln1 <- character(0); aln2 <- character(0)
  path <- list()
  while (i > 1 || j > 1) {
    path[[length(path) + 1]] <- c(i, j)
    move <- ptr[i, j]
    if (move == "D") {
      aln1 <- c(a[i - 1], aln1); aln2 <- c(b[j - 1], aln2); i <- i - 1; j <- j - 1
    } else if (move == "U") {
      aln1 <- c(a[i - 1], aln1); aln2 <- c("-", aln2); i <- i - 1
    } else if (move == "L") {
      aln1 <- c("-", aln1); aln2 <- c(b[j - 1], aln2); j <- j - 1
    } else break
  }
  path[[length(path) + 1]] <- c(1, 1)
  path <- rev(path)

  list(
    F = F, ptr = ptr, alignment = list(top = paste0(aln1, collapse = ""),
                                        bottom = paste0(aln2, collapse = "")),
    score = F[m + 1, n + 1], path = path
  )
}

# ---------------------- Gotoh (Affine gaps: M, Ix, Iy) ---------------------
NW_affine <- function(s1, s2, submat, gap_open, gap_extend) {
  go <- abs(gap_open); ge <- abs(gap_extend)
  a <- strsplit(s1, "")[[1]]
  b <- strsplit(s2, "")[[1]]
  m <- length(a); n <- length(b)

  negInf <- -1e9
  M  <- matrix(negInf, nrow = m + 1, ncol = n + 1)
  Ix <- matrix(negInf, nrow = m + 1, ncol = n + 1)  # gap in s2 (vertical)
  Iy <- matrix(negInf, nrow = m + 1, ncol = n + 1)  # gap in s1 (horizontal)

  ptrM  <- matrix("", nrow = m + 1, ncol = n + 1)
  ptrIx <- matrix("", nrow = m + 1, ncol = n + 1)
  ptrIy <- matrix("", nrow = m + 1, ncol = n + 1)

  M[1, 1] <- 0
  # Initialize first column (j = 0)
  for (i in 2:(m + 1)) {
    Ix[i, 1] <- -(go + (i - 2 + 1) * ge)  # -(go + (i-1)*ge)
    ptrIx[i, 1] <- if (i == 2) "M" else "Ix"
  }
  # Initialize first row (i = 0)
  for (j in 2:(n + 1)) {
    Iy[1, j] <- -(go + (j - 2 + 1) * ge)
    ptrIy[1, j] <- if (j == 2) "M" else "Iy"
  }

  for (i in 2:(m + 1)) {
    for (j in 2:(n + 1)) {
      s <- score_pair(a[i - 1], b[j - 1], submat)
      # M state
      prevs <- c(M[i - 1, j - 1], Ix[i - 1, j - 1], Iy[i - 1, j - 1])
      arg <- which.max(prevs)
      M[i, j] <- prevs[arg] + s
      ptrM[i, j] <- c("M", "Ix", "Iy")[arg]

      # Ix state (gap in s2, i advances)
      prevsIx <- c(M[i - 1, j] - (go + ge), Ix[i - 1, j] - ge)
      argIx <- which.max(prevsIx)
      Ix[i, j] <- prevsIx[argIx]
      ptrIx[i, j] <- c("M", "Ix")[argIx]

      # Iy state (gap in s1, j advances)
      prevsIy <- c(M[i, j - 1] - (go + ge), Iy[i, j - 1] - ge)
      argIy <- which.max(prevsIy)
      Iy[i, j] <- prevsIy[argIy]
      ptrIy[i, j] <- c("M", "Iy")[argIy]
    }
  }

  # Choose best final state
  finals <- c(M[m + 1, n + 1], Ix[m + 1, n + 1], Iy[m + 1, n + 1])
  st <- c("M", "Ix", "Iy")[which.max(finals)]
  score <- max(finals)

  # Traceback
  i <- m + 1; j <- n + 1; state <- st
  aln1 <- character(0); aln2 <- character(0)
  path <- list()
  while (i > 1 || j > 1) {
    path[[length(path) + 1]] <- c(i, j)
    if (state == "M") {
      prev <- ptrM[i, j]
      aln1 <- c(a[i - 1], aln1); aln2 <- c(b[j - 1], aln2)
      i <- i - 1; j <- j - 1; state <- prev
    } else if (state == "Ix") { # gap in s2
      prev <- ptrIx[i, j]
      aln1 <- c(a[i - 1], aln1); aln2 <- c("-", aln2)
      i <- i - 1; state <- prev
    } else { # Iy: gap in s1
      prev <- ptrIy[i, j]
      aln1 <- c("-", aln1); aln2 <- c(b[j - 1], aln2)
      j <- j - 1; state <- prev
    }
  }
  path[[length(path) + 1]] <- c(1, 1)
  path <- rev(path)

  list(M = M, Ix = Ix, Iy = Iy,
       alignment = list(top = paste0(aln1, collapse = ""),
                        bottom = paste0(aln2, collapse = "")),
       score = score, end_state = st, path = path)
}

# ---------------------- Table rendering with highlights ---------------------
make_matrix_table <- function(mat, s1, s2, path = NULL, title = NULL) {
  a <- c("-", strsplit(s1, "")[[1]])
  b <- c("-", strsplit(s2, "")[[1]])
  df <- as.data.frame(mat)
  colnames(df) <- b
#  rownames(df) <- a

  # Convert to character and format
  fmt <- matrix(sprintf("%d", round(mat)), nrow = nrow(mat), ncol = ncol(mat))

  # Build a data frame with rownames as first column for DT
  df_char <- data.frame(AA = a, fmt, check.names = FALSE)

  dt <- datatable(df_char, rownames = FALSE, escape = FALSE,
                  options = list(dom = 't', ordering = FALSE, pageLength = nrow(df_char)))

  if (!is.null(path)) {
    # path is list of (i,j) with 1-based indices in matrix space
    # Convert to DT cell styling (rows +1 due to AA column)
    for (p in path) {
      i <- p[1]; j <- p[2]
      dt <- formatStyle(dt, columns = j + 1, target = 'row',
                        backgroundColor = styleEqual(df_char$AA, c(rep(NA, nrow(df_char))),
                                                      default = NULL))
    }
    # Simpler: highlight by setting a specific background for those cells
    # We'll compute a style matrix mask
    mask <- matrix(FALSE, nrow = nrow(mat), ncol = ncol(mat))
    for (p in path) mask[p[1], p[2]] <- TRUE

    # Apply per-cell style using JS — simpler workaround via HTML span
    df_disp <- fmt
    for (i in seq_len(nrow(mat))) {
      for (j in seq_len(ncol(mat))) {
        if (mask[i, j]) df_disp[i, j] <- sprintf('<span style="background:#ffecb3;padding:2px 4px;border-radius:4px;">%s</span>', df_disp[i, j])
      }
    }
      df_char <- data.frame(AA = a, df_disp, check.names = FALSE)
      colnames(df_char) <- b
    dt <- datatable(df_char, rownames = FALSE, escape = FALSE,
                    options = list(dom = 't', ordering = FALSE, pageLength = nrow(df_char)))
  }

  if (!is.null(title)) dt <- dt %>% formatStyle(0, target = 'row')
  dt
}

# ------------------------------- UI ----------------------------------------
ui <- fluidPage(
  titlePanel("Needleman–Wunsch (Global Alignment) • BLOSUM50"),
  sidebarLayout(
    sidebarPanel(
      textInput("seq1", "Protein sequence A (1-letter)", value = "PAWHEAE"),
      textInput("seq2", "Protein sequence B (1-letter)", value = "HEAGAWGHEE"),
      radioButtons("gap_model", "Gap model", choices = c("Linear" = "linear", "Affine (Gotoh)" = "affine"), selected = "linear"),
      conditionalPanel(
        condition = "input.gap_model == 'linear'",
        numericInput("gap", "Gap penalty (linear)", value = 8, min = 0, step = 1)
      ),
      conditionalPanel(
        condition = "input.gap_model == 'affine'",
        numericInput("gap_open", "Gap open penalty", value = 10, min = 0, step = 1),
        numericInput("gap_ext", "Gap extension penalty", value = 1, min = 0, step = 1)
      ),
      actionButton("go", "Align"),
      br(),
      helpText("Notes: Uses Biostrings::BLOSUM50. Non-standard letters are dropped. Gap penalties are costs (positive numbers).")
    ),
    mainPanel(
      h4("Alignment Result"),
      verbatimTextOutput("aln_out"),
      h4("Score Matrix"),
      DTOutput("score_dt"),
      h4("Blosum50"),
      DTOutput("b50_dt"),
      
      conditionalPanel(
        condition = "input.gap_model == 'affine'",
        h4("Affine State Matrices (M, Ix, Iy)"),
        p("M: match/mismatch, Ix: gap in sequence B (vertical), Iy: gap in sequence A (horizontal)"),
        DTOutput("M_dt"),
        DTOutput("Ix_dt"),
        DTOutput("Iy_dt")
      )
    )
  )
)

# -------------------------------- Server -----------------------------------
server <- function(input, output, session) {
  observeEvent(input$go, {
    s1 <- clean_seq(input$seq1)
    s2 <- clean_seq(input$seq2)

    if (nchar(s1) == 0 || nchar(s2) == 0) {
      output$aln_out <- renderText("Please enter two valid protein sequences (AAs A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V).")
      output$score_dt <- renderDT(NULL)
      output$b50 <- renderDT(NULL)
      return()
    }

    if (input$gap_model == "linear") {
      res <- NW_linear(s1, s2, SUBMAT, gap = input$gap)
      output$aln_out <- renderText({
        paste0(
          res$alignment$top, "\n",
          res$alignment$bottom, "\n\n",
          "Score: ", round(res$score)
        )
      })
      #      output$b50_dt <-  renderDT(make_matrix_table(BLOSUM50,s1=rownames(BLOSUM50),s2=rownames(BLOSUM50),path = res$path))
            output$b50_dt <-  renderDT(BLOSUM50) 
 
      output$score_dt <- renderDT(make_matrix_table(res$F, s1, s2, path = res$path))
      
      output$M_dt <- renderDT(NULL); output$Ix_dt <- renderDT(NULL); output$Iy_dt <- renderDT(NULL)
    } else {
      res <- NW_affine(s1, s2, SUBMAT, gap_open = input$gap_open, gap_extend = input$gap_ext)
      output$aln_out <- renderText({
        paste0(
          res$alignment$top, "\n",
          res$alignment$bottom, "\n\n",
          "Score (", res$end_state, "): ", round(res$score)
        )
      })
      # For highlighting, show the max over states as the main matrix
      Fstar <- pmax(res$M, res$Ix, res$Iy)
      output$score_dt <- renderDT(make_matrix_table(Fstar, s1, s2, path = res$path))
      output$M_dt <- renderDT(make_matrix_table(res$M, s1, s2))
      output$Ix_dt <- renderDT(make_matrix_table(res$Ix, s1, s2))
      output$Iy_dt <- renderDT(make_matrix_table(res$Iy, s1, s2))
    }
  }, ignoreInit = TRUE)
}


