
#source("https://raw.githubusercontent.com/popgenDK/courses/refs/heads/main/BSA/dotplotShiny.R")
#shinyApp(ui, server)

# app.R
# Shiny dot plot app for DNA / RNA / Protein (seqinr) with:
# - .fsa support
# - plot updates only on button
# - self-compare option (FASTA 1 vs itself)
# - no point size/shape controls
# - show (truncated) sequences under the plot

# First-time setup:
# install.packages(c("shiny", "tools"))
# install.packages("seqinr")

library(shiny)
library(seqinr)
library(tools)

# ---------- helpers ----------

read_fasta_choices <- function(path, seqtype) {
  if (is.null(path) || !file.exists(path)) return(character(0))
  fa <- read.fasta(file = path, as.string = TRUE, seqtype = seqtype, forceDNAtolower = FALSE)
  recs <- names(fa)
  setNames(recs, recs)
}

get_seq_string <- function(path, record, seqtype) {
  req(path, record)
  fa <- read.fasta(file = path, as.string = TRUE, seqtype = seqtype, forceDNAtolower = FALSE)
  if (!(record %in% names(fa))) return(NULL)
  s <- getSequence(fa[[record]], as.string = TRUE)[[1]]
  toupper(s)
}

sanitize_seq <- function(chars, mode = c("DNA","RNA","Protein")) {
  mode <- match.arg(mode)
  chars <- toupper(chars)
  if (mode == "DNA") {
    keep <- c("A","C","G","T","R","Y","S","W","K","M","B","D","H","V","N")
    chars[!chars %in% keep] <- "N"
  } else if (mode == "RNA") {
    keep <- c("A","C","G","U","R","Y","S","W","K","M","B","D","H","V","N")
    chars[!chars %in% keep] <- "N"
  } else { # Protein
    keep <- c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V",
              "B","Z","X","*","-")
    chars[!chars %in% keep] <- "X"
  }
  chars
}

rev_comp_dna_string <- function(s) {
  vec <- s2c(s)
  rc <- rev(comp(vec))  # DNA complement
  c2s(rc)
}

truncate_100 <- function(s) {
  if (nchar(s) <= 100) return(s)
  paste0(substr(s, 1, 100), "...")
}

# Safe example sequences
exDNA1 <- "ATGCGTACGTTGACCTGATGATCGTAGCTAGCTGGATGATCGTGACCTGAT"
exDNA2 <- "ATGCGTTCGTTGTCCTGCTGATCGTGGCTAGATGCATGATCGTGACCTGAT"
exRNA1 <- chartr("T","U", exDNA1)
exRNA2 <- chartr("T","U", exDNA2)
exAA1  <- "MKTFFVLLLCTFTAFSSALA"
exAA2  <- "MKTLFVLLLCTFSALSSALA"

# ---------- UI ----------

ui <- fluidPage(
  titlePanel("Sequence Dot Plot (DNA / RNA / Protein) — seqinr"),
  sidebarLayout(
    sidebarPanel(
      h4("Inputs"),
      selectInput("seqtype", "Sequence type",
                  choices = c("DNA","RNA","Protein"), selected = "DNA"),

      fileInput("f1", "FASTA 1", accept = c(".fa",".fasta",".fas",".fsa")),
      uiOutput("rec1"),
      tags$hr(),

      checkboxInput("selfCompare", "Compare FASTA 1 with itself (ignore FASTA 2)", FALSE),

      conditionalPanel(
        condition = "!input.selfCompare",
        tagList(
          fileInput("f2", "FASTA 2", accept = c(".fa",".fasta",".fas",".fsa")),
          uiOutput("rec2")
        )
      ),

      conditionalPanel(
        condition = "input.seqtype == 'DNA' && !input.selfCompare",
        checkboxInput("rc2", "Reverse-complement FASTA 2 (DNA only)", FALSE)
      ),

      tags$hr(),
      numericInput("wsize", "Window size (wsize)", value = 1, min = 1, step = 1),
      numericInput("nmatch", "Match threshold (nmatch)", value = 1, min = 1, step = 1),

      checkboxInput("maskUnknown", "Mask unknowns (DNA/RNA: N; Protein: X/*/-)", TRUE),

      tags$hr(),
      actionButton("plotBtn", "Make dot plot", class = "btn-primary"),
      br(), br(),
      downloadButton("dl_png", "Download PNG"),
      width = 4
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Plot",
                 plotOutput("dotplot", height = "700px"),
                 tags$hr(),
                 h4("Sequences (first 100 characters)"),
                 verbatimTextOutput("seq_preview_1"),
                 verbatimTextOutput("seq_preview_2")
        ),
        tabPanel("Summary", verbatimTextOutput("summary")),
        tabPanel("Help",
          div(
            h4("How it works"),
            p("Compares two sequences and places a dot when short windows match exactly."),
            tags$ul(
              tags$li(strong("wsize:"), " sliding window length (k-mer or k-peptide)"),
              tags$li(strong("nmatch:"), " number of matching characters within each window"),
              tags$li(strong("Reverse-complement:"), " available for DNA (when not self-comparing)")
            ),
            h5("Tips"),
            tags$ul(
              tags$li("For long sequences, start with wsize ≥ 11; set nmatch close to wsize."),
              tags$li("For proteins, start with wsize 2–4; increase to reduce noise."),
              tags$li("Masking unknowns avoids spurious matches from ambiguous symbols.")
            )
          )
        )
      )
    )
  )
)

# ---------- server ----------

server <- function(input, output, session) {

  # Map UI selection to seqinr::read.fasta seqtype
  seqtype_to_seqinr <- reactive({
    switch(input$seqtype,
           "DNA"     = "DNA",
           "RNA"     = "RNA",
           "Protein" = "AA")
  })

  # Record selectors update with file + type
  observe({
    choices <- read_fasta_choices(input$f1$datapath, seqtype_to_seqinr())
    output$rec1 <- renderUI({
      if (length(choices)) {
        selectInput("rec1_sel", "Record (FASTA 1)", choices = choices, selected = choices[[1]])
      } else {
        helpText("Upload FASTA 1 to choose a record. Using example sequence if none provided.")
      }
    })
  })
  observe({
    choices <- read_fasta_choices(input$f2$datapath, seqtype_to_seqinr())
    output$rec2 <- renderUI({
      if (length(choices)) {
        selectInput("rec2_sel", "Record (FASTA 2)", choices = choices, selected = choices[[1]])
      } else {
        helpText("Upload FASTA 2 to choose a record. Using example sequence if none provided or if self-compare is on.")
      }
    })
  })

  # Example fallbacks by type
  example_pair <- reactive({
    switch(input$seqtype,
      "DNA"     = list(s1 = exDNA1, s2 = exDNA2),
      "RNA"     = list(s1 = exRNA1, s2 = exRNA2),
      "Protein" = list(s1 = exAA1,  s2 = exAA2)
    )
  })

  seq1_string <- reactive({
    if (!is.null(input$f1) && !is.null(input$rec1_sel)) {
      get_seq_string(input$f1$datapath, input$rec1_sel, seqtype_to_seqinr())
    } else {
      example_pair()$s1
    }
  })

  seq2_string_raw <- reactive({
    # Raw second sequence BEFORE self-compare / reverse complement handling
    if (!isTRUE(input$selfCompare)) {
      if (!is.null(input$f2) && !is.null(input$rec2_sel)) {
        get_seq_string(input$f2$datapath, input$rec2_sel, seqtype_to_seqinr())
      } else {
        example_pair()$s2
      }
    } else {
      # Self-compare: use FASTA 1 as source
      seq1_string()
    }
  })

  seq2_string <- reactive({
    s <- seq2_string_raw()
    if (identical(input$seqtype, "DNA") && isTRUE(input$rc2) && !isTRUE(input$selfCompare)) {
      s <- rev_comp_dna_string(s)
    }
    s
  })

  seq1_vec <- reactive({
    s <- seq1_string()
    v <- s2c(toupper(s))
    sanitize_seq(v, input$seqtype)
  })

  seq2_vec <- reactive({
    s <- seq2_string()
    v <- s2c(toupper(s))
    sanitize_seq(v, input$seqtype)
  })

  # Summary panel (updates live; plotting waits for button)
  output$summary <- renderPrint({
    s1 <- seq1_string(); s2 <- seq2_string()
    cat("Sequence type:", input$seqtype, "\n\n")
    cat("FASTA 1 record:", ifelse(is.null(input$rec1_sel), "(example)", input$rec1_sel), "\n")
    cat("Length:", nchar(s1), "\n\n")

    if (isTRUE(input$selfCompare)) {
      cat("Comparing FASTA 1 with itself\n\n")
    } else {
      cat("FASTA 2 record:", ifelse(is.null(input$rec2_sel), "(example)", input$rec2_sel), "\n")
      cat("Length:", nchar(s2), "\n\n")
    }

    cat("Parameters:\n")
    cat("  wsize   :", input$wsize, "\n")
    cat("  nmatch  :", input$nmatch, "\n")
    if (identical(input$seqtype, "DNA") && !isTRUE(input$selfCompare))
      cat("  revComp2:", isTRUE(input$rc2), "\n")
    cat("  maskUnknown:", isTRUE(input$maskUnknown), "\n\n")
  })

  # --- Plot updates only when the button is clicked ---
  plot_trigger <- reactiveVal(0)
  observeEvent(input$plotBtn, { plot_trigger(isolate(plot_trigger()) + 1) })

  maybe_mask_unknowns <- function(v) {
    if (!isTRUE(input$maskUnknown)) return(v)
    if (input$seqtype %in% c("DNA","RNA")) {
      v[v == "N"] <- "-"         # change to non-matching symbol
    } else {
      v[v %in% c("X","*","-")] <- "#"
    }
    v
  }

  make_plot <- function(filename = NULL, width = 1200, height = 1000, res = 120) {
    if (!is.null(filename)) png(filename, width = width, height = height, res = res)
    on.exit({ if (!is.null(filename)) dev.off() }, add = TRUE)

    s1 <- maybe_mask_unknowns(seq1_vec())
    s2 <- maybe_mask_unknowns(seq2_vec())

    # Guard for huge matrices with tiny windows
    L1 <- length(s1); L2 <- length(s2)
    if (L1 * L2 > 3e7 && input$wsize < 10) {
      par(mar = c(5, 5, 4, 2) + 0.1)
      plot.new()
      title(main = "Sequences too large for small wsize",
            sub = "Increase 'wsize' (>=10) or subset sequences.",
            cex.main = 1.2)
      return(invisible())
    }

    xlab <- ifelse(is.null(input$rec1_sel), paste0("Seq1 (", input$seqtype, ", example)"), input$rec1_sel)
    ylab <- if (isTRUE(input$selfCompare)) {
      paste0("Seq1 (", input$seqtype, ", self-compare)")
    } else {
      ifelse(is.null(input$rec2_sel), paste0("Seq2 (", input$seqtype, ", example)"), input$rec2_sel)
    }

    main <- sprintf("Dot plot (%s, wsize=%d, nmatch=%d)%s",
                    input$seqtype,
                    as.integer(input$wsize),
                    as.integer(input$nmatch),
                    ifelse(identical(input$seqtype,"DNA") && isTRUE(input$rc2) && !isTRUE(input$selfCompare),
                           "  [Seq2 rev-comp]", ""))

    par(mar = c(5, 5, 4, 2) + 0.1)
    dotPlot(
      s1, s2,
      wsize  = max(1L, as.integer(input$wsize)),
      nmatch = max(1L, as.integer(input$nmatch)),
      xlab   = xlab,
      ylab   = ylab,
      main   = main
    )
    box()
  }

  # Plot only when button clicked
  output$dotplot <- renderPlot({
    plot_trigger() # dependency
    make_plot()
  })

  # Sequence previews (also tied to button press)
  output$seq_preview_1 <- renderText({
    plot_trigger() # dependency to show snapshot at plotting time
    paste0("Seq1: ", truncate_100(seq1_string()))
  })

  output$seq_preview_2 <- renderText({
    plot_trigger() # dependency
    label <- if (isTRUE(input$selfCompare)) "Seq2 (self): " else "Seq2: "
    paste0(label, truncate_100(seq2_string()))
  })

  # Download handler uses the same plotting function
  output$dl_png <- downloadHandler(
    filename = function() {
      rec1 <- ifelse(is.null(input$rec1_sel), "seq1", input$rec1_sel)
      rec2 <- if (isTRUE(input$selfCompare)) rec1 else ifelse(is.null(input$rec2_sel), "seq2", input$rec2_sel)
      suffix <- if (identical(input$seqtype,"DNA") && isTRUE(input$rc2) && !isTRUE(input$selfCompare)) "_revcomp" else ""
      sprintf("dotplot_%s_%s_%s_w%d_n%d%s.png",
              input$seqtype,
              make.names(rec1),
              make.names(rec2),
              as.integer(input$wsize),
              as.integer(input$nmatch),
              suffix)
    },
    content = function(file) {
      make_plot(filename = file, width = 2000, height = 1600, res = 180)
    }
  )
}


