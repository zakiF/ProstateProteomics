library(shiny)
library(tidyverse)
library(survival)
library(survminer)
library(timeROC)
library(rms)
library(patchwork)
library(ggpubr)
library(shinycssloaders)

app_dir <- if (basename(getwd()) == "scripts") normalizePath(getwd()) else normalizePath(file.path(getwd(), "scripts"))
if (!file.exists(file.path(app_dir, "app.R"))) {
  app_dir <- normalizePath(getwd())
}

source(file.path(app_dir, "functions_shiny.R"))

shiny_data_dir <- file.path(app_dir, "ShinyData")
prep_file_rds <- file.path(shiny_data_dir, "11_shiny_prep.rds")
prep_file_rdata <- file.path(shiny_data_dir, "11_shiny_prep.RData")

if (file.exists(prep_file_rds)) {
  prep <- readRDS(prep_file_rds)
} else if (file.exists(prep_file_rdata)) {
  load(prep_file_rdata)
  prep <- prep_shiny
} else {
  stop(
    "Shiny data bundle not found in: ", shiny_data_dir,
    "\nRun scripts/11_ShinyPrep.R to generate scripts/ShinyData first."
  )
}

build_nomogram_obj <- function(df_complete) {
  dd <- rms::datadist(df_complete)
  assign("dd", dd, envir = .GlobalEnv)
  options(datadist = "dd")

  nom_fit <- rms::cph(
    survival::Surv(time_to_death, death) ~ log10_psa + log10_alp + log10_ldh + CPPS,
    data = df_complete, x = TRUE, y = TRUE, surv = TRUE
  )
  s24 <- function(lp) rms::survest(nom_fit, linear.predictors = lp, times = 24)$surv
  s36 <- function(lp) rms::survest(nom_fit, linear.predictors = lp, times = 36)$surv
  rms::nomogram(
    nom_fit,
    fun = list(s24, s36),
    funlabel = c("24m Survival", "36m Survival"),
    lp = FALSE
  )
}

ui <- fluidPage(
  titlePanel("CPPS Interactive Builder"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      tags$p("Paste protein symbols (Assay names), separated by comma, newline, or spaces."),
      textAreaInput("proteins", "Protein Set", rows = 8,
                    value = "DPY30\nNFU1\nLETM1, GDF15, TPR, CST7, FKBP5, LMNB2, SDC1, KRT19"),
      actionButton("calc", "Calculate", class = "btn-primary"),
      uiOutput("input_status"),
      br(), br(),
      downloadButton("download_pdf", "Download All Plots (PDF)")
    ),
    mainPanel(
      width = 9,
      tags$h4("Instructions"),
      tags$ul(
        tags$li("Enter one or more protein symbols, then click Calculate."),
        tags$li("You can separate proteins using commas, spaces, or new lines."),
        tags$li("Review results in the tabs: Summary, CPPS KM, Protein KMs, AUC, and Nomogram."),
        tags$li("Use Download All Plots (PDF) to export all figures into one report.")
      ),
      tabsetPanel(
        tabPanel("Summary", shinycssloaders::withSpinner(plotOutput("summary_plot", height = "700px"), type = 4)),
        tabPanel("CPPS KM", shinycssloaders::withSpinner(plotOutput("km_plot", height = "800px"), type = 4)),
        tabPanel("Protein KMs", shinycssloaders::withSpinner(plotOutput("protein_km_plot", height = "900px"), type = 4)),
        tabPanel("AUC", shinycssloaders::withSpinner(plotOutput("auc_plot", width = "100%", height = "380px"), type = 4)),
        tabPanel("Nomogram", shinycssloaders::withSpinner(plotOutput("nom_plot", width = "100%", height = "380px"), type = 4))
      )
    )
  )
)

server <- function(input, output, session) {
  status_text <- reactiveVal(NULL)
  status_type <- reactiveVal("message")
  full_cache <- reactiveVal(list(key = NULL, data = NULL))

  selected_proteins <- eventReactive(input$calc, {
    prot <- parse_proteins(input$proteins)
    shiny::validate(shiny::need(length(prot) > 0, "Please enter at least one protein."))

    available_assays <- prep$protein_stats %>%
      filter(!is.na(Assay), Assay != "") %>%
      pull(Assay) %>%
      unique()

    prot_upper <- toupper(prot)
    available_upper <- toupper(available_assays)
    missing_proteins <- prot[!(prot_upper %in% available_upper)]

    if (length(missing_proteins) > 0) {
      status_text(
        paste0(
          "Not found in dataset: ",
          paste(missing_proteins, collapse = ", ")
        )
      )
      status_type("warning")
    } else {
      status_text("All pasted proteins were found in the dataset.")
      status_type("message")
    }

    prot
  })

  summary_result <- eventReactive(input$calc, {
    prot <- selected_proteins()
    stats_sub <- prep$protein_stats %>%
      filter(Assay %in% prot) %>%
      distinct(Assay, .keep_all = TRUE)
    shiny::validate(shiny::need(nrow(stats_sub) > 0, "No selected proteins with available coefficients."))
    build_integration_panels_gg(stats_sub$Assay, prep = prep, cpps_method = "multivariate")
  })

  full_analysis <- reactive({
    prot <- selected_proteins()
    key <- paste(sort(unique(prot)), collapse = "|")
    cache <- full_cache()
    if (is.null(cache$key) || !identical(cache$key, key) || is.null(cache$data)) {
      cache <- list(key = key, data = build_analysis(prot, prep = prep))
      full_cache(cache)
    }
    cache$data
  })

  output$input_status <- renderUI({
    txt <- status_text()
    if (is.null(txt) || !nzchar(txt)) return(NULL)

    color <- if (identical(status_type(), "warning")) "#B23A00" else "#1B5E20"
    tags$div(
      style = paste0(
        "margin-top:8px; padding:8px 10px; border-radius:4px; ",
        "background:#F7F7F7; border:1px solid #D9D9D9; color:", color, "; ",
        "font-size:13px; white-space:normal;"
      ),
      txt
    )
  })

  output$summary_plot <- renderPlot({
    req(summary_result())
    print(summary_result())
  })

  output$km_plot <- renderPlot({
    req(full_analysis())
    print(full_analysis()$km_plot)
  })

  output$protein_km_plot <- renderPlot({
    req(full_analysis())
    km_plots <- full_analysis()$protein_km_plots
    shiny::validate(shiny::need(length(km_plots) > 0, "No individual protein KM plots available for selected proteins."))
    pg <- ggarrange(plotlist = km_plots[1:min(length(km_plots), 3)], ncol = 1, nrow = 3)
    print(pg)
  })

  output$auc_plot <- renderPlot({
    req(full_analysis())
    a <- full_analysis()
    auc24_clin <- round(a$auc_tbl$AUC_Clinical[a$auc_tbl$time_months == 24], 2)
    auc24_comb <- round(a$auc_tbl$AUC_Clinical_CPPS[a$auc_tbl$time_months == 24], 2)
    auc36_clin <- round(a$auc_tbl$AUC_Clinical[a$auc_tbl$time_months == 36], 2)
    auc36_comb <- round(a$auc_tbl$AUC_Clinical_CPPS[a$auc_tbl$time_months == 36], 2)

    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    plot(a$roc_clinical, time = 24, col = "#E74C3C", lwd = 3, main = "AUC at 24 months")
    plot(a$roc_combined, time = 24, col = "#2ECC71", lwd = 3, add = TRUE)
    legend(
      "bottomright",
      c(paste0("Clinical: ", auc24_clin), paste0("Clinical + CPPS: ", auc24_comb)),
      col = c("#E74C3C", "#2ECC71"),
      lwd = 3
    )

    plot(a$roc_clinical, time = 36, col = "#E74C3C", lwd = 3, main = "AUC at 36 months")
    plot(a$roc_combined, time = 36, col = "#2ECC71", lwd = 3, add = TRUE)
    legend(
      "bottomright",
      c(paste0("Clinical: ", auc36_clin), paste0("Clinical + CPPS: ", auc36_comb)),
      col = c("#E74C3C", "#2ECC71"),
      lwd = 3
    )
  })

  output$nom_plot <- renderPlot({
    req(full_analysis())
    nom <- build_nomogram_obj(full_analysis()$df_complete)
    par(mar = c(4, 4, 2, 1))
    plot(nom, main = "Nomogram: Clinical + CPPS")
  })

  output$download_pdf <- downloadHandler(
    filename = function() paste0("CPPS_Shiny_Report_", Sys.Date(), ".pdf"),
    content = function(file) {
      req(full_analysis())
      a <- full_analysis()

      pdf(file, width = 13.5, height = 8.2, onefile = TRUE)

      print(summary_result())
      print(a$km_plot)

      if (length(a$protein_km_plots) > 0) {
        idx <- split(seq_along(a$protein_km_plots), ceiling(seq_along(a$protein_km_plots) / 3))
        for (ii in idx) {
          pg <- ggarrange(plotlist = a$protein_km_plots[ii], ncol = 1, nrow = 3)
          print(pg)
        }
      }

      par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
      auc24_clin <- round(a$auc_tbl$AUC_Clinical[a$auc_tbl$time_months == 24], 2)
      auc24_comb <- round(a$auc_tbl$AUC_Clinical_CPPS[a$auc_tbl$time_months == 24], 2)
      auc36_clin <- round(a$auc_tbl$AUC_Clinical[a$auc_tbl$time_months == 36], 2)
      auc36_comb <- round(a$auc_tbl$AUC_Clinical_CPPS[a$auc_tbl$time_months == 36], 2)

      plot(a$roc_clinical, time = 24, col = "#E74C3C", lwd = 3, main = "AUC at 24 months")
      plot(a$roc_combined, time = 24, col = "#2ECC71", lwd = 3, add = TRUE)
      legend(
        "bottomright",
        c(paste0("Clinical: ", auc24_clin), paste0("Clinical + CPPS: ", auc24_comb)),
        col = c("#E74C3C", "#2ECC71"),
        lwd = 3
      )

      plot(a$roc_clinical, time = 36, col = "#E74C3C", lwd = 3, main = "AUC at 36 months")
      plot(a$roc_combined, time = 36, col = "#2ECC71", lwd = 3, add = TRUE)
      legend(
        "bottomright",
        c(paste0("Clinical: ", auc36_clin), paste0("Clinical + CPPS: ", auc36_comb)),
        col = c("#E74C3C", "#2ECC71"),
        lwd = 3
      )

      nom <- build_nomogram_obj(a$df_complete)
      par(mar = c(4, 4, 2, 1))
      plot(nom, main = "Nomogram: Clinical + CPPS")

      dev.off()
    }
  )
}

shinyApp(ui, server)
