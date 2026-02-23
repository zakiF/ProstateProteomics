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
  tags$head(
    tags$style(HTML("
      #app-loading {
        position: fixed;
        z-index: 9999;
        inset: 0;
        background: #ffffff;
        display: flex;
        align-items: center;
        justify-content: center;
        padding: 24px;
      }
      .app-loading-card {
        max-width: 760px;
        border: 1px solid #dddddd;
        border-radius: 8px;
        padding: 20px 24px;
        background: #fafafa;
      }
      .app-loading-title {
        font-size: 22px;
        font-weight: 700;
        margin-bottom: 10px;
      }
      .app-loading-text {
        font-size: 14px;
        line-height: 1.45;
        color: #333333;
        margin-bottom: 14px;
      }
      .app-loading-note {
        font-size: 12px;
        color: #666666;
      }
      .app-loading-spinner {
        width: 28px;
        height: 28px;
        border: 3px solid #d9d9d9;
        border-top-color: #4a4a4a;
        border-radius: 50%;
        animation: app-spin 0.9s linear infinite;
        margin-bottom: 10px;
      }
      @keyframes app-spin {
        to { transform: rotate(360deg); }
      }
    ")),
    tags$script(HTML("
      function hideAppLoading() {
        var el = document.getElementById('app-loading');
        if (el) el.style.display = 'none';
      }
      document.addEventListener('shiny:connected', function() {
        hideAppLoading();
      });
      // Fallback: never block UI if event timing is missed
      window.setTimeout(hideAppLoading, 6000);
    "))
  ),
  tags$div(
    id = "app-loading",
    tags$div(
      class = "app-loading-card",
      tags$div(class = "app-loading-title", "Please wait while the application is loading"),
      tags$div(class = "app-loading-note", "this should take <1 minute."),
      tags$div(class = "app-loading-spinner"),
      tags$div(class = "app-loading-text", "")
    )
  ),
  titlePanel("mCRPC Prognostic Protein Explorer"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      conditionalPanel(
        condition = "input.main_tab && input.main_tab != 'Welcome'",
        tags$p("Paste protein symbols, separated by comma, newline, or spaces."),
        textAreaInput("proteins", "Protein Set", rows = 8,
                      value = "DPY30\nNFU1\nLETM1, GDF15, TPR, CST7, FKBP5, LMNB2, SDC1, KRT19"),
        actionButton("calc", "Calculate", class = "btn-primary"),
        uiOutput("input_status"),
        br(), br(),
        # downloadButton("download_pdf", "Download All Plots (PDF)")
      )
    ),
    mainPanel(
      width = 9,
      conditionalPanel(
        condition = "input.main_tab && input.main_tab != 'Welcome'",
        tags$h4("Instructions"),
        tags$ul(
          tags$li("Enter one or more protein symbols, then click Calculate."),
          tags$li("You can separate proteins using commas, spaces, or new lines.")
        )
      ),
      tabsetPanel(
        id = "main_tab",
        selected = "Welcome",
        tabPanel(
          "Welcome",
          tags$h3("Welcome!"),
          tags$h4("Experimental overview"),
          tags$p(
            "To characterize the progressive proteomic landscape of prostate cancer, we performed high-throughput ",
            "plasma proteomic profiling across localized prostate cancer (",
            tags$b("local PC"),
            "), metastatic hormone-sensitive prostate cancer (",
            tags$b("mHSPC"),
            "), and metastatic castration-resistant prostate cancer (",
            tags$b("mCRPC"),
            ")."
          ),
          tags$p(
            "We identified proteins overexpressed in mCRPC that were associated with poor survival and had potential ",
            "as therapeutic targets. These findings were then validated in tissue using a tissue microarray and an ",
            "independent tissue proteomics dataset."
          ),
          tags$h4("Using this resource"),
          tags$p(
            "This dataset is designed to support hypothesis-driven discovery and rapid clinical translation in mCRPC ",
            "proteomics. It enables users to test whether a custom protein signature carries prognostic signal and to ",
            "contextualize proteins by biological and therapeutic relevance."
          ),
          tags$p("Explore the tabs to navigate the available analyses and outputs."),
          uiOutput("welcome_overview_figure"),
          NULL
        ),
        tabPanel(
          "CPPS KM",
          tags$p(
            "A Composite Proteomic Prognostic Score (CPPS) was developed as a weighted linear combination of ",
            "protein expression values derived from multivariable Cox regression coefficients."
          ),
          tags$p(
            "Patients were stratified by the median CPPS into high- and low-risk groups, and higher CPPS was ",
            "associated with significantly worse overall survival."
          ),
          shinycssloaders::withSpinner(plotOutput("km_plot", height = "800px"), type = 4)
        ),
        tabPanel(
          "AUC",
          tags$p(
            "Time-dependent ROC analysis evaluates survival discrimination at fixed follow-up times. ",
            "Here, AUC is shown at 24 and 36 months for a clinical-only model (PSA, LDH, and ALP) and a combined clinical + CPPS model."
          ),
          tags$p(
            "In the manuscript, integrating CPPS with clinical factors improved survival discrimination compared with ",
            "clinical variables alone."
          ),
          shinycssloaders::withSpinner(plotOutput("auc_plot", width = "100%", height = "380px"), type = 4)
        ),
        # tabPanel("Protein KMs", shinycssloaders::withSpinner(plotOutput("protein_km_plot", height = "900px"), type = 4)),
        tabPanel(
          "Nomogram",
          tags$p(
            "The nomogram integrates CPPS with established clinical prognostic factors (PSA, ALP, and LDH) to ",
            "estimate individualized survival probability."
          ),
          tags$p(
            "This provides patient-level 24- and 36-month risk estimates beyond broad clinical risk-grouping alone."
          ),
          shinycssloaders::withSpinner(plotOutput("nom_plot", width = "100%", height = "380px"), type = 4)
        ),
        tabPanel(
          "Summary",
          tags$p(
            "This panel links plasma mCRPC signals to external tissue proteomics context and key protein annotations ",
            "to support biological interpretation of your selected signature."
          ),
          tags$p(
            "Tissue proteomics expression was taken from the public mCRPC tissue dataset reported by Latonen et al.",
            " (",
            tags$a(
              href = "https://www.nature.com/articles/s41467-018-03573-6",
              target = "_blank",
              "Nature Communications 2018"
            ),
            ")."
          ),
          shinycssloaders::withSpinner(plotOutput("summary_plot", height = "700px"), type = 4),
          tags$hr(),
          tags$h5("How to read this panel"),
          tags$p(tags$b("Proteins:"), " list of proteins queried."),
          tags$p(tags$b("In Tissue:"), " protein detected in Latonen et al. tissue proteomics dataset."),
          tags$p(tags$b("Biomarker:"), " prior evidence that the protein is a biomarker in prostate cancer."),
          tags$p(tags$b("Surface Protein:"), " protein annotated as a cell-surface candidate."),
          tags$p(tags$b("Druggable:"), " protein with a documented pharmacologically active drug-target interaction in DrugBank or supported by published literature."),
          tags$p(tags$b("Tissue FC:"), " log2 fold-change of mCRPC vs local PC in tissue (Latonen et al.)."),
          tags$p(tags$b("Plasma FC:"), " log2 fold-change of mCRPC vs local PC in plasma (this study)."),
          tags$p(tags$b("Scaled Exp:"), " within-platform relative expression on a 0-1 scale, used for comparison across proteins."),
          tags$p(tags$b("log2HR (OS):"), " multivariable hazard ratio for overall survival."),
          tags$p(
            tags$b("Color guide:"),
            " in expression panels, ",
            tags$span("Local PC", style = "color:#AEC7E8; font-weight:600;"),
            " and ",
            tags$span("mCRPC", style = "color:#98DF8A; font-weight:600;"),
            ". In FC panels, ",
            tags$span("blue", style = "color:#4575B4; font-weight:600;"),
            " indicates lower in mCRPC, ",
            tags$span("white", style = "color:#808080; font-weight:600;"),
            " indicates no change, and ",
            tags$span("red", style = "color:#D73027; font-weight:600;"),
            " indicates higher in mCRPC."
          )
        )
      )
    )
  )
)

server <- function(input, output, session) {
  status_text <- reactiveVal(NULL)
  status_type <- reactiveVal("message")
  full_cache <- reactiveVal(list(key = NULL, data = NULL))
  overview_fig_path <- file.path(app_dir, "www", "Fig_Shiny.png")

  output$welcome_overview_figure <- renderUI({
    if (!file.exists(overview_fig_path)) {
      return(
        tags$div(
          style = "margin: 10px 0 16px 0; padding: 10px; border: 1px solid #d9d9d9; background: #fafafa; color: #555;",
          "Figure file not found: scripts/www/Fig_Shiny.png"
        )
      )
    }
    tags$div(
      style = "margin: 10px 0 16px 0;",
      tags$img(
        src = "Fig_Shiny.png",
        style = "max-width: 100%; height: auto; border: 1px solid #d9d9d9;"
      )
    )
  })

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
  }, ignoreNULL = FALSE)

  summary_result <- eventReactive(input$calc, {
    prot <- selected_proteins()
    stats_sub <- prep$protein_stats %>%
      filter(Assay %in% prot) %>%
      distinct(Assay, .keep_all = TRUE)
    shiny::validate(shiny::need(nrow(stats_sub) > 0, "No selected proteins with available coefficients."))
    build_integration_panels_gg(stats_sub$Assay, prep = prep, cpps_method = "multivariate")
  }, ignoreNULL = FALSE)

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
    plot(a$roc_clinical, time = 24, col = "#E74C3C", lwd = 3, main = "", title = FALSE)
    plot(a$roc_combined, time = 24, col = "#2ECC71", lwd = 3, add = TRUE, title = FALSE)
    title(main = "24 months")
    legend(
      "bottomright",
      c(paste0("Clinical: ", auc24_clin), paste0("Clinical + CPPS: ", auc24_comb)),
      col = c("#E74C3C", "#2ECC71"),
      lwd = 3
    )

    plot(a$roc_clinical, time = 36, col = "#E74C3C", lwd = 3, main = "", title = FALSE)
    plot(a$roc_combined, time = 36, col = "#2ECC71", lwd = 3, add = TRUE, title = FALSE)
    title(main = "36 months")
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

      plot(a$roc_clinical, time = 24, col = "#E74C3C", lwd = 3, main = "", title = FALSE)
      plot(a$roc_combined, time = 24, col = "#2ECC71", lwd = 3, add = TRUE, title = FALSE)
      title(main = "24 months")
      legend(
        "bottomright",
        c(paste0("Clinical: ", auc24_clin), paste0("Clinical + CPPS: ", auc24_comb)),
        col = c("#E74C3C", "#2ECC71"),
        lwd = 3
      )

      plot(a$roc_clinical, time = 36, col = "#E74C3C", lwd = 3, main = "", title = FALSE)
      plot(a$roc_combined, time = 36, col = "#2ECC71", lwd = 3, add = TRUE, title = FALSE)
      title(main = "36 months")
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
