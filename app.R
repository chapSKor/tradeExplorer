# app.R
# Pharma Trade Explorer (t-SNE / UMAP)
# - Local CSV picker (no upload dialog)
# - Robust feature construction for mixed numeric + categorical data
#   * Handles spaces in column names (e.g., `Unit of measure`)
#   * Drops categorical columns with <2 levels in the current subset
# - Optional filter: drop rows with Quantity == 0 before embedding
# - UMAP/t-SNE with progress + error notifications
# - Memory management: frees large matrices after embedding + "Clear embedding"
# - Fast plotting via plotly scattergl (WebGL)
# - Plot controls: point size, density contours, convex hulls by commodity,
#   medoid "centroids", and zoom subplot (when coloring by Commodity)
#
# FIXES INCLUDED (per your last message):
# 1) Zoom plot color now matches master plot (stable global color map).
# 2) Convex hulls won’t nuke the plot (guards against degenerate polygons).
# 3) "Centroids" now use a medoid-like center and are drawn on top of points.

# ---- Packages ----
library(shiny)
library(dplyr)
library(readr)
library(Matrix)
library(irlba)
library(plotly)
library(tidyr)
library(MASS)        # kde2d for density contours
library(Rtsne)       # t-SNE
library(uwot)        # UMAP

# ---------- Helpers ----------

# ---- Comparison helpers ----

safe_num <- function(x) x[is.finite(x) & !is.na(x)]

# Cliff's delta: P(X>Y) - P(X<Y)
cliffs_delta <- function(x, y, max_n = 5000) {
  x <- safe_num(x); y <- safe_num(y)
  if (length(x) < 2 || length(y) < 2) return(NA_real_)
  
  # downsample to keep outer() from blowing up memory
  if (length(x) > max_n) x <- sample(x, max_n)
  if (length(y) > max_n) y <- sample(y, max_n)
  
  mat <- outer(x, y, "-")
  (sum(mat > 0) - sum(mat < 0)) / (length(x) * length(y))
}


# 1D Wasserstein distance (empirical) via quantiles
wasserstein_1d <- function(x, y, ngrid = 200) {
  x <- sort(safe_num(x)); y <- sort(safe_num(y))
  if (length(x) < 2 || length(y) < 2) return(NA_real_)
  p <- seq(0, 1, length.out = ngrid)
  qx <- stats::quantile(x, probs = p, type = 7, names = FALSE)
  qy <- stats::quantile(y, probs = p, type = 7, names = FALSE)
  mean(abs(qx - qy))
}

# Build shared histogram probabilities for binned metrics (JSD + overlap)
hist_probs <- function(x, breaks) {
  x <- safe_num(x)
  if (length(x) < 2) return(rep(0, length(breaks) - 1))
  h <- hist(x, breaks = breaks, plot = FALSE)
  p <- h$counts / sum(h$counts)
  p[is.na(p)] <- 0
  p
}

# Jensen–Shannon divergence (base-2, bounded [0,1] for discrete distributions)
jsd <- function(p, q, eps = 1e-12) {
  p <- p + eps; q <- q + eps
  p <- p / sum(p); q <- q / sum(q)
  m <- 0.5 * (p + q)
  kl <- function(a, b) sum(a * log2(a / b))
  0.5 * kl(p, m) + 0.5 * kl(q, m)
}

# Overlap coefficient: sum(min(p, q)) in [0,1]
overlap_coeff <- function(p, q) {
  sum(pmin(p, q))
}

# Summary stats for a vector
summ_stats <- function(x) {
  x <- safe_num(x)
  if (length(x) == 0) {
    return(list(n = 0, mean = NA, median = NA, sd = NA, iqr = NA,
                p10 = NA, p25 = NA, p75 = NA, p90 = NA))
  }
  qs <- stats::quantile(x, probs = c(0.10, 0.25, 0.75, 0.90), names = FALSE)
  list(
    n = length(x),
    mean = mean(x),
    median = stats::median(x),
    sd = stats::sd(x),
    iqr = stats::IQR(x),
    p10 = qs[1], p25 = qs[2], p75 = qs[3], p90 = qs[4]
  )
}

if ("package:MASS" %in% search()) {
  # no-op, just a reminder — we explicitly namespace dplyr verbs in helpers
}

load_trade_data <- function(path) {
  # File has an extra first line "Imports"
  df <- readr::read_csv(path, skip = 1, show_col_types = FALSE)
  
  df <- df %>%
    mutate(
      Period = as.Date(Period),
      Year = as.integer(format(Period, "%Y")),
      State = na_if(State, "N/A"),
      `Unit of measure` = na_if(`Unit of measure`, "N/A")
    )
  
  df
}

make_feature_matrix <- function(df, numeric_cols, categorical_cols, pca_k = 50, seed = 1) {
  set.seed(seed)
  
  # ---- numeric: log1p + scale ----
  num <- df %>%
    dplyr::select(dplyr::all_of(numeric_cols)) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ tidyr::replace_na(.x, 0))) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ log1p(.x))) %>%
    as.matrix()
  
  num <- scale(num)
  
  # ---- categorical: safe one-hot ----
  safe_names <- function(x) {
    ifelse(make.names(x) != x, paste0("`", x, "`"), x)
  }
  
  # Keep only categorical columns that have >=2 non-NA unique values in this df
  keep_cat <- categorical_cols[
    vapply(categorical_cols, function(cn) {
      if (!cn %in% names(df)) return(FALSE)
      vals <- df[[cn]]
      nlv  <- length(unique(vals[!is.na(vals)]))
      nlv >= 2
    }, logical(1))
  ]
  
  if (length(keep_cat) > 0) {
    df_cat <- df
    for (cn in keep_cat) {
      df_cat[[cn]] <- droplevels(factor(df_cat[[cn]]))
    }
    
    cat_terms <- safe_names(keep_cat)
    form <- as.formula(paste0("~ 0 + ", paste(cat_terms, collapse = " + ")))
    cat_sparse <- Matrix::sparse.model.matrix(form, data = df_cat)
  } else {
    cat_sparse <- Matrix::Matrix(0, nrow = nrow(df), ncol = 0, sparse = TRUE)
  }
  
  X <- cbind(num, cat_sparse)
  X <- as(X, "dgCMatrix")
  
  # ---- PCA via truncated SVD ----
  if (!is.null(pca_k) && pca_k > 0) {
    k <- min(pca_k, ncol(X) - 1, nrow(X) - 1)
    k <- max(k, 2)
    sv <- irlba::irlba(X, nv = k)
    X_pca <- sv$u %*% diag(sv$d)
    colnames(X_pca) <- paste0("PC", seq_len(ncol(X_pca)))
    return(X_pca)
  } else {
    return(as.matrix(X))
  }
}

run_tsne <- function(X, perplexity = 30, theta = 0.5, seed = 1) {
  set.seed(seed)
  X <- as.matrix(X)
  out <- Rtsne::Rtsne(
    X,
    perplexity = perplexity,
    theta = theta,
    dims = 2,
    check_duplicates = FALSE,
    verbose = FALSE
  )
  out$Y
}

run_umap <- function(X, n_neighbors = 15, min_dist = 0.1, metric = "cosine", seed = 1) {
  set.seed(seed)
  uwot::umap(
    X,
    n_neighbors = n_neighbors,
    min_dist = min_dist,
    metric = metric,
    n_components = 2,
    verbose = FALSE,
    ret_model = FALSE,
    n_threads = 1  # stability (esp. on some Mac setups)
  )
}

# Stable, high-contrast color map for categorical labels (used across main + zoom)
make_color_map <- function(levels_vec) {
  base <- c(
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
    "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"
  )
  cols <- grDevices::colorRampPalette(base)(length(levels_vec))
  setNames(cols, levels_vec)
}

# Convex hulls per group with guards (avoid degenerates)
compute_hulls <- function(df, group_col = "Commodity", xcol = "Dim1", ycol = "Dim2", min_n = 10) {
  if (!(group_col %in% names(df))) return(dplyr::tibble())
  
  df %>%
    dplyr::filter(!is.na(.data[[group_col]])) %>%
    dplyr::group_by(.data[[group_col]]) %>%
    dplyr::filter(dplyr::n() >= min_n) %>%
    dplyr::group_modify(~{
      pts <- .x %>% dplyr::distinct(.data[[xcol]], .data[[ycol]])
      if (nrow(pts) < 3) return(dplyr::tibble())
      
      idx <- chull(pts[[xcol]], pts[[ycol]])
      pts[idx, , drop = FALSE] %>%
        dplyr::mutate(!!group_col := unique(.x[[group_col]])[1])
    }) %>%
    dplyr::ungroup() %>%
    dplyr::select(dplyr::all_of(c(group_col, xcol, ycol)))
}


# "Centroids" as medoid-ish representative point (robust for UMAP curved clusters)
compute_medoids <- function(df, group_col = "Commodity", xcol = "Dim1", ycol = "Dim2", min_n = 10) {
  if (!(group_col %in% names(df))) return(dplyr::tibble())
  
  df %>%
    dplyr::filter(!is.na(.data[[group_col]])) %>%
    dplyr::group_by(.data[[group_col]]) %>%
    dplyr::filter(dplyr::n() >= min_n) %>%
    dplyr::group_modify(~{
      # Convert to matrix explicitly
      X <- as.matrix(.x[, c(xcol, ycol)])
      
      # FIX: Use base:: to avoid Matrix package conflicts
      ctr <- base::colMeans(X) 
      
      # Calculate distances to the center
      d <- (X[, 1] - ctr[1])^2 + (X[, 2] - ctr[2])^2
      
      # Return the single row that is closest to the mean (the medoid)
      .x[which.min(d), c(group_col, xcol, ycol), drop = FALSE]
    }) %>%
    dplyr::ungroup()
}

# ---------- UI ----------
ui <- fluidPage(
  
  tags$head(
    tags$style(HTML("
      .sidebar-section-1 {
        background-color: #ceeaec;
        padding: 12px;
        border-radius: 6px;
        margin-bottom: 15px;
      }
      .sidebar-section-2 {
        background-color: #dcdcdc;
        padding: 12px;
        border-radius: 6px;
        margin-bottom: 15px;
      }
    "))
  ),
  
  titlePanel("Canadian Pharma Trade Imports Explorer"),
  
  sidebarLayout(
    
    sidebarPanel(
      selectInput("file_local", "Pick a CSV in this folder", choices = character(0)),
      actionButton("refresh_files", "Refresh file list"),
      
      # --- 1) Value Histogram ---
      div(class = "sidebar-section-1",
          h3("Value Histogram"),
          
          checkboxInput("hist_enable", "Enable histogram panel", value = TRUE),
          
          uiOutput("hist_province_ui"),
          uiOutput("hist_commodity_ui"),
          uiOutput("hist_year_ui"),
          
          hr(),
          h4("Histogram grouping"),
          
          radioButtons(
            "hist_group_by",
            "Group histogram traces by",
            choices = c("Province", "Country", "State (US only)"),
            selected = "Province",
            inline = TRUE
          ),
          
          uiOutput("hist_origin_country_ui"),
          uiOutput("hist_origin_state_ui"),
          
          checkboxInput("hist_drop_zero_qty", "Drop rows with Quantity = 0 (histogram)", value = FALSE),
          
          sliderInput("hist_bins", "Bins", min = 10, max = 120, value = 40, step = 5),
          
          checkboxInput("hist_log_value", "Log10(Value) for histogram", value = TRUE),
          checkboxInput("hist_density", "Plot density instead of counts", value = FALSE),
          
          sliderInput("hist_alpha", "Transparency (alpha)", min = 0.05, max = 1, value = 0.35, step = 0.05)
      ),
      
      # --- 2) Comparisons ---
      div(class = "sidebar-section-2",
          h3("Comparisons"),
          
          selectInput(
            "cmp_mode", "Comparison mode",
            choices = c("Province vs Province", "Year vs Year (within Province)")
          ),
          
          uiOutput("cmp_prov_ui"),
          uiOutput("cmp_year_ui"),
          uiOutput("cmp_commodity_ui"),
          
          checkboxInput("cmp_use_hist_filters", "Use histogram commodity/year filters", value = TRUE),
          
          hr(),
          h4("Origin filter (optional)"),
          
          checkboxInput("cmp_origin_enable", "Filter by origin country/state", value = FALSE),
          uiOutput("cmp_origin_country_ui"),
          uiOutput("cmp_origin_state_ui"),
          
          actionButton("cmp_run", "Run comparison")
      ),
      
      # --- 3) Custom trend plot ---
      div(class = "sidebar-section-1",
          h3("Custom trend plot"),
          checkboxInput("cmp_trend_enable", "Show trend plot", value = TRUE),
          
          radioButtons(
            "cmp_trend_mode",
            "Lines represent",
            choices = c("Provinces (select commodities)" = "province_lines",
                        "Commodities (select provinces)" = "commodity_lines"),
            selected = "province_lines"
          ),
          
          uiOutput("cmp_trend_prov_ui"),
          uiOutput("cmp_trend_comm_ui"),
          checkboxInput("cmp_trend_drop_zero_qty", "Drop rows with Quantity = 0 (trend)", value = FALSE),
          checkboxInput("cmp_trend_log_value", "Log10(Value) for trend", value = TRUE),
          actionButton("cmp_trend_run", "Update trend plot")
      ),
      
      
      # --- 5) Trade Clustering Map ---
      div(class = "sidebar-section-2",
          h3("Trade Clustering Map"),
          
          # --- 4) Filters ---
          div(class = "sidebar-section-2",
              h4("Clustering Filters"),
              uiOutput("year_ui"),
              uiOutput("province_ui"),
              uiOutput("country_ui"),
              uiOutput("commodity_ui")
          ),
          
          
          radioButtons("method", "Cluster Method", choices = c("UMAP", "t-SNE"), inline = TRUE),
          
          sliderInput("sample_n", "Sample size (for speed)", min = 500, max = 50000, value = 8000, step = 500),
          numericInput("seed", "Random seed", value = 1, min = 1),
          
          hr(),
          h4("Feature construction"),
          checkboxInput("drop_zero_qty", "Drop rows with Quantity = 0 before embedding", value = FALSE),
          
          checkboxGroupInput(
            "cat_cols", "Categorical columns",
            choices = c("Commodity", "Province", "Country", "State", "Unit of measure"),
            selected = c("Commodity", "Province", "Country")
          ),
          checkboxGroupInput(
            "num_cols", "Numeric columns",
            choices = c("Value ($)", "Quantity"),
            selected = c("Value ($)", "Quantity")
          ),
          sliderInput("pca_k", "PCA dims before embedding", min = 10, max = 200, value = 50, step = 10),
          
          hr(),
          h4("Method params"),
          conditionalPanel(
            condition = "input.method == 't-SNE'",
            sliderInput("perplexity", "Perplexity", min = 5, max = 80, value = 30, step = 1),
            sliderInput("theta", "Theta (speed/accuracy)", min = 0, max = 0.9, value = 0.5, step = 0.05)
          ),
          conditionalPanel(
            condition = "input.method == 'UMAP'",
            sliderInput("n_neighbors", "n_neighbors", min = 5, max = 100, value = 15, step = 1),
            sliderInput("min_dist", "min_dist", min = 0, max = 1, value = 0.1, step = 0.05),
            selectInput("umap_metric", "Metric", choices = c("cosine", "euclidean", "manhattan"), selected = "cosine")
          ),
          
          hr(),
          h4("Plot"),
          selectInput(
            "color_by", "Color by",
            choices = c("Province", "Year", "Country", "Commodity", "State", "Unit of measure"),
            selected = "Province"
          ),
          
          fluidRow(
            column(6, actionButton("go", "Run embedding", class = "btn-primary")),
            column(6, actionButton("clear", "Clear embedding"))
          ),
          
          hr(),
          h4("Plot overlays"),
          sliderInput("pt_size", "Point size", min = 2, max = 12, value = 5, step = 1),
          
          checkboxInput("show_density", "Density contours", value = FALSE),
          sliderInput("density_bins", "Density grid resolution", min = 40, max = 200, value = 80, step = 10),
          
          checkboxInput("show_hulls", "Convex hulls (by Commodity)", value = FALSE),
          checkboxInput("show_centroids", "Centroids (by Commodity)", value = FALSE),
          
          conditionalPanel(
            condition = "input.color_by == 'Commodity'",
            selectizeInput(
              "focus_commodity",
              "Zoom panel: select commodity",
              choices = NULL,
              multiple = FALSE,
              options = list(placeholder = "Pick one")
            ),
            checkboxInput("show_zoom_panel", "Show zoomed panel", value = TRUE)
          )
      )
    ),## END sidebarPanel
    
    mainPanel(
      
      hr(),
      plotlyOutput("hist_plot", height = "500px"),
      verbatimTextOutput("hist_status"),
      
      hr(),
      h3("Comparison Metrics"),
      tableOutput("cmp_table"),
      verbatimTextOutput("cmp_note"), 
      
      hr(),
      h3("Province Trend (Median ± IQR over Years)"),
      plotlyOutput("cmp_trend_plot", height = "450px"),
      verbatimTextOutput("cmp_trend_status"),
      
      plotlyOutput("plt", height = "900px"),
      hr(),
      verbatimTextOutput("status")
      
      
    )#END mainPanel
  )#END sidebar layout
)#END fluidpage

# ---------- Server ----------
server <- function(input, output, session) {
  
  raw_data <- reactiveVal(NULL)
  embedding_store <- reactiveVal(NULL)
  
  # Populate dropdown with CSVs from working directory
  observe({
    csvs <- list.files(getwd(), pattern = "\\.csv$", full.names = FALSE)
    sel  <- if (length(csvs) > 0) csvs[1] else character(0)
    updateSelectInput(session, "file_local", choices = csvs, selected = sel)
  })
  
  observeEvent(input$refresh_files, {
    csvs <- list.files(getwd(), pattern = "\\.csv$", full.names = FALSE)
    sel  <- if (length(csvs) > 0) csvs[1] else character(0)
    updateSelectInput(session, "file_local", choices = csvs, selected = sel)
  })
  
  # Load selected CSV
  observeEvent(input$file_local, {
    req(input$file_local)
    if (!file.exists(input$file_local)) return()
    df <- load_trade_data(input$file_local)
    raw_data(df)
    embedding_store(NULL)
  }, ignoreInit = FALSE)
  
  # Filter UI components
  output$year_ui <- renderUI({
    df <- raw_data(); req(df)
    yr <- sort(unique(df$Year))
    sliderInput("year_range", "Year range",
                min = min(yr), max = max(yr),
                value = c(min(yr), max(yr)), step = 1)
  })
  
  output$province_ui <- renderUI({
    df <- raw_data(); req(df)
    choices <- sort(unique(df$Province))
    selectizeInput("province", "Province",
                   choices = choices, selected = NULL, multiple = TRUE,
                   options = list(placeholder = "All"))
  })
  
  output$country_ui <- renderUI({
    df <- raw_data(); req(df)
    choices <- sort(unique(df$Country))
    selectizeInput("country", "Country",
                   choices = choices, selected = NULL, multiple = TRUE,
                   options = list(placeholder = "All"))
  })
  
  output$commodity_ui <- renderUI({
    df <- raw_data(); req(df)
    choices <- sort(unique(df$Commodity))
    selectizeInput("commodity", "Commodity",
                   choices = choices, selected = NULL, multiple = TRUE,
                   options = list(placeholder = "All", maxOptions = 2000))
  })
  
  # Base filtering + optional drop Quantity==0
  filtered_data <- reactive({
    df <- raw_data(); req(df)
    req(input$year_range)
    
    df2 <- df %>%
      filter(Year >= input$year_range[1], Year <= input$year_range[2])
    
    if (!is.null(input$province) && length(input$province) > 0) {
      df2 <- df2 %>% filter(Province %in% input$province)
    }
    if (!is.null(input$country) && length(input$country) > 0) {
      df2 <- df2 %>% filter(Country %in% input$country)
    }
    if (!is.null(input$commodity) && length(input$commodity) > 0) {
      df2 <- df2 %>% filter(Commodity %in% input$commodity)
    }
    
    if (isTRUE(input$drop_zero_qty) && "Quantity" %in% names(df2)) {
      df2 <- df2 %>% filter(!is.na(Quantity) & Quantity != 0)
    }
    
    df2
  })
  
  # Compute + store embedding
  observeEvent(input$go, {
    
    df <- filtered_data()
    req(nrow(df) >= 200)
    
    withProgress(message = "Computing embedding…", value = 0, {
      
      incProgress(0.10, detail = "Sampling rows")
      set.seed(input$seed)
      n_take <- min(input$sample_n, nrow(df))
      idx <- sample(seq_len(nrow(df)), n_take)
      dsub <- df[idx, , drop = FALSE]
      
      req(length(input$num_cols) > 0)
      req(length(input$cat_cols) > 0)
      
      incProgress(0.40, detail = "Building features (one-hot + PCA)")
      X <- make_feature_matrix(
        dsub,
        numeric_cols = input$num_cols,
        categorical_cols = input$cat_cols,
        pca_k = input$pca_k,
        seed = input$seed
      )
      
      incProgress(0.75, detail = paste("Running", input$method))
      coords <- tryCatch({
        if (input$method == "t-SNE") {
          max_perp <- max(5, floor((nrow(X) - 1) / 3))
          perp <- min(input$perplexity, max_perp)
          run_tsne(X, perplexity = perp, theta = input$theta, seed = input$seed)
        } else {
          run_umap(X, n_neighbors = input$n_neighbors, min_dist = input$min_dist,
                   metric = input$umap_metric, seed = input$seed)
        }
      }, error = function(e) {
        showNotification(paste("Embedding failed:", e$message), type = "error", duration = 10)
        return(NULL)
      })
      
      rm(X)
      gc()
      
      req(!is.null(coords))
      
      incProgress(0.95, detail = "Finalizing")
      coords <- as.data.frame(coords)
      colnames(coords) <- c("Dim1", "Dim2")
      
      keep_cols <- intersect(
        c("Province", "Year", "Country", "State", "Unit of measure",
          "Commodity", "Value ($)", "Quantity"),
        names(dsub)
      )
      
      out <- dsub %>%
        dplyr::select(dplyr::all_of(keep_cols)) %>%
        dplyr::mutate(Dim1 = coords$Dim1, Dim2 = coords$Dim2)
      
      embedding_store(out)
      showNotification("Embedding complete.", type = "message", duration = 2)
      gc()
    })
  })
  
  # Clear embedding
  observeEvent(input$clear, {
    embedding_store(NULL)
    gc()
    showNotification("Cleared embedding and ran garbage collection.", type = "message", duration = 3)
  })
  
  # Populate zoom commodity choices when embedding exists
  observe({
    res <- embedding_store()
    if (is.null(res)) return()
    if (!("Commodity" %in% names(res))) return()
    comms <- sort(unique(res$Commodity))
    updateSelectizeInput(session, "focus_commodity", choices = comms, server = TRUE)
  })
  
  output$status <- renderPrint({
    cat("Shiny getwd():", getwd(), "\n")
    cat("CSVs here:", paste(list.files(getwd(), pattern="\\.csv$"), collapse=", "), "\n\n")
    
    df <- raw_data()
    if (is.null(df)) {
      cat("Pick a CSV from the dropdown to begin.\n")
      return()
    }
    
    cat("Loaded file:", input$file_local, "\n")
    cat("Total rows:", nrow(df), " | Columns:", ncol(df), "\n")
    
    df_f <- filtered_data()
    cat("Filtered rows (pre-sample):", nrow(df_f), "\n")
    cat("Zero-Quantity rows dropped:", ifelse(isTRUE(input$drop_zero_qty), "ON", "OFF"), "\n")
    
    if (!is.null(embedding_store())) {
      cat("Embedding stored rows:", nrow(embedding_store()), "\n")
    } else {
      cat("Embedding stored: none (click 'Run embedding')\n")
    }
    
    cat("\nTip: If it slows down, reduce Sample size, reduce PCA dims, and/or uncheck high-cardinality categories like Commodity.\n")
  })
  
  # ---- Plot (with stable colors + overlays + zoom) ----
  output$plt <- renderPlotly({
    res <- embedding_store()
    validate(need(!is.null(res), "Choose a CSV, adjust filters/features, then click 'Run embedding'."))
    
    color_col <- input$color_by
    validate(need(color_col %in% names(res), "Chosen color column not available in plotted data."))
    
    # Build hover text (full info)
    hover_text <- apply(res, 1, function(row) {
      paste(
        if ("Province" %in% names(res)) paste0("Province: ", row["Province"], "<br>") else "",
        if ("Year" %in% names(res)) paste0("Year: ", row["Year"], "<br>") else "",
        if ("Country" %in% names(res)) paste0("Country: ", row["Country"], "<br>") else "",
        if ("State" %in% names(res)) paste0("State: ", row["State"], "<br>") else "",
        if ("Unit of measure" %in% names(res)) paste0("Unit: ", row["Unit of measure"], "<br>") else "",
        if ("Commodity" %in% names(res)) paste0("Commodity: ", row["Commodity"], "<br>") else "",
        if ("Value ($)" %in% names(res)) paste0("Value ($): ", row["Value ($)"], "<br>") else "",
        if ("Quantity" %in% names(res)) paste0("Quantity: ", row["Quantity"]) else ""
      )
    })
    
    # Stable factor levels + stable color map across all traces (main + zoom + centroids)
    levels_all <- sort(unique(res[[color_col]]))
    res[[color_col]] <- factor(res[[color_col]], levels = levels_all)
    col_map <- make_color_map(levels_all)
    
    # Per-point colors (fast + stable)
    col_vec <- unname(col_map[as.character(res[[color_col]])])
    
    # MAIN scatter (single WebGL trace, with explicit per-point colors)
    p_main <- plotly::plot_ly(
      data = res,
      x = ~Dim1, y = ~Dim2,
      type = "scattergl", mode = "markers",
      marker = list(size = input$pt_size, opacity = 0.85, color = col_vec),
      text = hover_text,
      hoverinfo = "text",
      showlegend = TRUE
    )
    
    # Add legend-only traces so the legend shows correct category ↔ color mapping
    legend_df <- data.frame(lbl = levels_all, x = 0, y = 0, stringsAsFactors = FALSE)
    legend_cols <- unname(col_map[legend_df$lbl])
    
    p_main <- p_main %>%
      add_trace(
        data = legend_df,
        x = ~x, y = ~y,
        type = "scatter",
        mode = "markers",
        marker = list(size = 9, color = legend_cols),
        name = ~lbl,
        hoverinfo = "skip",
        visible = "legendonly",
        inherit = FALSE,
        showlegend = TRUE
      ) %>%
      layout(
        title = paste(input$method, "embedding"),
        width = 900,
        height = 900,
        xaxis = list(title = "Dim1"),
        yaxis = list(title = "Dim2", scaleanchor = "x"),
        legend = list(
          orientation = "h",
          x = 0, y = -0.25,
          xanchor = "left",
          yanchor = "top",
          font = list(size = 10)
        ),
        margin = list(b = 260)
      )
    
    # Density contours (added as separate contour trace)
    if (isTRUE(input$show_density)) {
      ngrid <- input$density_bins
      kd <- MASS::kde2d(res$Dim1, res$Dim2, n = ngrid)
      
      p_main <- p_main %>%
        add_trace(
          x = kd$x, y = kd$y, z = kd$z,
          type = "contour",
          showscale = FALSE,
          contours = list(coloring = "lines"),
          line = list(width = 1),
          hoverinfo = "skip",
          inherit = FALSE,
          showlegend = FALSE
        )
    }
    
    # Convex hulls (by Commodity) — safe guards so plot won't disappear
    if (isTRUE(input$show_hulls) && "Commodity" %in% names(res)) {
      hulls <- compute_hulls(res, group_col = "Commodity", xcol = "Dim1", ycol = "Dim2", min_n = 10)
      
      if (nrow(hulls) > 0) {
        for (comm in unique(hulls$Commodity)) {
          poly <- hulls %>% filter(Commodity == comm) %>% distinct(Dim1, Dim2, .keep_all = TRUE)
          if (nrow(poly) < 3) next
          poly <- bind_rows(poly, poly[1, ]) # close polygon
          
          p_main <- p_main %>%
            add_trace(
              data = poly,
              x = ~Dim1, y = ~Dim2,
              type = "scatter",
              mode = "lines",
              fill = "toself",
              opacity = 0.12,
              line = list(width = 1),
              hoverinfo = "skip",
              showlegend = FALSE,
              inherit = FALSE
            )
        }
      }
    }

    # Centroids (by Commodity) — robust layout shapes
    if (isTRUE(input$show_centroids) && "Commodity" %in% names(res)) {
      cents <- compute_medoids(res, group_col = "Commodity", xcol = "Dim1", ycol = "Dim2", min_n = 10) %>%
        dplyr::filter(is.finite(Dim1), is.finite(Dim2))
      
      if (nrow(cents) > 0) {
        # Limit clutter
        topN <- 80
        counts <- res %>%
          dplyr::filter(!is.na(Commodity)) %>%
          dplyr::count(Commodity, name = "n") %>%
          dplyr::arrange(dplyr::desc(n)) %>%
          dplyr::slice_head(n = topN)
        
        cents <- cents %>% dplyr::inner_join(counts, by = "Commodity")
        
        if (nrow(cents) > 0) {
          xr <- range(res$Dim1, na.rm = TRUE)
          yr <- range(res$Dim2, na.rm = TRUE)
          r  <- 0.008 * max(diff(xr), diff(yr))
          
          # FIX: Reference axes dynamically. 
          # If a subplot exists, these need to point to the first subplot axes.
          target_x <- if (isTRUE(input$show_zoom_panel) && identical(input$color_by, "Commodity")) "x" else "x"
          target_y <- if (isTRUE(input$show_zoom_panel) && identical(input$color_by, "Commodity")) "y" else "y"
          
          centroid_shapes <- lapply(seq_len(nrow(cents)), function(i) {
            list(
              type = "circle",
              xref = target_x, yref = target_y,
              x0 = cents$Dim1[i] - r, x1 = cents$Dim1[i] + r,
              y0 = cents$Dim2[i] - r, y1 = cents$Dim2[i] + r,
              line = list(color = "black", width = 1),
              fillcolor = "black",
              opacity = 0.5,
              layer = "above" # Changed to "above" so they don't hide under the points
            )
          })
          
          p_main <- p_main %>% layout(shapes = centroid_shapes)
        }
      }
    }
    
    
    # Zoom subplot when coloring by Commodity (preserve same colors)
    if (isTRUE(input$show_zoom_panel) &&
        identical(input$color_by, "Commodity") &&
        !is.null(input$focus_commodity) &&
        "Commodity" %in% names(res)) {
      
      zdf <- res %>% filter(Commodity == input$focus_commodity)
      
      if (nrow(zdf) >= 5) {
        # zoom ranges
        xr <- range(zdf$Dim1, na.rm = TRUE)
        yr <- range(zdf$Dim2, na.rm = TRUE)
        pad_x <- 0.1 * diff(xr); pad_y <- 0.1 * diff(yr)
        xr <- xr + c(-pad_x, pad_x)
        yr <- yr + c(-pad_y, pad_y)
        
        hover_zoom <- apply(zdf, 1, function(row) {
          paste(
            paste0("Commodity: ", row["Commodity"], "<br>"),
            if ("Value ($)" %in% names(zdf)) paste0("Value ($): ", row["Value ($)"], "<br>") else "",
            if ("Quantity" %in% names(zdf)) paste0("Quantity: ", row["Quantity"]) else ""
          )
        })
        
        # Match color exactly to master map
        zdf$Commodity <- factor(zdf$Commodity, levels = levels_all)
        zoom_col <- unname(col_map[as.character(zdf$Commodity)])
        
        p_zoom <- plotly::plot_ly(
          data = zdf,
          x = ~Dim1, y = ~Dim2,
          type = "scattergl", mode = "markers",
          marker = list(size = max(3, input$pt_size), opacity = 0.9, color = zoom_col),
          text = hover_zoom,
          hoverinfo = "text",
          showlegend = FALSE
        ) %>%
          layout(
            title = paste0("Zoom: ", input$focus_commodity),
            xaxis = list(title = "Dim1", range = xr),
            yaxis = list(title = "Dim2", range = yr, scaleanchor = "x"),
            margin = list(b = 60)
          )
        
        return(plotly::subplot(p_main, p_zoom, nrows = 1, widths = c(0.65, 0.35), shareY = FALSE))
      }
    }
    
    p_main
  })
  
  output$hist_province_ui <- renderUI({
    df <- raw_data(); req(df)
    provs <- sort(unique(df$Province))
    selectizeInput(
      "hist_province", "Province(s)",
      choices = provs,
      selected = NULL,
      multiple = TRUE,
      options = list(placeholder = "All provinces")
    )
  })
  
  output$hist_commodity_ui <- renderUI({
    df <- raw_data(); req(df)
    comms <- sort(unique(df$Commodity))
    selectizeInput(
      "hist_commodity", "Commodity(ies)",
      choices = comms,
      selected = NULL,
      multiple = TRUE,
      options = list(placeholder = "All commodities", maxOptions = 2000)
    )
  })
  
  output$hist_year_ui <- renderUI({
    df <- raw_data(); req(df)
    yrs <- sort(unique(df$Year))
    selectizeInput(
      "hist_year", "Year(s)",
      choices = yrs,
      selected = NULL,
      multiple = TRUE,
      options = list(placeholder = "All years")
    )
  })
  
  output$hist_origin_country_ui <- renderUI({
    df <- raw_data(); req(df)
    cn <- sort(unique(df$Country))
    selectizeInput(
      "hist_origin_country", "Origin country filter (optional)",
      choices = cn, selected = NULL, multiple = TRUE,
      options = list(placeholder = "All origin countries")
    )
  })
  
  output$hist_origin_state_ui <- renderUI({
    df <- raw_data(); req(df)
    # only meaningful when Country is US + state exists
    st <- sort(unique(df$State[!is.na(df$State)]))
    selectizeInput(
      "hist_origin_state", "Origin state filter (US only, optional)",
      choices = st, selected = NULL, multiple = TRUE,
      options = list(placeholder = "All US states")
    )
  })
  
  hist_data <- reactive({
    req(input$hist_enable)
    df <- raw_data(); req(df)
    
    df2 <- df
    
    # Year selection (exact years)
    if (!is.null(input$hist_year) && length(input$hist_year) > 0) {
      df2 <- df2 %>% dplyr::filter(Year %in% as.integer(input$hist_year))
    } else {
      # fallback to your existing year_range if present
      if (!is.null(input$year_range)) {
        df2 <- df2 %>% dplyr::filter(Year >= input$year_range[1], Year <= input$year_range[2])
      }
    }
    
    # Province selection
    if (!is.null(input$hist_province) && length(input$hist_province) > 0) {
      df2 <- df2 %>% dplyr::filter(Province %in% input$hist_province)
    }
    
    # Commodity selection
    if (!is.null(input$hist_commodity) && length(input$hist_commodity) > 0) {
      df2 <- df2 %>% dplyr::filter(Commodity %in% input$hist_commodity)
    }
    
    # Origin country filter (optional)
    if (!is.null(input$hist_origin_country) && length(input$hist_origin_country) > 0) {
      df2 <- df2 %>% dplyr::filter(Country %in% input$hist_origin_country)
    }
    
    # Origin state filter (optional) — typically only makes sense for US rows
    if (!is.null(input$hist_origin_state) && length(input$hist_origin_state) > 0) {
      df2 <- df2 %>% dplyr::filter(State %in% input$hist_origin_state)
    }
    
    # Drop Quantity == 0 if requested
    if (isTRUE(input$hist_drop_zero_qty) && "Quantity" %in% names(df2)) {
      df2 <- df2 %>% dplyr::filter(!is.na(Quantity) & Quantity != 0)
    }
    
    # Keep only usable rows for value
    if ("Value ($)" %in% names(df2)) {
      df2 <- df2 %>% dplyr::filter(!is.na(`Value ($)`), is.finite(`Value ($)`))
    }
    
    df2
  })
  
  output$hist_plot <- renderPlotly({
    req(input$hist_enable)
    
    df2 <- hist_data()
    validate(need(nrow(df2) > 0, "Histogram: no rows after filtering."))
    
    validate(need("Value ($)" %in% names(df2), "Histogram requires a column named 'Value ($)'."))
    
    # Value transform (optional)
    x <- df2$`Value ($)`
    if (isTRUE(input$hist_log_value)) {
      x <- log10(pmax(x, 1e-6))  # avoid log10(0)
      xlab <- "log10(Value ($))"
    } else {
      xlab <- "Value ($)"
    }
    df2 <- df2 %>% dplyr::mutate(.hist_x = x)
    
    # Choose grouping variable (default remains Province)
    group_var <- input$hist_group_by
    grp <- dplyr::case_when(
      group_var == "Province" ~ as.character(df2$Province),
      group_var == "Country"  ~ as.character(df2$Country),
      group_var == "State (US only)" ~ as.character(df2$State),
      TRUE ~ as.character(df2$Province)
    )
    
    df2 <- df2 %>% dplyr::mutate(.grp = grp)
    df2$.grp[is.na(df2$.grp) | df2$.grp == ""] <- "Unknown"
    
    # If grouping by state, restrict to US rows (prevents nonsense states)
    if (identical(group_var, "State (US only)")) {
      df2 <- df2 %>% dplyr::filter(Country %in% c("United States", "USA", "United States of America"))
    }
    
    # Cap to top-N groups for legibility
    max_groups <- 12
    grp_counts <- df2 %>% dplyr::count(.grp, name = "n") %>% dplyr::arrange(dplyr::desc(n))
    keep_grp <- head(grp_counts$.grp, max_groups)
    df_plot <- df2 %>% dplyr::filter(.grp %in% keep_grp)
    
    p <- plot_ly()
    groups <- sort(unique(df_plot$.grp))
    for (g in groups) {
      dd <- df_plot %>% dplyr::filter(.grp == g)
      p <- p %>%
        add_histogram(
          data = dd,
          x = ~.hist_x,
          name = g,
          opacity = input$hist_alpha,
          nbinsx = input$hist_bins,
          histnorm = if (isTRUE(input$hist_density)) "probability density" else ""
        )
    }

    p %>%
      layout(
        barmode = "overlay",
        xaxis = list(title = xlab),
        yaxis = list(title = if (isTRUE(input$hist_density)) "Density" else "Count"),
        legend = list(orientation = "h", x = 0, y = -0.2, xanchor = "left"),
        margin = list(b = 120)
      )
  })
  
  output$hist_status <- renderPrint({
    req(input$hist_enable)
    df2 <- hist_data()
    cat("Histogram rows:", nrow(df2), "\n")
    cat("Provinces (unique):", length(unique(df2$Province)), "\n")
    cat("Commodities (unique):", length(unique(df2$Commodity)), "\n")
    cat("Years (unique):", length(unique(df2$Year)), "\n")
  })
  
  output$cmp_prov_ui <- renderUI({
    df <- raw_data(); req(df)
    provs <- sort(unique(df$Province))
    
    if (input$cmp_mode == "Province vs Province") {
      tagList(
        selectizeInput("cmp_prov_a", "Province A", choices = provs, selected = provs[1], multiple = FALSE),
        selectizeInput("cmp_prov_b", "Province B", choices = provs, selected = provs[min(2, length(provs))], multiple = FALSE)
      )
    } else {
      # Year vs Year (within Province)
      selectizeInput("cmp_prov_one", "Province", choices = provs, selected = provs[1], multiple = FALSE)
    }
  })
  
  
  output$cmp_year_ui <- renderUI({
    df <- raw_data(); req(df)
    yrs <- sort(unique(df$Year))
    
    if (input$cmp_mode == "Province vs Province") {
      selectizeInput("cmp_years", "Years (optional filter)", choices = yrs, selected = NULL, multiple = TRUE,
                     options = list(placeholder = "All years (or use year_range/hist_year)"))
    } else {
      tagList(
        selectizeInput("cmp_year_a", "Year A", choices = yrs, selected = yrs[1], multiple = FALSE),
        selectizeInput("cmp_year_b", "Year B", choices = yrs, selected = yrs[min(2, length(yrs))], multiple = FALSE)
      )
    }
  })
  
  output$cmp_commodity_ui <- renderUI({
    df <- raw_data(); req(df)
    comms <- sort(unique(df$Commodity))
    selectizeInput(
      "cmp_commodity", "Commodity(ies) (comparison filter)",
      choices = comms,
      selected = NULL,
      multiple = TRUE,
      options = list(placeholder = "All commodities", maxOptions = 2000)
    )
  })
  
  output$cmp_origin_country_ui <- renderUI({
    req(input$cmp_origin_enable)
    df <- raw_data(); req(df)
    
    countries <- sort(unique(df$Country))
    selectizeInput(
      "cmp_origin_country",
      "Origin country (optional)",
      choices = countries,
      selected = NULL,
      multiple = FALSE,
      options = list(placeholder = "All countries")
    )
  })
  
  output$cmp_origin_state_ui <- renderUI({
    req(input$cmp_origin_enable)
    df <- raw_data(); req(df)
    
    # Only show state filter if the chosen country is US-like
    us_vals <- c("United States", "USA", "United States of America")
    if (is.null(input$cmp_origin_country) || !(input$cmp_origin_country %in% us_vals)) {
      return(NULL)
    }
    
    states <- sort(unique(df$State[!is.na(df$State)]))
    selectizeInput(
      "cmp_origin_state",
      "Origin state (US only, optional)",
      choices = states,
      selected = NULL,
      multiple = TRUE,
      options = list(placeholder = "All US states")
    )
  })
  
  
  cmp_vectors <- eventReactive(input$cmp_run, {
    df <- raw_data(); req(df)
    validate(need("Value ($)" %in% names(df), "Comparison requires 'Value ($)' column."))
    
    df2 <- df
    
    # Optionally use histogram filters
    if (isTRUE(input$cmp_use_hist_filters)) {
      # Year selection from hist_year if set; otherwise fall back to year_range if available
      if (!is.null(input$hist_year) && length(input$hist_year) > 0) {
        df2 <- df2 %>% dplyr::filter(Year %in% as.integer(input$hist_year))
      } else if (!is.null(input$year_range)) {
        df2 <- df2 %>% dplyr::filter(Year >= input$year_range[1], Year <= input$year_range[2])
      }
      
      # Commodity selection from hist_commodity if set
      if (!is.null(input$hist_commodity) && length(input$hist_commodity) > 0) {
        df2 <- df2 %>% dplyr::filter(Commodity %in% input$hist_commodity)
      }
      
      # Drop Quantity==0 using histogram toggle
      if (isTRUE(input$hist_drop_zero_qty) && "Quantity" %in% names(df2)) {
        df2 <- df2 %>% dplyr::filter(!is.na(Quantity) & Quantity != 0)
      }
    } else {
      # If not using hist filters, still respect your global year_range if present
      if (!is.null(input$year_range)) {
        df2 <- df2 %>% dplyr::filter(Year >= input$year_range[1], Year <= input$year_range[2])
      }
    }
    
    # Additional year filter for province-vs-province mode (optional)
    if (input$cmp_mode == "Province vs Province") {
      if (!is.null(input$cmp_years) && length(input$cmp_years) > 0) {
        df2 <- df2 %>% dplyr::filter(Year %in% as.integer(input$cmp_years))
      }
    }
    
    # ---- OPTIONAL origin filters (single country, then optional US state) ----
    if (isTRUE(input$cmp_origin_enable)) {
      
      # Country filter (single)
      if (!is.null(input$cmp_origin_country) && nzchar(input$cmp_origin_country)) {
        df2 <- df2 %>% dplyr::filter(Country == input$cmp_origin_country)
      }
      
      # State filter only if country is US-like and user selected states
      us_vals <- c("United States", "USA", "United States of America")
      if (!is.null(input$cmp_origin_country) &&
          input$cmp_origin_country %in% us_vals &&
          !is.null(input$cmp_origin_state) &&
          length(input$cmp_origin_state) > 0) {
        
        df2 <- df2 %>% dplyr::filter(State %in% input$cmp_origin_state)
      }
    }
    
    
    # Comparison-specific commodity filter (applies to both modes)
    if (!is.null(input$cmp_commodity) && length(input$cmp_commodity) > 0) {
      df2 <- df2 %>% dplyr::filter(Commodity %in% input$cmp_commodity)
    }
    
    
    # Value transform consistent with histogram
    xval <- df2$`Value ($)`
    if (isTRUE(input$hist_log_value)) {
      xval <- log10(pmax(xval, 1e-6))
    }
    df2 <- df2 %>% dplyr::mutate(.val = xval)
    
    # Split into A and B according to mode
    if (input$cmp_mode == "Province vs Province") {
      req(input$cmp_prov_a, input$cmp_prov_b)
      
      a <- df2 %>% dplyr::filter(Province == input$cmp_prov_a) %>% dplyr::pull(.val)
      b <- df2 %>% dplyr::filter(Province == input$cmp_prov_b) %>% dplyr::pull(.val)
      
      list(x = a, y = b,
           label_x = paste0(input$cmp_prov_a),
           label_y = paste0(input$cmp_prov_b))
    } else {
      req(input$cmp_prov_one, input$cmp_year_a, input$cmp_year_b)
      
      a <- df2 %>% dplyr::filter(Province == input$cmp_prov_one, Year == as.integer(input$cmp_year_a)) %>% dplyr::pull(.val)
      b <- df2 %>% dplyr::filter(Province == input$cmp_prov_one, Year == as.integer(input$cmp_year_b)) %>% dplyr::pull(.val)
      
      list(x = a, y = b,
           label_x = paste0(input$cmp_prov_one, " ", input$cmp_year_a),
           label_y = paste0(input$cmp_prov_one, " ", input$cmp_year_b))
    }
    
  })
  
  output$cmp_table <- renderTable({
    v <- cmp_vectors()
    x <- v$x; y <- v$y
    
    validate(need(length(safe_num(x)) >= 10, paste("Too few rows for:", v$label_x)))
    validate(need(length(safe_num(y)) >= 10, paste("Too few rows for:", v$label_y)))
    
    sx <- summ_stats(x)
    sy <- summ_stats(y)
    
    # shared breaks for binned metrics (use hist_bins, consistent with histogram)
    allv <- c(safe_num(x), safe_num(y))
    # ensure enough spread
    if (length(unique(allv)) < 5) {
      jsd_val <- NA_real_
      ovl_val <- NA_real_
    } else {
      breaks <- pretty(range(allv), n = input$hist_bins)
      p <- hist_probs(x, breaks)
      q <- hist_probs(y, breaks)
      jsd_val <- jsd(p, q)
      ovl_val <- overlap_coeff(p, q)
    }
    
    # KS statistic (guard against constant vectors)
    ksD <- NA_real_
    ksP <- NA_real_
    if (length(unique(safe_num(x))) > 1 && length(unique(safe_num(y))) > 1) {
      kt <- suppressWarnings(stats::ks.test(safe_num(x), safe_num(y)))
      ksD <- as.numeric(kt$statistic)
      ksP <- as.numeric(kt$p.value)
    }
    
    out <- data.frame(
      Metric = c(
        "n",
        "mean",
        "median",
        "sd",
        "IQR",
        "p10",
        "p25",
        "p75",
        "p90",
        "Δ median (X - Y)",
        "ratio median (X / Y)",
        "Δ mean (X - Y)",
        "Cliff's delta",
        "Wasserstein (1D)",
        "KS D",
        "KS p-value",
        "Jensen–Shannon divergence",
        "Overlap coefficient"
      ),
      X = c(
        sx$n, sx$mean, sx$median, sx$sd, sx$iqr, sx$p10, sx$p25, sx$p75, sx$p90,
        sx$median - sy$median,
        sx$median / sy$median,
        sx$mean - sy$mean,
        cliffs_delta(x, y),
        wasserstein_1d(x, y),
        ksD,
        ksP,
        jsd_val,
        ovl_val
      ),
      Y = c(
        sy$n, sy$mean, sy$median, sy$sd, sy$iqr, sy$p10, sy$p25, sy$p75, sy$p90,
        NA, NA, NA, NA, NA, NA, NA, NA, NA
      ),
      check.names = FALSE
    )
    
    names(out)[2] <- v$label_x
    names(out)[3] <- v$label_y
    
    out
  }, digits = 4, striped = TRUE, bordered = TRUE, spacing = "s")
  output$cmp_note <- renderPrint({
    v <- cmp_vectors()
    x_n <- length(safe_num(v$x))
    y_n <- length(safe_num(v$y))
    
    cat("Scale:", ifelse(isTRUE(input$hist_log_value), "log10(Value ($))", "Value ($)"), "\n")
    cat("X =", v$label_x, "  (n =", x_n, "rows)\n")
    cat("Y =", v$label_y, "  (n =", y_n, "rows)\n")
    cat("n = number of trade records used for each group after the active comparison filters.\n\n")
    
    if (x_n < 30 || y_n < 30) {
      cat("⚠ Note: One or both groups have small n; metrics may be noisy.\n\n")
    }
    
    cat("Interpretation hints:\n")
    cat("- Wasserstein: average shift between distributions (same units as the scale above).\n")
    cat("- KS D: max CDF separation (0 to 1). Larger = more different.\n")
    cat("- JSD: bounded divergence (0 to ~1). Larger = more different.\n")
    cat("- Overlap: 0–1, larger = more overlap (more similar).\n")
    cat("- Cliff’s delta: -1..1, positive means X tends to be larger than Y.\n")
    
    cat("Origin filter enabled:", ifelse(isTRUE(input$cmp_origin_enable), "YES", "NO"), "\n")
    if (isTRUE(input$cmp_origin_enable)) {
      cat("Origin country:", ifelse(!is.null(input$cmp_origin_country) && nzchar(input$cmp_origin_country),
                                    input$cmp_origin_country, "All"), "\n")
      if (!is.null(input$cmp_origin_state) && length(input$cmp_origin_state) > 0) {
        cat("US states:", paste(input$cmp_origin_state, collapse = ", "), "\n")
      } else {
        cat("US states: All (or not applicable)\n")
      }
    }
    cat("\n")
    
  })
  
  
  # ---- Trend plot: Median ± IQR over years (Province or Commodity lines) ----
  
  # Trend selectors (independent of histogram)
  output$cmp_trend_prov_ui <- renderUI({
    req(input$cmp_trend_enable)
    df <- raw_data(); req(df)
    provs <- sort(unique(df$Province))
    
    # If user wants commodity lines, they must select province(s)
    # If user wants province lines, province selection is optional (default = all)
    if (identical(input$cmp_trend_mode, "commodity_lines")) {
      selectizeInput(
        "cmp_trend_provinces",
        "Province(s) to include",
        choices = provs,
        selected = head(provs, 3),
        multiple = TRUE,
        options = list(placeholder = "Select 1+ provinces")
      )
    } else {
      selectizeInput(
        "cmp_trend_provinces",
        "Province(s) (optional)",
        choices = provs,
        selected = NULL,
        multiple = TRUE,
        options = list(placeholder = "All provinces")
      )
    }
  })
  
  output$cmp_trend_comm_ui <- renderUI({
    req(input$cmp_trend_enable)
    df <- raw_data(); req(df)
    comms <- sort(unique(df$Commodity))
    
    # If user wants province lines, they must select commodity(ies)
    # If user wants commodity lines, commodity selection is optional (default = all)
    if (identical(input$cmp_trend_mode, "province_lines")) {
      selectizeInput(
        "cmp_trend_commodities",
        "Commodity(ies) to include",
        choices = comms,
        selected = head(comms, 5),
        multiple = TRUE,
        options = list(placeholder = "Select 1+ commodities", maxOptions = 2000)
      )
    } else {
      selectizeInput(
        "cmp_trend_commodities",
        "Commodity(ies) (optional)",
        choices = comms,
        selected = NULL,
        multiple = TRUE,
        options = list(placeholder = "All commodities", maxOptions = 2000)
      )
    }
  })
  
  # Event-driven compute so it doesn't re-run constantly
  cmp_trend_data <- eventReactive(input$cmp_trend_run, {
    req(input$cmp_trend_enable)
    df <- raw_data(); req(df)
    
    validate(need(all(c("Province","Year","Commodity","Value ($)") %in% names(df)),
                  "Trend plot requires Province, Year, Commodity, and Value ($) columns."))
    
    # force the trend window (as requested)
    year_min <- 2020L
    year_max <- 2025L
    
    df2 <- df %>%
      dplyr::filter(
        !is.na(Province), !is.na(Year), !is.na(Commodity),
        !is.na(`Value ($)`), is.finite(`Value ($)`),
        Year >= year_min, Year <= year_max
      )
    
    # optional drop Quantity==0 for trend only
    if (isTRUE(input$cmp_trend_drop_zero_qty) && "Quantity" %in% names(df2)) {
      df2 <- df2 %>% dplyr::filter(!is.na(Quantity) & Quantity != 0)
    }
    
    # independent filters (NOT histogram)
    if (!is.null(input$cmp_trend_provinces) && length(input$cmp_trend_provinces) > 0) {
      df2 <- df2 %>% dplyr::filter(Province %in% input$cmp_trend_provinces)
    }
    if (!is.null(input$cmp_trend_commodities) && length(input$cmp_trend_commodities) > 0) {
      df2 <- df2 %>% dplyr::filter(Commodity %in% input$cmp_trend_commodities)
    }
    
    # enforce required selection depending on mode
    if (identical(input$cmp_trend_mode, "province_lines")) {
      validate(need(!is.null(input$cmp_trend_commodities) && length(input$cmp_trend_commodities) > 0,
                    "Select 1+ commodities to plot province trends."))
    } else {
      validate(need(!is.null(input$cmp_trend_provinces) && length(input$cmp_trend_provinces) > 0,
                    "Select 1+ provinces to plot commodity trends."))
    }
    
    # scale for trend only (independent)
    if (isTRUE(input$cmp_trend_log_value)) {
      df2 <- df2 %>% dplyr::mutate(.val = log10(pmax(`Value ($)`, 1e-6)))
      ylab <- "log10(Value ($))"
    } else {
      df2 <- df2 %>% dplyr::mutate(.val = `Value ($)`)
      ylab <- "Value ($)"
    }
    
    # summarise
    if (identical(input$cmp_trend_mode, "province_lines")) {
      # lines = Province, pooled across selected commodities
      summ <- df2 %>%
        dplyr::group_by(Province, Year) %>%
        dplyr::summarise(
          n = dplyr::n(),
          q25 = stats::quantile(.val, 0.25, na.rm = TRUE, names = FALSE),
          med = stats::median(.val, na.rm = TRUE),
          q75 = stats::quantile(.val, 0.75, na.rm = TRUE, names = FALSE),
          .groups = "drop"
        ) %>%
        dplyr::arrange(Province, Year)
      
      list(summ = summ, line_dim = "Province", ylab = ylab,
           subtitle = paste0("Commodities: ", paste(input$cmp_trend_commodities, collapse = ", ")))
    } else {
      # lines = Commodity, pooled across selected provinces
      summ <- df2 %>%
        dplyr::group_by(Commodity, Year) %>%
        dplyr::summarise(
          n = dplyr::n(),
          q25 = stats::quantile(.val, 0.25, na.rm = TRUE, names = FALSE),
          med = stats::median(.val, na.rm = TRUE),
          q75 = stats::quantile(.val, 0.75, na.rm = TRUE, names = FALSE),
          .groups = "drop"
        ) %>%
        dplyr::arrange(Commodity, Year)
      
      list(summ = summ, line_dim = "Commodity", ylab = ylab,
           subtitle = paste0("Provinces: ", paste(input$cmp_trend_provinces, collapse = ", ")))
    }
  })
  
  output$cmp_trend_plot <- renderPlotly({
    req(input$cmp_trend_enable)
    td <- cmp_trend_data()
    summ <- td$summ
    validate(need(nrow(summ) > 0, "Trend plot: no rows after filtering (2020–2025)."))
    
    line_dim <- td$line_dim
    ylab <- td$ylab
    
    # stable colors for whichever thing is lines
    lvls <- sort(unique(as.character(summ[[line_dim]])))
    cols <- make_color_map(lvls)
    
    p <- plotly::plot_ly()
    
    for (g in lvls) {
      dd <- summ %>%
        dplyr::filter(.data[[line_dim]] == g) %>%
        dplyr::arrange(Year)
      
      if (nrow(dd) < 1) next
      
      err_plus  <- dd$q75 - dd$med
      err_minus <- dd$med - dd$q25
      
      p <- p %>%
        plotly::add_trace(
          data = dd,
          x = ~Year,
          y = ~med,
          type = "scatter",
          mode = "lines+markers",
          name = g,
          line = list(width = 2, color = unname(cols[g])),
          marker = list(size = 6, color = unname(cols[g])),
          error_y = list(
            type = "data",
            array = err_plus,
            arrayminus = err_minus,
            visible = TRUE,
            thickness = 1.2,
            width = 3
          ),
          text = ~paste0(
            line_dim, ": ", .data[[line_dim]],
            "<br>Year: ", Year,
            "<br>n: ", n,
            "<br>Median: ", round(med, 3),
            "<br>Q25: ", round(q25, 3),
            "<br>Q75: ", round(q75, 3)
          ),
          hoverinfo = "text"
        )
    }
    
    p %>%
      plotly::layout(
        title = list(text = paste0("Trend (Median ± IQR), 2020–2025<br><sup>", td$subtitle, "</sup>")),
        xaxis = list(title = "Year", tickmode = "linear"),
        yaxis = list(title = ylab),
        legend = list(orientation = "h", x = 0, y = -0.25, xanchor = "left"),
        margin = list(b = 120)
      )
  })
  
  output$cmp_trend_status <- renderPrint({
    req(input$cmp_trend_enable)
    
    cat("Trend window: 2020–2025\n")
    cat("Mode:", ifelse(identical(input$cmp_trend_mode, "province_lines"),
                        "Lines = Provinces (filter by commodities)",
                        "Lines = Commodities (filter by provinces)"), "\n")
    
    cat("Trend scale:", ifelse(isTRUE(input$cmp_trend_log_value), "log10(Value ($))", "Value ($)"), "\n")
    
    if (!is.null(input$cmp_trend_commodities) && length(input$cmp_trend_commodities) > 0) {
      cat("Selected commodities:", length(input$cmp_trend_commodities), "\n")
    } else {
      cat("Selected commodities: (none / all)\n")
    }
    
    if (!is.null(input$cmp_trend_provinces) && length(input$cmp_trend_provinces) > 0) {
      cat("Selected provinces:", length(input$cmp_trend_provinces), "\n")
    } else {
      cat("Selected provinces: (none / all)\n")
    }
    
    cat("Click 'Update trend plot' after changing selections.\n")
  })
  
  
  
}

shinyApp(ui, server)

