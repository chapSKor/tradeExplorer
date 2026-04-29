
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
  df <- readr::read_csv(path, skip = 0, show_col_types = FALSE, trim_ws = TRUE)
  
  # Back-compat: old files had a junk first row like "Imports"
  if (!("Period" %in% names(df))) {
    df <- readr::read_csv(path, skip = 1, show_col_types = FALSE, trim_ws = TRUE)
  }
  
  # Drop fully blank Period rows
  df <- df %>%
    dplyr::filter(!(is.na(Period) | trimws(as.character(Period)) == ""))
  
  # Robust Period parse: try YYYY-mm-dd first, then m/d/YYYY
  p <- as.character(df$Period)
  d1 <- as.Date(p, format = "%Y-%m-%d")
  if (all(is.na(d1))) {
    d1 <- as.Date(p, format = "%m/%d/%Y")
  }
  df$Period <- d1
  
  # If still failing, stop with a useful message (prevents silent grey-out)
  if (all(is.na(df$Period))) {
    stop("Period parsing failed. Example Period values: ",
         paste(head(p, 5), collapse = ", "),
         ". Expected YYYY-mm-dd or m/d/YYYY (e.g., 2/1/1992).")
  }
  
  df <- df %>%
    mutate(
      Year = as.integer(format(Period, "%Y")),
      State = na_if(State, "N/A"),
      `unit_of_measure` = na_if(`unit_of_measure`, "N/A")
    )
  
  df
}

load_trade_data_df <- function(df) {
  # Drop fully blank Period rows
  df <- df %>%
    dplyr::filter(!(is.na(Period) | trimws(as.character(Period)) == ""))
  
  # Robust Period parse: try YYYY-mm-dd first, then m/d/YYYY
  p <- as.character(df$Period)
  d1 <- as.Date(p, format = "%Y-%m-%d")
  if (all(is.na(d1))) {
    d1 <- as.Date(p, format = "%m/%d/%Y")
  }
  df$Period <- d1
  
  if (all(is.na(df$Period))) {
    stop("Period parsing failed. Example Period values: ",
         paste(head(p, 5), collapse = ", "),
         ". Expected YYYY-mm-dd or m/d/YYYY (e.g., 2/1/1992).")
  }
  
  df <- df %>%
    dplyr::mutate(
      Year = as.integer(format(Period, "%Y")),
      State = dplyr::na_if(State, "N/A"),
      `unit_of_measure` = dplyr::na_if(`unit_of_measure`, "N/A")
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
make_color_map <- function(levels_vec, seed = 123) {
  levels_vec <- sort(unique(as.character(levels_vec)))
  n <- length(levels_vec)
  
  set.seed(seed)
  
  # Randomize order to avoid similar adjacent colors
  shuffled_levels <- sample(levels_vec, size = n, replace = FALSE)
  
  # Spectral palette (much better than default or jet)
  pal <- grDevices::colorRampPalette(
    rev(RColorBrewer::brewer.pal(11, "Spectral"))
  )(n)
  
  # Assign colors
  col_map <- setNames(pal, shuffled_levels)
  
  # Return in original order
  col_map[levels_vec]
}

# ---- ADD THIS RIGHT BELOW ----
make_topn_color_map <- function(levels_vec, seed = 123) {
  levels_vec <- sort(unique(as.character(levels_vec)))
  n <- length(levels_vec)
  
  set.seed(seed)
  
  qual_col_pals <- RColorBrewer::brewer.pal.info[
    RColorBrewer::brewer.pal.info$category == "qual", ,
    drop = FALSE
  ]
  
  col_vector <- unlist(
    mapply(
      RColorBrewer::brewer.pal,
      qual_col_pals$maxcolors,
      rownames(qual_col_pals)
    )
  )
  
  col_vector <- unique(col_vector)
  
  if (n > length(col_vector)) {
    stop("Not enough distinct qualitative colors available for this number of categories.")
  }
  
  # Shuffle available colors, then assign first n to the levels
  chosen_cols <- sample(col_vector, n, replace = FALSE)
  
  setNames(chosen_cols, levels_vec)
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


app_dir <- function() {
  d <- shiny::getShinyOption("appDir")
  if (!is.null(d)) return(d)
  getwd()
}


######## ~~~~~ UI ~~~~~ ########
tradeDashboardUI <- function(id, title = "Trade Dashboard") {
  ns <- shiny::NS(id)
  
  fluidPage(
    h2(title),
    
    sidebarLayout(
      sidebarPanel(
        
        div(class = "sidebar-section-1",
            h3("Qualitative Analyses"),
            
            checkboxInput(ns("hist_enable"), "Enable histogram panel", value = TRUE),
            
            uiOutput(ns("hist_province_ui")),
            uiOutput(ns("hist_commodity_ui")),
            uiOutput(ns("hist_year_ui")),
            
            hr(),
            h4("Histogram heatmap grouping"),
            
            radioButtons(
              ns("hist_group_by"),
              "Group histogram traces by",
              choices = c("Province", "Country", "State (US only)"),
              selected = "Province",
              inline = TRUE
            ),
            
            uiOutput(ns("hist_origin_country_ui")),
            uiOutput(ns("hist_origin_state_ui")),
            
            checkboxInput(ns("hist_drop_zero_qty"), "Drop rows with Quantity = 0 (histogram)", value = TRUE),
            
            sliderInput(ns("hist_bins"), "Bins", min = 10, max = 120, value = 40, step = 5),
            
            checkboxInput(ns("hist_log_value"), "Log10(Value) for histogram", value = TRUE),
            checkboxInput(ns("hist_density"), "Plot density instead of counts", value = FALSE),
            
            sliderInput(ns("hist_alpha"), "Transparency (alpha)", min = 0.05, max = 1, value = 0.35, step = 0.05),
            actionButton(ns("hist_run"), "Run histogram", class = "btn-primary")
        ),
        
        div(class = "sidebar-section-2",
            h3("Comparisons"),
            
            selectInput(
              ns("cmp_mode"), "Comparison mode",
              choices = c("Province vs Province", "Year vs Year (within Province)")
            ),
            
            uiOutput(ns("cmp_prov_ui")),
            uiOutput(ns("cmp_year_ui")),
            uiOutput(ns("cmp_commodity_ui")),
            
            checkboxInput(ns("cmp_use_hist_filters"), "Use histogram commodity/year filters", value = FALSE),
            
            hr(),
            h4("Origin filter (optional)"),
            
            checkboxInput(ns("cmp_origin_enable"), "Filter by origin country/state", value = FALSE),
            uiOutput(ns("cmp_origin_country_ui")),
            uiOutput(ns("cmp_origin_state_ui")),
            
            actionButton(ns("cmp_run"), "Run comparison")
        ),
        
        div(class = "sidebar-section-1",
            h3("Custom trend plot"),
            checkboxInput(ns("cmp_trend_enable"), "Show trend plot", value = TRUE),
            uiOutput(ns("cmp_trend_year_ui")),
            
            radioButtons(
              ns("cmp_trend_mode"),
              "Lines represent",
              choices = c(
                "Provinces (select commodities)" = "province_lines",
                "Commodities (select provinces)" = "commodity_lines"
              ),
              selected = "province_lines"
            ),
            
            uiOutput(ns("cmp_trend_prov_ui")),
            uiOutput(ns("cmp_trend_comm_ui")),
            checkboxInput(ns("cmp_trend_drop_zero_qty"), "Drop rows with Quantity = 0 (trend)", value = TRUE),
            checkboxInput(ns("cmp_trend_log_value"), "Log10(Value) for trend", value = TRUE),
            actionButton(ns("cmp_trend_run"), "Update trend plot"),
            checkboxInput(ns("cmp_trend_forecast_enable"), "Overlay 5-year trend forecast", value = FALSE),
            
            checkboxGroupInput(
              ns("cmp_trend_forecast_methods"),
              "Forecast methods",
              choices = c("Linear" = "linear", "Nonlinear" = "nonlinear"),
              selected = c("linear", "nonlinear")
            ),
            
            checkboxInput(ns("cmp_trend_forecast_ci"), "Show forecast confidence intervals", value = TRUE),
            
            hr(),
            h4("Forecasting (next 5 years)"),
            actionButton(ns("forecast_run"), "Run forecast", class = "btn-success"),
            helpText("Uses all available years up to 2025. Drops Quantity=0 rows automatically."),
            
            h4("Forecast plot controls"),
            
            checkboxGroupInput(
              ns("fc_methods"),
              "Methods to show",
              choices = c(
                "ARIMA (univariate)" = "arima_uni",
                "ARIMA (pooled aggregate)" = "arima_pooled",
                "Correlation (factor model)" = "factor",
                "Random Forest (lag features)" = "rf"
              ),
              selected = c("arima_uni", "rf")
            ),
            
            checkboxInput(ns("fc_log_y"), "Plot on log10 scale (recommended)", value = TRUE),
            
            sliderInput(
              ns("fc_y_clip"),
              "Y-axis clipping (quantile) to avoid one method flattening the plot",
              min = 0.90, max = 1.00, value = 0.995, step = 0.001
            ),
            
            checkboxInput(ns("fc_show_ci"), "Show confidence interval ribbons", value = TRUE)
        ),
        
        div(class = "sidebar-section-2",
            h3("Trade Clustering Map"),
            
            div(class = "sidebar-section-2",
                h4("Clustering Filters"),
                uiOutput(ns("year_ui")),
                uiOutput(ns("province_ui")),
                uiOutput(ns("country_ui")),
                uiOutput(ns("commodity_ui"))
            ),
            
            radioButtons(ns("method"), "Cluster Method", choices = c("UMAP", "t-SNE"), inline = TRUE),
            
            sliderInput(ns("sample_n"), "Sample size (for speed)", min = 500, max = 50000, value = 8000, step = 500),
            numericInput(ns("seed"), "Random seed", value = 1, min = 1),
            
            hr(),
            h4("Feature construction"),
            checkboxInput(ns("drop_zero_qty"), "Drop rows with Quantity = 0 before embedding", value = TRUE),
            
            checkboxGroupInput(
              ns("cat_cols"), "Categorical columns",
              choices = c("Commodity", "Province", "Country", "State", "unit_of_measure"),
              selected = c("Commodity", "Province", "Country")
            ),
            
            checkboxGroupInput(
              ns("num_cols"), "Numeric columns",
              choices = c("value_dollars", "Quantity"),
              selected = c("value_dollars", "Quantity")
            ),
            
            sliderInput(ns("pca_k"), "PCA dims before embedding", min = 10, max = 200, value = 50, step = 10),
            
            hr(),
            h4("Method params"),
            
            conditionalPanel(
              condition = sprintf("input['%s'] == 't-SNE'", ns("method")),
              sliderInput(ns("perplexity"), "Perplexity", min = 5, max = 80, value = 30, step = 1),
              sliderInput(ns("theta"), "Theta (speed/accuracy)", min = 0, max = 0.9, value = 0.5, step = 0.05)
            ),
            
            conditionalPanel(
              condition = sprintf("input['%s'] == 'UMAP'", ns("method")),
              sliderInput(ns("n_neighbors"), "n_neighbors", min = 5, max = 100, value = 15, step = 1),
              sliderInput(ns("min_dist"), "min_dist", min = 0, max = 1, value = 0.1, step = 0.05),
              selectInput(ns("umap_metric"), "Metric", choices = c("cosine", "euclidean", "manhattan"), selected = "cosine")
            ),
            
            hr(),
            h4("Plot"),
            
            selectInput(
              ns("color_by"), "Color by",
              choices = c("Province", "Year", "Country", "Commodity", "State", "unit_of_measure"),
              selected = "Province"
            ),
            
            fluidRow(
              column(6, actionButton(ns("go"), "Run embedding", class = "btn-primary")),
              column(6, actionButton(ns("clear"), "Clear embedding"))
            ),
            
            hr(),
            h4("Plot overlays"),
            sliderInput(ns("pt_size"), "Point size", min = 2, max = 12, value = 5, step = 1),
            checkboxInput(ns("top_n_only"), "Only plot top 15 commodities", value = FALSE),
            checkboxInput(ns("recolor_top_n"), "Re-interpolate colors for top 15 plot", value = FALSE),
            checkboxInput(ns("show_embedding_legend"), "Show legend", value = TRUE),
            checkboxInput(ns("show_density"), "Density contours", value = FALSE),
            sliderInput(ns("density_bins"), "Density grid resolution", min = 40, max = 200, value = 80, step = 10),
            
            checkboxInput(ns("show_hulls"), "Convex hulls (by Commodity)", value = FALSE),
            checkboxInput(ns("show_centroids"), "Centroids (by Commodity)", value = FALSE),
            
            conditionalPanel(
              condition = sprintf("input['%s'] == 'Commodity'", ns("color_by")),
              selectizeInput(
                ns("focus_commodity"),
                "Zoom panel: select commodity",
                choices = NULL,
                multiple = FALSE,
                options = list(placeholder = "Pick one")
              ),
              checkboxInput(ns("show_zoom_panel"), "Show zoomed panel", value = TRUE)
            )
        )
      ),
      
      mainPanel(
        hr(),
        plotlyOutput(ns("hist_plot"), height = "500px"),
        verbatimTextOutput(ns("hist_status")),
        
        conditionalPanel(
          condition = sprintf(
            "input['%s'] == 'Country' || input['%s'] == 'State (US only)'",
            ns("hist_group_by"), ns("hist_group_by")
          ),
          uiOutput(ns("hist_summary_text")),
          tableOutput(ns("hist_summary_table"))
        ),
        
        conditionalPanel(
          condition = sprintf("input['%s'] == 'State (US only)'", ns("hist_group_by")),
          h4("U.S. State Heatmap"),
          plotlyOutput(ns("hist_state_map"), height = "550px")
        ),
        
        conditionalPanel(
          condition = sprintf(
            "input['%s'] == 'Province' && (!input['%s'] || input['%s'].length == 0)",
            ns("hist_group_by"), ns("hist_province"), ns("hist_province")
          ),
          h4("Canadian Province Heatmap"),
          plotOutput(ns("hist_province_map"), height = "550px")
        ),
        
        hr(),
        h3("Comparison Metrics"),
        tableOutput(ns("cmp_table")),
        verbatimTextOutput(ns("cmp_note")),
        
        hr(),
        h3("Province Trend (Median ± IQR over Years)"),
        plotlyOutput(ns("cmp_trend_plot"), height = "900px"),
        verbatimTextOutput(ns("cmp_trend_status")),
        
        hr(),
        h3("Forecast (next 5 years)"),
        plotlyOutput(ns("forecast_plot"), height = "520px"),
        tableOutput(ns("forecast_table")),
        verbatimTextOutput(ns("forecast_status")),
        
        div(
          style = "max-width: 950px; margin:auto;",
          plotlyOutput(ns("plt"), height = "950px", width = "950px")
        ),
        
        hr(),
        verbatimTextOutput(ns("status"))
      )
    )
  )
}

######## ~~~~~ Server ~~~~~ ########
######## ~~~~~ SERVER ~~~~~ ########
tradeDashboardServer <- function(
    id,
    parquet_file,
    sorted_vars_file,
    table_name = "trade",
    monthly = FALSE,
    data_label = "trade",
    forecast_min_year = NULL,
    forecast_last_year = NULL,
    forecast_horizon_years = 5L
) {
  shiny::moduleServer(id, function(input, output, session) {
    
    # ---- Forecasting config ----
    FORECAST_MIN_YEAR <- if (!is.null(forecast_min_year)) {
      as.integer(forecast_min_year)
    } else if (isTRUE(monthly)) {
      2000L
    } else {
      1988L
    }
    
    FORECAST_LAST_YEAR <- if (!is.null(forecast_last_year)) {
      as.integer(forecast_last_year)
    } else if (isTRUE(monthly)) {
      2026L
    } else {
      2025L
    }
    
    FORECAST_HORIZON_YEARS <- as.integer(forecast_horizon_years)
    
    con_current <- NULL
    session$onSessionEnded(function() {
      if (!is.null(con_current)) {
        try(DBI::dbDisconnect(con_current, shutdown = TRUE), silent = TRUE)
      }
    })
    
    # session$onSessionEnded(function() {
    #   con <- con_rv()
    #   if (!is.null(con)) {
    #     try(DBI::dbDisconnect(con, shutdown = TRUE), silent = TRUE)
    #   }
    # })
    
    #raw_data <- reactiveVal(NULL)
    embedding_store <- reactiveVal(NULL)
    
    choices_cache <- reactiveValues(
      years = NULL,
      year_months = NULL,
      provinces = NULL,
      countries = NULL,
      commodities = NULL,
      states = NULL,
      units = NULL
    )
    
    con_rv <- reactiveVal(NULL)     # holds DBI connection
    trade_con <- reactive({ req(con_rv()); con_rv() })
    trade_tbl_rv <- reactiveVal(NULL)  # holds dbplyr tbl() pointing to the view
    
    observeEvent(TRUE, {
      tryCatch({
        message("Startup [", data_label, "]: begin init()")
        
        pq_path <- file.path(app_dir(), "data", parquet_file)
        if (!file.exists(pq_path)) {
          showNotification(paste("Missing parquet file:", pq_path,
                                 "(run the appropriate build_parquet script locally and redeploy)"),
                           type = "error", duration = 15)
          stop("Missing parquet: ", pq_path)
        }
        
        cache_path <- file.path(app_dir(), "sortedVars", sorted_vars_file)
        if (!file.exists(cache_path)) {
          showNotification(paste("Missing cache file:", cache_path,
                                 "(run the appropriate build_sorted_vars script locally and redeploy)"),
                           type = "error", duration = 15)
          stop("Missing sortedVars cache: ", cache_path)
        }
        
        # DuckDB: use an on-disk db file in tempdir (safe on shinyapps.io)
        dbfile <- file.path(tempdir(), paste0(table_name, ".duckdb"))
        con <- DBI::dbConnect(duckdb::duckdb(), dbdir = dbfile)
        
        # Make sure we clean up
        session$onSessionEnded(function() {
          try(DBI::dbDisconnect(con, shutdown = TRUE), silent = TRUE)
        })
        
        DBI::dbExecute(con, sprintf(
          "CREATE OR REPLACE VIEW trade AS
   SELECT
     Period,
     try_cast(Period AS DATE) AS Period_date,
     year(try_cast(Period AS DATE)) AS Year,
     Province, Country, State, Commodity,
     unit_of_measure,
     value_dollars,
     Quantity
   FROM parquet_scan('%s')",
          gsub("'", "''", pq_path)
        ))
        
        # ingle, correct place to set everything
        con_current <<- con
        con_rv(con)
        trade_tbl_rv(dplyr::tbl(con, "trade"))
        
        # Load cached UI choices
        sv <- readRDS(cache_path)
        needed <- c("years","provinces","countries","commodities","states","units")
        miss <- setdiff(needed, names(sv))
        if (length(miss) > 0) stop("sorted_vars.rds missing: ", paste(miss, collapse = ", "))
        
        choices_cache$years       <- sv$years
        choices_cache$year_months <- if ("year_months" %in% names(sv)) sv$year_months else NULL
        choices_cache$provinces   <- sv$provinces
        choices_cache$countries   <- sv$countries
        choices_cache$commodities <- sv$commodities
        choices_cache$states      <- sv$states
        choices_cache$units       <- sv$units
        
        message("Startup [", data_label, "]: init() complete")
        
      }, error = function(e) {
        message("Startup [", data_label, "] ERROR: ", conditionMessage(e))
        showNotification(paste("Startup error:", conditionMessage(e)),
                         type = "error", duration = 15)
      })
    }, once = TRUE)
    
    observeEvent(input$cmp_origin_country, {
      us_vals <- c("United States", "USA", "United States of America")
      if (is.null(input$cmp_origin_country) || !(input$cmp_origin_country %in% us_vals)) {
        updateSelectizeInput(session, "cmp_origin_state", selected = character(0))
      }
    }, ignoreInit = TRUE)
    
    
    
    # --- Bonus: same idea for histogram (because state filter only makes sense for US rows) ---
    observeEvent(input$hist_origin_country, {
      us_vals <- c("United States", "USA", "United States of America")
      # hist_origin_country is multiple=TRUE
      if (is.null(input$hist_origin_country) || !any(input$hist_origin_country %in% us_vals)) {
        updateSelectizeInput(session, "hist_origin_state", selected = character(0))
      }
    }, ignoreInit = TRUE)
    
    # Filter UI components
    output$year_ui <- renderUI({
      yr <- choices_cache$years
      req(yr)
      sliderInput(session$ns("year_range"), "Year range",
                  min = min(yr), max = max(yr),
                  value = c(min(yr), max(yr)), step = 1)
    })
    
    output$province_ui <- renderUI({
      provs <- choices_cache$provinces
      req(provs)
      selectizeInput(session$ns("province"), "Province",
                     choices = provs, selected = NULL, multiple = TRUE,
                     options = list(placeholder = "All"))
    })
    
    output$country_ui <- renderUI({
      cn <- choices_cache$countries
      req(cn)
      selectizeInput(session$ns("country"), "Country",
                     choices = cn, selected = NULL, multiple = TRUE,
                     options = list(placeholder = "All"))
    })
    
    output$commodity_ui <- renderUI({
      comms <- choices_cache$commodities
      req(comms)
      selectizeInput(session$ns("commodity"), "Commodity",
                     choices = comms, selected = NULL, multiple = TRUE,
                     options = list(placeholder = "All", maxOptions = 2000))
    })
    
    
    filtered_tbl <- reactive({
      req(input$year_range)
      
      tt <- trade_tbl_rv()
      req(tt)
      
      tt <- tt %>%
        filter(Year >= input$year_range[1], Year <= input$year_range[2])
      
      if (!is.null(input$province) && length(input$province) > 0) {
        tt <- tt %>% filter(Province %in% input$province)
      }
      if (!is.null(input$country) && length(input$country) > 0) {
        tt <- tt %>% filter(Country %in% input$country)
      }
      if (!is.null(input$commodity) && length(input$commodity) > 0) {
        tt <- tt %>% filter(Commodity %in% input$commodity)
      }
      
      if (isTRUE(input$drop_zero_qty)) {
        tt <- tt %>% filter(!is.na(Quantity), Quantity != 0)
      }
      
      tt
    })
    
    # --- NEW: map full U.S. state names to USPS abbreviations for plotly choropleth ---
    us_state_lookup <- c(
      "Alabama" = "AL", "Alaska" = "AK", "Arizona" = "AZ", "Arkansas" = "AR",
      "California" = "CA", "Colorado" = "CO", "Connecticut" = "CT", "Delaware" = "DE",
      "District of Columbia" = "DC", "Florida" = "FL", "Georgia" = "GA", "Hawaii" = "HI",
      "Idaho" = "ID", "Illinois" = "IL", "Indiana" = "IN", "Iowa" = "IA",
      "Kansas" = "KS", "Kentucky" = "KY", "Louisiana" = "LA", "Maine" = "ME",
      "Maryland" = "MD", "Massachusetts" = "MA", "Michigan" = "MI", "Minnesota" = "MN",
      "Mississippi" = "MS", "Missouri" = "MO", "Montana" = "MT", "Nebraska" = "NE",
      "Nevada" = "NV", "New Hampshire" = "NH", "New Jersey" = "NJ", "New Mexico" = "NM",
      "New York" = "NY", "North Carolina" = "NC", "North Dakota" = "ND", "Ohio" = "OH",
      "Oklahoma" = "OK", "Oregon" = "OR", "Pennsylvania" = "PA", "Rhode Island" = "RI",
      "South Carolina" = "SC", "South Dakota" = "SD", "Tennessee" = "TN", "Texas" = "TX",
      "Utah" = "UT", "Vermont" = "VT", "Virginia" = "VA", "Washington" = "WA",
      "West Virginia" = "WV", "Wisconsin" = "WI", "Wyoming" = "WY"
    )
    # --- NEW: U.S. state choropleth for histogram summary ---
    output$hist_state_map <- renderPlotly({
      req(input$hist_enable)
      req(input$hist_run > 0)
      req(identical(input$hist_group_by, "State (US only)"))
      
      df <- hist_importer_summary_tbl()
      shiny::validate(shiny::need(nrow(df) > 0, "No U.S. state rows available for mapping."))
      
      # Keep only valid state names and map to abbreviations
      df <- df %>%
        dplyr::filter(!is.na(importer), nzchar(importer)) %>%
        dplyr::mutate(
          state_abbr = unname(us_state_lookup[importer]),
          log10_value = log10(pmax(import_value_dollars, 1))
        ) %>%
        dplyr::filter(!is.na(state_abbr))
      
      shiny::validate(shiny::need(nrow(df) > 0, "No mappable U.S. states found after name matching."))
      
      # White -> warm orange -> complementary blue
      pharma_colorscale <- list(
        list(0.00, "#FFFFFF"),
        list(0.35, "#F4A261"),
        list(1.00, "#3A86C8")
      )
      
      plotly::plot_ly(
        data = df,
        type = "choropleth",
        locationmode = "USA-states",
        locations = ~state_abbr,
        z = ~log10_value,
        text = ~paste0(
          "State: ", importer,
          "<br>Total imports ($): ", format(round(import_value_dollars, 2), big.mark = ","),
          "<br>log10(Value): ", round(log10_value, 3),
          "<br>Rows: ", n_rows
        ),
        hoverinfo = "text",
        colorscale = pharma_colorscale,
        marker = list(line = list(color = "white", width = 0.8)),
        colorbar = list(title = "Value (log dollars)")
      ) %>%
        plotly::layout(
          geo = list(
            scope = "usa",
            projection = list(type = "albers usa"),
            showlakes = TRUE,
            lakecolor = toRGB("white")
          ),
          margin = list(l = 20, r = 20, t = 40, b = 20),
          title = "Filtered U.S. state exports to Canada"
        )
    })
    
    # NEW: Provinces 
    
    # --- NEW: Canadian province choropleth for histogram summary ---
    output$hist_province_map <- renderPlot({
      req(input$hist_enable)
      req(input$hist_run > 0)
      req(identical(input$hist_group_by, "Province"))
      req(is.null(input$hist_province) || length(input$hist_province) == 0)
      
      df <- hist_province_summary_tbl()
      
      df <- df %>%
        dplyr::filter(!is.na(importer), nzchar(importer)) %>%
        dplyr::mutate(
          map_name = normalize_ca_province(importer),
          log10_value = log10(pmax(import_value_dollars, 1))
        )
      
      plot_df <- canada_provinces_sf %>%
        dplyr::left_join(df, by = "map_name")
      
      shiny::validate(
        shiny::need(any(is.finite(plot_df$log10_value), na.rm = TRUE),
                    "No province totals available for mapping.")
      )
      
      ggplot2::ggplot(plot_df) +
        ggplot2::geom_sf(ggplot2::aes(fill = log10_value), color = "white", linewidth = 0.6) +
        ggplot2::scale_fill_gradientn(
          colours = c("#FFFFFF", "#F4A261", "#3A86C8"),
          na.value = "grey90",
          name = "Value (Log CAD)"
        ) +
        ggplot2::labs(
          title = "Filtered Canadian imports by province",
          subtitle = "Current histogram filters: commodity, year, and all provinces"
        ) +
        ggplot2::theme_void() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 16),
          plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 11),
          legend.title = ggplot2::element_text(size = 11),
          legend.text = ggplot2::element_text(size = 10),
          plot.margin = ggplot2::margin(10, 10, 10, 10)
        )
    })
    
    #### Forceast Trends based directly off of trend plot
    
    trend_projection_tbl <- function(df, year_horizon = 5L) {
      # df must contain: Year, med, q25, q75, n
      df <- df %>%
        dplyr::arrange(Year) %>%
        dplyr::filter(is.finite(Year), is.finite(med), !is.na(med))
      
      if (nrow(df) < 5) return(NULL)
      
      future_years <- seq.int(max(df$Year) + 1L, max(df$Year) + year_horizon)
      newdata <- data.frame(Year = future_years)
      
      out_list <- list()
      
      # historical spread proxy
      hist_half_iqr <- (df$q75 - df$q25) / 2
      spread_ref <- stats::median(hist_half_iqr[is.finite(hist_half_iqr)], na.rm = TRUE)
      if (!is.finite(spread_ref)) spread_ref <- 0
      
      # ---------- Linear ----------
      fit_lm <- stats::lm(med ~ Year, data = df)
      pr_lm <- stats::predict(fit_lm, newdata = newdata, interval = "confidence", level = 0.95)
      
      out_list$linear <- dplyr::tibble(
        Year = future_years,
        method = "Linear",
        point = as.numeric(pr_lm[, "fit"]),
        lo = as.numeric(pr_lm[, "lwr"]) - spread_ref,
        hi = as.numeric(pr_lm[, "upr"]) + spread_ref
      )
      
      # ---------- Nonlinear ----------
      fit_ns <- stats::lm(med ~ splines::ns(Year, df = min(3, max(1, nrow(df) - 2))), data = df)
      pr_ns <- stats::predict(fit_ns, newdata = newdata, interval = "confidence", level = 0.95)
      
      out_list$nonlinear <- dplyr::tibble(
        Year = future_years,
        method = "Nonlinear",
        point = as.numeric(pr_ns[, "fit"]),
        lo = as.numeric(pr_ns[, "lwr"]) - spread_ref,
        hi = as.numeric(pr_ns[, "upr"]) + spread_ref
      )
      
      dplyr::bind_rows(out_list)
    }
    ########
    
    observeEvent(input$go, {
      
      req(con_rv())
      con <- con_rv()
      req(input$year_range)
      
      # Build WHERE clause from clustering filters
      where <- c(sprintf("Year BETWEEN %d AND %d",
                         as.integer(input$year_range[1]),
                         as.integer(input$year_range[2])))
      
      if (!is.null(input$province) && length(input$province) > 0) {
        where <- c(where, sprintf("Province IN (%s)",
                                  paste(DBI::dbQuoteString(con, input$province), collapse = ",")))
      }
      if (!is.null(input$country) && length(input$country) > 0) {
        where <- c(where, sprintf("Country IN (%s)",
                                  paste(DBI::dbQuoteString(con, input$country), collapse = ",")))
      }
      if (!is.null(input$commodity) && length(input$commodity) > 0) {
        where <- c(where, sprintf("Commodity IN (%s)",
                                  paste(DBI::dbQuoteString(con, input$commodity), collapse = ",")))
      }
      if (isTRUE(input$drop_zero_qty)) {
        where <- c(where, "Quantity IS NOT NULL AND Quantity != 0")
      }
      
      where_sql <- paste(where, collapse = " AND ")
      
      #max_pull <- 200000L
      max_pull <- max(200000L, as.integer(input$sample_n) * 5L)
      max_pull <- min(max_pull, 300000L)
      seed <- if (is.null(input$seed) || is.na(input$seed)) 1L else as.integer(input$seed)
      
      sql <- sprintf("
  SELECT
    Province, Country, State, Commodity, unit_of_measure,
    Year, value_dollars, Quantity
  FROM (
    SELECT
      Province, Country, State, Commodity, unit_of_measure,
      Year, value_dollars, Quantity
    FROM trade
    WHERE %s
  ) t
  USING SAMPLE RESERVOIR(%d ROWS)
  REPEATABLE (%d)
", where_sql, max_pull, seed)
      
      message("EMBED SQL:\n", sql)  # <-- remove after confirmed
      
      tryCatch({
        withProgress(message = "Computing embedding…", value = 0, {
          
          incProgress(0.05, detail = "Querying DuckDB sample")
          df <- DBI::dbGetQuery(con, sql)
          
          shiny::validate(shiny::need(nrow(df) >= 200,
                                      "Too few rows after filters for embedding."))
          
          incProgress(0.20, detail = "Sampling rows in R")
          set.seed(seed)
          n_take <- min(as.integer(input$sample_n), nrow(df))
          idx <- sample.int(nrow(df), n_take)
          dsub <- df[idx, , drop = FALSE]
          
          req(length(input$num_cols) > 0)
          req(length(input$cat_cols) > 0)
          
          incProgress(0.45, detail = "Building features (one-hot + PCA)")
          X <- make_feature_matrix(
            dsub,
            numeric_cols = input$num_cols,
            categorical_cols = input$cat_cols,
            pca_k = input$pca_k,
            seed = seed
          )
          
          incProgress(0.80, detail = paste("Running", input$method))
          coords <- if (identical(input$method, "t-SNE")) {
            max_perp <- max(5, floor((nrow(X) - 1) / 3))
            perp <- min(as.integer(input$perplexity), max_perp)
            run_tsne(X, perplexity = perp, theta = input$theta, seed = seed)
          } else {
            run_umap(X, n_neighbors = input$n_neighbors, min_dist = input$min_dist,
                     metric = input$umap_metric, seed = seed)
          }
          
          rm(X); gc()
          
          incProgress(0.95, detail = "Finalizing")
          coords <- as.data.frame(coords)
          colnames(coords) <- c("Dim1", "Dim2")
          
          out <- dsub %>%
            dplyr::mutate(Dim1 = coords$Dim1, Dim2 = coords$Dim2)
          
          embedding_store(out)
          showNotification("Embedding complete.", type = "message", duration = 2)
        })
      }, error = function(e) {
        showNotification(paste("Run embedding failed:", conditionMessage(e)),
                         type = "error", duration = 12)
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
      req(con_rv())
      con <- con_rv()
      
      cat("Shiny getwd():", getwd(), "\n")
      cat("Data folder:", file.path(app_dir(), "data"), "\n")
      cat("Files in data/:", paste(list.files(file.path(app_dir(), "data")), collapse=", "), "\n\n")
      
      # total rows in trade
      n_total <- DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM trade")$n
      cat("Total rows (", data_label, "):", n_total, "\n", sep = "")
      
      # filtered rows estimate (after current clustering filters)
      req(input$year_range)
      where <- c(sprintf("Year BETWEEN %d AND %d",
                         as.integer(input$year_range[1]),
                         as.integer(input$year_range[2])))
      
      if (!is.null(input$province) && length(input$province) > 0) {
        where <- c(where, sprintf("Province IN (%s)",
                                  paste(DBI::dbQuoteString(con, input$province), collapse=",")))
      }
      if (!is.null(input$country) && length(input$country) > 0) {
        where <- c(where, sprintf("Country IN (%s)",
                                  paste(DBI::dbQuoteString(con, input$country), collapse=",")))
      }
      if (!is.null(input$commodity) && length(input$commodity) > 0) {
        where <- c(where, sprintf("Commodity IN (%s)",
                                  paste(DBI::dbQuoteString(con, input$commodity), collapse=",")))
      }
      if (isTRUE(input$drop_zero_qty)) {
        where <- c(where, "Quantity IS NOT NULL AND Quantity != 0")
      }
      
      where_sql <- paste(where, collapse=" AND ")
      n_filt <- DBI::dbGetQuery(con, sprintf("SELECT COUNT(*) AS n FROM trade WHERE %s", where_sql))$n
      cat("Filtered rows (clustering filters):", n_filt, "\n")
      
      if (!is.null(embedding_store())) {
        cat("Embedding stored rows:", nrow(embedding_store()), "\n")
      } else {
        cat("Embedding stored: none (click 'Run embedding')\n")
      }
    })
    
    # ---- Plot (with stable colors + overlays + zoom) ----
    output$plt <- renderPlotly({
      res <- embedding_store()
      shiny::validate(shiny::need(!is.null(res), "Choose a CSV, adjust filters/features, then click 'Run embedding'."))
      
      color_col <- input$color_by
      shiny::validate(shiny::need(color_col %in% names(res), "Chosen color column not available in plotted data."))
      
      # Optional: restrict to top 15 categories in the currently displayed subset
      if (isTRUE(input$top_n_only) && identical(color_col, "Commodity")) {
        keep_levels <- res %>%
          dplyr::filter(!is.na(Commodity)) %>%
          dplyr::count(Commodity, name = "n", sort = TRUE) %>%
          dplyr::slice_head(n = 15) %>%
          dplyr::pull(Commodity) %>%
          as.character()
        
        res <- res %>%
          dplyr::filter(as.character(Commodity) %in% keep_levels)
      }
      
      if (identical(color_col, "Commodity")) {
        if (isTRUE(input$top_n_only) && isTRUE(input$recolor_top_n)) {
          message("TOP-15 RECOLOR BRANCH ACTIVE")
          visible_comms <- sort(unique(as.character(res$Commodity[!is.na(res$Commodity)])))
          col_map <- make_topn_color_map(visible_comms, seed = 123)
          levels_all <- visible_comms
          res[[color_col]] <- factor(as.character(res[[color_col]]), levels = visible_comms)
        } else {
          col_map <- commodity_color_map()
          levels_all <- sort(unique(as.character(res[[color_col]])))
          res[[color_col]] <- factor(as.character(res[[color_col]]), levels = names(col_map))
        }
      } else {
        levels_all <- sort(unique(as.character(res[[color_col]])))
        res[[color_col]] <- factor(as.character(res[[color_col]]), levels = levels_all)
        col_map <- make_color_map(levels_all, seed = 123)
      }
      
      # Build hover text (full info)
      hover_text <- apply(res, 1, function(row) {
        paste(
          if ("Province" %in% names(res)) paste0("Province: ", row["Province"], "<br>") else "",
          if ("Year" %in% names(res)) paste0("Year: ", row["Year"], "<br>") else "",
          if ("Country" %in% names(res)) paste0("Country: ", row["Country"], "<br>") else "",
          if ("State" %in% names(res)) paste0("State: ", row["State"], "<br>") else "",
          if ("unit_of_measure" %in% names(res)) paste0("Unit: ", row["unit_of_measure"], "<br>") else "",
          if ("Commodity" %in% names(res)) paste0("Commodity: ", row["Commodity"], "<br>") else "",
          if ("value_dollars" %in% names(res)) paste0("value_dollars: ", row["value_dollars"], "<br>") else "",
          if ("Quantity" %in% names(res)) paste0("Quantity: ", row["Quantity"]) else ""
        )
      })
      
      # Stable factor levels + stable color map across all traces (main + zoom + centroids)
      # levels_all <- sort(unique(res[[color_col]]))
      # res[[color_col]] <- factor(res[[color_col]], levels = levels_all)
      # col_map <- make_color_map(levels_all)
      
      if (identical(color_col, "Commodity")) {
        col_map <- commodity_color_map()
        levels_all <- intersect(names(col_map), sort(unique(as.character(res[[color_col]]))))
        res[[color_col]] <- factor(as.character(res[[color_col]]), levels = names(col_map))
      } else {
        levels_all <- sort(unique(as.character(res[[color_col]])))
        res[[color_col]] <- factor(as.character(res[[color_col]]), levels = levels_all)
        col_map <- make_color_map(levels_all)
      }
      
      # # Per-point colors (fast + stable)
      # col_vec <- unname(col_map[as.character(res[[color_col]])])
      # 
      # # MAIN scatter (single WebGL trace, with explicit per-point colors)
      # p_main <- plotly::plot_ly(
      #   data = res,
      #   x = ~Dim1, y = ~Dim2,
      #   type = "scattergl", mode = "markers",
      #   marker = list(size = input$pt_size, opacity = 0.85, color = col_vec),
      #   text = hover_text,
      #   hoverinfo = "text",
      #   showlegend = TRUE
      # )
      # 
      # # Add legend-only traces so the legend shows correct category ↔ color mapping
      # # legend_df <- data.frame(lbl = levels_all, x = 0, y = 0, stringsAsFactors = FALSE)
      # # legend_cols <- unname(col_map[legend_df$lbl])
      # 
      # legend_levels <- if (identical(color_col, "Commodity")) {
      #   sort(unique(as.character(res$Commodity[!is.na(res$Commodity)])))
      # } else {
      #   levels_all
      # }
      # 
      # legend_df <- data.frame(lbl = legend_levels, x = 0, y = 0, stringsAsFactors = FALSE)
      # legend_cols <- unname(col_map[legend_df$lbl])
      # 
      # legend_spec <- if (isTRUE(input$show_embedding_legend)) {
      #   list(
      #     orientation = "h",
      #     x = 0, y = -0.25,
      #     xanchor = "left",
      #     yanchor = "top",
      #     font = list(size = 10)
      #   )
      # } else {
      #   list()
      # }
      # 
      # bottom_margin <- if (isTRUE(input$show_embedding_legend)) 260 else 60
      # 
      # p_main <- p_main %>%
      #   add_trace(
      #     data = legend_df,
      #     x = ~x, y = ~y,
      #     type = "scatter",
      #     mode = "markers",
      #     marker = list(size = 9, color = legend_cols),
      #     name = ~lbl,
      #     hoverinfo = "skip",
      #     visible = "legendonly",
      #     inherit = FALSE,
      #     showlegend = isTRUE(input$show_embedding_legend)
      #   ) %>%
      #   layout(
      #     title = paste(input$method, "embedding"),
      #     width = 900,
      #     height = 900,
      #     xaxis = if (isTRUE(input$show_tsne_axes)) {
      #       list(
      #         title = list(text = "Dim1", font = list(size = 18)),
      #         tickfont = list(size = 16),
      #         zeroline = FALSE
      #       )
      #     } else {
      #       list(
      #         title = "",
      #         showticklabels = FALSE,
      #         showgrid = FALSE,
      #         zeroline = FALSE
      #       )
      #     },
      #     yaxis = if (isTRUE(input$show_tsne_axes)) {
      #       list(
      #         title = list(text = "Dim2", font = list(size = 18)),
      #         tickfont = list(size = 16),
      #         zeroline = FALSE,
      #         scaleanchor = "x"
      #       )
      #     } else {
      #       list(
      #         title = "",
      #         showticklabels = FALSE,
      #         showgrid = FALSE,
      #         zeroline = FALSE,
      #         scaleanchor = "x"
      #       )
      #     },
      #     showlegend = isTRUE(input$show_embedding_legend),
      #     legend = legend_spec,
      #     margin = list(b = bottom_margin)
      #   )
      # 
      legend_spec <- if (isTRUE(input$show_embedding_legend)) {
        list(
          orientation = "h",
          x = 0, y = -0.25,
          xanchor = "left",
          yanchor = "top",
          font = list(size = 10)
        )
      } else {
        list()
      }
      
      bottom_margin <- if (isTRUE(input$show_embedding_legend)) 260 else 60
      
      legend_levels <- res %>%
        dplyr::filter(!is.na(.data[[color_col]])) %>%
        dplyr::count(.data[[color_col]], name = "n", sort = TRUE) %>%
        dplyr::pull(.data[[color_col]]) %>%
        as.character()
      
      p_main <- plotly::plot_ly()
      
      for (g in legend_levels) {
        dd <- res[as.character(res[[color_col]]) == g, , drop = FALSE]
        if (nrow(dd) == 0) next
        
        dd_hover <- hover_text[as.character(res[[color_col]]) == g]
        
        p_main <- p_main %>%
          add_trace(
            data = dd,
            x = ~Dim1, y = ~Dim2,
            type = "scattergl",
            mode = "markers",
            name = g,
            marker = list(
              size = input$pt_size,
              opacity = 0.85,
              color = unname(col_map[g])
            ),
            text = dd_hover,
            hoverinfo = "text",
            showlegend = isTRUE(input$show_embedding_legend)
          )
      }
      
      
      
      p_main <- p_main %>%
        layout(
          title = paste(input$method, "embedding"),
          width = 900,
          height = 900,
          xaxis = if (isTRUE(input$show_tsne_axes)) {
            list(
              title = list(text = "Dim1", font = list(size = 18)),
              tickfont = list(size = 16),
              zeroline = FALSE
            )
          } else {
            list(
              title = "",
              showticklabels = FALSE,
              showgrid = FALSE,
              zeroline = FALSE
            )
          },
          yaxis = if (isTRUE(input$show_tsne_axes)) {
            list(
              title = list(text = "Dim2", font = list(size = 18)),
              tickfont = list(size = 16),
              zeroline = FALSE,
              scaleanchor = "x"
            )
          } else {
            list(
              title = "",
              showticklabels = FALSE,
              showgrid = FALSE,
              zeroline = FALSE,
              scaleanchor = "x"
            )
          },
          showlegend = isTRUE(input$show_embedding_legend),
          legend = legend_spec,
          margin = list(b = bottom_margin)
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
              if ("value_dollars" %in% names(zdf)) paste0("value_dollars: ", row["value_dollars"], "<br>") else "",
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
      provs <- choices_cache$provinces
      req(provs)
      selectizeInput(session$ns("hist_province"), "Province(s)",
                     choices = provs, selected = NULL, multiple = TRUE,
                     options = list(placeholder = "All provinces"))
    })
    
    output$hist_commodity_ui <- renderUI({
      comms <- choices_cache$commodities
      req(comms)
      selectizeInput(session$ns("hist_commodity"), "Commodity(ies)",
                     choices = comms, selected = NULL, multiple = TRUE,
                     options = list(placeholder = "All commodities", maxOptions = 2000))
    })
    
    output$hist_year_ui <- renderUI({
      yrs <- choices_cache$years
      req(yrs)
      selectizeInput(session$ns("hist_year"), "Year(s)",
                     choices = yrs, selected = NULL, multiple = TRUE,
                     options = list(placeholder = "All years"))
    })
    
    commodity_color_map <- reactive({
      req(choices_cache$commodities)
      comms <- sort(unique(choices_cache$commodities))
      make_color_map(comms, seed = 123)
    })
    
    output$hist_origin_country_ui <- renderUI({
      cn <- choices_cache$countries
      req(cn)
      selectizeInput(session$ns("hist_origin_country"), "Origin country filter (optional)",
                     choices = cn, selected = NULL, multiple = TRUE,
                     options = list(placeholder = "All origin countries"))
    })
    
    output$hist_origin_state_ui <- renderUI({
      st <- choices_cache$states
      req(st)
      selectizeInput(session$ns("hist_origin_state"), "Origin state filter (US only, optional)",
                     choices = st, selected = NULL, multiple = TRUE,
                     options = list(placeholder = "All US states"))
    })
    
    
    ## NEW CODE ##
    
    # --- NEW: helper to build histogram WHERE clause (same logic as hist_bins_tbl) ---
    build_hist_where_sql <- function(con, input) {
      where <- c("value_dollars IS NOT NULL")
      
      # Year selection (exact hist_year OR fallback to year_range)
      if (!is.null(input$hist_year) && length(input$hist_year) > 0) {
        yrs <- as.integer(input$hist_year)
        where <- c(where, sprintf("Year IN (%s)", paste(yrs, collapse = ",")))
      } else if (!is.null(input$year_range)) {
        where <- c(where, sprintf("Year BETWEEN %d AND %d",
                                  as.integer(input$year_range[1]),
                                  as.integer(input$year_range[2])))
      }
      
      if (!is.null(input$hist_province) && length(input$hist_province) > 0) {
        where <- c(where, sprintf("Province IN (%s)",
                                  paste(DBI::dbQuoteString(con, input$hist_province), collapse = ",")))
      }
      
      if (!is.null(input$hist_commodity) && length(input$hist_commodity) > 0) {
        where <- c(where, sprintf("Commodity IN (%s)",
                                  paste(DBI::dbQuoteString(con, input$hist_commodity), collapse = ",")))
      }
      
      if (!is.null(input$hist_origin_country) && length(input$hist_origin_country) > 0) {
        where <- c(where, sprintf("Country IN (%s)",
                                  paste(DBI::dbQuoteString(con, input$hist_origin_country), collapse = ",")))
      }
      
      if (!is.null(input$hist_origin_state) && length(input$hist_origin_state) > 0) {
        where <- c(where, sprintf("State IN (%s)",
                                  paste(DBI::dbQuoteString(con, input$hist_origin_state), collapse = ",")))
      }
      
      if (isTRUE(input$hist_drop_zero_qty)) {
        where <- c(where, "Quantity IS NOT NULL AND Quantity != 0")
      }
      
      # If grouping by state, restrict to US
      if (identical(input$hist_group_by, "State (US only)")) {
        where <- c(where, "Country IN ('United States','USA','United States of America')")
      }
      
      paste(where, collapse = " AND ")
    }
    
    # NEW: Provinces 
    # --- NEW: normalize province names for joining to map geometry ---
    normalize_ca_province <- function(x) {
      x <- trimws(as.character(x))
      dplyr::case_when(
        x %in% c("Newfoundland", "Newfoundland and Labrador") ~ "Newfoundland and Labrador",
        x %in% c("Prince Edward Island", "PEI", "P.E.I.") ~ "Prince Edward Island",
        x %in% c("Nova Scotia") ~ "Nova Scotia",
        x %in% c("New Brunswick") ~ "New Brunswick",
        x %in% c("Quebec", "Québec") ~ "Quebec",
        x %in% c("Ontario") ~ "Ontario",
        x %in% c("Manitoba") ~ "Manitoba",
        x %in% c("Saskatchewan") ~ "Saskatchewan",
        x %in% c("Alberta") ~ "Alberta",
        x %in% c("British Columbia") ~ "British Columbia",
        x %in% c("Yukon", "Yukon Territory") ~ "Yukon",
        x %in% c("Northwest Territories", "NWT") ~ "Northwest Territories",
        x %in% c("Nunavut") ~ "Nunavut",
        TRUE ~ x
      )
    }
    
    # --- NEW: Canada provinces sf object (built once) ---
    canada_provinces_sf <- local({
      shp <- rnaturalearth::ne_states(country = "canada", returnclass = "sf")
      
      shp %>%
        dplyr::mutate(map_name = normalize_ca_province(name_en)) %>%
        dplyr::group_by(map_name) %>%
        dplyr::summarise(geometry = sf::st_union(geometry), .groups = "drop") %>%
        sf::st_make_valid()
    })
    
    # --- NEW: province summary table for Canada heatmap ---
    hist_province_summary_tbl <- eventReactive(input$hist_run, {
      req(input$hist_enable)
      req(con_rv())
      con <- con_rv()
      
      req(identical(input$hist_group_by, "Province"))
      req(is.null(input$hist_province) || length(input$hist_province) == 0)
      
      where_sql <- build_hist_where_sql(con, input)
      
      sql <- sprintf("
    SELECT
      Province AS importer,
      COUNT(*) AS n_rows,
      SUM(CAST(value_dollars AS DOUBLE)) AS import_value_dollars,
      SUM(CAST(Quantity AS DOUBLE)) AS import_quantity
    FROM trade
    WHERE %s
    GROUP BY 1
    ORDER BY import_value_dollars DESC
  ", where_sql)
      
      df <- DBI::dbGetQuery(con, sql)
      shiny::validate(shiny::need(nrow(df) > 0, "No province rows available for mapping."))
      df
    })
    
    # --- NEW: summary table for border-policy collaborators (Country or US State) ---
    hist_importer_summary_tbl <- eventReactive(input$hist_run, {
      req(input$hist_enable)
      req(con_rv())
      con <- con_rv()
      
      req(identical(input$hist_group_by, "Country") || identical(input$hist_group_by, "State (US only)"))
      
      where_sql <- build_hist_where_sql(con, input)
      group_col <- if (identical(input$hist_group_by, "Country")) "Country" else "State"
      
      sql <- sprintf("
    SELECT
      %s AS importer,
      COUNT(*) AS n_rows,
      SUM(CAST(value_dollars AS DOUBLE)) AS import_value_dollars,
      SUM(CAST(Quantity AS DOUBLE)) AS import_quantity
    FROM trade
    WHERE %s
    GROUP BY 1
    ORDER BY import_value_dollars DESC
    LIMIT 50
  ", group_col, where_sql)
      
      df <- DBI::dbGetQuery(con, sql)
      shiny::validate(shiny::need(nrow(df) > 0, "No rows after filtering for the summary table."))
      df
    })
    
    ####### Forcasting Helpers ######
    
    # ---- Forecasting helpers (additive; does not touch existing features) ----
    
    # Safer log transform for annual totals; avoids eps decisions
    f_log1p <- function(x) log1p(pmax(x, 0))
    f_expm1 <- function(x) pmax(expm1(x), 0)
    
    # Convert a vector to a "contiguous enough" annual ts by filling missing years with NA
    make_year_grid <- function(df, year_min, year_max) {
      grid <- data.frame(Year = seq.int(year_min, year_max), stringsAsFactors = FALSE)
      dplyr::left_join(grid, df, by = "Year")
    }
    
    # ARIMA fallback that does not require the forecast package
    # Returns list(mean, lo, hi) on the modeled scale
    arima_forecast_base <- function(y, h, level = 0.95) {
      # y is numeric with possible NAs
      yy <- y
      ok <- is.finite(yy) & !is.na(yy)
      if (sum(ok) < 6) stop("Too few non-missing years for ARIMA (need at least ~6).")
      
      # Fit on non-missing segment only (but keep order)
      # We'll fit to the observed sequence after dropping NAs
      y_obs <- yy[ok]
      
      # Simple auto-ish: try a couple candidates; pick by AIC
      cand <- list(
        c(0,1,1),
        c(1,1,1),
        c(1,1,0),
        c(0,1,0),
        c(2,1,2)
      )
      fits <- lapply(cand, function(ord) {
        try(stats::arima(y_obs, order = ord), silent = TRUE)
      })
      aics <- vapply(fits, function(f) if (inherits(f, "try-error")) Inf else AIC(f), numeric(1))
      best <- fits[[which.min(aics)]]
      if (inherits(best, "try-error")) stop("ARIMA fitting failed for this series.")
      
      pr <- stats::predict(best, n.ahead = h)
      mu <- as.numeric(pr$pred)
      se <- as.numeric(pr$se)
      
      z <- stats::qnorm((1 + level) / 2)
      lo <- mu - z * se
      hi <- mu + z * se
      list(mean = mu, lo = lo, hi = hi)
    }
    
    # Factor-based "correlation method" (Method 2B):
    # - build panel matrix (years x series), log1p
    # - standardize per series
    # - SVD -> k factors
    # - forecast each factor with ARIMA
    # - reconstruct each series' forecasts (on standardized scale), then unstandardize
    factor_forecast <- function(Ymat, h, k = 2, level = 0.95) {
      # Ymat: (T x M), may contain NA
      Tn <- nrow(Ymat); M <- ncol(Ymat)
      if (M < 2) stop("Need at least 2 series for correlation/factor forecast.")
      
      # log1p + standardize per series with NA-safe mean/sd
      Z <- apply(Ymat, 2, function(v) f_log1p(v))
      mu <- apply(Z, 2, function(v) mean(v, na.rm = TRUE))
      sdv <- apply(Z, 2, function(v) stats::sd(v, na.rm = TRUE))
      sdv[!is.finite(sdv) | sdv == 0] <- 1
      
      Zs <- sweep(Z, 2, mu, "-")
      Zs <- sweep(Zs, 2, sdv, "/")
      
      # Fill NAs with 0 for factor extraction (on standardized scale)
      Zs_fill <- Zs
      Zs_fill[!is.finite(Zs_fill)] <- 0
      
      k <- max(1, min(k, min(Tn - 1, M)))
      sv <- base::svd(Zs_fill, nu = k, nv = k)
      F <- sv$u[, 1:k, drop = FALSE] %*% diag(sv$d[1:k], nrow = k)  # factors over time
      L <- sv$v[, 1:k, drop = FALSE]                                 # loadings (M x k)
      
      # Forecast each factor with ARIMA (base)
      Fc_mean <- matrix(NA_real_, nrow = h, ncol = k)
      Fc_lo   <- matrix(NA_real_, nrow = h, ncol = k)
      Fc_hi   <- matrix(NA_real_, nrow = h, ncol = k)
      
      for (j in seq_len(k)) {
        fr <- arima_forecast_base(F[, j], h = h, level = level)
        Fc_mean[, j] <- fr$mean
        Fc_lo[, j]   <- fr$lo
        Fc_hi[, j]   <- fr$hi
      }
      
      # Reconstruct series forecasts on standardized scale:
      # Zhat = Fhat %*% t(L)
      Zhat_mean <- Fc_mean %*% t(L)  # (h x M)
      
      # For intervals: conservative approach by propagating factor bands (not exact)
      Zhat_lo <- Fc_lo %*% t(L)
      Zhat_hi <- Fc_hi %*% t(L)
      
      # Unstandardize and inverse-transform back to dollars
      # z = zhat * sd + mu, y = expm1(z)
      unstd <- function(Zhat) {
        Z_un <- sweep(Zhat, 2, sdv, "*")
        Z_un <- sweep(Z_un, 2, mu, "+")
        apply(Z_un, 2, f_expm1)
      }
      
      list(
        mean = unstd(Zhat_mean),
        lo   = unstd(Zhat_lo),
        hi   = unstd(Zhat_hi)
      )
    }
    
    # --- NEW: paragraph describing selections ---
    output$hist_summary_text <- renderUI({
      req(input$hist_enable)
      req(input$hist_run > 0)
      req(identical(input$hist_group_by, "Country") || identical(input$hist_group_by, "State (US only)"))
      
      # Years
      years_txt <- if (!is.null(input$hist_year) && length(input$hist_year) > 0) {
        paste(sort(as.integer(input$hist_year)), collapse = ", ")
      } else if (!is.null(input$year_range)) {
        paste0(as.integer(input$year_range[1]), "–", as.integer(input$year_range[2]))
      } else {
        "All years"
      }
      
      # Provinces
      prov_txt <- if (!is.null(input$hist_province) && length(input$hist_province) > 0) {
        paste(input$hist_province, collapse = ", ")
      } else {
        "All provinces"
      }
      
      # Commodities
      comm_txt <- if (!is.null(input$hist_commodity) && length(input$hist_commodity) > 0) {
        paste(input$hist_commodity, collapse = ", ")
      } else {
        "All commodities"
      }
      
      # Mode label
      mode_txt <- if (identical(input$hist_group_by, "Country")) "Country" else "US State"
      
      tagList(
        tags$p(
          tags$b("Summary (", mode_txt, "): "),
          "Years: ", years_txt, " | ",
          "Provinces: ", prov_txt, " | ",
          "Commodities: ", comm_txt
        )
      )
    })
    
    # --- NEW: table output ---
    output$hist_summary_table <- renderTable({
      req(input$hist_enable)
      req(input$hist_run > 0)
      req(identical(input$hist_group_by, "Country") || identical(input$hist_group_by, "State (US only)"))
      
      df <- hist_importer_summary_tbl()
      df$import_value_dollars <- round(df$import_value_dollars, 2)
      df$import_quantity <- round(df$import_quantity, 2)
      df
    }, striped = TRUE, hover = TRUE, spacing = "s", width = "100%")
    
    ##  ##  ##  ##  ##  ##  ##  ##
    
    hist_bins_tbl <- eventReactive(input$hist_run, {
      req(input$hist_enable)
      req(con_rv())
      con <- con_rv()
      
      # Build WHERE clauses
      where <- c("value_dollars IS NOT NULL")
      
      # Year selection (exact hist_year OR fallback to year_range)
      if (!is.null(input$hist_year) && length(input$hist_year) > 0) {
        yrs <- as.integer(input$hist_year)
        where <- c(where, sprintf("Year IN (%s)", paste(yrs, collapse = ",")))
      } else if (!is.null(input$year_range)) {
        where <- c(where, sprintf("Year BETWEEN %d AND %d",
                                  as.integer(input$year_range[1]),
                                  as.integer(input$year_range[2])))
      }
      
      if (!is.null(input$hist_province) && length(input$hist_province) > 0) {
        where <- c(where, sprintf("Province IN (%s)",
                                  paste(DBI::dbQuoteString(con, input$hist_province), collapse = ",")))
      }
      
      if (!is.null(input$hist_commodity) && length(input$hist_commodity) > 0) {
        where <- c(where, sprintf("Commodity IN (%s)",
                                  paste(DBI::dbQuoteString(con, input$hist_commodity), collapse = ",")))
      }
      
      if (!is.null(input$hist_origin_country) && length(input$hist_origin_country) > 0) {
        where <- c(where, sprintf("Country IN (%s)",
                                  paste(DBI::dbQuoteString(con, input$hist_origin_country), collapse = ",")))
      }
      
      if (!is.null(input$hist_origin_state) && length(input$hist_origin_state) > 0) {
        where <- c(where, sprintf("State IN (%s)",
                                  paste(DBI::dbQuoteString(con, input$hist_origin_state), collapse = ",")))
      }
      
      if (isTRUE(input$hist_drop_zero_qty)) {
        where <- c(where, "Quantity IS NOT NULL AND Quantity != 0")
      }
      
      # If grouping by state, restrict to US
      if (identical(input$hist_group_by, "State (US only)")) {
        where <- c(where, "Country IN ('United States','USA','United States of America')")
      }
      
      where_sql <- paste(where, collapse = " AND ")
      
      # Choose grouping column
      group_col <- if (identical(input$hist_group_by, "Province")) {
        "Province"
      } else if (identical(input$hist_group_by, "Country")) {
        "Country"
      } else {
        "State"
      }
      
      # Choose x expression
      if (isTRUE(input$hist_log_value)) {
        x_expr <- "log10(GREATEST(value_dollars, 1e-6))"
      } else {
        x_expr <- "value_dollars"
      }
      
      nb <- as.integer(input$hist_bins)
      
      # 1) Get top groups (limit legend clutter, same as your old logic)
      top_sql <- sprintf("
    SELECT %s AS grp, COUNT(*) AS n
    FROM trade
    WHERE %s
    GROUP BY 1
    ORDER BY n DESC
    LIMIT 12
  ", group_col, where_sql)
      
      top_groups <- DBI::dbGetQuery(con, top_sql)
      shiny::validate(shiny::need(nrow(top_groups) > 0, "Histogram: no rows after filtering."))
      
      grp_in <- paste(DBI::dbQuoteString(con, top_groups$grp), collapse = ",")
      
      # 2) Get min/max of x over kept groups (for stable binning)
      rng_sql <- sprintf("
    SELECT MIN(%s) AS xmin, MAX(%s) AS xmax
    FROM trade
    WHERE %s AND %s IN (%s)
  ", x_expr, x_expr, where_sql, group_col, grp_in)
      
      rng <- DBI::dbGetQuery(con, rng_sql)
      xmin <- rng$xmin[1]
      xmax <- rng$xmax[1]
      shiny::validate(shiny::need(is.finite(xmin) && is.finite(xmax) && xmax > xmin, "Histogram: insufficient value range after filtering."))    
      # 3) Bin + count in SQL (tiny result)
      # width_bucket returns 1..nb (we’ll keep that)
      bin_width <- (xmax - xmin) / nb
      
      bin_sql <- sprintf("
  SELECT
    %s AS grp,
    LEAST(%d, GREATEST(1,
      1 + CAST(FLOOR((%s - %f) / %f) AS INTEGER)
    )) AS bin,
    COUNT(*) AS n
  FROM trade
  WHERE %s
    AND %s IN (%s)
    AND %s IS NOT NULL
  GROUP BY 1, 2
  ORDER BY 1, 2
", group_col, nb, x_expr, xmin, bin_width, where_sql, group_col, grp_in, x_expr)
      
      counts <- DBI::dbGetQuery(con, bin_sql)
      
      # Add bin centers for plotting
      bin_width <- (xmax - xmin) / nb
      counts$bin_center <- xmin + (counts$bin - 0.5) * bin_width
      counts$bin_width  <- bin_width
      
      # Optional density normalization per-group
      if (isTRUE(input$hist_density)) {
        counts <- counts |>
          dplyr::group_by(grp) |>
          dplyr::mutate(density = n / (sum(n) * bin_width)) |>
          dplyr::ungroup()
      }
      
      list(
        counts = counts,
        xmin = xmin,
        xmax = xmax,
        xlab = if (isTRUE(input$hist_log_value)) "Value (log CAD)" else "Value (CAD)",
        ylab = if (isTRUE(input$hist_density)) "Density" else "Count"
      )
    })
    
    output$hist_plot <- renderPlotly({
      req(input$hist_enable)
      
      h <- hist_bins_tbl()
      counts <- h$counts
      shiny::validate(shiny::need(nrow(counts) > 0, "Histogram: no rows after filtering."))
      
      p <- plotly::plot_ly()
      
      groups <- sort(unique(counts$grp))
      for (g in groups) {
        dd <- counts[counts$grp == g, , drop = FALSE]
        
        yvals <- if (isTRUE(input$hist_density)) dd$density else dd$n
        
        p <- p %>%
          plotly::add_bars(
            data = dd,
            x = ~bin_center,
            y = yvals,
            name = g,
            opacity = input$hist_alpha,
            hoverinfo = "text",
            text = ~paste0(
              "Group: ", grp,
              "<br>x: ", round(bin_center, 4),
              "<br>n: ", n,
              if (isTRUE(input$hist_density)) paste0("<br>density: ", signif(density, 4)) else ""
            )
          )
      }
      
      p %>%
        plotly::layout(
          barmode = "overlay",
          xaxis = list(title = h$xlab),
          yaxis = list(title = h$ylab),
          legend = list(orientation = "h", x = 0, y = -0.2, xanchor = "left"),
          margin = list(b = 120)
        )
    })
    
    output$hist_status <- renderPrint({
      req(input$hist_enable)
      
      if (is.null(input$hist_run) || input$hist_run == 0) {
        cat("Click 'Run histogram' to compute.\n")
        return()
      }
      
      h <- hist_bins_tbl()
      counts <- h$counts
      
      cat("Histogram groups:", length(unique(counts$grp)), "\n")
      cat("Histogram bins:", length(unique(counts$bin)), "\n")
      cat("x range:", signif(h$xmin, 4), "to", signif(h$xmax, 4), "\n")
    })
    
    output$cmp_prov_ui <- renderUI({
      provs <- choices_cache$provinces
      req(provs)
      
      if (input$cmp_mode == "Province vs Province") {
        tagList(
          selectizeInput(session$ns("cmp_prov_a"), "Province A", choices = provs, selected = provs[1], multiple = FALSE),
          selectizeInput(session$ns("cmp_prov_b"), "Province B", choices = provs, selected = provs[min(2, length(provs))], multiple = FALSE)
        )
      } else {
        selectizeInput(session$ns("cmp_prov_one"), "Province", choices = provs, selected = provs[1], multiple = FALSE)
      }
    })
    
    
    output$cmp_year_ui <- renderUI({
      yrs <- choices_cache$years
      req(yrs)
      
      if (input$cmp_mode == "Province vs Province") {
        selectizeInput(session$ns("cmp_years"), "Years (optional filter)",
                       choices = yrs, selected = NULL, multiple = TRUE,
                       options = list(placeholder = "All years (or use year_range/hist_year)"))
      } else {
        tagList(
          selectizeInput(session$ns("cmp_year_a"), "Year A", choices = yrs, selected = yrs[1], multiple = FALSE),
          selectizeInput(session$ns("cmp_year_b"), "Year B", choices = yrs, selected = yrs[min(2, length(yrs))], multiple = FALSE)
        )
      }
    })
    
    output$cmp_commodity_ui <- renderUI({
      comms <- choices_cache$commodities
      req(comms)
      selectizeInput(session$ns("cmp_commodity"), "Commodity(ies) (comparison filter)",
                     choices = comms, selected = NULL, multiple = TRUE,
                     options = list(placeholder = "All commodities", maxOptions = 2000))
    })
    
    output$cmp_origin_country_ui <- renderUI({
      req(input$cmp_origin_enable)
      cn <- choices_cache$countries
      req(cn)
      
      selectizeInput(session$ns("cmp_origin_country"), "Origin country (optional)",
                     choices = cn, selected = NULL, multiple = FALSE,
                     options = list(placeholder = "All countries"))
    })
    
    output$cmp_origin_state_ui <- renderUI({
      req(input$cmp_origin_enable)
      
      us_vals <- c("United States", "USA", "United States of America")
      if (is.null(input$cmp_origin_country) || !(input$cmp_origin_country %in% us_vals)) return(NULL)
      
      states <- choices_cache$states
      req(states)
      
      selectizeInput(session$ns("cmp_origin_state"), "Origin state (US only, optional)",
                     choices = states, selected = NULL, multiple = TRUE,
                     options = list(placeholder = "All US states"))
    })
    
    cmp_vectors <- eventReactive(input$cmp_run, {
      
      req(con_rv())
      con <- con_rv()
      
      # optional sanity check once (DuckDB)
      cols <- DBI::dbGetQuery(con, "DESCRIBE trade")$column_name
      shiny::validate(shiny::need("value_dollars" %in% cols, "DuckDB view 'trade' is missing value_dollars"))
      
      # --------- Build WHERE clause from UI ---------
      where <- c("value_dollars IS NOT NULL")
      
      # Use histogram filters if requested
      if (isTRUE(input$cmp_use_hist_filters)) {
        
        if (!is.null(input$hist_year) && length(input$hist_year) > 0) {
          yrs <- as.integer(input$hist_year)
          where <- c(where, sprintf("Year IN (%s)", paste(yrs, collapse = ",")))
        } else if (!is.null(input$year_range)) {
          where <- c(where, sprintf("Year BETWEEN %d AND %d",
                                    as.integer(input$year_range[1]),
                                    as.integer(input$year_range[2])))
        }
        
        if (!is.null(input$hist_commodity) && length(input$hist_commodity) > 0) {
          where <- c(where, sprintf("Commodity IN (%s)",
                                    paste(DBI::dbQuoteString(con, input$hist_commodity), collapse = ",")))
        }
        
        if (isTRUE(input$hist_drop_zero_qty)) {
          where <- c(where, "Quantity IS NOT NULL AND Quantity != 0")
        }
        
      } else {
        # Otherwise use global year_range (if present)
        if (!is.null(input$year_range)) {
          where <- c(where, sprintf("Year BETWEEN %d AND %d",
                                    as.integer(input$year_range[1]),
                                    as.integer(input$year_range[2])))
        }
      }
      
      # Comparison commodity filter
      if (!is.null(input$cmp_commodity) && length(input$cmp_commodity) > 0) {
        where <- c(where, sprintf("Commodity IN (%s)",
                                  paste(DBI::dbQuoteString(con, input$cmp_commodity), collapse = ",")))
      }
      
      # Optional origin filters
      if (isTRUE(input$cmp_origin_enable)) {
        
        if (!is.null(input$cmp_origin_country) && nzchar(input$cmp_origin_country)) {
          where <- c(where, sprintf("Country = %s", DBI::dbQuoteString(con, input$cmp_origin_country)))
        }
        
        us_vals <- c("United States", "USA", "United States of America")
        if (!is.null(input$cmp_origin_country) &&
            input$cmp_origin_country %in% us_vals &&
            !is.null(input$cmp_origin_state) &&
            length(input$cmp_origin_state) > 0) {
          where <- c(where, sprintf("State IN (%s)",
                                    paste(DBI::dbQuoteString(con, input$cmp_origin_state), collapse = ",")))
        }
      }
      
      where_sql <- paste(where, collapse = " AND ")
      
      # --------- Value expression (log or raw) ---------
      val_expr <- if (isTRUE(input$hist_log_value)) {
        "log10(GREATEST(value_dollars, 1e-6))"
      } else {
        "value_dollars"
      }
      
      # --------- Deterministic bounded sampling (NO SORT) ---------
      # cap how many rows we ever pull per group
      max_pull <- 150000L
      #seed <- as.integer(input$seed %||% 1)
      seed <- if (is.null(input$seed) || is.na(input$seed)) 1L else as.integer(input$seed)
      
      pull_group <- function(extra_where) {
        # Filter first, then reservoir sample from the filtered relation.
        # REPEATABLE(seed) makes it deterministic for a given seed.
        sql <- sprintf("
      SELECT val
      FROM (
        SELECT %s AS val
        FROM trade
        WHERE %s AND (%s) AND %s IS NOT NULL
      ) t
      TABLESAMPLE RESERVOIR(%d ROWS) REPEATABLE(%d)
    ", val_expr, where_sql, extra_where, val_expr, max_pull, seed)
        
        v <- DBI::dbGetQuery(con, sql)$val
        v <- v[is.finite(v) & !is.na(v)]
        v
      }
      
      # --------- Define groups A/B based on mode ---------
      if (identical(input$cmp_mode, "Province vs Province")) {
        req(input$cmp_prov_a, input$cmp_prov_b)
        
        # Optional cmp_years filter (only in this mode)
        if (!is.null(input$cmp_years) && length(input$cmp_years) > 0) {
          yrs2 <- as.integer(input$cmp_years)
          where_sql <- paste(where_sql, sprintf("AND Year IN (%s)", paste(yrs2, collapse = ",")))
        }
        
        a_where <- sprintf("Province = %s", DBI::dbQuoteString(con, input$cmp_prov_a))
        b_where <- sprintf("Province = %s", DBI::dbQuoteString(con, input$cmp_prov_b))
        
        a <- pull_group(a_where)
        b <- pull_group(b_where)
        
        list(
          x = a, y = b,
          label_x = input$cmp_prov_a,
          label_y = input$cmp_prov_b
        )
        
      } else {
        # Year vs Year within Province
        req(input$cmp_prov_one, input$cmp_year_a, input$cmp_year_b)
        
        a_where <- sprintf("Province = %s AND Year = %d",
                           DBI::dbQuoteString(con, input$cmp_prov_one),
                           as.integer(input$cmp_year_a))
        
        b_where <- sprintf("Province = %s AND Year = %d",
                           DBI::dbQuoteString(con, input$cmp_prov_one),
                           as.integer(input$cmp_year_b))
        
        a <- pull_group(a_where)
        b <- pull_group(b_where)
        
        list(
          x = a, y = b,
          label_x = paste0(input$cmp_prov_one, " ", input$cmp_year_a),
          label_y = paste0(input$cmp_prov_one, " ", input$cmp_year_b)
        )
      }
    })
    
    output$cmp_table <- renderTable({
      v <- cmp_vectors()
      x <- v$x; y <- v$y
      
      shiny::validate(shiny::need(length(safe_num(x)) >= 10, paste("Too few rows for:", v$label_x)))
      shiny::validate(shiny::need(length(safe_num(y)) >= 10, paste("Too few rows for:", v$label_y)))
      
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
      
      cat("Scale:", ifelse(isTRUE(input$hist_log_value), "log10(value_dollars)", "value_dollars"), "\n")
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
      
      provs <- choices_cache$provinces
      req(provs)
      
      if (identical(input$cmp_trend_mode, "commodity_lines")) {
        selectizeInput(
          session$ns("cmp_trend_provinces"),
          "Province(s) to include",
          choices = provs,
          selected = head(provs, 3),
          multiple = TRUE,
          options = list(placeholder = "Select 1+ provinces")
        )
      } else {
        selectizeInput(
          session$ns("cmp_trend_provinces"),
          "Province(s) (optional)",
          choices = provs,
          selected = NULL,
          multiple = TRUE,
          options = list(placeholder = "All provinces")
        )
      }
    })
    
    output$cmp_trend_year_ui <- renderUI({
      req(input$cmp_trend_enable)
      
      yrs <- choices_cache$years
      req(yrs)
      
      sliderInput(
        session$ns("cmp_trend_year_range"),
        "Trend year range",
        min = min(yrs), max = max(yrs),
        value = c(max(min(yrs), 2020), min(max(yrs), 2025)),
        step = 1, sep = ""
      )
    })
    
    output$cmp_trend_comm_ui <- renderUI({
      req(input$cmp_trend_enable)
      
      comms <- choices_cache$commodities
      req(comms)
      
      selectizeInput(
        session$ns("cmp_trend_commodities"),
        if (identical(input$cmp_trend_mode, "province_lines")) {
          "Commodity(ies) to include"
        } else {
          "Commodity(ies) (optional)"
        },
        choices = comms,
        selected = NULL,
        multiple = TRUE,
        options = list(placeholder = "All commodities", maxOptions = 2000)
      )
    })
    
    # Event-driven compute so it doesn't re-run constantly
    # helper for SQL IN (...) clauses
    sql_in <- function(x) {
      x <- x[!is.na(x) & nzchar(x)]
      if (length(x) == 0) return(NULL)
      paste0("('", gsub("'", "''", x), "')", collapse = ", ")
    }
    
    cmp_trend_data <- eventReactive(input$cmp_trend_run, {
      req(input$cmp_trend_enable)
      req(con_rv())
      con <- con_rv()
      
      req(input$cmp_trend_year_range)
      year_min <- as.integer(input$cmp_trend_year_range[1])
      year_max <- as.integer(input$cmp_trend_year_range[2])
      
      # Required selections depending on mode
      # if (identical(input$cmp_trend_mode, "province_lines")) {
      #   shiny::validate(shiny::need(!is.null(input$cmp_trend_commodities) && length(input$cmp_trend_commodities) > 0,
      #                 "Select 1+ commodities to plot province trends."))
      # } else {
      #   shiny::validate(shiny::need(!is.null(input$cmp_trend_provinces) && length(input$cmp_trend_provinces) > 0,
      #                 "Select 1+ provinces to plot commodity trends."))
      # }
      
      if (identical(input$cmp_trend_mode, "commodity_lines")) {
        shiny::validate(shiny::need(!is.null(input$cmp_trend_provinces) && length(input$cmp_trend_provinces) > 0,
                                    "Select 1+ provinces to plot commodity trends."))
      }
      
      # Choose value expression inside SQL (log or raw)
      if (isTRUE(input$cmp_trend_log_value)) {
        ylab <- "log10(Value ($))"
        val_expr <- "log10(GREATEST(value_dollars, 1e-6))"
      } else {
        ylab <- "Value ($)"
        val_expr <- "value_dollars"
      }
      
      # Optional filters -> SQL WHERE fragments
      where <- c(
        sprintf("Year BETWEEN %d AND %d", year_min, year_max),
        "value_dollars IS NOT NULL"
      )
      
      if (isTRUE(input$cmp_trend_drop_zero_qty)) {
        where <- c(where, "Quantity IS NOT NULL AND Quantity != 0")
      }
      
      # Optional Province filter
      if (!is.null(input$cmp_trend_provinces) && length(input$cmp_trend_provinces) > 0) {
        where <- c(where, sprintf("Province IN (%s)",
                                  paste(DBI::dbQuoteString(con, input$cmp_trend_provinces), collapse = ",")))
      }
      
      # Optional Commodity filter
      if (!is.null(input$cmp_trend_commodities) && length(input$cmp_trend_commodities) > 0) {
        where <- c(where, sprintf("Commodity IN (%s)",
                                  paste(DBI::dbQuoteString(con, input$cmp_trend_commodities), collapse = ",")))
      }
      
      where_sql <- paste(where, collapse = " AND ")
      
      # Grouping depends on mode
      if (identical(input$cmp_trend_mode, "province_lines")) {
        group_dim <- "Province"
        subtitle <- paste0("Commodities: ", paste(input$cmp_trend_commodities, collapse = ", "))
      } else {
        group_dim <- "Commodity"
        subtitle <- paste0("Provinces: ", paste(input$cmp_trend_provinces, collapse = ", "))
      }
      
      # DuckDB quantiles (exact). If you want faster, we can swap to approx_quantile().
      sql <- sprintf("
    SELECT
      %s AS group_dim,
      Year,
      COUNT(*) AS n,
      quantile_cont(%s, 0.25) AS q25,
      quantile_cont(%s, 0.50) AS med,
      quantile_cont(%s, 0.75) AS q75
    FROM trade
    WHERE %s
    GROUP BY 1, 2
    ORDER BY 1, 2
  ", group_dim, val_expr, val_expr, val_expr, where_sql)
      
      summ <- DBI::dbGetQuery(con, sql)
      
      shiny::validate(shiny::need(nrow(summ) > 0,
                                  paste0("Trend plot: no rows after filtering (", year_min, "–", year_max, ").")))
      
      # Return in the same format your plotting code expects
      list(
        summ = summ %>% dplyr::rename(!!group_dim := group_dim),
        line_dim = group_dim,
        ylab = ylab,
        subtitle = subtitle
      )
    })
    
    cmp_trend_forecast_data <- reactive({
      req(input$cmp_trend_enable)
      req(isTRUE(input$cmp_trend_forecast_enable))
      
      td <- cmp_trend_data()
      summ <- td$summ
      line_dim <- td$line_dim
      
      lvls <- sort(unique(as.character(summ[[line_dim]])))
      
      fc <- lapply(lvls, function(g) {
        dd <- summ %>%
          dplyr::filter(.data[[line_dim]] == g) %>%
          dplyr::arrange(Year)
        
        pr <- trend_projection_tbl(dd, year_horizon = 5L)
        if (is.null(pr)) return(NULL)
        
        pr[[line_dim]] <- g
        pr
      })
      
      dplyr::bind_rows(fc)
    })
    
    # ---- Forecasting compute (new button; uses ALL years up to 2025) ----
    forecast_results <- eventReactive(input$forecast_run, {
      req(con_rv())
      con <- con_rv()
      
      req(input$cmp_trend_enable)
      req(input$cmp_trend_mode)
      
      year_min <- as.integer(FORECAST_MIN_YEAR)
      year_max <- as.integer(FORECAST_LAST_YEAR)
      h <- as.integer(FORECAST_HORIZON_YEARS)
      future_years <- seq.int(year_max + 1L, year_max + h)
      
      # Determine what "series" means based on trend mode + selections
      mode <- input$cmp_trend_mode
      
      # Filters: ALWAYS drop Quantity=0 for forecasting
      where <- c(
        sprintf("Year BETWEEN %d AND %d", year_min, year_max),
        "value_dollars IS NOT NULL",
        "Quantity IS NOT NULL AND Quantity != 0"
      )
      
      # In your trend UI:
      # - province_lines means: lines=Province, user selects commodities
      # - commodity_lines means: lines=Commodity, user selects provinces
      if (identical(mode, "province_lines")) {
        req(!is.null(input$cmp_trend_commodities) && length(input$cmp_trend_commodities) > 0)
        where <- c(where, sprintf("Commodity IN (%s)",
                                  paste(DBI::dbQuoteString(con, input$cmp_trend_commodities), collapse = ",")))
        series_dim <- "Province"
        series_vals <- input$cmp_trend_provinces
        # If provinces not specified, allow “all provinces” BUT it can be huge;
        # keep v1 safe: require at least one province if none selected
        if (is.null(series_vals) || length(series_vals) == 0) {
          stop("Forecasting: please select 1+ provinces (so we don’t forecast every province).")
        }
        where <- c(where, sprintf("Province IN (%s)",
                                  paste(DBI::dbQuoteString(con, series_vals), collapse = ",")))
        subtitle <- paste0("Commodities: ", paste(input$cmp_trend_commodities, collapse = ", "))
        
      } else {
        req(!is.null(input$cmp_trend_provinces) && length(input$cmp_trend_provinces) > 0)
        where <- c(where, sprintf("Province IN (%s)",
                                  paste(DBI::dbQuoteString(con, input$cmp_trend_provinces), collapse = ",")))
        series_dim <- "Commodity"
        series_vals <- input$cmp_trend_commodities
        if (is.null(series_vals) || length(series_vals) == 0) {
          stop("Forecasting: please select 1+ commodities (so we don’t forecast every commodity).")
        }
        where <- c(where, sprintf("Commodity IN (%s)",
                                  paste(DBI::dbQuoteString(con, series_vals), collapse = ",")))
        subtitle <- paste0("Provinces: ", paste(input$cmp_trend_provinces, collapse = ", "))
      }
      
      where_sql <- paste(where, collapse = " AND ")
      
      # Pull annual totals for each series (in DuckDB, small result)
      sql <- sprintf("
      SELECT
        %s AS series,
        Year,
        SUM(value_dollars) AS total_value
      FROM trade
      WHERE %s
      GROUP BY 1, 2
      ORDER BY 1, 2
    ", series_dim, where_sql)
      
      annual <- DBI::dbGetQuery(con, sql)
      shiny::validate(shiny::need(nrow(annual) > 0, "Forecasting: no rows after filtering."))
      
      # Fill missing years with NA per series (to handle gaps)
      year_grid <- seq.int(year_min, year_max)
      
      # Build wide matrix: rows=years, cols=series
      series_levels <- sort(unique(annual$series))
      wide <- tidyr::pivot_wider(annual, names_from = series, values_from = total_value)
      wide <- make_year_grid(wide, year_min = year_min, year_max = year_max)
      
      # Ensure ordering
      wide <- wide[order(wide$Year), , drop = FALSE]
      Ymat <- as.matrix(wide[, setdiff(names(wide), "Year"), drop = FALSE])
      colnames(Ymat) <- setdiff(names(wide), "Year")
      
      # METHOD 1: ARIMA per series (on log1p(total_value))
      method1 <- lapply(colnames(Ymat), function(s) {
        y <- Ymat[, s]
        fr <- arima_forecast_base(f_log1p(y), h = h, level = 0.95)
        
        # Back-transform to dollars (median-like)
        tibble::tibble(
          series = s,
          Year = future_years,
          method = "ARIMA (univariate)",
          point = f_expm1(fr$mean),
          lo = f_expm1(fr$lo),
          hi = f_expm1(fr$hi)
        )
      })
      m1_tbl <- dplyr::bind_rows(method1)
      
      # METHOD 2A: pooled aggregate over selected series (single series)
      pooled <- rowSums(Ymat, na.rm = TRUE)
      fr_pool <- arima_forecast_base(f_log1p(pooled), h = h, level = 0.95)
      m2a_tbl <- tibble::tibble(
        series = "ALL_SELECTED (pooled)",
        Year = future_years,
        method = "ARIMA (pooled aggregate)",
        point = f_expm1(fr_pool$mean),
        lo = f_expm1(fr_pool$lo),
        hi = f_expm1(fr_pool$hi)
      )
      
      # METHOD 2B: correlation/factor method (only if >=2 series)
      m2b_tbl <- NULL
      if (ncol(Ymat) >= 2) {
        ff <- factor_forecast(Ymat, h = h, k = min(2, ncol(Ymat)), level = 0.95)
        m2b_tbl <- tibble::as_tibble(ff$mean)
        m2b_lo  <- tibble::as_tibble(ff$lo)
        m2b_hi  <- tibble::as_tibble(ff$hi)
        
        m2b_tbl$Year <- future_years
        m2b_lo$Year  <- future_years
        m2b_hi$Year  <- future_years
        
        m2b_long <- tidyr::pivot_longer(m2b_tbl, -Year, names_to = "series", values_to = "point")
        m2b_long$method <- "Correlation (factor model)"
        
        lo_long <- tidyr::pivot_longer(m2b_lo, -Year, names_to = "series", values_to = "lo")
        hi_long <- tidyr::pivot_longer(m2b_hi, -Year, names_to = "series", values_to = "hi")
        
        m2b_tbl <- dplyr::left_join(m2b_long, lo_long, by = c("Year","series")) |>
          dplyr::left_join(hi_long, by = c("Year","series"))
      }
      
      # METHOD 3: Random Forest regression (optional; only if ranger is installed)
      m3_tbl <- NULL
      if (requireNamespace("ranger", quietly = TRUE)) {
        # Simple supervised setup per series:
        # features: lag1, lag2, lag3 of log1p(total_value)
        make_rf_forecast <- function(y_raw) {
          y <- f_log1p(y_raw)
          # build training df
          df <- data.frame(
            Year = year_grid,
            y = y,
            lag1 = dplyr::lag(y, 1),
            lag2 = dplyr::lag(y, 2),
            lag3 = dplyr::lag(y, 3)
          )
          df <- df[is.finite(df$y) & is.finite(df$lag1) & is.finite(df$lag2), , drop = FALSE]
          if (nrow(df) < 8) stop("Too few rows for RF (need >=8 after lagging).")
          
          fit <- ranger::ranger(y ~ lag1 + lag2 + lag3, data = df, num.trees = 400)
          
          # recursive forecasting
          preds <- numeric(h)
          # start from last observed (use last non-NA y)
          y_last <- y
          last_idx <- max(which(is.finite(y_last) & !is.na(y_last)))
          cur <- y_last[last_idx]
          cur_l1 <- y_last[last_idx]
          cur_l2 <- y_last[last_idx - 1]
          cur_l3 <- y_last[last_idx - 2]
          
          for (i in seq_len(h)) {
            nd <- data.frame(lag1 = cur_l1, lag2 = cur_l2, lag3 = cur_l3)
            pr <- predict(fit, data = nd)$predictions
            preds[i] <- as.numeric(pr)
            # shift lags
            cur_l3 <- cur_l2
            cur_l2 <- cur_l1
            cur_l1 <- preds[i]
          }
          
          # crude interval: use OOB prediction error sd as constant band (v1)
          # (keeps things simple; better intervals later)
          resid <- df$y - predict(fit, data = df)$predictions
          s <- stats::sd(resid, na.rm = TRUE)
          z <- stats::qnorm(0.975)
          lo <- preds - z * s
          hi <- preds + z * s
          
          list(point = preds, lo = lo, hi = hi)
        }
        
        m3_list <- lapply(colnames(Ymat), function(s) {
          y <- Ymat[, s]
          rf <- make_rf_forecast(y)
          tibble::tibble(
            series = s,
            Year = future_years,
            method = "Random Forest (lag features)",
            point = f_expm1(rf$point),
            lo = f_expm1(rf$lo),
            hi = f_expm1(rf$hi)
          )
        })
        m3_tbl <- dplyr::bind_rows(m3_list)
      }
      
      # Build "actual" long for plotting (dollars)
      actual_long <- data.frame(Year = wide$Year, stringsAsFactors = FALSE)
      actual_long <- cbind(actual_long, as.data.frame(Ymat))
      actual_long <- tidyr::pivot_longer(actual_long, -Year, names_to = "series", values_to = "total_value") |>
        dplyr::mutate(method = "Actual")
      
      # Combine forecast tables
      fc_tbl <- dplyr::bind_rows(m1_tbl, m2a_tbl, m2b_tbl, m3_tbl) |>
        dplyr::mutate(
          point = as.numeric(point),
          lo = as.numeric(lo),
          hi = as.numeric(hi)
        )
      
      list(
        subtitle = subtitle,
        actual_long = actual_long,
        forecast_tbl = fc_tbl,
        year_min = year_min,
        year_max = year_max,
        future_years = future_years
      )
    })
    
    output$cmp_trend_plot <- renderPlotly({
      req(input$cmp_trend_enable)
      td <- cmp_trend_data()
      summ <- td$summ
      #shiny::validate(shiny::need(nrow(summ) > 0, "Trend plot: no rows after filtering (2020–2025)."))
      yr <- input$cmp_trend_year_range
      shiny::validate(shiny::need(nrow(summ) > 0,
                                  paste0("Trend plot: no rows after filtering (", yr[1], "–", yr[2], ").")))
      
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
              width = 3,
              color = unname(cols[g])
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
      
      # ---- Optional overlay forecast tied to the currently displayed trend plot ----
      if (isTRUE(input$cmp_trend_forecast_enable)) {
        fc <- cmp_trend_forecast_data()
        
        if (!is.null(fc) && nrow(fc) > 0) {
          show_methods <- input$cmp_trend_forecast_methods
          if (is.null(show_methods) || length(show_methods) == 0) {
            show_methods <- character(0)
          }
          
          fc <- fc %>%
            dplyr::filter(
              (method == "Linear" & "linear" %in% show_methods) |
                (method == "Nonlinear" & "nonlinear" %in% show_methods)
            )
          
          if (nrow(fc) > 0) {
            for (g in lvls) {
              dfg <- fc %>%
                dplyr::filter(.data[[line_dim]] == g) %>%
                dplyr::arrange(Year)
              
              if (nrow(dfg) == 0) next
              
              # color matched to the historical series
              base_col <- unname(cols[g])
              
              # CI ribbons
              if (isTRUE(input$cmp_trend_forecast_ci)) {
                for (m in unique(dfg$method)) {
                  ddr <- dfg %>% dplyr::filter(method == m)
                  if (nrow(ddr) == 0) next
                  
                  p <- p %>%
                    plotly::add_ribbons(
                      data = ddr,
                      x = ~Year,
                      ymin = ~lo,
                      ymax = ~hi,
                      name = paste0(g, " ", m, " CI"),
                      showlegend = FALSE,
                      hoverinfo = "skip",
                      opacity = if (m == "Linear") 0.10 else 0.16,
                      line = list(width = 0)
                    )
                }
              }
              
              # forecast lines
              if ("linear" %in% show_methods) {
                ddl <- dfg %>% dplyr::filter(method == "Linear")
                if (nrow(ddl) > 0) {
                  p <- p %>%
                    plotly::add_trace(
                      data = ddl,
                      x = ~Year,
                      y = ~point,
                      type = "scatter",
                      mode = "lines+markers",
                      name = paste0(g, " forecast (Linear)"),
                      line = list(width = 2, dash = "dash", color = base_col),
                      marker = list(size = 6, color = base_col),
                      hoverinfo = "text",
                      text = ~paste0(
                        line_dim, ": ", .data[[line_dim]],
                        "<br>Method: Linear",
                        "<br>Year: ", Year,
                        "<br>Forecast: ", round(point, 3),
                        "<br>CI: [", round(lo, 3), ", ", round(hi, 3), "]"
                      )
                    )
                }
              }
              
              if ("nonlinear" %in% show_methods) {
                ddn <- dfg %>% dplyr::filter(method == "Nonlinear")
                if (nrow(ddn) > 0) {
                  p <- p %>%
                    plotly::add_trace(
                      data = ddn,
                      x = ~Year,
                      y = ~point,
                      type = "scatter",
                      mode = "lines+markers",
                      name = paste0(g, " forecast (Nonlinear)"),
                      line = list(width = 2, dash = "dot", color = base_col),
                      marker = list(size = 6, color = base_col),
                      hoverinfo = "text",
                      text = ~paste0(
                        line_dim, ": ", .data[[line_dim]],
                        "<br>Method: Nonlinear",
                        "<br>Year: ", Year,
                        "<br>Forecast: ", round(point, 3),
                        "<br>CI: [", round(lo, 3), ", ", round(hi, 3), "]"
                      )
                    )
                }
              }
            }
          }
        }
      }
      
      p %>%
        plotly::layout(
          title = list(
            text = paste0(
              "Trend (Median \u00B1 IQR), ",
              min(summ$Year, na.rm = TRUE), "\u2013", max(summ$Year, na.rm = TRUE),
              if (isTRUE(input$cmp_trend_forecast_enable)) " + 5-year forecast" else "",
              "<br><sup>", td$subtitle, "</sup>"
            )
          ),
          xaxis = list(
            title = list(text = "Year", font = list(size = 18)),
            tickfont = list(size = 18)
          ),
          yaxis = list(
            title = list(
              text = "Value",
              font = list(size = 18)   # <-- increase title size
            ),
            tickfont = list(size = 18) # <-- increase tick label size
          ),
          legend = list(orientation = "h", x = 0, y = -0.25, xanchor = "left"),
          margin = list(b = 120)
        )
    })
    
    output$forecast_table <- renderTable({
      fr <- forecast_results()
      tbl <- fr$forecast_tbl
      
      # Keep it readable: show best 3 methods by default? For now show all rows.
      # Round dollars
      tbl <- tbl |>
        dplyr::mutate(
          point = round(point, 2),
          lo = round(lo, 2),
          hi = round(hi, 2)
        ) |>
        dplyr::arrange(method, series, Year)
      
      tbl
    }, striped = TRUE, bordered = TRUE, spacing = "s", width = "100%")
    
    output$forecast_status <- renderPrint({
      if (is.null(input$forecast_run) || input$forecast_run == 0) {
        cat("Click 'Run forecast' to compute 5-year projections.\n")
        return()
      }
      fr <- forecast_results()
      cat("Forecast horizon:", FORECAST_HORIZON_YEARS, "years\n")
      cat("Year window used:", fr$year_min, "–", fr$year_max, "\n")
      cat("Selection:", fr$subtitle, "\n")
      cat("Actual rows:", nrow(fr$actual_long), "\n")
      cat("Forecast rows:", nrow(fr$forecast_tbl), "\n")
      cat("Methods included:", paste(sort(unique(fr$forecast_tbl$method)), collapse = " | "), "\n")
    })
    
    output$forecast_plot <- renderPlotly({
      fr <- forecast_results()
      actual <- fr$actual_long
      fc <- fr$forecast_tbl
      
      shiny::validate(shiny::need(nrow(actual) > 0, "No actual data to plot."))
      shiny::validate(shiny::need(nrow(fc) > 0, "No forecast results to plot."))
      
      # ---- Map method labels -> method keys used by the UI ----
      fc <- fc %>%
        dplyr::mutate(
          method_key = dplyr::case_when(
            method == "ARIMA (univariate)" ~ "arima_uni",
            method == "ARIMA (pooled aggregate)" ~ "arima_pooled",
            method == "Correlation (factor model)" ~ "factor",
            method == "Random Forest (lag features)" ~ "rf",
            TRUE ~ "other"
          )
        )
      
      # ---- Filter to selected methods ----
      req(input$fc_methods)
      fc <- fc %>% dplyr::filter(method_key %in% input$fc_methods)
      shiny::validate(shiny::need(nrow(fc) > 0, "No forecasts left after method filtering."))
      
      # ---- Y transform for plotting only (keeps forecast math unchanged) ----
      eps <- 1e-6
      
      actual <- actual %>%
        dplyr::mutate(
          y_plot = if (isTRUE(input$fc_log_y)) log10(pmax(total_value, eps)) else total_value
        )
      
      fc <- fc %>%
        dplyr::mutate(
          point_plot = if (isTRUE(input$fc_log_y)) log10(pmax(point, eps)) else point,
          lo_plot    = if (isTRUE(input$fc_log_y)) log10(pmax(lo, eps)) else lo,
          hi_plot    = if (isTRUE(input$fc_log_y)) log10(pmax(hi, eps)) else hi
        )
      
      # ---- Robust y-axis range (clip upper tail) ----
      y_candidates <- c(actual$y_plot, fc$point_plot)
      if (isTRUE(input$fc_show_ci)) {
        y_candidates <- c(y_candidates, fc$hi_plot)
      }
      y_candidates <- y_candidates[is.finite(y_candidates) & !is.na(y_candidates)]
      
      shiny::validate(shiny::need(length(y_candidates) >= 2, "Not enough values to set y-axis range."))
      
      ymin <- min(y_candidates, na.rm = TRUE)
      ymax <- as.numeric(stats::quantile(y_candidates, probs = input$fc_y_clip, na.rm = TRUE))
      
      if (!is.finite(ymax) || ymax <= ymin) {
        ymax <- max(y_candidates, na.rm = TRUE)
      }
      
      y_title <- if (isTRUE(input$fc_log_y)) "Annual Total Import value (log CAD)" else "Annual total import value (CAD)"
      
      # ---- Color map across ALL plotted series (actual + forecast + pooled) ----
      all_series <- sort(unique(c(actual$series, fc$series)))
      cols <- make_color_map(all_series)
      
      p <- plotly::plot_ly()
      
      # ---- Actual: solid ----
      for (s in sort(unique(actual$series))) {
        dd <- actual %>% dplyr::filter(series == s) %>% dplyr::arrange(Year)
        
        p <- p %>% plotly::add_trace(
          data = dd,
          x = ~Year, y = ~y_plot,
          type = "scatter", mode = "lines+markers",
          name = paste0(s, " (Actual)"),
          line = list(width = 2, color = unname(cols[s])),
          marker = list(size = 5, color = unname(cols[s])),
          hoverinfo = "text",
          text = ~paste0(
            "Series: ", series,
            "<br>Year: ", Year,
            "<br>Total $: ", round(total_value, 2),
            if (isTRUE(input$fc_log_y)) paste0("<br>log10: ", round(y_plot, 4)) else ""
          )
        )
      }
      
      # ---- Forecast: dashed, optional ribbons ----
      method_order <- c(
        "ARIMA (univariate)",
        "ARIMA (pooled aggregate)",
        "Correlation (factor model)",
        "Random Forest (lag features)"
      )
      plot_methods <- intersect(method_order, unique(fc$method))
      
      for (m in plot_methods) {
        ddm <- fc %>% dplyr::filter(method == m)
        
        for (s in sort(unique(ddm$series))) {
          dds <- ddm %>% dplyr::filter(series == s) %>% dplyr::arrange(Year)
          
          # CI ribbon (optional)
          if (isTRUE(input$fc_show_ci)) {
            p <- p %>% plotly::add_ribbons(
              data = dds,
              x = ~Year,
              ymin = ~lo_plot,
              ymax = ~hi_plot,
              name = paste0(s, " (", m, " CI)"),
              showlegend = FALSE,
              hoverinfo = "skip",
              opacity = 0.12,
              line = list(width = 0)
            )
          }
          
          # Forecast line (dashed)
          p <- p %>% plotly::add_trace(
            data = dds,
            x = ~Year, y = ~point_plot,
            type = "scatter", mode = "lines+markers",
            name = paste0(s, " (", m, ")"),
            line = list(width = 2, dash = "dash", color = unname(cols[s])),
            marker = list(size = 5, color = unname(cols[s])),
            hoverinfo = "text",
            text = ~paste0(
              "Series: ", series,
              "<br>Method: ", method,
              "<br>Year: ", Year,
              "<br>Point $: ", round(point, 2),
              "<br>CI $: [", round(lo, 2), ", ", round(hi, 2), "]",
              if (isTRUE(input$fc_log_y)) paste0(
                "<br>Point log10: ", round(point_plot, 4),
                "<br>CI log10: [", round(lo_plot, 4), ", ", round(hi_plot, 4), "]"
              ) else ""
            )
          )
        }
      }
      
      p %>% plotly::layout(
        title = list(text = paste0("Forecast (next 5 years)<br><sup>", fr$subtitle, "</sup>")),
        xaxis = list(title = "Year", tickmode = "linear"),
        yaxis = list(title = y_title, range = c(ymin, ymax)),
        legend = list(orientation = "h", x = 0, y = -0.25, xanchor = "left"),
        margin = list(b = 160)
      )
    })
    
    output$cmp_trend_status <- renderPrint({
      req(input$cmp_trend_enable)
      
      req(input$cmp_trend_year_range)
      cat("Trend window:", input$cmp_trend_year_range[1], "–", input$cmp_trend_year_range[2], "\n")
      
      cat("Mode:", ifelse(identical(input$cmp_trend_mode, "province_lines"),
                          "Lines = Provinces (filter by commodities)",
                          "Lines = Commodities (filter by provinces)"), "\n")
      
      cat("Trend scale:", ifelse(isTRUE(input$cmp_trend_log_value), "log10(value_dollars)", "value_dollars"), "\n")
      
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
    
  })
}

