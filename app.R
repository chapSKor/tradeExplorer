# app.R
# Pharma Trade Explorer (t-SNE / UMAP)
# - Local CSV picker (no upload dialog)
# - Robust feature construction for mixed numeric + categorical data
#   * Handles spaces in column names (e.g., `unit_of_measure`)
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
library(RColorBrewer)
library(dplyr)
library(dbplyr)
library(readr)
library(Matrix)
library(irlba)
library(plotly)
library(tidyr)
library(MASS)        # kde2d for density contours
library(Rtsne)       # t-SNE
library(uwot)        # UMAP

library(DBI)
library(duckdb)
library(ranger)

#Nature earth data is a large package, can run slow. 
library(sf)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(pak)
library(splines)


# ---------- Helpers ----------
source("modules/mod_trade_dashboard.R")

app_dir <- function() {
  # Most reliable in Shiny contexts
  d <- shiny::getShinyOption("appDir")
  if (!is.null(d)) return(d)
  # Fallback
  getwd()
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


# ---------- UI ----------

ui <- navbarPage(
  title = "PharmaTrace",
  
  header = tags$head(
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
  
  tabPanel(
    title = "Imports",
    tradeDashboardUI("imports", title = "Imports Dashboard")
  ),
  
  tabPanel(
    title = "Exports",
    tradeDashboardUI("exports", title = "Exports Dashboard")
  )
)

server <- function(input, output, session) {
  
  tradeDashboardServer(
    id = "imports",
    parquet_file = "all_CADpharma_imports_1988_2025.parquet",
    sorted_vars_file = "sorted_vars.rds",
    table_name = "trade_imports",
    monthly = FALSE,
    data_label = "imports"
  )
  
  tradeDashboardServer(
    id = "exports",
    parquet_file = "all_CADpharma_exports_2000_2026_byMonth.parquet",
    sorted_vars_file = "sorted_vars_EXPORTS_2026.rds",
    table_name = "trade_exports",
    monthly = TRUE,
    data_label = "exports"
  )
}

shinyApp(ui, server)

