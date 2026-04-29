# scripts/build_parquet.R
# scripts/build_parquet.R
# Convert CSV -> Parquet with standardized column names for the app.

library(readr)
library(dplyr)
library(arrow)

in_csv  <- file.path("data", "all_CADpharma_imports_1988_2025.csv")
#out_pq  <- file.path("data", "all_CADpharma_imports_1988_2025.parquet")
out_pq  <- file.path("data", "all_CADpharma_imports_BROWSER_DATA_2020_2025.parquet")

df <- readr::read_csv(in_csv, show_col_types = FALSE, trim_ws = TRUE)

# Keep only needed cols (using original CSV names), then rename to app-friendly names
df <- df %>%
  dplyr::select(
    Period, Province, Country, State, Commodity,
    `Unit of measure`, `Value ($)`, Quantity
  ) %>%
  dplyr::rename(
    unit_of_measure = `Unit of measure`,
    value_dollars   = `Value ($)`
  )

arrow::write_parquet(df, out_pq, compression = "zstd")
cat("Wrote:", out_pq, "\n")
cat("Columns:", paste(names(df), collapse = ", "), "\n")

pq <- arrow::read_parquet("data/all_CADpharma_imports_BROWSER_DATA_2020_2025.parquet")
names(pq)