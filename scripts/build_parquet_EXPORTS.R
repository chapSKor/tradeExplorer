# scripts/build_parquet_EXPORTS.R
# Convert the monthly Canadian pharmaceutical EXPORT CSV -> Parquet with
# standardized column names matching the PharmaTrace import app schema.
#
# Input:
#   data_raw/EXPORT_2000-2026_byMonth.csv
# Output:
#   data/all_CADpharma_exports_2000_2026_byMonth.parquet
#
# Notes:
# - Export records are monthly, not annual. Period is preserved as a Date.
# - Year and Month are also added for the future Exports tab controls.
# - The core app-facing names are kept parallel to imports:
#     Period, Province, Country, State, Commodity,
#     unit_of_measure, value_dollars, Quantity

library(readr)
library(dplyr)
library(arrow)
library(lubridate)

in_csv <- file.path("data_raw", "EXPORT_2000-2026_byMonth.csv")
out_dir <- "data"
out_pq <- file.path(out_dir, "all_CADpharma_exports_2000_2026_byMonth.parquet")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(in_csv)) {
  stop("Missing export CSV: ", in_csv)
}

read_trade_csv <- function(path) {
  df <- readr::read_csv(path, show_col_types = FALSE, trim_ws = TRUE)

  # Back-compat with files whose first row is a title such as "Exports".
  if (!("Period" %in% names(df))) {
    df <- readr::read_csv(path, skip = 1, show_col_types = FALSE, trim_ws = TRUE)
  }

  if (!("Period" %in% names(df))) {
    stop("Could not find a Period column after reading: ", path)
  }

  df
}

parse_period_monthly <- function(x) {
  x_chr <- trimws(as.character(x))

  # Try common date/month encodings seen in Statistics Canada-style exports.
  parsed <- suppressWarnings(as.Date(x_chr, format = "%Y-%m-%d"))

  if (all(is.na(parsed))) parsed <- suppressWarnings(as.Date(x_chr, format = "%m/%d/%Y"))
  if (all(is.na(parsed))) parsed <- suppressWarnings(as.Date(x_chr, format = "%Y/%m/%d"))
  if (all(is.na(parsed))) parsed <- suppressWarnings(as.Date(paste0(x_chr, "-01"), format = "%Y-%m-%d"))
  if (all(is.na(parsed))) parsed <- suppressWarnings(as.Date(paste0(x_chr, "/01"), format = "%Y/%m/%d"))
  if (all(is.na(parsed))) parsed <- suppressWarnings(lubridate::ymd(x_chr, quiet = TRUE))
  if (all(is.na(parsed))) parsed <- suppressWarnings(lubridate::mdy(x_chr, quiet = TRUE))
  if (all(is.na(parsed))) parsed <- suppressWarnings(lubridate::my(x_chr, quiet = TRUE))
  if (all(is.na(parsed))) parsed <- suppressWarnings(lubridate::ym(x_chr, quiet = TRUE))

  parsed
}

standardize_trade_cols <- function(df) {
  required_cols <- c("Period", "Province", "Country", "State", "Commodity", "Unit of measure", "Value ($)", "Quantity")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop(
      "Export CSV is missing required columns: ", paste(missing_cols, collapse = ", "),
      "\nAvailable columns: ", paste(names(df), collapse = ", ")
    )
  }

  out <- df %>%
    dplyr::select(
      Period, Province, Country, State, Commodity,
      `Unit of measure`, `Value ($)`, Quantity
    ) %>%
    dplyr::rename(
      unit_of_measure = `Unit of measure`,
      value_dollars   = `Value ($)`
    ) %>%
    dplyr::filter(!(is.na(Period) | trimws(as.character(Period)) == "")) %>%
    dplyr::mutate(
      Period = parse_period_monthly(Period)
    )

  if (all(is.na(out$Period))) {
    stop(
      "Period parsing failed. Example Period values: ",
      paste(head(as.character(df$Period), 5), collapse = ", "),
      ". Expected values like YYYY-mm-dd, m/d/YYYY, YYYY-mm, YYYY/mm, Mon YYYY, or Month YYYY."
    )
  }

  out %>%
    dplyr::mutate(
      Year = as.integer(format(Period, "%Y")),
      Month = as.integer(format(Period, "%m")),
      YearMonth = format(Period, "%Y-%m"),
      State = dplyr::na_if(State, "N/A"),
      unit_of_measure = dplyr::na_if(unit_of_measure, "N/A"),
      value_dollars = suppressWarnings(as.numeric(value_dollars)),
      Quantity = suppressWarnings(as.numeric(Quantity))
    ) %>%
    dplyr::select(
      Period, Year, Month, YearMonth,
      Province, Country, State, Commodity,
      unit_of_measure, value_dollars, Quantity
    )
}

df_raw <- read_trade_csv(in_csv)
df <- standardize_trade_cols(df_raw)
df <- df %>%
  mutate(
    YearMonthDate = as.Date(paste0(YearMonth, "-01"))
  )

arrow::write_parquet(df, out_pq, compression = "zstd")

cat("Wrote: ", out_pq, "\n", sep = "")
cat("Rows: ", nrow(df), "\n", sep = "")
cat("Columns: ", paste(names(df), collapse = ", "), "\n", sep = "")
cat("Date range: ", min(df$Period, na.rm = TRUE), " to ", max(df$Period, na.rm = TRUE), "\n", sep = "")
cat("YearMonth range: ", min(df$YearMonth, na.rm = TRUE), " to ", max(df$YearMonth, na.rm = TRUE), "\n", sep = "")

range(as.Date(df$Period), na.rm = TRUE)
