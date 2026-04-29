# scripts/build_sorted_vars_EXPORTS.R
# Build UI choice caches for the monthly PharmaTrace EXPORT data.
#
# Input:
#   data_raw/EXPORT_2000-2026_byMonth.csv
# Output:
#   sortedVars/sorted_vars_EXPORTS_2026.rds
#
# The exports cache mirrors the import sorted_vars.rds structure and adds
# Month / YearMonth fields for monthly export-tab controls.

library(readr)
library(dplyr)
library(lubridate)

in_path <- file.path("data_raw", "EXPORT_2000-2026_byMonth.csv")
out_dir <- "sortedVars"
out_file <- file.path(out_dir, "sorted_vars_EXPORTS_2026.rds")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(in_path)) {
  stop("Missing export CSV: ", in_path)
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

required_cols <- c("Period", "Province", "Country", "State", "Commodity", "Unit of measure")

meta <- read_trade_csv(in_path)
missing_cols <- setdiff(required_cols, names(meta))
if (length(missing_cols) > 0) {
  stop(
    "Export CSV is missing required columns: ", paste(missing_cols, collapse = ", "),
    "\nAvailable columns: ", paste(names(meta), collapse = ", ")
  )
}

meta <- meta %>%
  dplyr::select(Period, Province, Country, State, Commodity, `Unit of measure`) %>%
  dplyr::filter(!(is.na(Period) | trimws(as.character(Period)) == "")) %>%
  dplyr::mutate(
    Period = parse_period_monthly(Period),
    Year = as.integer(format(Period, "%Y")),
    Month = as.integer(format(Period, "%m")),
    YearMonth = format(Period, "%Y-%m"),
    State = dplyr::na_if(State, "N/A"),
    `Unit of measure` = dplyr::na_if(`Unit of measure`, "N/A")
  )

if (all(is.na(meta$Period))) {
  stop(
    "Period parsing failed. Example Period values: ",
    paste(head(as.character(meta$Period), 5), collapse = ", "),
    ". Expected values like YYYY-mm-dd, m/d/YYYY, YYYY-mm, YYYY/mm, Mon YYYY, or Month YYYY."
  )
}

month_labels <- stats::setNames(month.name, seq_len(12))

sorted_vars_exports <- list(
  years       = sort(unique(meta$Year[!is.na(meta$Year)])),
  months      = sort(unique(meta$Month[!is.na(meta$Month)])),
  month_labels = month_labels[as.character(sort(unique(meta$Month[!is.na(meta$Month)])))],
  year_months = sort(unique(meta$YearMonth[!is.na(meta$YearMonth)])),
  provinces   = sort(unique(meta$Province[!is.na(meta$Province)])),
  countries   = sort(unique(meta$Country[!is.na(meta$Country)])),
  commodities = sort(unique(meta$Commodity[!is.na(meta$Commodity)])),
  states      = sort(unique(meta$State[!is.na(meta$State)])),
  units       = sort(unique(meta$`Unit of measure`[!is.na(meta$`Unit of measure`)]))
)

saveRDS(sorted_vars_exports, out_file)

cat("Wrote: ", out_file, "\n", sep = "")
cat("Years: ", min(sorted_vars_exports$years), "-", max(sorted_vars_exports$years), "\n", sep = "")
cat("Year-months: ", min(sorted_vars_exports$year_months), "-", max(sorted_vars_exports$year_months), "\n", sep = "")
cat("Provinces: ", length(sorted_vars_exports$provinces), "\n", sep = "")
cat("Countries: ", length(sorted_vars_exports$countries), "\n", sep = "")
cat("Commodities: ", length(sorted_vars_exports$commodities), "\n", sep = "")
cat("States: ", length(sorted_vars_exports$states), "\n", sep = "")
cat("Units: ", length(sorted_vars_exports$units), "\n", sep = "")
