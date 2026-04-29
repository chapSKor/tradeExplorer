library(readr)
library(dplyr)

in_path  <- "data/all_CADpharma_imports_1988_2025.csv"
out_dir  <- "sortedVars"
out_file <- file.path(out_dir, "sorted_vars.rds")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Read only the columns needed for UI choices
meta <- readr::read_csv(
  in_path,
  show_col_types = FALSE,
  trim_ws = TRUE,
  col_select = c(Period, Province, Country, State, Commodity, `Unit of measure`)
)

# Parse Year from Period (your same robust logic)
p <- as.character(meta$Period)
d1 <- as.Date(p, format = "%Y-%m-%d")
if (all(is.na(d1))) d1 <- as.Date(p, format = "%m/%d/%Y")
meta$Year <- as.integer(format(d1, "%Y"))

sorted_vars <- list(
  years       = sort(unique(meta$Year[!is.na(meta$Year)])),
  provinces   = sort(unique(meta$Province)),
  countries   = sort(unique(meta$Country)),
  commodities = sort(unique(meta$Commodity)),
  states      = sort(unique(meta$State[!is.na(meta$State)])),
  units       = sort(unique(meta$`Unit of measure`))
)

saveRDS(sorted_vars, out_file)
cat("Wrote:", out_file, "\n")

library(arrow)
pq <- arrow::read_parquet("data/all_CADpharma_imports_1988_2025.parquet", as_data_frame = FALSE)
pq$schema

#############################################
#SCAN data #############################################
#############################################

library(readr)
library(dplyr)
library(stringr)
library(tidyr)

in_path  <- "data/all_CADpharma_imports_1988_2025.csv"
out_dir  <- "technical_validation"
out_csv  <- file.path(out_dir, "province_import_coverage_summary.csv")
out_rds  <- file.path(out_dir, "province_import_coverage_summary.rds")
meta_rds <- file.path(out_dir, "province_import_coverage_metadata.rds")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Helper: parse Year from Period
# -----------------------------
parse_year_from_period <- function(x) {
  x <- as.character(x)
  
  d1 <- as.Date(x, format = "%Y-%m-%d")
  if (all(is.na(d1))) {
    d1 <- as.Date(x, format = "%m/%d/%Y")
  }
  
  yr <- as.integer(format(d1, "%Y"))
  
  # Fallback: extract first 4-digit year if Date parsing fails
  bad <- is.na(yr)
  if (any(bad)) {
    y2 <- stringr::str_extract(x[bad], "(19|20)\\d{2}")
    yr[bad] <- suppressWarnings(as.integer(y2))
  }
  
  yr
}

# --------------------------------------------
# Helper: standardize and screen province field
# --------------------------------------------
clean_province <- function(x) {
  x <- str_squish(as.character(x))
  
  x <- case_when(
    x %in% c("Nfld.", "Newfoundland") ~ "Newfoundland and Labrador",
    x %in% c("P.E.I.", "PEI") ~ "Prince Edward Island",
    x %in% c("Que.", "Québec") ~ "Quebec",
    x %in% c("Alta.") ~ "Alberta",
    x %in% c("B.C.") ~ "British Columbia",
    x %in% c("N.W.T.", "Northwest Terr.") ~ "Northwest Territories",
    TRUE ~ x
  )
  
  x
}

is_valid_province <- function(x) {
  !is.na(x) &
    x != "" &
    !tolower(x) %in% c(
      "unknown", "n/a", "na", "not available", "other", "all provinces"
    )
}

# -----------------------------
# Read only needed columns
# -----------------------------
raw <- readr::read_csv(
  in_path,
  show_col_types = FALSE,
  trim_ws = TRUE,
  col_select = c(Period, Province, Country, State, Commodity, `Unit of measure`)
)

# -----------------------------
# Basic cleaning
# -----------------------------
dat <- raw %>%
  mutate(
    Year = parse_year_from_period(Period),
    Province = clean_province(Province),
    Commodity = str_squish(as.character(Commodity))
  )

# Preserve unscreened province values for audit
all_province_values <- sort(unique(dat$Province))

# ------------------------------------------
# Screen to valid Canadian destination fields
# ------------------------------------------
dat_valid <- dat %>%
  filter(is_valid_province(Province), !is.na(Year), !is.na(Commodity), Commodity != "")

# ------------------------------------------
# Province-level coverage summary
# ------------------------------------------
province_summary <- dat_valid %>%
  group_by(Province) %>%
  summarise(
    n_records = n(),
    n_unique_commodities = n_distinct(Commodity),
    first_year = min(Year, na.rm = TRUE),
    last_year  = max(Year, na.rm = TRUE),
    n_years_with_data = n_distinct(Year),
    year_span = paste0(first_year, "\u2013", last_year),
    .groups = "drop"
  ) %>%
  arrange(Province)

# ------------------------------------------
# National summary objects for manuscript use
# ------------------------------------------
national_summary <- list(
  total_tracked_destination_provinces = n_distinct(dat_valid$Province),
  tracked_destination_provinces = sort(unique(dat_valid$Province)),
  total_unique_commodities_across_all_tracked_provinces = n_distinct(dat_valid$Commodity),
  full_year_range = c(min(dat_valid$Year, na.rm = TRUE), max(dat_valid$Year, na.rm = TRUE))
)

# ------------------------------------------
# Metadata / interpretive note
# ------------------------------------------
summary_metadata <- list(
  note = paste(
    "Province denotes the Canadian destination / province-of-entry field recorded in the source trade data.",
    "These records indicate where pharmaceutical imports arrived into Canada for reporting purposes,",
    "but do not necessarily identify the final end destination of products within Canada.",
    "For example, shipments recorded as entering through Ontario or British Columbia may ultimately",
    "have been redistributed to northern or other downstream jurisdictions."
  ),
  unscreened_unique_province_values = all_province_values,
  screened_in_provinces = sort(unique(dat_valid$Province)),
  screened_out_province_values = setdiff(all_province_values, sort(unique(dat_valid$Province))),
  national_summary = national_summary
)

# ------------------------------------------
# Write outputs
# ------------------------------------------
readr::write_csv(province_summary, out_csv)
saveRDS(province_summary, out_rds)
saveRDS(summary_metadata, meta_rds)

# ------------------------------------------
# Console output
# ------------------------------------------
cat("Wrote province summary CSV: ", out_csv, "\n", sep = "")
cat("Wrote province summary RDS: ", out_rds, "\n", sep = "")
cat("Wrote metadata RDS:         ", meta_rds, "\n\n", sep = "")

cat("Tracked destination provinces:\n")
print(national_summary$tracked_destination_provinces)

cat("\nTotal tracked destination provinces: ",
    national_summary$total_tracked_destination_provinces, "\n", sep = "")
cat("Total unique commodities across all tracked provinces: ",
    national_summary$total_unique_commodities_across_all_tracked_provinces, "\n", sep = "")
cat("Overall year coverage: ",
    national_summary$full_year_range[1], "\u2013", national_summary$full_year_range[2], "\n\n", sep = "")

cat("Province-level summary:\n")
print(province_summary, n = Inf)

cat("\nInterpretive note:\n")
cat(summary_metadata$note, "\n")

