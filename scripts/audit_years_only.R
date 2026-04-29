suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

f <- "data/all_CADpharma_imports_1988_2025.csv"  # change if needed

df <- readr::read_csv(f, skip = 0, show_col_types = FALSE, trim_ws = TRUE)
if (!("Period" %in% names(df))) {
  df <- readr::read_csv(f, skip = 1, show_col_types = FALSE, trim_ws = TRUE)
}

df <- df %>% filter(!(is.na(Period) | trimws(as.character(Period)) == ""))

p <- as.character(df$Period)

# Try ISO then MDY; keep whichever parses more rows
d_iso <- as.Date(p, format = "%Y-%m-%d")
d_mdy <- as.Date(p, format = "%m/%d/%Y")
d <- if (sum(!is.na(d_mdy)) > sum(!is.na(d_iso))) d_mdy else d_iso

df$Period_parsed <- d
df$Year_parsed <- as.integer(format(df$Period_parsed, "%Y"))

cat("Rows:", nrow(df), "\n")
cat("Parsed Period:", sum(!is.na(df$Period_parsed)), "rows\n")

if (all(is.na(df$Year_parsed))) {
  cat("ERROR: Could not parse Period into dates.\n")
  cat("First 10 raw Period values:\n")
  print(head(p, 10))
  stop("Fix Period parsing (MDY vs DMY vs extra text).")
}

cat("Year range:", min(df$Year_parsed, na.rm = TRUE), "to", max(df$Year_parsed, na.rm = TRUE), "\n")
cat("Period range:", min(df$Period_parsed, na.rm = TRUE), "to", max(df$Period_parsed, na.rm = TRUE), "\n")

# Write years to a file
yrs <- sort(unique(df$Year_parsed))
writeLines(as.character(yrs), "audit_unique_years.txt")
cat("Wrote: audit_unique_years.txt (n=", length(yrs), ")\n", sep = "")


f <- "data/all_CADpharma_imports_1988_2025.csv"
df <- readr::read_csv(f, show_col_types = FALSE)
head(df$Period, 3)
tail(df$Period, 3)


