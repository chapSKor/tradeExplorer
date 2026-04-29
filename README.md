# PharmaTrace

---

# PharmaTrace

**PharmaTrace** is an interactive R Shiny application for exploring, visualizing, comparing, clustering, and forecasting Canadian pharmaceutical trade flows.

The application was developed to support exploratory analysis of Canadian pharmaceutical supply chains, with emphasis on pharmaceutical imports and exports across provinces, origin/destination countries, U.S. states, commodity classes, units of measurement, trade value, and quantity.

PharmaTrace is designed as both a research tool and a decision-support dashboard for users interested in pharmaceutical trade resilience, supply chain structure, border-policy analysis, and data-driven hypothesis generation.

---

## 🚀 Live Features

- 📊 Interactive trade value histograms
- 🗺️ Geographic heatmaps (Canada & U.S.)
- ⚖️ Distributional comparison metrics
- 📈 Custom trend plots
- 🔮 Forecasting (ARIMA, RF, factor models)
- 🧠 UMAP / t-SNE clustering of trade structure

---

## 📦 Data Coverage

| Dataset  | Resolution | Time Range        |
|----------|------------|------------------|
| Imports  | Annual     | 1988–2025        |
| Exports  | Monthly    | 2000–2026        |

---

## 🧪 Quick Start

```r
setwd("path/to/PharmaTrace")
shiny::runApp()

PharmaTrace is organized into two main dashboard tabs:

### 1. Imports

The **Imports** tab contains the original fully functional PharmaTrace dashboard for Canadian pharmaceutical import records.

Current import data structure:

- Annual Canadian pharmaceutical import records
- Time range: **1988–2025**
- Province-level Canadian import destinations
- Country-level trade origins
- U.S. state-level origin data where available
- Commodity descriptions
- Unit of measurement
- Trade value in Canadian dollars
- Quantity

The Imports tab currently supports:

- Interactive value histograms
- Province, country, and U.S. state filtering
- Commodity filtering
- Year filtering
- Canadian province heatmaps
- U.S. state heatmaps
- Distributional comparisons
- Custom trend plots
- Forecasting
- UMAP and t-SNE clustering maps

---

### 2. Exports

The **Exports** tab mirrors the Imports dashboard structure but is designed for Canadian pharmaceutical export records.

Current export data structure:

- Monthly Canadian pharmaceutical export records
- Time range: **2000-01 to 2026-02**
- Province-level Canadian export origins
- Country-level export destinations
- U.S. state-level destination data where available
- Commodity descriptions
- Unit of measurement
- Trade value in Canadian dollars
- Quantity
- Monthly fields:
  - `Month`
  - `YearMonth`
  - `YearMonthDate`

The Exports tab currently loads from the export Parquet/RDS pipeline and mirrors the Imports dashboard architecture. Some features may still behave in yearly mode until monthly-specific plotting and forecasting logic is fully finalized.

---

## Key Features

### Interactive Histograms

The histogram panel allows users to examine trade-value distributions under flexible filters.

Users can filter by:

- Province
- Commodity
- Year
- Origin/destination country
- U.S. state
- Quantity-zero exclusion

Histogram options include:

- Group traces by province
- Group traces by country
- Group traces by U.S. state
- Log-scale trade value
- Density or count view
- Adjustable bin number
- Adjustable transparency

---

### Geographic Heatmaps

PharmaTrace includes geographic summaries linked to the active histogram filters.

Current map outputs include:

- **Canadian province heatmap**
  - Displays filtered import/export value by province
- **U.S. state heatmap**
  - Displays filtered U.S. state-level trade totals where state-level data are available

These maps help identify regional concentration, trade dependence, and potential geographic vulnerabilities.

---

### Distributional Comparisons

The comparison module allows users to compare trade-value distributions across selected groups.

Supported comparison modes include:

- Province vs province
- Year vs year within a province

Computed comparison metrics include:

- Sample size
- Mean
- Median
- Standard deviation
- Interquartile range
- Percentiles
- Difference in median
- Ratio of medians
- Difference in mean
- Cliff’s delta
- Wasserstein distance
- Kolmogorov–Smirnov statistic
- Kolmogorov–Smirnov p-value
- Jensen–Shannon divergence
- Distributional overlap coefficient

These metrics are intended for exploratory comparison rather than formal causal inference.

---

### Custom Trend Plots

The custom trend module allows users to examine trade-value trends over time.

Users can plot:

- Province-level lines for selected commodities
- Commodity-level lines for selected provinces

Trend summaries are based on:

- Median trade value
- Interquartile range
- Number of records

The trend plot supports:

- Year-range selection
- Province filtering
- Commodity filtering
- Optional quantity-zero exclusion
- Optional log-scale plotting
- Optional five-year trend projection overlay

For export data, monthly plotting support is being incorporated using `YearMonthDate`.

---

### Forecasting

The forecasting module provides exploratory five-year projections based on selected trend-series inputs.

Current forecasting methods include:

- ARIMA per series
- ARIMA on pooled aggregate series
- Correlation/factor-model forecasting
- Random forest lag-feature forecasting

Forecast outputs include:

- Interactive forecast plots
- Point estimates
- Confidence interval-style bands
- Forecast summary table
- Forecast status output

Forecasting is currently structured around annual series and is best developed for the Imports tab. Monthly-aware export forecasting will require additional handling of monthly time steps, seasonality, and forecast horizon definitions.

---

### UMAP and t-SNE Clustering

PharmaTrace includes an interactive clustering module for exploring high-dimensional trade structure.

Supported embedding methods:

- UMAP
- t-SNE

Users can configure:

- Year range
- Province
- Country
- Commodity
- Sample size
- Random seed
- Numeric features
- Categorical features
- PCA dimensions before embedding
- UMAP parameters
- t-SNE parameters

Available plot features include:

- Point size control
- Coloring by province, year, country, commodity, state, or unit
- Optional density contours
- Optional convex hulls by commodity
- Optional medoid-like commodity centroids
- Optional commodity zoom panel
- Stable color mapping across plots

The clustering tool is intended for exploratory visualization and hypothesis generation.

---

## Data Sources and Processing

The application uses preprocessed Parquet files and sorted-variable RDS caches.

### Expected Folder Structure

```text
project-root/
├── app.R
├── appRun.R
├── modules/
│   └── mod_trade_dashboard.R
├── data/
│   ├── all_CADpharma_imports_1988_2025.parquet
│   └── all_CADpharma_exports_2000_2026_byMonth.parquet
├── data_raw/
│   └── EXPORT_2000-2026_byMonth.csv
├── sortedVars/
│   ├── sorted_vars.rds
│   └── sorted_vars_EXPORTS_2026.rds
└── scripts/
    ├── build_parquet.R
    ├── build_parquet_EXPORTS.R
    ├── build_sorted_vars.R
    └── build_sorted_vars_EXPORTS.R
