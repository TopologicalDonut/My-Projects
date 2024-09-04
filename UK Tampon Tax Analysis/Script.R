library(tidyverse)
library(lubridate)
library(ggplot2)
library(tsibble)
library(fable)
library(feasts)
library(tseries)
library(patchwork)
library(kableExtra)

# ---- Merge and Clean Raw Data ----

data_list <- list.files("ONS_data/Raw", pattern = "*.csv", full.names = TRUE) %>%
  setdiff("ONS_data/Raw/CPI.csv")

product_data <- data_list %>% 
  map_df(~read_csv(.x) %>% 
           rename_all(toupper))

product_data_clean <- product_data %>%
  select(INDEX_DATE, ITEM_ID, ITEM_DESC, ALL_GM_INDEX) %>%
  mutate(INDEX_DATE = ym(INDEX_DATE))

write_csv(product_data_clean, "ONS_data/Processed/merged_product_data_clean.csv")

# ---- Table ----

product_data <- read_csv("ONS_data/Processed/merged_product_data_clean.csv", show_col_types = FALSE)

example_data <- product_data %>%
  filter(ITEM_ID == 610310) %>%
  filter(month(INDEX_DATE) %in% c(12, 1, 2, 3) & year(INDEX_DATE) < 2010) %>%
  slice(3:10)

kbl(example_data,
    col.names = (c("Index Date", "Item ID", "Item",
                   "Price Index")),
    align = c('l', 'c', 'c', 'c'),
    booktabs = T,
    linesep = "",
    digits = 2,
    caption = "First 10 Rows of Cleaned Data") %>%
  kable_styling(latex_options = c("striped", "hold_position"))

# ---- Analysis Prep ----

create_item_data <- function(item_id) {
  product_data %>%
    filter(ITEM_ID == item_id) %>%
    mutate(Month = tsibble::yearmonth(INDEX_DATE),
           tax_dummy = as.integer(INDEX_DATE >= as.Date("2021-01-01"))) %>%
    select(-INDEX_DATE) %>%
    as_tsibble(index = Month)
}

rebase_cpi <- function(data){
  january_indices <- data %>%
    filter(month(Month) == 1) %>%
    mutate(year = year(Month),
           jan_index = ALL_GM_INDEX) %>%
    as_tibble() %>%
    select(year, jan_index) %>%
    mutate(
      cumulative_factor = cumprod(jan_index / 100)
    )
  
  data_rebased_cpi <- data %>%
    mutate(year = year(Month)) %>%
    left_join(january_indices, by = "year") %>%
    mutate(
      cumulative_factor = if_else(month(Month) == 1, lag(cumulative_factor), cumulative_factor) %>%
        replace_na(1),
      rebased_index = ALL_GM_INDEX * cumulative_factor,
      log_rebased_index = log(rebased_index)
    ) %>%
    select(-year, -ALL_GM_INDEX, -jan_index, -cumulative_factor)
  return(data_rebased_cpi)
}

analyze_item <- function(item_id) {
  item_data <- create_item_data(item_id)
  item_name <- item_data$ITEM_DESC[1]
  item_data_rebased <- rebase_cpi(item_data)
  
  arima_model <- item_data_rebased %>%
    model(ARIMA(log_rebased_index ~ tax_dummy, stepwise = FALSE))
  
  acf_plot <- augment(arima_model) %>%
    ACF(.resid, lag_max = 24) %>%
    autoplot() +
    ggtitle(paste("ACF of Residuals for", item_name))
  
  resid_vs_fitted <- ggplot(augment(arima_model), aes(x = .fitted, y = .resid)) +
    geom_point(alpha = 0.5) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    labs(title = paste("Residuals vs Fitted for", item_name),
         x = "Fitted Values",
         y = "Residuals")
  
  estimates_info <- tidy(arima_model) %>%
    filter(term == "tax_dummy") %>%
    select(estimate, std.error, p.value)
  
  # Return results
  list(
    item_id = item_id,
    item_name = item_name,
    data = item_data_rebased,
    model = arima_model,
    acf_plot =acf_plot,
    resid_vs_fitted = resid_vs_fitted,
    estimates_info = estimates_info
  )
}

# ---- ACF Plot ----

tampon_data <- create_item_data(520206) %>%
  rebase_cpi()

tampon_data %>% 
  ACF(log_rebased_index, lag_max = 48) %>%
  autoplot()

# ---- Tampon Analysis ----

tampon_analysis <- analyze_item(520206)

kbl(tampon_analysis$estimates_info,
    col.names = c("Estimate", "Std. Error", "p-value"),
    align = c('c', 'c', 'c'),
    booktabs = T,
    linesep = "",
    digits = 4,
    caption = "Estimate of Effect of Tampon Tax Abolition") %>%
  kable_styling(latex_options = c("striped", "hold_position"))

# ---- Tampon Resid Graphs ----

tampon_resid_and_ACF <- tampon_analysis$resid_vs_fitted / tampon_analysis$acf_plot +
  plot_layout(heights = c(1,1))
print(tampon_resid_and_ACF)

# ----

autoplot(tampon_data_rebased_cpi, log_rebased_index) +
  geom_vline(xintercept = as.Date("2021-01-01"), color = "red", linetype = "dashed")

testplots <- test$model %>% gg_tsresiduals()

tidy(test)

# Extract residuals and fitted values
testplotadd <- 
  ggplot(augment(test$model),aes(x = .fitted, y = .resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  labs(title = paste("Residuals vs Fitted Values"),
       x = "Fitted Values",
       y = "Residuals")

item_ids <- c(520213, 430536, 520249, 520241)

results <- lapply(item_ids, analyze_item)

# Access results
results[[1]]$model %>% report()  # Report for first item
results[[1]]$ts_plot  # Time series plot for first item
results[[1]]$residual_plots  # Residual plots for first item
results[[1]]$bp_test  # BP test results for first item

summary_table <- do.call(rbind, lapply(results, function(x) {
  tibble(
    item_id = x$item_id,
    item_name = x$item_name,
    tax_dummy_coef = x$
    tax_dummy_pvalue = x$tax_dummy_pvalue
  )
}))

print(summary_table)

residual_plots_list <- lapply(results, function(x) x$residual_plots)

residual_plots <- flatten(lapply(residual_plots_list, function(x) x[1:2]))

# Arrange these plots
do.call(grid.arrange, c(residual_plots, ncol = 4))