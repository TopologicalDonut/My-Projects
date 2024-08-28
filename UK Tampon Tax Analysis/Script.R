library(tidyverse)
library(lubridate)
library(ggplot2)
library(fable)
library(feasts)
library(tseries)

# ---- Merge Raw Data ----
  
data_list <- list.files("ONS_data/raw_data", pattern = "*.csv", full.names = TRUE) %>%
  setdiff("ONS_data/raw_data/CPI.csv")

# Read all CSV files, standardize column names, and combine
product_data <- data_list %>% 
  map_df(~read_csv(.x) %>% 
           rename_all(toupper))

write_csv(product_data, "ONS_data/merged_product_data.csv")

# ---- Data Cleaning ----

product_data <- read_csv("ONS_data/merged_product_data.csv")

tampon_data <- product_data %>%
  filter(ITEM_ID == 520206) %>%
  mutate(INDEX_DATE = ym(INDEX_DATE)) %>%
  mutate(tax_dummy = as.integer(INDEX_DATE >= as.Date("2021-01-01"))) %>%
  mutate(Month = tsibble::yearmonth(INDEX_DATE)) %>%
  select(INDEX_DATE, Month, ITEM_ID, ITEM_DESC, ALL_GM_INDEX, tax_dummy) %>%
  as_tsibble(index = Month)

oil_data <- product_data %>%
  filter(ITEM_ID == 610310) %>%
  mutate(INDEX_DATE = ym(INDEX_DATE)) %>%
  mutate(tax_dummy = as.integer(INDEX_DATE >= as.Date("2021-01-01"))) %>%
  mutate(Month = tsibble::yearmonth(INDEX_DATE)) %>%
  select(INDEX_DATE, Month, ITEM_ID, ITEM_DESC, ALL_GM_INDEX, tax_dummy) %>%
  as_tsibble(index = Month)

january_indices_oil <- oil_data %>%
  filter(month(INDEX_DATE) == 1) %>%
  mutate(year = year(INDEX_DATE),
         jan_index = ALL_GM_INDEX) %>%
  as_tibble() %>%
  select(year, jan_index) %>%
  mutate(
    cumulative_factor = cumprod(jan_index / 100)
  )

oil_data_rebased_cpi <-oil_data %>%
  mutate(year = year(INDEX_DATE)) %>%
  left_join(january_indices_oil, by = "year") %>%
  mutate(
    jan_index = if_else(month(INDEX_DATE) == 1, lag(jan_index), jan_index) %>%
      replace_na(100),
    cumulative_factor = if_else(month(INDEX_DATE) == 1, lag(cumulative_factor), cumulative_factor) %>%
      replace_na(1),
    rebased_index = ALL_GM_INDEX * cumulative_factor
  )

# Extract January indices
january_indices <- tampon_data %>%
  filter(month(INDEX_DATE) == 1) %>%
  mutate(year = year(INDEX_DATE),
         jan_index = ALL_GM_INDEX) %>%
  as_tibble() %>%
  select(year, jan_index) %>%
  mutate(
    cumulative_factor = cumprod(jan_index / 100)
  )

tampon_data_rebased_cpi <-tampon_data %>%
  mutate(year = year(INDEX_DATE)) %>%
  left_join(january_indices, by = "year") %>%
  mutate(
    jan_index = if_else(month(INDEX_DATE) == 1, lag(jan_index), jan_index) %>%
      replace_na(100),
    cumulative_factor = if_else(month(INDEX_DATE) == 1, lag(cumulative_factor), cumulative_factor) %>%
      replace_na(1),
    rebased_index = ALL_GM_INDEX * cumulative_factor
    )

tampon_data_rebased_cpi %>% 
  ACF(tampon_data_rebased_cpi$rebased_index, lag_max = 48) %>%
  autoplot()

Box.test(tampon_data_rebased_cpi$rebased_index, lag = 10, type = "Ljung-Box")

tampon_data_rebased_cpi %>%
  gg_season(y = rebased_index)

autoplot(tampon_data_rebased_cpi, rebased_index)+
  geom_vline(xintercept = as.Date("2021-01-01"), color = "red", linetype = "dashed")

test <- tampon_data_rebased_cpi %>%
  model(ARIMA(log(rebased_index) ~ tax_dummy, stepwise = FALSE)) %>%
  report()

test %>% gg_tsresiduals()

tidy(test)

bp_result <- bptest(log(rebased_index) ~ INDEX_DATE, data = tampon_data_rebased_cpi)

print(bp_result)

summary(lm(log(rebased_index) ~ INDEX_DATE, data = tampon_data_rebased_cpi))
