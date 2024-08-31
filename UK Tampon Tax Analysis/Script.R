library(tidyverse)
library(lubridate)
library(ggplot2)
library(fable)
library(feasts)
library(tseries)

# ---- Merge and Clean Raw Data ----

data_list <- list.files("ONS_data/Raw", pattern = "*.csv", full.names = TRUE) %>%
  setdiff("ONS_data/Raw/CPI.csv")

product_data <- data_list %>% 
  map_df(~read_csv(.x) %>% 
           rename_all(toupper))

product_data_clean <- product_data %>%
  select(INDEX_DATE, ITEM_ID, ITEM_DESC, ALL_GM_INDEX) %>%
  mutate(INDEX_DATE = ym(INDEX_DATE)) %>%
  mutate(tax_dummy = as.integer(INDEX_DATE >= as.Date("2021-01-01"))) %>%
  mutate(Month = tsibble::yearmonth(INDEX_DATE)) %>%
  as_tibble(index = Month)

write_csv(product_data_clean, "ONS_data/Processed/merged_product_data_clean.csv")

# ---- Analysis Prep ----

product_data <- read_csv("ONS_data/Processed/merged_product_data_clean.csv") 

create_item_data <- function(item_id) {
  product_data %>%
    filter(ITEM_ID == item_id)
}

rebase_cpi <- function(data){
  january_indices <- data %>%
    filter(month(INDEX_DATE) == 1) %>%
    mutate(year = year(INDEX_DATE),
           jan_index = ALL_GM_INDEX) %>%
    as_tibble() %>%
    select(year, jan_index) %>%
    mutate(
      cumulative_factor = cumprod(jan_index / 100)
    )
  
  data_rebased_cpi <- data %>%
    mutate(year = year(INDEX_DATE)) %>%
    left_join(january_indices, by = "year") %>%
    mutate(
      jan_index = if_else(month(INDEX_DATE) == 1, lag(jan_index), jan_index) %>%
        replace_na(100),
      cumulative_factor = if_else(month(INDEX_DATE) == 1, lag(cumulative_factor), cumulative_factor) %>%
        replace_na(1),
      rebased_index = ALL_GM_INDEX * cumulative_factor
    ) %>%
    select(-year, -jan_index, -cumulative_factor, -ALL_GM_INDEX)
  return(data_rebased_cpi)
}

tampon_data <- create_item_data(520206)
tampon_data_rebased_cpi <- rebase_cpi(tampon_data)

tampon_data_rebased_cpi %>% 
  ACF(tampon_data_rebased_cpi$rebased_index, lag_max = 48) %>%
  autoplot()

tampon_data_rebased_cpi %>%
  gg_season(y = rebased_index)

autoplot(tampon_data_rebased_cpi, rebased_index)+
  geom_vline(xintercept = as.Date("2021-01-01"), color = "red", linetype = "dashed")

test <- tampon_data_rebased_cpi %>%
  model(ARIMA(log(rebased_index) ~ tax_dummy, stepwise = FALSE)) %>%
  report()

test %>% gg_tsresiduals()

tidy(test)