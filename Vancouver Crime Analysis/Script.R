library(tidyverse)
library(lubridate)
library(ggplot2)
library(patchwork)
library(kableExtra)

# ---- DataCleaning ----

data_crime <- read.csv("Data/crimedata_csv_AllNeighbourhoods_AllYears.csv")
data_pop <- read.csv("Data/VancouverPop.csv")
data_weather <- read.csv("Data/weatherstats_vancouver_daily.csv")

data_weather_relevant <- data_weather %>%
  select(date, avg_temperature, rain) %>%
  mutate(date = as_date(date))

data_crime_with_date <- data_crime %>%
  mutate(date = make_date(YEAR,MONTH,DAY))

merged_data <- data_crime_with_date %>%
  left_join(data_pop, by = c("YEAR" = "Year")) %>%
  left_join(data_weather_relevant, by = "date")

get_dst_start <- function(date) {
  as.Date(sapply(date, function(date) {
    if (as.integer(format(date, "%Y")) < 2007) {
      # Before 2007, DST started on the first Sunday in April
      dst_date <- as.Date(paste0(format(date, "%Y"), "-04-01"))
      dst_date + (7 - as.integer(format(dst_date, "%u")))
    } else {
      # From 2007 onwards, DST starts on the second Sunday in March
      dst_date <- as.Date(paste0(format(date, "%Y"), "-03-01"))
      dst_date + (7 - as.integer(format(dst_date, "%u"))) + 7
    }
  }))
}

merged_data_clean <- merged_data %>%
  group_by(date, Population, avg_temperature, rain) %>%
  summarise(
    num_property_crime = sum(TYPE %in% c("Break and Enter Commercial", "Break and Enter Residential/Other", 
                  "Other Theft", "Theft from Vehicle", "Theft of Bicycle", "Theft of Vehicle")),
    num_violent_crime = sum(TYPE %in% c("Homicide", "Offence Against a Person")),
    .groups = "drop"
    ) %>%
  mutate(num_property_crime_perht = num_property_crime / Population * 100000,
         num_violent_crime_perht = num_violent_crime / Population * 100000,
         dst_start_date = get_dst_start(date),
         days_from_dst = as.integer(date - dst_start_date),
         dst_dummy = as.integer(days_from_dst >= 0),
         day_of_week = as.factor(wday(date))) %>%
  filter(days_from_dst >= -60 & days_from_dst <= 60) %>%
  drop_na()

# ---- Plots ----

avg_crime_by_day <- merged_data_clean %>%
  group_by(days_from_dst, dst_dummy) %>%
  summarise(
    avg_property_crime = mean(num_property_crime_perht),
    avg_violent_crime = mean(num_violent_crime_perht),
    .groups = "drop"
  )

property_crime_plot <- ggplot(avg_crime_by_day, aes(x = days_from_dst, avg_property_crime)) +
  geom_point(alpha = 0.5) +
  geom_smooth(data = filter(avg_crime_by_day, dst_dummy == 0), 
              aes(x = days_from_dst, y = avg_property_crime), se = FALSE, colour = "blue") +
  geom_smooth(data = filter(avg_crime_by_day, dst_dummy == 1), 
              aes(x = days_from_dst, y = avg_property_crime), se = FALSE, colour = "red") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(title = "Property Crimes Around DST Start", x = "Days from DST Start", y = "Average Number of Property Crimes") + 
  theme_minimal()

violent_crime_plot <- ggplot(avg_crime_by_day, aes(x = days_from_dst, avg_violent_crime)) +
  geom_point(alpha = 0.5) +
  geom_smooth(data = filter(avg_crime_by_day, dst_dummy == 0), 
              aes(x = days_from_dst, y = avg_violent_crime), se = FALSE, colour = "blue") +
  geom_smooth(data = filter(avg_crime_by_day, dst_dummy == 1), 
              aes(x = days_from_dst, y = avg_violent_crime), se = FALSE, colour = "red") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(title = "Violent Crimes Around DST Start", x = "Days from DST Start", y = "Average Number of Violent Crimes") + 
  theme_minimal()

combined_plot <- property_crime_plot + violent_crime_plot
print(combined_plot)

run_model <- function(outcome, bandwidth) {
  formula <- as.formula(paste(outcome, "~ days_from_dst*dst_dummy + day_of_week + rain + avg_temperature"))
  lm(formula, data = merged_data_clean, subset = abs(days_from_dst) <= bandwidth)
}

model_property <- run_model("num_property_crime_perht", 60)

model_property_cubic <- lm(
  num_property_crime_perht ~ dst_dummy*(days_from_dst + I(days_from_dst^2) + I(days_from_dst^3)) 
  + day_of_week + rain + avg_temperature, data = merged_data_clean)

model_violent <- lm(num_violent_crime_perht ~ days_from_dst*dst_dummy + day_of_week, data=data_clean, subset = abs(days_from_dst) <= 21)

