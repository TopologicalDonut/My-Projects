library(tidyverse)
library(lubridate)

data <- read.csv("Data/crimedata_csv_AllNeighbourhoods_AllYears.csv")

dataFiltered <- data %>%
  group_by(YEAR, MONTH, DAY) %>%
  summarise(
    numPropertyCrime = sum(TYPE %in% c("Break and Enter Commercial", "Break and Enter Residential/Other", 
                  "Other Theft", "Theft from Vehicle", "Theft of Bicycle", "Theft of Vehicle")),
    numViolentCrime = sum(TYPE %in% c("Homicide", "Offence Against a Person"))
    ) %>%
  mutate(Date = make_date(YEAR,MONTH,DAY))