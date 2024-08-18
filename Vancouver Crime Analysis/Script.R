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
  mutate(Date = make_date(YEAR,MONTH,DAY),
         dstStartDate = case_when(
           YEAR == 2003 ~ as.Date("2003-04-06"),
           YEAR == 2004 ~ as.Date("2004-04-04"),
           YEAR == 2005 ~ as.Date("2005-04-03"),
           YEAR == 2006 ~ as.Date("2006-04-02"),
           YEAR == 2007 ~ as.Date("2007-03-11"),
           YEAR == 2008 ~ as.Date("2008-03-09"),
           YEAR == 2009 ~ as.Date("2009-03-08"),
           YEAR == 2010 ~ as.Date("2010-03-14"),
           YEAR == 2011 ~ as.Date("2011-03-13"),
           YEAR == 2012 ~ as.Date("2012-03-11"),
           YEAR == 2013 ~ as.Date("2013-03-10"),
           YEAR == 2014 ~ as.Date("2014-03-09"),
           YEAR == 2015 ~ as.Date("2015-03-08"),
           YEAR == 2016 ~ as.Date("2016-03-13"),
           YEAR == 2017 ~ as.Date("2017-03-12"),
           YEAR == 2018 ~ as.Date("2018-03-11"),
           YEAR == 2019 ~ as.Date("2019-03-10"),
           YEAR == 2020 ~ as.Date("2020-03-08"),
           YEAR == 2021 ~ as.Date("2021-03-14"),
           YEAR == 2022 ~ as.Date("2022-03-13"),
           YEAR == 2023 ~ as.Date("2023-03-12")
         ),
         daysFromDST = as.integer(Date - dstStartDate),
         dstDummy = as.integer(daysFromDST >= 0))
