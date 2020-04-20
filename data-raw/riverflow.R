## code to prepare `riverflow` dataset goes here
RiverFlow = read.csv('data-raw/river_flow_data.csv', header=TRUE)
usethis::use_data(RiverFlow, overwrite = TRUE)
