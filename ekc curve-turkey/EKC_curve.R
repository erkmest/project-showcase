library(readr)
library(pastecs)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(zoo)
library(reshape2)
library(urca)
library(tseries)
library(lmtest)
library(sandwich)
library(nortest)
library(car)
library(tibble)

#Read the dataset 
turkey_CO2 <- read_delim("turkey_CO2 dataset.csv", delim = ";")

# Rename the columns' name
turkey_CO2 <- turkey_CO2 %>%
  rename(
    gdp_capita = `GDP_per_capita_2015_constant_$`,
    total_fuel_CO2 = `CO2_emissions_from_fossil_fuels_and_cement production(thousand metric tons of C)`,
    solid_CO2 = `Emissions_from_solid_fuel_consumption`,
    liquid_CO2 = `Emissions_from_liquid_fuel_consumption`,
    gas_CO2 = `Emissions_from_gas_fuel_consumption`,
    cement_CO2 = `Emissions from cement production`,
    CO2_fuel_capita = `CO2_emission_per_capita (metric tons of carbon)`,
    renewable_energy = `renewable_twh`,
    pop_growth = `Population growth (annual %)`
  )

# Removing quotation marks and replace commas
turkey_CO2$`pop_growth` <- gsub('"', '', turkey_CO2$`pop_growth`)  
turkey_CO2$`pop_growth` <- gsub(',', '.', turkey_CO2$`pop_growth`) 
# Convert string to numeric
turkey_CO2$`pop_growth` <- as.numeric(turkey_CO2$`pop_growth`)
# Replace 0 with NA 
turkey_CO2$gas_CO2[turkey_CO2$gas_CO2 == 0] <- NA


turkey_CO2$gas_CO2 <- zoo::na.approx(turkey_CO2$gas_CO2, na.rm = FALSE)


# Fitting total CO2 a log model
turkey_CO2$Log_gdp_capita <- log(turkey_CO2$gdp_capita)
turkey_CO2$Log_CO2_fuel_capita <- log(turkey_CO2$CO2_fuel_capita)

# Plot the observed data log gdp per capita and CO2 per capita
ggplot(turkey_CO2, aes(x = Log_gdp_capita)) +
  geom_point(aes(y = Log_CO2_fuel_capita)) +
  labs(title = "EKC: CO2 Emissions Per Capita vs. Log GDP per Capita", 
       x = "Log GDP per Capita (2015 constant $)", 
       y = "Log Total CO2 Emissions Per Capita (metric tons of carbon)")



#Analysis for EKC curve

CO2_sector <- read_delim("/Users/mesuterk/kuznet curve/EKC Curve for Turkey/CO2_emissions_by_sector.csv", delim = ";")

#Removing empty column
CO2_sector <- CO2_sector[, -ncol(CO2_sector)]

#Rename the column
CO2_sector <- CO2_sector %>%
  rename(
    Year = Category, 
    agriculture = `Agriculture(t)`,
    buildings = `Buildings(t)`,
    fuel_exp = `Fuel Exploitation(t)`,
    industrial_comb = `Industrial Combustion(t)`,
    power = `Power Industry(t)`,
    processes = `Processes(t)`,
    transport = `Transport(t)`,
    waste = `Waste(t)`,
    CO2_sector_capita = `Total CO2/cap`,
    total_sector_CO2 = `Total CO2 Emission`
  )
# Convert character columns to numeric
CO2_sector <- CO2_sector %>%
  mutate(
    agriculture = as.numeric(gsub("\\.", "", agriculture)),
    buildings = as.numeric(gsub("\\.", "", buildings)),
    fuel_exp = as.numeric(gsub("\\.", "", fuel_exp)),
    industrial_comb = as.numeric(gsub("\\.", "", industrial_comb)),
    power = as.numeric(gsub("\\.", "", power)),
    processes = as.numeric(gsub("\\.", "", processes)),
    transport = as.numeric(gsub("\\.", "", transport)),
    total_sector_CO2 = as.numeric(gsub("\\.", "", total_sector_CO2))
  )

CO2_sector$Year <- as.numeric(as.character(CO2_sector$Year))


# Filter CO2_sector
CO2_sector <- subset(CO2_sector, Year <= 2017)

# Merging the datasets on the year
merged_data <- merge(CO2_sector, turkey_CO2, by = "Year")



merged_data <- merged_data %>%
  mutate(Log_total_sector_CO2 = log(total_sector_CO2),
         Log_CO2_sector_capita = log(CO2_sector_capita),
         Log_power = log(power))

colnames(merged_data)

# ADF Test
adf_result_log_gdp_capita <- adf.test(merged_data$Log_gdp_capita, alternative = "stationary")
adf_result_log_CO2_sector_capita <- adf.test(merged_data$Log_CO2_sector_capita, alternative = "stationary")
adf_result_log_power <- adf.test(merged_data$Log_power, alternative = "stationary")
print(adf_result_log_gdp_capita)
print(adf_result_log_CO2_sector_capita)
print(adf_result_log_power)

merged_data$diff_log_gdp_capita <- c(NA, diff(merged_data$Log_gdp_capita))
merged_data$diff_log_CO2_sector_capita <- c(NA, diff(merged_data$Log_CO2_sector_capita))
merged_data$diff_log_power <- c(NA, diff(merged_data$Log_power))

# ADF Test on the Differenced Data
adf_result_diff_log_gdp_capita <- adf.test(na.omit(merged_data$diff_log_gdp_capita), alternative = "stationary")
adf_result_diff_log_CO2_capita <- adf.test(na.omit(merged_data$diff_log_CO2_sector_capita), alternative = "stationary")
adf_result_diff_log_power <- adf.test(na.omit(merged_data$diff_log_power), alternative = "stationary")
print(adf_result_diff_log_gdp_capita)
print(adf_result_diff_log_CO2_capita)
print(adf_result_diff_log_power)

#diff log -log regression
model_ekc <- lm(diff_log_CO2_sector_capita ~ diff_log_gdp_capita + I(diff_log_gdp_capita^2) +
                  diff_log_power,
                data = merged_data)
summary(model_ekc)


# Basic plots
par(mfrow = c(2, 2))
plot(model_ekc)

#autocorrelation of residuals with the Durbin-Watson test
dwtest(model_ekc)

#check for normality of residuals
qqPlot(model_ekc, main="Q-Q Plot for Model Residuals")  # This function is part of the car package

# Check for homoscedasticity
ncvTest(model_ekc)

# Resetting plot settings
par(mfrow = c(1, 1))


# Scatter plot
ggplot(merged_data, aes(x = Log_gdp_capita, y = Log_CO2_sector_capita)) +
  geom_point(alpha = 0.6, color = "blue", size = 3) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), color = "red", size = 1) +
  labs(title = "Environmental Kuznets Curve",
       x = "Log-GDP-per-Capita(2015 constant $)",
       y = "Log-CO2-Emissions-Per-Capita(metric tons)") +
  theme_minimal(base_size =12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10), 
    panel.grid.major = element_line(color = "gray", linetype = "solid",  linewidth = 0.5),
    panel.grid.minor = element_line(color = "gray", linetype = "solid",  linewidth = 0.5)
  )

# Plotting code
ggplot(merged_data, aes(x = gdp_capita, y = Log_CO2_sector_capita)) +
  geom_point(alpha = 0.6, color = "blue", size = 3) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), color = "red", size = 1) +
  labs(title = "Environmental Kuznets Curve Including Power",
       x = "GDP per Capita (2015 constant $)",
       y = "Log CO2 Emissions Per Capita (metric tons)") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        panel.grid.major = element_line(color = "gray", linetype = "solid", linewidth = 0.5),
        panel.grid.minor = element_line(color = "gray", linetype = "solid", linewidth = 0.5)
  )


##Find the turning point!!
model_ekc_point <- lm(Log_CO2_sector_capita ~ gdp_capita + I(gdp_capita^2),
                data = merged_data)
# Get the coefficients from the model
coefficients <- coef(model_ekc_point)

# Calculate the turning point
a <- coefficients[3]  
b <- coefficients[2]
turning_point <- -b / (2 * a)

# Output the turning point
print(turning_point)
###GDP per capita at turning point is $10.379 ###












