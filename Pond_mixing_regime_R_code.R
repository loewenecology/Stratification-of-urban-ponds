
##### Salinization, warming, and loss of water clarity inhibit mixing of small urban ponds #####

##################
# Clear any variables from the R environment
##################

rm(list=ls())

##################
# Load R packages
##################

#install.packages("devtools")
#devtools::install_github("acgold/compensateR")

library(ggplot2)
library(tidyverse)
library(lubridate)
library(grid) 
library(gridExtra) 
library(ggmap)
library(rLakeAnalyzer)
library(compensateR)
library(interactions) 
library(ggeffects)
library(scales)
library(ggsn)
library(gtable)

##################
# Load data in R environment
##################

setwd("C://Users//charl//OneDrive - Iowa State University//Postdoc//Alliance//profiles//") #set working directory

profiles = read.table("Final_profile_data.csv", header = T, sep = ",") #read in depth profile data
sites = read.table("Final_site_data.csv", header = T, sep = ",") #read in depth profile data

##################
# Create theme for plotting
##################

workingtheme.depth <- theme(strip.background = element_blank(),
                            
                            panel.border = element_rect(colour = "black", fill = NA),
                            panel.grid.major = element_line(colour = "grey"),
                            panel.grid.minor = element_blank(),
                            panel.background = element_rect(colour = "black", fill = "white", linewidth = 0.75),
                            panel.spacing = unit(0.12, "cm"),
                            
                            axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", size = 12),
                            axis.text.y = element_text(colour = "black", size = 12),
                            axis.ticks = element_line(colour = "black", linewidth = 0.4),
                            axis.title.y = element_text(colour = "black", size = 16),
                            axis.title.x = element_text(colour = "black", size = 16))

workingtheme.regressions <- theme(strip.background = element_blank(),
                                  
                                  panel.border = element_rect(colour = "black", fill = NA),
                                  panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank(),
                                  panel.background = element_rect(colour = "black", fill = "white", linewidth = 0.75),
                                  panel.spacing = unit(0.12, "cm"),
                                  
                                  axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", size = 12),
                                  axis.text.y = element_text(colour = "black", size = 12),
                                  axis.ticks = element_line(colour = "black", linewidth = 0.4),
                                  axis.title.y = element_text(colour = "black", size = 16),
                                  axis.title.x = element_text(colour = "black", size = 16),
                                  legend.position = c(0.3, 0.7))

##################
# Drop outlier
##################

profiles <- subset(profiles, Site != "C45")
sites <- subset(sites, Site != "C45")

##################
# Calculate physical parameters
##################

## Summarize site characteristics
# Create variable for maximum depth at each location during each sampling event
max.d.jul <- subset(profiles, Parameter == "Water Temperature" & Event == "July") %>% 
  group_by(Site) %>% 
  summarise(Max.depth = max(Depth, na.rm = T))

max.d.aug <- subset(profiles, Parameter == "Water Temperature" & Event == "August") %>% 
  group_by(Site) %>% 
  summarise(Max.depth = max(Depth, na.rm = T))

max.d.sep <- subset(profiles, Parameter == "Water Temperature" & Event == "September") %>% 
  group_by(Site) %>% 
  summarise(Max.depth = max(Depth, na.rm = T))

# Create variable for surface water temperature at 0.3-m depth at each location during each sampling event
surface.temp.jul <- subset(profiles, Parameter == "Water Temperature" & Event == "July" & Depth == "0.3")  %>%
group_by(Site)

surface.temp.jul <- rename(surface.temp.jul, Surface.temperature = Value)
surface.temp.jul$Surface.temperature <-  as.numeric(surface.temp.jul$Surface.temperature)

surface.temp.aug <- subset(profiles, Parameter == "Water Temperature" & Event == "August" & Depth == "0.3")  %>%
  group_by(Site)

surface.temp.aug <- rename(surface.temp.aug, Surface.temperature = Value)
surface.temp.aug$Surface.temperature <-  as.numeric(surface.temp.aug$Surface.temperature)

surface.temp.sep <- subset(profiles, Parameter == "Water Temperature" & Event == "September" & Depth == "0.3")  %>%
  group_by(Site)

surface.temp.sep <- rename(surface.temp.sep, Surface.temperature = Value)
surface.temp.sep$Surface.temperature <-  as.numeric(surface.temp.sep$Surface.temperature)

# Add variables to seasonal data frames
physical.jul <- left_join(sites, max.d.jul, by = c("Site"))
physical.jul <- left_join(physical.jul, surface.temp.jul, by = c("Site"))

physical.aug <- left_join(sites, max.d.aug, by = c("Site"))
physical.aug <- left_join(physical.aug, surface.temp.aug, by = c("Site"))

physical.sep <- left_join(sites, max.d.sep, by = c("Site"))
physical.sep <- left_join(physical.sep, surface.temp.sep, by = c("Site"))

# Add Secchi depth variables
secchi.jul <- subset(profiles, Parameter == "Secchi" & Event == "July")  %>%
  group_by(Site)

secchi.jul <- secchi.jul %>%
  select(Site, Value)

secchi.jul <- rename(secchi.jul, Secchi.depth = Value)
secchi.jul$Secchi.depth <-  as.numeric(secchi.jul$Secchi.depth)

physical.jul <- left_join(physical.jul, secchi.jul, by = c("Site"))

secchi.aug <- subset(profiles, Parameter == "Secchi" & Event == "August")  %>%
  group_by(Site)

secchi.aug <- secchi.aug %>%
  select(Site, Value)

secchi.aug <- rename(secchi.aug, Secchi.depth = Value)
secchi.aug$Secchi.depth <-  as.numeric(secchi.aug$Secchi.depth)

physical.aug <- left_join(physical.aug, secchi.aug, by = c("Site"))

secchi.sep <- subset(profiles, Parameter == "Secchi" & Event == "September")  %>%
  group_by(Site)

secchi.sep <- secchi.sep %>%
  select(Site, Value)

secchi.sep <- rename(secchi.sep, Secchi.depth = Value)
secchi.sep$Secchi.depth <-  as.numeric(secchi.sep$Secchi.depth)

physical.sep <- left_join(physical.sep, secchi.sep, by = c("Site"))

# Calculate water densities NOT accounting for salinity
wd.a <- function(data) {
  water.density(as.numeric(data$Value))
}

profiles.nested.a <- subset(profiles, Parameter == "Water Temperature") %>%
  group_by(Site) %>%
  nest()

wd.calc.a <- profiles.nested.a %>%
  mutate(fit = map(data, wd.a))

unnested1.a <- unnest(wd.calc.a, data)
unnested2.a <- unnest(wd.calc.a, fit)
unnested1.a$wd.a.value <- unnested2.a$fit

# Calculate salinity from specific conductance measured at each location during each sampling event
sal.data <- subset(profiles, Parameter == "Specific Conductance")
sal.data <- mutate(sal.data, Salinity = salinity(as.numeric(sal.data$Value)))

select.sal.data <- sal.data %>%
  select(Site, Depth, Event, Salinity)

# Calculate water densities accounting for salinity
wd.s <- function(data) {
  water.density(as.numeric(data$Value), data$Salinity)
}

wd.data.s <- left_join(profiles, select.sal.data, by = c("Site", "Event", "Depth"))

profiles.nested.s <- subset(wd.data.s, Parameter == "Water Temperature") %>%
  group_by(Site) %>%
  nest()

wd.s.calc <- profiles.nested.s %>%
  mutate(fit = map(data, wd.s))

unnested1.s <- unnest(wd.s.calc, data)
unnested2.s <- unnest(wd.s.calc, fit)
unnested1.s$wd.s.value <- unnested2.s$fit

# Calculate differences in water densities calculated with and without accounting for salinity
unnested1.s$wd.a.value <- unnested1.a$wd.a.value
unnested1.s$wd.diff.value <- abs(unnested1.s$wd.s.value - unnested1.s$wd.a.value)

# Calculate change in water densities with depth (intervals) NOT accounting for salinity
int.wd.a <- function(data) {
  data$wd.a.value - lag(data$wd.a.value)
}

unnested1.a.nested <- unnested1.a %>%
  group_by(Site, Event) %>%
  nest()

int.wd.a.nested <- unnested1.a.nested %>%
  mutate(wd.a.int.value = map(data, int.wd.a))

int.wd.a.unnested.1 <- unnest(int.wd.a.nested, data)
int.wd.a.unnested.2 <- unnest(int.wd.a.nested, wd.a.int.value)
int.wd.a.unnested.1$wd.a.int.value <- abs(int.wd.a.unnested.2$wd.a.int.value)

# Calculate change in water densities with depth (intervals) accounting for salinity
int.wd.s <- function(data) {
  data$wd.s.value - lag(data$wd.s.value)
}

unnested1.s.nested <- unnested1.s %>%
  group_by(Site, Event) %>%
  nest()

int.wd.s.nested <- unnested1.s.nested %>%
  mutate(wd.s.int.value = map(data, int.wd.s))

int.wd.s.unnested.1 <- unnest(int.wd.s.nested, data)
int.wd.s.unnested.2 <- unnest(int.wd.s.nested, wd.s.int.value)
int.wd.s.unnested.1$wd.s.int.value <- abs(int.wd.s.unnested.2$wd.s.int.value)

# Calculate proportion of water density change with depth (intervals) related to salinity
int.wd.s.unnested.1$wd.a.int.value <- int.wd.a.unnested.1$wd.a.int.value
int.wd.s.unnested.1$wd.prop.int.value <- 1 - (int.wd.s.unnested.1$wd.a.int.value / int.wd.s.unnested.1$wd.s.int.value)
int.wd.s.unnested.1$wd.prop.int.value[int.wd.s.unnested.1$wd.prop.int.value < 0] <- 0

# Calculate mean water density change with depth (intervals) at each location in during July sampling event
int.wd.mean.jul.s <- subset(int.wd.s.unnested.1, Event == "July") %>%
  group_by(Site) %>%
  summarise(int.wd.mean.jul.s = mean(abs(wd.s.int.value), na.rm = TRUE))

# Calculate max water density change with depth (intervals) at each location in during July sampling event
int.wd.max.jul.s <- subset(int.wd.s.unnested.1, Event == "July") %>%
  group_by(Site) %>%
  summarise(int.wd.max.jul.s = max(abs(wd.s.int.value), na.rm = TRUE))

# Calculate mean proportion of water density change with depth (intervals) related to salinity at each location in during July sampling event
int.wd.prop.mean.jul.s <- subset(int.wd.s.unnested.1, Event == "July") %>%
  group_by(Site) %>%
  summarise(int.wd.prop.mean.jul.s = mean(abs(wd.prop.int.value), na.rm = TRUE))

# Calculate minimum proportion of water density change with depth (intervals) related to salinity at each location in during July sampling event
int.wd.prop.min.jul.s <- subset(int.wd.s.unnested.1, Event == "July") %>%
  group_by(Site) %>%
  summarise(int.wd.prop.min.jul.s = min(wd.prop.int.value, na.rm = TRUE))

# Calculate maximum proportion of water density change with depth (intervals) related to salinity at each location in during July sampling event
int.wd.prop.max.jul.s <- subset(int.wd.s.unnested.1, Event == "July") %>%
  group_by(Site) %>%
  summarise(int.wd.prop.max.jul.s = max(wd.prop.int.value, na.rm = TRUE))

# Calculate maximum salinity at each location in during July sampling event
sal.max.jul.s <- subset(int.wd.s.unnested.1, Event == "July") %>%
  group_by(Site) %>%
  summarise(Sal.max = max(abs(Salinity), na.rm = TRUE))

# Add variables to seasonal data frames
physical.jul <- left_join(physical.jul, int.wd.mean.jul.s, by = c("Site"))
physical.jul <- left_join(physical.jul, int.wd.max.jul.s, by = c("Site"))
physical.jul <- left_join(physical.jul, int.wd.prop.mean.jul.s, by = c("Site"))
physical.jul <- left_join(physical.jul, int.wd.prop.min.jul.s, by = c("Site"))
physical.jul <- left_join(physical.jul, int.wd.prop.max.jul.s, by = c("Site"))
physical.jul <- left_join(physical.jul, sal.max.jul.s, by = c("Site"))

# Calculate mean water density change with depth (intervals) at each location in during August sampling event
int.wd.mean.aug.s <- subset(int.wd.s.unnested.1, Event == "August") %>%
  group_by(Site) %>%
  summarise(int.wd.mean.aug.s = mean(abs(wd.s.int.value), na.rm = TRUE))

# Calculate max water density change with depth (intervals) at each location in during August sampling event
int.wd.max.aug.s <- subset(int.wd.s.unnested.1, Event == "August") %>%
  group_by(Site) %>%
  summarise(int.wd.max.aug.s = max(abs(wd.s.int.value), na.rm = TRUE))

# Calculate mean proportion of water density change with depth (intervals) related to salinity at each location in during August sampling event
int.wd.prop.mean.aug.s <- subset(int.wd.s.unnested.1, Event == "August") %>%
  group_by(Site) %>%
  summarise(int.wd.prop.mean.aug.s = mean(abs(wd.prop.int.value), na.rm = TRUE))

# Calculate minimum proportion of water density change with depth (intervals) related to salinity at each location in during August sampling event
int.wd.prop.min.aug.s <- subset(int.wd.s.unnested.1, Event == "August") %>%
  group_by(Site) %>%
  summarise(int.wd.prop.min.aug.s = min(wd.prop.int.value, na.rm = TRUE))

# Calculate maximum proportion of water density change with depth (intervals) related to salinity at each location in during August sampling event
int.wd.prop.max.aug.s <- subset(int.wd.s.unnested.1, Event == "August") %>%
  group_by(Site) %>%
  summarise(int.wd.prop.max.aug.s = max(wd.prop.int.value, na.rm = TRUE))

# Calculate maximum salinity at each location in during August sampling event
sal.max.aug.s <- subset(int.wd.s.unnested.1, Event == "August") %>%
  group_by(Site) %>%
  summarise(Sal.max = max(abs(Salinity), na.rm = TRUE))

# Add variables to seasonal data frames
physical.aug <- left_join(physical.aug, int.wd.mean.aug.s, by = c("Site"))
physical.aug <- left_join(physical.aug, int.wd.max.aug.s, by = c("Site"))
physical.aug <- left_join(physical.aug, int.wd.prop.mean.aug.s, by = c("Site"))
physical.aug <- left_join(physical.aug, int.wd.prop.min.aug.s, by = c("Site"))
physical.aug <- left_join(physical.aug, int.wd.prop.max.aug.s, by = c("Site"))
physical.aug <- left_join(physical.aug, sal.max.aug.s, by = c("Site"))

# Calculate mean water density change with depth (intervals) at each location in during September sampling event
int.wd.mean.sep.s <- subset(int.wd.s.unnested.1, Event == "September") %>%
  group_by(Site) %>%
  summarise(int.wd.mean.sep.s = mean(abs(wd.s.int.value), na.rm = TRUE))

# Calculate max water density change with depth (intervals) at each location in during September sampling event
int.wd.max.sep.s <- subset(int.wd.s.unnested.1, Event == "September") %>%
  group_by(Site) %>%
  summarise(int.wd.max.sep.s = max(abs(wd.s.int.value), na.rm = TRUE))

# Calculate mean proportion of water density change with depth (intervals) related to salinity at each location in during September sampling event
int.wd.prop.mean.sep.s <- subset(int.wd.s.unnested.1, Event == "September") %>%
  group_by(Site) %>%
  summarise(int.wd.prop.mean.sep.s = mean(abs(wd.prop.int.value), na.rm = TRUE))

# Calculate minimum proportion of water density change with depth (intervals) related to salinity at each location in during September sampling event
int.wd.prop.min.sep.s <- subset(int.wd.s.unnested.1, Event == "September") %>%
  group_by(Site) %>%
  summarise(int.wd.prop.min.sep.s = min(wd.prop.int.value, na.rm = TRUE))

# Calculate maximum proportion of water density change with depth (intervals) related to salinity at each location in during September sampling event
int.wd.prop.max.sep.s <- subset(int.wd.s.unnested.1, Event == "September") %>%
  group_by(Site) %>%
  summarise(int.wd.prop.max.sep.s = max(wd.prop.int.value, na.rm = TRUE))

# Calculate maximum salinity at each location in during September sampling event
sal.max.sep.s <- subset(int.wd.s.unnested.1, Event == "September") %>%
  group_by(Site) %>%
  summarise(Sal.max = max(abs(Salinity), na.rm = TRUE))

# Add variables to seasonal data frames
physical.sep <- left_join(physical.sep, int.wd.mean.sep.s, by = c("Site"))
physical.sep <- left_join(physical.sep, int.wd.max.sep.s, by = c("Site"))
physical.sep <- left_join(physical.sep, int.wd.prop.mean.sep.s, by = c("Site"))
physical.sep <- left_join(physical.sep, int.wd.prop.min.sep.s, by = c("Site"))
physical.sep <- left_join(physical.sep, int.wd.prop.max.sep.s, by = c("Site"))
physical.sep <- left_join(physical.sep, sal.max.sep.s, by = c("Site"))

# Calculate number of sites experiencing anoxia in August
nrow(subset(profiles, Parameter == "Dissolved Oxygen" & Units == "mg/L" & Value < 0.5 & Event == "August")  %>%
       group_by(Site) %>%
       summarise(count = sum(Value > 0.5)))

# Summary statistics for site characteristics
sum(int.wd.mean.jul.s$int.wd.mean.jul.s > 0.0287)
sum(int.wd.mean.aug.s$int.wd.mean.aug.s > 0.0287)
sum(int.wd.mean.sep.s$int.wd.mean.sep.s > 0.0287)

min(int.wd.s.unnested.1$wd.prop.int.value, na.rm = TRUE)
max(int.wd.s.unnested.1$wd.prop.int.value, na.rm = TRUE)
mean(int.wd.s.unnested.1$wd.prop.int.value, na.rm = TRUE)

min(physical.jul$Sal.max, na.rm = TRUE)
max(physical.jul$Sal.max, na.rm = TRUE)
mean(physical.jul$Sal.max, na.rm = TRUE)

min(physical.jul$Surface.temperature, na.rm = TRUE)
max(physical.jul$Surface.temperature, na.rm = TRUE)
mean(physical.jul$Surface.temperature, na.rm = TRUE)

min(physical.jul$Secchi.depth, na.rm = TRUE)
max(physical.jul$Secchi.depth, na.rm = TRUE)
mean(physical.jul$Secchi.depth, na.rm = TRUE)

min(physical.jul$Max.depth, na.rm = TRUE)
max(physical.jul$Max.depth, na.rm = TRUE)
mean(physical.jul$Max.depth, na.rm = TRUE)

min(physical.aug$Sal.max, na.rm = TRUE)
max(physical.aug$Sal.max, na.rm = TRUE)
mean(physical.aug$Sal.max, na.rm = TRUE)

min(physical.aug$Surface.temperature, na.rm = TRUE)
max(physical.aug$Surface.temperature, na.rm = TRUE)
mean(physical.aug$Surface.temperature, na.rm = TRUE)

min(physical.aug$Secchi.depth, na.rm = TRUE)
max(physical.aug$Secchi.depth, na.rm = TRUE)
mean(physical.aug$Secchi.depth, na.rm = TRUE)

min(physical.aug$Max.depth, na.rm = TRUE)
max(physical.aug$Max.depth, na.rm = TRUE)
mean(physical.aug$Max.depth, na.rm = TRUE)

min(physical.sep$Sal.max, na.rm = TRUE)
max(physical.sep$Sal.max, na.rm = TRUE)
mean(physical.sep$Sal.max, na.rm = TRUE)

min(physical.sep$Secchi.depth, na.rm = TRUE)
max(physical.sep$Secchi.depth, na.rm = TRUE)
mean(physical.sep$Secchi.depth, na.rm = TRUE)

min(physical.sep$Surface.temperature, na.rm = TRUE)
max(physical.sep$Surface.temperature, na.rm = TRUE)
mean(physical.sep$Surface.temperature, na.rm = TRUE)

min(physical.sep$Max.depth, na.rm = TRUE)
max(physical.sep$Max.depth, na.rm = TRUE)
mean(physical.sep$Max.depth, na.rm = TRUE)

min(sites$Surface.area, na.rm = TRUE)
max(sites$Surface.area, na.rm = TRUE)
mean(sites$Surface.area, na.rm = TRUE)

min(sites$T.Phosphorus, na.rm = TRUE)
max(sites$T.Phosphorus, na.rm = TRUE)
mean(sites$T.Phosphorus, na.rm = TRUE)

min(sites$T.Sodium, na.rm = TRUE)
max(sites$T.Sodium, na.rm = TRUE)
mean(sites$T.Sodium, na.rm = TRUE)

min(sites$D.Chloride, na.rm = TRUE)
max(sites$D.Chloride, na.rm = TRUE)
mean(sites$D.Chloride, na.rm = TRUE)

##################
# Mapping
##################

## Figure 1. requires combining plots and minor graphic edits
# Obtain primary basemap
basemap <- get_stamenmap( bbox = c(left = -79.86, bottom = 43.61, right = -79.64, top = 43.8), zoom = 12, maptype = "terrain")

# Generate map showing salinity and maximum depth across sampling locations
(Salinity_Map <- ggmap(basemap) +
    geom_point(data = physical.aug, alpha = 1, aes(Longitude, Latitude, color = Sal.max, size = Max.depth)) +
    scale_colour_gradient(name = "Max salinity (PSU)", low = "blue4", high = "red") +
    xlab("Longitude") + ylab("Latitude") +
    scalebar(x.min = attr(basemap, "bb")[[2]],
             y.min = attr(basemap, "bb")[[1]],
             x.max = attr(basemap, "bb")[[4]],
             y.max = attr(basemap, "bb")[[3]],
             dist = 2, anchor = c(x = -79.85, y = 43.79),
             transform = T, location = "bottomleft", st.size = 5, st.dist = 0.032, dist_unit = "km") +
    theme(legend.position = "bottom"))

# Generate map showing surface water temperature and maximum depth across sampling locations
(Temperature_Map <- ggmap(basemap) +
    geom_point(data = physical.aug, alpha = 1, aes(Longitude, Latitude, color = Surface.temperature, size = Max.depth)) +
    scale_colour_gradient(name = "Surface temperature (°C)", low = "blue4", high = "red") +
    xlab("Longitude") + ylab("Latitude") +
    scalebar(x.min = attr(basemap, "bb")[[2]],
             y.min = attr(basemap, "bb")[[1]],
             x.max = attr(basemap, "bb")[[4]],
             y.max = attr(basemap, "bb")[[3]],
             dist = 2, anchor = c(x = -79.85, y = 43.79),
             transform = T, location = "bottomleft", st.size = 5, st.dist = 0.032, dist_unit = "km") +
    theme(legend.position = "bottom"))

# Obtain secondary basemap
basemap2 <- get_stamenmap( bbox = c(left = -81.75, bottom = 41.80, right = -77.75, top = 45.69), zoom = 8,
                             maptype = "terrain-background")

# Generate secondary map showing location of study area
(samplingmap2 <- ggmap(basemap2) + 
    geom_point(data = physical.aug, alpha = 1, aes(Longitude, Latitude, color = "brown")) +
    xlab("Longitude") + ylab("Latitude") + labs(color = "Sampling location") +
    scalebar(x.min = attr(basemap, "bb")[[2]],
             y.min = attr(basemap, "bb")[[1]],
             x.max = attr(basemap, "bb")[[4]],
             y.max = attr(basemap, "bb")[[3]],
             dist = 50, anchor = c(x = -81.65, y = 45.58),
             transform = T, location = "bottomleft", st.size = 5, st.dist = 0.5, dist_unit = "km", height = 0.4) +
    theme(legend.position = "bottom"))

##################
# Vertical density  profiles
##################

## Figure 2. requires minor edits with a graphics editor
# Plot vertical density profiles 
profiles$Site = with(profiles, reorder(Site, Depth, FUN = max))
unnested1.s$Site = with(unnested1.s, reorder(Site, Depth, FUN = max))
unnested1.a$Site = with(unnested1.a, reorder(Site, Depth, FUN = max))

# Plot water density profiles NOT taking salinity into account for July
(plot.dens.a.combo.jul <- ggplot(data = subset(unnested1.a, Event == "July"),
                                 aes(x = wd.a.value, colour = Site, y = Depth)) +
    geom_path(alpha = 0.4, linewidth = 4) +
    labs(title = "July",
         x = bquote(paste("Density without salinity (kg/", m^3, ")")),
         y = "Depth from surface (m)") +
    scale_y_reverse(limits = c(4.1,0.1), breaks = c(0, 1, 2, 3, 4)) +
    scale_x_continuous(limits = c(996, 1005), breaks = c(996, 998, 1000, 1002, 1004)) +
    geom_point(alpha = 0.4) +
    scale_colour_discrete(name = "Sites") +
    workingtheme.depth +
    theme(legend.position = "none"))

# Plot water density profiles NOT taking salinity into account for August
(plot.dens.a.combo.aug <- ggplot(data = subset(unnested1.a, Event == "August"),
                                 aes(x = wd.a.value, colour = Site, y = Depth)) +
    geom_path(alpha = 0.4, linewidth = 4) +
    labs(title = "August",
         x = bquote(paste("Density without salinity (kg/", m^3, ")")),
         y = "Depth from surface (m)") +
    scale_y_reverse(limits = c(4.1,0.1), breaks = c(0, 1, 2, 3, 4)) +
    scale_x_continuous(limits = c(996, 1005), breaks = c(996, 998, 1000, 1002, 1004)) +
    geom_point(alpha = 0.4) +
    scale_colour_discrete(name = "Sites") +
    workingtheme.depth +
    theme(legend.position = "none"))

# Plot water density profiles taking salinity into account for July
(plot.dens.s.combo.jul <- ggplot(data = subset(unnested1.s, Event == "July"),
                                 aes(x = wd.s.value, colour = reorder(factor(as.factor(Site)), Value), y = Depth)) +
    geom_path(alpha = 0.4, linewidth = 4) +
    labs(title = "July",
         x = bquote(paste("Density with salinity (kg/", m^3, ")")),
         y = "Depth from surface (m)") +
    scale_y_reverse(limits = c(4.1,0.1), breaks = c(0, 1, 2, 3, 4)) +
    scale_x_continuous(limits = c(996, 1005), breaks = c(996, 998, 1000, 1002, 1004)) +
    geom_point(alpha = 0.4) +
    scale_colour_discrete(name = "Sites") +
    workingtheme.depth +
    theme(legend.position = "none"))

# Plot water density profiles taking salinity into account for August
(plot.dens.s.combo.aug <- ggplot(data = subset(unnested1.s, Event == "August"),
                                 aes(x = wd.s.value, colour = Site, y = Depth)) +
    geom_path(alpha = 0.4, linewidth = 4) +
    labs(title = "August",
         x = bquote(paste("Density with salinity (kg/", m^3, ")")),
         y = "Depth from surface (m)") +
    scale_y_reverse(limits = c(4.1,0.1), breaks = c(0, 1, 2, 3, 4)) +
    scale_x_continuous(limits = c(996, 1005), breaks = c(996, 998, 1000, 1002, 1004)) +
    geom_point(alpha = 0.4) +
    scale_colour_discrete(name = "Sites") +
    workingtheme.depth +
    theme(legend.position = "none"))

# Plot profiles of differences in water densities with and without taking salinity into account for July 
(plot.dens.diff.combo.jul <- ggplot(data = subset(unnested1.s, Event == "July"),
                                    aes(x = wd.diff.value, colour = Site, y = Depth)) +
    geom_path(alpha = 0.4, linewidth = 4) +
    labs(title = "July",
         x = bquote(paste("Contribution of salinity (kg/", m^3, ")")),
         y = "Depth from surface (m)") +
    scale_y_reverse(limits = c(4.1,0.1), breaks = c(0, 1, 2, 3, 4)) +
    scale_x_continuous(limits = c(0, 6.4), breaks = c(0, 2, 4, 6)) +
    geom_point(alpha = 0.4) +
    scale_colour_discrete(name = "Sites") +
    workingtheme.depth +
    theme(legend.position = "none"))

# Plot profiles of differences in water densities with and without taking salinity into account for August 
(plot.dens.diff.combo.aug <- ggplot(data = subset(unnested1.s, Event == "August"),
                                    aes(x = wd.diff.value, colour = Site, y = Depth)) +
    geom_path(alpha = 0.4, linewidth = 4) +
    labs(title = "August",
         x = bquote(paste("Contribution of salinity (kg/", m^3, ")")),
         y = "Depth from surface (m)") +
    scale_y_reverse(limits = c(4.1,0.1), breaks = c(0, 1, 2, 3, 4)) +
    scale_x_continuous(limits = c(0, 6.4), breaks = c(0, 2, 4, 6)) +
    geom_point(alpha = 0.4) +
    scale_colour_discrete(name = "Sites") +
    workingtheme.depth +
    theme(legend.position = "none"))

# Arrange profile plots in vertical orientation
plot.dens.a.combo.jul <- ggplotGrob(plot.dens.a.combo.jul)
plot.dens.a.combo.aug <- ggplotGrob(plot.dens.a.combo.aug)
plot.dens.s.combo.jul <- ggplotGrob(plot.dens.s.combo.jul)
plot.dens.s.combo.aug <- ggplotGrob(plot.dens.s.combo.aug)
plot.dens.diff.combo.jul <- ggplotGrob(plot.dens.diff.combo.jul)
plot.dens.diff.combo.aug <- ggplotGrob(plot.dens.diff.combo.aug)

fig2.c1 <- rbind(plot.dens.a.combo.jul,
                 plot.dens.s.combo.jul,
                 plot.dens.diff.combo.jul)

fig2.c1$widths <- unit.pmax(plot.dens.a.combo.jul$widths,
                            plot.dens.s.combo.jul$widths,
                            plot.dens.diff.combo.jul$widths)

fig2.c2 <- rbind(plot.dens.a.combo.aug,
                 plot.dens.s.combo.aug,
                 plot.dens.diff.combo.aug)

fig2.c2$widths <- unit.pmax(plot.dens.a.combo.aug$widths,
                            plot.dens.s.combo.aug$widths,
                            plot.dens.diff.combo.aug$widths)

grid.arrange(fig2.c1, fig2.c2, ncol = 2)

##################
# Stratification thresholds
##################

## Figure 3. requires combining plots and minor graphic edits
# Plot scatterplots identifying salinity, temperature, and depth thresholds for July
physical.jul <- physical.jul %>%
  mutate(Strat = if_else(int.wd.mean.jul.s > 0.0287, 'Stratified', 'Mixed'))

physical.jul %>% filter(Strat == "Mixed") %>% summarise(Mixed.Sal.max = max(Sal.max))
physical.jul %>% filter(Strat == "Mixed") %>% summarise(Mixed.Temp.max = max(Surface.temperature))
physical.jul %>% filter(Strat == "Mixed") %>% summarise(Mixed.Depth.max = max(Max.depth))

(plot.strat.temp.sal.jul <- ggplot(data = physical.jul) +
    geom_point(aes(x = Surface.temperature, y = Sal.max, color = Strat), size = 4) +
    geom_hline(aes(yintercept = Mixed.Sal.max), alpha = 0.1, linewidth = 3, data = . %>% filter(Strat == "Mixed") %>% summarise(Mixed.Sal.max = max(Sal.max))) +
    geom_vline(aes(xintercept = Mixed.Temp.max), alpha = 0.1, linewidth = 3, data = . %>% filter(Strat == "Mixed") %>% summarise(Mixed.Temp.max = max(Surface.temperature))) +
    labs(title = "July",
         x = "Water temperature (°C)",
         y = "Max Salinity (PSU)") +
    scale_x_continuous(limits = c(21.1, 28.3), breaks = c(21, 22, 23, 24, 25, 26, 27, 28)) +
    scale_y_continuous(limits = c(0, 8.4), breaks = c(0, 2, 4, 6, 8)) +
    workingtheme.regressions)

(plot.strat.depth.sal.jul <- ggplot(data = physical.jul) +
    geom_point(aes(x = Max.depth, y = Sal.max, color = Strat), size = 4) +
    geom_hline(aes(yintercept = Mixed.Sal.max), alpha = 0.1, linewidth = 3, data = . %>% filter(Strat == "Mixed") %>% summarise(Mixed.Sal.max = max(Sal.max))) +
    geom_vline(aes(xintercept = Mixed.Depth.max), alpha = 0.1, linewidth = 3, data = . %>% filter(Strat == "Mixed") %>% summarise(Mixed.Depth.max = max(Max.depth))) +
    labs(x = "Max depth (m)",
         y = "Max Salinity (PSU)") +
    scale_x_continuous(limits = c(0, 4.1), breaks = c(0, 1, 2, 3, 4)) +
    scale_y_continuous(limits = c(0, 8.4), breaks = c(0, 2, 4, 6, 8)) +
    workingtheme.regressions)

(plot.strat.depth.temp.jul <- ggplot(data = physical.jul) +
    geom_point(aes(x = Max.depth, y = Surface.temperature, color = Strat), size = 4) +
    geom_hline(aes(yintercept = Mixed.Temp.max), alpha = 0.1, linewidth = 3, data = . %>% filter(Strat == "Mixed") %>% summarise(Mixed.Temp.max = max(Surface.temperature))) +
    geom_vline(aes(xintercept = Mixed.Depth.max), alpha = 0.1, linewidth = 3, data = . %>% filter(Strat == "Mixed") %>% summarise(Mixed.Depth.max = max(Max.depth))) +
    labs(x = "Max depth (m)",
         y = "Water temperature (°C)") +
    scale_x_continuous(limits = c(0, 4.1), breaks = c(0, 1, 2, 3, 4)) +
    scale_y_continuous(limits = c(21.1, 28.3), breaks = c(21, 22, 23, 24, 25, 26, 27, 28)) +
    workingtheme.regressions)

# Plot scatterplots identifying salinity, temperature, and depth thresholds for August
physical.aug <- physical.aug %>%
  mutate(Strat = if_else(int.wd.mean.aug.s > 0.0287, 'Stratified', 'Mixed'))

physical.aug %>% filter(Strat == "Mixed") %>% summarise(Mixed.Sal.max = max(Sal.max))
physical.aug %>% filter(Strat == "Mixed") %>% summarise(Mixed.Temp.max = max(Surface.temperature))
physical.aug %>% filter(Strat == "Mixed") %>% summarise(Mixed.Depth.max = max(Max.depth))

(plot.strat.temp.sal.aug <- ggplot(data = physical.aug) +
    geom_point(aes(x = Surface.temperature, y = Sal.max, color = Strat), size = 4) +
    geom_hline(aes(yintercept = Mixed.Sal.max), alpha = 0.1, linewidth = 3, data = . %>% filter(Strat == "Mixed") %>% summarise(Mixed.Sal.max = max(Sal.max))) +
    geom_vline(aes(xintercept = Mixed.Temp.max), alpha = 0.1, linewidth = 3, data = . %>% filter(Strat == "Mixed") %>% summarise(Mixed.Temp.max = max(Surface.temperature))) +
    labs(title = "August",
         x = "Water temperature (°C)",
         y = "Max Salinity (PSU)") +
    scale_x_continuous(limits = c(21.1, 28.3), breaks = c(21, 22, 23, 24, 25, 26, 27, 28)) +
    scale_y_continuous(limits = c(0, 8.4), breaks = c(0, 2, 4, 6, 8)) +
    workingtheme.regressions)

(plot.strat.depth.sal.aug <- ggplot(data = physical.aug) +
    geom_point(aes(x = Max.depth, y = Sal.max, color = Strat), size = 4) +
    geom_hline(aes(yintercept = Mixed.Sal.max), alpha = 0.1, linewidth = 3, data = . %>% filter(Strat == "Mixed") %>% summarise(Mixed.Sal.max = max(Sal.max))) +
    geom_vline(aes(xintercept = Mixed.Depth.max), alpha = 0.1, linewidth = 3, data = . %>% filter(Strat == "Mixed") %>% summarise(Mixed.Depth.max = max(Max.depth))) +
    labs(x = "Max depth (m)",
         y = "Max Salinity (PSU)") +
    scale_x_continuous(limits = c(0, 4.1), breaks = c(0, 1, 2, 3, 4)) +
    scale_y_continuous(limits = c(0, 8.4), breaks = c(0, 2, 4, 6, 8)) +
    workingtheme.regressions)

(plot.strat.depth.temp.aug <- ggplot(data = physical.aug) +
    geom_point(aes(x = Max.depth, y = Surface.temperature, color = Strat), size = 4) +
    geom_hline(aes(yintercept = Mixed.Temp.max), alpha = 0.1, linewidth = 3, data = . %>% filter(Strat == "Mixed") %>% summarise(Mixed.Temp.max = max(Surface.temperature))) +
    geom_vline(aes(xintercept = Mixed.Depth.max), alpha = 0.1, linewidth = 3, data = . %>% filter(Strat == "Mixed") %>% summarise(Mixed.Depth.max = max(Max.depth))) +
    labs(x = "Max depth (m)",
         y = "Water temperature (°C)") +
    scale_x_continuous(limits = c(0, 4.1), breaks = c(0, 1, 2, 3, 4)) +
    scale_y_continuous(limits = c(21.1, 28.3), breaks = c(21, 22, 23, 24, 25, 26, 27, 28)) +
    workingtheme.regressions)

# Arrange profile plots for vertical facet
plot.strat.temp.sal.jul <- ggplotGrob(plot.strat.temp.sal.jul)
plot.strat.depth.sal.jul <- ggplotGrob(plot.strat.depth.sal.jul)
plot.strat.depth.temp.jul <- ggplotGrob(plot.strat.depth.temp.jul)
plot.strat.temp.sal.aug <- ggplotGrob(plot.strat.temp.sal.aug)
plot.strat.depth.sal.aug <- ggplotGrob(plot.strat.depth.sal.aug)
plot.strat.depth.temp.aug <- ggplotGrob(plot.strat.depth.temp.aug)

fig3.c1 <- rbind(plot.strat.temp.sal.jul,
                 plot.strat.depth.sal.jul,
                 plot.strat.depth.temp.jul)

fig3.c1$widths <- unit.pmax(plot.strat.depth.temp.jul$widths,
                            plot.strat.depth.sal.jul$widths,
                            plot.strat.depth.temp.jul$widths)

fig3.c2 <- rbind(plot.strat.temp.sal.aug,
                 plot.strat.depth.sal.aug,
                 plot.strat.depth.temp.aug)

fig3.c2$widths <- unit.pmax(plot.strat.temp.sal.aug$widths,
                            plot.strat.depth.sal.aug$widths,
                            plot.strat.depth.temp.aug$widths)

grid.arrange(fig3.c1, fig3.c2, ncol = 2)

##################
# Regression models
##################

## Center/scale predictors
physical.jul$Sal.max.scale <- c(scale(physical.jul$Sal.max, center = TRUE, scale = TRUE))
physical.jul$Max.depth.scale <- c(scale(physical.jul$Max.depth, center = TRUE, scale = TRUE))
physical.jul$Surface.temperature.scale <- c(scale(physical.jul$Surface.temperature, center = TRUE, scale = TRUE))
physical.jul$Surface.area.scale <- c(scale(physical.jul$Surface.area, center = TRUE, scale = TRUE))
physical.jul$Secchi.depth.scale <- c(scale(physical.jul$Secchi.depth, center = TRUE, scale = TRUE))

physical.aug$Sal.max.scale <- c(scale(physical.aug$Sal.max, center = TRUE, scale = TRUE))
physical.aug$Max.depth.scale <- c(scale(physical.aug$Max.depth, center = TRUE, scale = TRUE))
physical.aug$Surface.temperature.scale <- c(scale(physical.aug$Surface.temperature, center = TRUE, scale = TRUE))
physical.aug$Surface.area.scale <- c(scale(physical.aug$Surface.area, center = TRUE, scale = TRUE))
physical.aug$Secchi.depth.scale <- c(scale(physical.aug$Secchi.depth, center = TRUE, scale = TRUE))

physical.sep$Sal.max.scale <- c(scale(physical.sep$Sal.max, center = TRUE, scale = TRUE))
physical.sep$Max.depth.scale <- c(scale(physical.sep$Max.depth, center = TRUE, scale = TRUE))
physical.sep$Surface.temperature.scale <- c(scale(physical.sep$Surface.temperature, center = TRUE, scale = TRUE))
physical.sep$Surface.area.scale <- c(scale(physical.sep$Surface.area, center = TRUE, scale = TRUE))
physical.sep$Secchi.depth.scale <- c(scale(physical.sep$Secchi.depth, center = TRUE, scale = TRUE))

## Examine predictor correlations
corr.jul <- physical.jul %>%
  select(Sal.max.scale, Surface.temperature.scale, Secchi.depth.scale, Max.depth.scale, Surface.area.scale)
round(cor(corr.jul, use ="complete.obs"), 2)

corr.aug <- physical.aug %>%
  select(Sal.max.scale, Surface.temperature.scale, Secchi.depth.scale, Max.depth.scale, Surface.area.scale)
round(cor(corr.aug, use ="complete.obs"), 2)

corr.sep <- physical.sep %>%
  select(Sal.max.scale, Surface.temperature.scale, Secchi.depth.scale, Max.depth.scale, Surface.area.scale)
round(cor(corr.sep, use ="complete.obs"), 2)

## Figure 4. requires minor edits with a graphical editor
# Plot July response metric histogram
hist(physical.jul$int.wd.mean.jul.s)

# Run full July regression model
model.int.wd.mean.jul.s <- glm(int.wd.mean.jul.s ~ Sal.max.scale + Surface.temperature.scale + Secchi.depth.scale + Max.depth.scale + Surface.area.scale +
                                 Sal.max.scale:Surface.temperature.scale + Sal.max.scale:Secchi.depth.scale + Sal.max.scale:Max.depth.scale + Sal.max.scale:Surface.area.scale, data = physical.jul, family = gaussian(link = "log"))
summary(model.int.wd.mean.jul.s)

# Run P-value adjustments
p.adjust(coef(summary(model.int.wd.mean.jul.s))[,4], "fdr")

# Plot July regression models
model.int.wd.mean.jul.s.pred <- ggpredict(model.int.wd.mean.jul.s, terms = "Sal.max.scale [all]")
model.int.wd.mean.jul.s.pred$x <- model.int.wd.mean.jul.s.pred$x * sd(physical.jul$Sal.max) + mean(physical.jul$Sal.max) 

(plot.model.int.wd.mean.sal.jul.s.pred <- ggplot() +
    geom_ribbon(data = model.int.wd.mean.jul.s.pred, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = physical.jul, alpha = 0.8, size = 4, aes(x = Sal.max, y = int.wd.mean.jul.s, color = Surface.temperature)) +
    geom_line(data = model.int.wd.mean.jul.s.pred, aes(x, predicted), linewidth = 1, linetype = "longdash") +
    labs(title = "July",
         y = bquote(paste("Density change (kg/", m^3, ")")),
         x = "Maximum salinity (PSU)",
         color = "Surface temperature (°C)") +
    scale_x_continuous(limits = c(0, 8.4), breaks = c(0, 2, 4, 6, 8)) +
    scale_y_continuous(limits = c(0, 0.65), breaks = c(0, 0.2, 0.4, 0.6)) +
    workingtheme.regressions)

model.int.wd.mean.jul.s.pred <- ggpredict(model.int.wd.mean.jul.s, terms = "Surface.temperature.scale [all]")
model.int.wd.mean.jul.s.pred$x <- model.int.wd.mean.jul.s.pred$x * sd(physical.jul$Surface.temperature) + mean(physical.jul$Surface.temperature) 

(plot.model.int.wd.mean.temp.jul.s.pred <- ggplot() +
    geom_ribbon(data = model.int.wd.mean.jul.s.pred, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = physical.jul, alpha = 0.8, size = 4, aes(x = Surface.temperature, y = int.wd.mean.jul.s, color = Sal.max)) +
    geom_line(data = model.int.wd.mean.jul.s.pred, aes(x, predicted), linewidth = 1, linetype = "longdash") +
    labs(title = "July",
         y = bquote(paste("Density change (kg/", m^3, ")")),
         x = "Surface temperature (°C)",
         color = "Max salinity (PSU)") +
    scale_x_continuous(limits = c(22.5, 28.3), breaks = c(23, 24, 25, 26, 27, 28)) +
    scale_y_continuous(limits = c(0, 0.43), breaks = c(0, 0.1, 0.2, 0.3, 0.4)) +
    workingtheme.regressions)

model.int.wd.mean.jul.s.pred <- ggpredict(model.int.wd.mean.jul.s, terms = "Secchi.depth.scale [all]")
model.int.wd.mean.jul.s.pred$x <- model.int.wd.mean.jul.s.pred$x * sd(physical.jul$Secchi.depth) + mean(physical.jul$Secchi.depth) 

(plot.model.int.wd.mean.secchi.jul.s.pred <- ggplot() +
    geom_ribbon(data = model.int.wd.mean.jul.s.pred, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = physical.jul, alpha = 0.8, size = 4, aes(x = Secchi.depth, y = int.wd.mean.jul.s, color = Sal.max)) +
    geom_line(data = model.int.wd.mean.jul.s.pred, aes(x, predicted), linewidth = 1, linetype = "longdash") +
    labs(title = "July",
         y = bquote(paste("Density change (kg/", m^3, ")")),
         x = "Secchi depth (m)",
         color = "Max salinity (PSU)") +
    scale_x_continuous(limits = c(0.2, 3), breaks = c(0, 1, 2, 3)) +
    scale_y_continuous(limits = c(0, 0.43), breaks = c(0, 0.1, 0.2, 0.3, 0.4)) +
    workingtheme.regressions)

model.int.wd.mean.jul.s.pred <- ggpredict(model.int.wd.mean.jul.s, terms = "Max.depth.scale [all]")
model.int.wd.mean.jul.s.pred$x <- model.int.wd.mean.jul.s.pred$x * sd(physical.jul$Max.depth) + mean(physical.jul$Max.depth) 

(plot.model.int.wd.mean.depth.jul.s.pred <- ggplot() +
    geom_ribbon(data = model.int.wd.mean.jul.s.pred, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = physical.jul, alpha = 0.8, size = 4, aes(x = Max.depth, y = int.wd.mean.jul.s, color = Sal.max)) +
    geom_line(data = model.int.wd.mean.jul.s.pred, aes(x, predicted), linewidth = 1, linetype = "longdash") +
    labs(title = "July",
         y = bquote(paste("Density change (kg/", m^3, ")")),
         x = "Max depth (m)",
         color = "Max salinity (PSU)") +
    scale_x_continuous(limits = c(0.5, 4.1), breaks = c(0, 1, 2, 3, 4)) +
    scale_y_continuous(limits = c(0, 0.43), breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5)) +
    workingtheme.regressions)

model.int.wd.mean.jul.s.pred <- ggpredict(model.int.wd.mean.jul.s, terms = "Surface.area.scale [all]")
model.int.wd.mean.jul.s.pred$x <- model.int.wd.mean.jul.s.pred$x * sd(physical.jul$Surface.area) + mean(physical.jul$Surface.area) 

(plot.model.int.wd.mean.area.jul.s.pred <- ggplot() +
    geom_ribbon(data = model.int.wd.mean.jul.s.pred, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = physical.jul, alpha = 0.8, size = 4, aes(x = Surface.area, y = int.wd.mean.jul.s, color = Sal.max)) +
    geom_line(data = model.int.wd.mean.jul.s.pred, aes(x, predicted), linewidth = 1, linetype = "longdash") +
    labs(title = "July",
         y = bquote(paste("Density change (kg/", m^3, ")")),
         x = bquote(paste("Surface area (", m^2, ")")),
         color = "Max salinity (PSU)") +
    scale_x_continuous(limits = c(0, 24100), breaks = c(0, 6000, 12000, 18000, 24000)) +
    scale_y_continuous(limits = c(0, 0.43), breaks = c(0, 0.1, 0.2, 0.3, 0.4 ,0.5)) +
    workingtheme.regressions)

# Plot August response metric histogram
hist(physical.aug$int.wd.mean.aug.s)

# Run full August regression model
model.int.wd.mean.aug.s <- glm(int.wd.mean.aug.s ~ Sal.max.scale + Surface.temperature.scale + Secchi.depth.scale + Max.depth.scale + Surface.area.scale +
                                 Sal.max.scale:Surface.temperature.scale + Sal.max.scale:Secchi.depth.scale + Sal.max.scale:Max.depth.scale + Sal.max.scale:Surface.area.scale, data = physical.aug, family = gaussian(link = "log"))
summary(model.int.wd.mean.aug.s)

# Run P-value adjustments
p.adjust(coef(summary(model.int.wd.mean.aug.s))[,4], "fdr")

# Plot August regression models
model.int.wd.mean.aug.s.pred <- ggpredict(model.int.wd.mean.aug.s, terms = "Sal.max.scale [all]")
model.int.wd.mean.aug.s.pred$x <- model.int.wd.mean.aug.s.pred$x * sd(physical.aug$Sal.max) + mean(physical.aug$Sal.max) 

(plot.model.int.wd.mean.sal.aug.s.pred <- ggplot() +
    geom_ribbon(data = model.int.wd.mean.aug.s.pred, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = physical.aug, alpha = 0.8, size = 4, aes(x = Sal.max, y = int.wd.mean.aug.s, color = Surface.temperature)) +
    geom_line(data = model.int.wd.mean.aug.s.pred, aes(x, predicted), linewidth = 1, linetype = "longdash") +
    labs(title = "August",
         y = bquote(paste("Density change (kg/", m^3, ")")),
         x = "Maximum salinity (PSU)",
         color = "Surface temperature (°C)") +
    scale_x_continuous(limits = c(0, 7.2), breaks = c(0, 2, 4, 6)) +
    scale_y_continuous(limits = c(0, 0.63), breaks = c(0, 0.1, 0.2, 0.4, 0.6)) +
    workingtheme.regressions)

model.int.wd.mean.aug.s.pred <- ggpredict(model.int.wd.mean.aug.s, terms = "Surface.temperature.scale [all]")
model.int.wd.mean.aug.s.pred$x <- model.int.wd.mean.aug.s.pred$x * sd(physical.aug$Surface.temperature) + mean(physical.aug$Surface.temperature) 

(plot.model.int.wd.mean.temp.aug.s.pred <- ggplot() +
    geom_ribbon(data = model.int.wd.mean.aug.s.pred, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = physical.aug, alpha = 0.8, size = 4, aes(x = Surface.temperature, y = int.wd.mean.aug.s, color = Sal.max)) +
    geom_line(data = model.int.wd.mean.aug.s.pred, aes(x, predicted), linewidth = 1, linetype = "longdash") +
    labs(title = "August",
         y = bquote(paste("Density change (kg/", m^3, ")")),
         x = "Surface temperature (°C)",
         color = "Max salinity (PSU)") +
    scale_x_continuous(limits = c(21, 25.7), breaks = c(21, 22, 23, 24, 25, 26)) +
    scale_y_continuous(limits = c(0, 0.3), breaks = c(0, 0.1, 0.2, 0.3)) +
    workingtheme.regressions)

model.int.wd.mean.aug.s.pred <- ggpredict(model.int.wd.mean.aug.s, terms = "Secchi.depth.scale [all]")
model.int.wd.mean.aug.s.pred$x <- model.int.wd.mean.aug.s.pred$x * sd(physical.aug$Secchi.depth) + mean(physical.aug$Secchi.depth) 

(plot.model.int.wd.mean.secchi.aug.s.pred <- ggplot() +
    geom_ribbon(data = model.int.wd.mean.aug.s.pred, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = physical.aug, alpha = 0.8, size = 4, aes(x = Secchi.depth, y = int.wd.mean.aug.s, color = Sal.max)) +
    geom_line(data = model.int.wd.mean.aug.s.pred, aes(x, predicted), linewidth = 1, linetype = "longdash") +
    labs(title = "August",
         y = bquote(paste("Density change (kg/", m^3, ")")),
         x = "Secchi depth (m)",
         color = "Max salinity (PSU)") +
    scale_x_continuous(limits = c(0, 2), breaks = c(0, 1, 2)) +
    scale_y_continuous(limits = c(0, 0.3), breaks = c(0, 0.1, 0.2, 0.3)) +
    workingtheme.regressions)

model.int.wd.mean.aug.s.pred <- ggpredict(model.int.wd.mean.aug.s, terms = "Max.depth.scale [all]")
model.int.wd.mean.aug.s.pred$x <- model.int.wd.mean.aug.s.pred$x * sd(physical.aug$Max.depth) + mean(physical.aug$Max.depth) 

(plot.model.int.wd.mean.depth.aug.s.pred <- ggplot() +
    geom_ribbon(data = model.int.wd.mean.aug.s.pred, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = physical.aug, alpha = 0.8, size = 4, aes(x = Max.depth, y = int.wd.mean.aug.s, color = Sal.max)) +
    geom_line(data = model.int.wd.mean.aug.s.pred, aes(x, predicted), linewidth = 1, linetype = "longdash") +
    labs(title = "August",
         y = bquote(paste("Density change (kg/", m^3, ")")),
         x = "Max depth (m)",
         color = "Max salinity (PSU)") +
    scale_x_continuous(limits = c(0.4, 4), breaks = c(0, 1, 2, 3, 4)) +
    scale_y_continuous(limits = c(0, 0.3), breaks = c(0, 0.1, 0.2, 0.3)) +
    workingtheme.regressions)

model.int.wd.mean.aug.s.pred <- ggpredict(model.int.wd.mean.aug.s, terms = "Surface.area.scale [all]")
model.int.wd.mean.aug.s.pred$x <- model.int.wd.mean.aug.s.pred$x * sd(physical.aug$Surface.area) + mean(physical.aug$Surface.area) 

(plot.model.int.wd.mean.area.aug.s.pred <- ggplot() +
    geom_ribbon(data = model.int.wd.mean.aug.s.pred, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = physical.aug, alpha = 0.8, size = 4, aes(x = Surface.area, y = int.wd.mean.aug.s, color = Sal.max)) +
    geom_line(data = model.int.wd.mean.aug.s.pred, aes(x, predicted), linewidth = 1, linetype = "longdash") +
    labs(title = "August",
         y = bquote(paste("Density change (kg/", m^3, ")")),
         x = bquote(paste("Surface area (", m^2, ")")),
         color = "Max salinity (PSU)") +
    scale_x_continuous(limits = c(0, 24100), breaks = c(0, 6000, 12000, 18000, 24000)) +
    scale_y_continuous(limits = c(0, 0.3), breaks = c(0.1, 0.2, 0.3, 0.4, 0.5)) +
    workingtheme.regressions)

# Arrange profile plots with vertical orientation
plot.model.int.wd.mean.sal.jul.s.pred <- ggplotGrob(plot.model.int.wd.mean.sal.jul.s.pred)
plot.model.int.wd.mean.sal.aug.s.pred <- ggplotGrob(plot.model.int.wd.mean.sal.aug.s.pred)
plot.model.int.wd.mean.temp.jul.s.pred <- ggplotGrob(plot.model.int.wd.mean.temp.jul.s.pred)
plot.model.int.wd.mean.temp.aug.s.pred <- ggplotGrob(plot.model.int.wd.mean.temp.aug.s.pred)
plot.model.int.wd.mean.depth.jul.s.pred <- ggplotGrob(plot.model.int.wd.mean.depth.jul.s.pred)
plot.model.int.wd.mean.depth.aug.s.pred <- ggplotGrob(plot.model.int.wd.mean.depth.aug.s.pred)
plot.model.int.wd.mean.area.jul.s.pred <- ggplotGrob(plot.model.int.wd.mean.area.jul.s.pred)
plot.model.int.wd.mean.area.aug.s.pred <- ggplotGrob(plot.model.int.wd.mean.area.aug.s.pred)
plot.model.int.wd.mean.secchi.jul.s.pred <- ggplotGrob(plot.model.int.wd.mean.secchi.jul.s.pred)
plot.model.int.wd.mean.secchi.aug.s.pred <- ggplotGrob(plot.model.int.wd.mean.secchi.aug.s.pred)

fig4.c1 <- rbind(plot.model.int.wd.mean.sal.jul.s.pred,
                 plot.model.int.wd.mean.temp.jul.s.pred,
                 plot.model.int.wd.mean.secchi.jul.s.pred,
                 plot.model.int.wd.mean.depth.jul.s.pred,
                 plot.model.int.wd.mean.area.jul.s.pred)

fig4.c1$widths <- unit.pmax(plot.model.int.wd.mean.sal.jul.s.pred$widths,
                            plot.model.int.wd.mean.temp.jul.s.pred$widths,
                            plot.model.int.wd.mean.secchi.jul.s.pred$widths,
                            plot.model.int.wd.mean.depth.jul.s.pred$widths,
                            plot.model.int.wd.mean.area.jul.s.pred$widths)

fig4.c2 <- rbind(plot.model.int.wd.mean.sal.aug.s.pred,
                 plot.model.int.wd.mean.temp.aug.s.pred,
                 plot.model.int.wd.mean.secchi.aug.s.pred,
                 plot.model.int.wd.mean.depth.aug.s.pred,
                 plot.model.int.wd.mean.area.aug.s.pred)

fig4.c2$widths <- unit.pmax(plot.model.int.wd.mean.sal.aug.s.pred$widths,
                            plot.model.int.wd.mean.temp.aug.s.pred$widths,
                            plot.model.int.wd.mean.secchi.aug.s.pred$widths,
                            plot.model.int.wd.mean.depth.aug.s.pred$widths,
                            plot.model.int.wd.mean.area.aug.s.pred$widths)

grid.arrange(fig4.c1, fig4.c2, ncol = 2)

##################
# Supplemental analyses
##################

## Figure S1. requires minor edits with a graphics editor
## Run and plot dissolved chloride regressed on max salinity
model.sal.cl.aug <- glm(Sal.max ~ D.Chloride, data = physical.aug, family = gaussian(link = "identity"))
summary(model.sal.cl.aug)

model.sal.cl.aug.pred <- ggpredict(model.sal.cl.aug, terms = "D.Chloride [all]")

(plot.model.sal.cl.aug.pred <- ggplot() +
    geom_ribbon(data = model.sal.cl.aug.pred, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = physical.aug, alpha = 0.8, size = 4, aes(x = D.Chloride, y = Sal.max, color = T.Sodium)) +
    geom_line(data = model.sal.cl.aug.pred, aes(x, predicted), linewidth = 1, linetype = "longdash") +
    labs(title = "August",
         x = "Dissolved chloride (mg/L)",
         y = "Maximum salinity (PSU)",
         color = "Total Sodium (mg/L)") +
    workingtheme.regressions)

## Figure S2. requires minor edits with a graphics editor
# Plot temperature, salinity, and dissolved oxygen profiles 
profiles$Site = with(profiles, reorder(Site, Depth, FUN = max))
unnested1.s$Site = with(unnested1.s, reorder(Site, Depth, FUN = max))
unnested1.a$Site = with(unnested1.a, reorder(Site, Depth, FUN = max))

# Plot temperature profiles for July
(plot.temp.combo.jul <- ggplot(data = subset(profiles, Parameter == "Water Temperature" & Event == "July"),
                               aes(x = as.numeric(Value), colour = Site, y = Depth)) +
    geom_path(alpha = 0.4, linewidth = 4) +
    labs(title = "July",
         x = "Water temperature (°C)",
         y = "Depth from surface (m)") +
    scale_y_reverse(limits = c(4.1,0.1), breaks = c(0, 1, 2, 3, 4)) +
    scale_x_continuous(limits = c(11.3, 28.7), breaks = c(12, 16, 20, 24, 28)) +
    geom_point(alpha = 0.4) +
    scale_colour_discrete(name = "Sites") +
    workingtheme.depth +
    theme(legend.position = "none"))

# Plot temperature profiles for August
(plot.temp.combo.aug <- ggplot(data = subset(profiles, Parameter == "Water Temperature" & Event == "August"),
                               aes(x = as.numeric(Value), colour = Site, y = Depth)) +
    geom_path(alpha = 0.4, linewidth = 4) +
    labs(title = "August",
         x = "Water temperature (°C)",
         y = "Depth from surface (m)") +
    scale_y_reverse(limits = c(4.1,0.1), breaks = c(0, 1, 2, 3, 4)) +
    scale_x_continuous(limits = c(11.3, 28.7), breaks = c(12, 16, 20, 24, 28)) +
    geom_point(alpha = 0.4) +
    scale_colour_discrete(name = "Sites") +
    workingtheme.depth +
    theme(legend.position = "none"))

# Plot salinity profiles for July
(plot.sal.combo.jul <- ggplot(data = subset(unnested1.s, Event == "July"), 
                              aes(x = Salinity, colour = Site, y = Depth)) +
    geom_path(alpha = 0.4, linewidth = 4) +
    labs(title = "July",
         x = "Salinity (PSU)",
         y = "Depth from surface (m)") +
    scale_y_reverse(limits = c(4.1,0.1), breaks = c(0, 1, 2, 3, 4)) +
    scale_x_continuous(limits = c(0, 8.4), breaks = c(0, 2, 4, 6 ,8)) +
    geom_point(alpha = 0.4) +
    scale_colour_discrete(name = "Sites") +
    workingtheme.depth +
    theme(legend.position = "none"))

# Plot salinity profiles for August
(plot.sal.combo.aug <- ggplot(data = subset(unnested1.s, Event == "August"), 
                              aes(x = Salinity, colour = Site, y = Depth)) +
    geom_path(alpha = 0.4, linewidth = 4) +
    labs(title = "August",
         x = "Salinity (PSU)",
         y = "Depth from surface (m)") +
    scale_y_reverse(limits = c(4.1,0.1), breaks = c(0, 1, 2, 3, 4)) +
    scale_x_continuous(limits = c(0, 8.4), breaks = c(0, 2, 4, 6 ,8)) +
    geom_point(alpha = 0.4) +
    scale_colour_discrete(name = "Sites") +
    workingtheme.depth +
    theme(legend.position = "none"))

# Plot dissolved oxygen concentration profiles for August
(plot.docon.combo.aug <- ggplot(data = subset(profiles, Parameter == "Dissolved Oxygen" & Units == "mg/L"),
                                aes(x = as.numeric(Value), colour = Site, y = Depth)) +
    geom_path(alpha = 0.4, linewidth = 4) +
    labs(title = "August",
         x = "Dissolved oxygen (mg/L)",
         y = "Depth from surface (m)") +
    scale_y_reverse(limits = c(4.1,0.1), breaks = c(0, 1, 2, 3, 4)) +
    scale_x_continuous(limits = c(0, 22), breaks = c(0, 5, 10, 15, 20)) +
    geom_point(alpha = 0.4) +
    scale_colour_discrete(name = "Sites") +
    workingtheme.depth +
    theme(legend.position = "none"))

# Arrange profile plots in vertical orientation
plot.temp.combo.jul <- ggplotGrob(plot.temp.combo.jul)
plot.temp.combo.aug <- ggplotGrob(plot.temp.combo.aug)
plot.sal.combo.jul <- ggplotGrob(plot.sal.combo.jul)
plot.sal.combo.aug <- ggplotGrob(plot.sal.combo.aug)
plot.docon.combo.aug <- ggplotGrob(plot.docon.combo.aug)

# Vertical orientation requires minor edits using a graphics editor
figS2.c1 <- rbind(plot.temp.combo.jul,
                  plot.sal.combo.jul,
                  plot.sal.combo.jul)

figS2.c1$widths <- unit.pmax(plot.temp.combo.jul$widths,
                             plot.sal.combo.jul$widths,
                             plot.sal.combo.jul$widths)

figS2.c2 <- rbind(plot.temp.combo.aug,
                  plot.sal.combo.aug,
                  plot.docon.combo.aug)

figS2.c2$widths <- unit.pmax(plot.temp.combo.aug$widths,
                             plot.sal.combo.aug$widths,
                             plot.docon.combo.aug$widths)

grid.arrange(figS2.c1, figS2.c2, ncol = 2)

## Figure S3-S5. requires minor edits with a graphics editor
# Plot temperature profile facet
(plot.temp.facet <- ggplot(data = subset(profiles, Parameter == "Water Temperature"), 
                           aes(x = as.numeric(Value), y = Depth, colour = factor(as.factor(Event)))) +
    geom_path(alpha = 0.4, linewidth = 2) +
    labs(x = "Water temperature (°C)",
         y = "Depth from surface (m)") +
    scale_y_reverse() +
    scale_colour_discrete(name = "Sampling event", breaks = c("July","August","September")) +
    facet_wrap(~reorder(Site, as.numeric(Depth))) +
    workingtheme.depth)

# Plot salinity profile facet
(plot.sal.facet <- ggplot(data = subset(unnested1.s),
                          aes(x = Salinity, y = Depth, colour = factor(as.factor(Event)))) +
    geom_path(alpha = 0.4, linewidth = 2) +
    labs(x = "Salinity (PSU)",
         y = "Depth from surface (m)") +
    scale_y_reverse() +
    scale_colour_discrete(name = "Sampling event", breaks = c("July","August","September")) +
    facet_wrap(~reorder(Site, as.numeric(Depth))) +
    workingtheme.depth)

# Plot dissolved oxygen concentration profile facet
(plot.docon.facet <- ggplot(data = subset(profiles, Parameter == "Dissolved Oxygen" & Units == "mg/L"),
                            aes(x = as.numeric(Value), y = Depth, colour = factor(as.factor(Event)))) +
    geom_path(alpha = 0.4, linewidth = 2) +
    labs(x = "Dissolved oxygen (mg/L)",
         y = "Depth from surface (m)") +
    scale_y_reverse() +
    scale_colour_discrete(name = "Sampling event", breaks = c("July","August","September")) +
    facet_wrap(~reorder(Site, as.numeric(Depth))) +
    workingtheme.depth)

## Figure S6-S8. requires minor edits with a graphics editor
# Prepare data
profiles.c45 = read.table("Final_profile_data.csv", header = T, sep = ",")
profiles.c45 <- subset(profiles.c45, Site == "C45")

sal.data.c45 <- subset(profiles.c45, Parameter == "Specific Conductance")
sal.data.c45 <- mutate(sal.data.c45, Salinity = salinity(as.numeric(sal.data.c45$Value)))

sal.data.c45 <- sal.data.c45 %>%
  select(Site, Depth, Event, Salinity)

profiles.c45 <- left_join(profiles.c45, sal.data.c45, by = c("Site", "Event", "Depth"))

# Plot temperature profile for outlier location
(plot.temp.c45 <- ggplot(data = subset(profiles.c45, Parameter == "Water Temperature"),
                         aes(x = as.numeric(Value), y = Depth, colour = factor(as.factor(Event)))) +
    geom_path(alpha = 0.4, linewidth = 4) +
    labs(title = "Outlier C45",
         x = "Water temperature (°C)",
         y = "Depth from surface (m)") +
    geom_point(alpha = 0.4) +
    scale_y_reverse() +
    scale_colour_discrete(name = "Sampling event", breaks = c("July","August","September")) +
    workingtheme.depth)

# Plot salinity profile for outlier location
(plot.sal.c45 <- ggplot(data = subset(profiles.c45, Parameter == "Water Temperature"),
                        aes(x = Salinity, y = Depth, colour = factor(as.factor(Event)))) +
    geom_path(alpha = 0.4, linewidth = 4) +
    labs(title = "Outlier C45",
         x = "Salinity (PSU)",
         y = "Depth from surface (m)") +
    scale_y_reverse() +
    geom_point(alpha = 0.4) +
    scale_colour_discrete(name = "Sampling event", breaks = c("July","August","September")) +
    workingtheme.depth)

# Plot dissolved oxygen profile for outlier location
(plot.docon.combo.aug <- ggplot(data = subset(profiles.c45, Parameter == "Dissolved Oxygen" & Units == "mg/L"), 
                                aes(x = as.numeric(Value), y = Depth)) +
    geom_path(alpha = 0.4, linewidth = 4, colour = "#00BA38") +
    labs(title = "Outlier C45 dissolved oxygen profile",
         x = "Dissolved oxygen (mg/L)",
         y = "Depth from surface (m)") +
    scale_y_reverse() +
    geom_point(alpha = 0.4, colour = "#00BA38") +
    scale_colour_discrete(name = "Sites") +
    workingtheme.depth)

## Figure S9-S11. requires minor edits with a graphics editor
# Plot water density profile facet NOT taking salinity into account
(plot.dens.a.facet <- ggplot(data = unnested1.a, aes(x = wd.a.value, y = Depth, colour = factor(as.factor(Event)))) +
    geom_path(alpha = 0.4, linewidth = 2) +
    labs(x = bquote(paste("Water density without salinity (kg/", m^3, ")")),
         y = "Depth from surface (m)") +
    scale_y_reverse() +
    scale_colour_discrete(name = "Event", breaks = c("July","August","September")) +
    facet_wrap(~reorder(Site, as.numeric(Depth))) +
    workingtheme.depth)

# Plot water density profile facets taking salinity into account
(plot.dens.s.facet <- ggplot(data = unnested1.s, aes(x = wd.s.value, y = Depth, colour = factor(as.factor(Event)))) +
    geom_path(alpha = 0.4, linewidth = 2) +
    labs(x = bquote(paste("Water density with salinity (kg/", m^3, ")")),
         y = "Depth from surface (m)") +
    scale_y_reverse() +
    scale_colour_discrete(name = "Event", breaks = c("July","August","September")) +
    facet_wrap(~reorder(Site, as.numeric(Depth))) +
    workingtheme.depth)

# Plot profile facet of differences in water densities with and without taking salinity into account
(plot.dens.diff.facet <- ggplot(data = unnested1.s, aes(x = wd.diff.value, y = Depth, colour = factor(as.factor(Event)))) +
    geom_path(alpha = 0.4, linewidth = 2) +
    labs(x = bquote(paste("Contribution of salinity to water density (kg/", m^3, ")")),
         y = "Depth from surface (m)") +
    scale_y_reverse(limits = c(4.1,0.1), breaks = c(0, 1, 2, 3, 4)) +
    scale_x_continuous(limits = c(0, 6.4), breaks = c(0, 2, 4, 6)) +
    scale_colour_discrete(name = "Event", breaks = c("July","August","September")) +
    facet_wrap(~reorder(Site, as.numeric(Depth))) +
    workingtheme.depth)

## Figure S12. requires minor edits with a graphical editor
# Run full July regression model
model.int.wd.mean.jul.s <- glm(int.wd.mean.jul.s ~ Sal.max.scale + Surface.temperature.scale + Secchi.depth.scale + Max.depth.scale + Surface.area.scale +
                                 Sal.max.scale:Surface.temperature.scale + Sal.max.scale:Secchi.depth.scale + Sal.max.scale:Max.depth.scale + Sal.max.scale:Surface.area.scale, data = physical.jul, family = gaussian(link = "log"))

# Plot July interaction plots
(model.int.wd.mean.jul.s.salxtemp <- interact_plot(model.int.wd.mean.jul.s, pred = Sal.max.scale, modx = Surface.temperature.scale, plot.points = TRUE, 
                                                   point.size = 4, point.alpha = 0.8,
                                                   y = bquote(paste("Density change (kg/", m^3, ")")),
                                                   x = "Max salinity",
                                                   legend.main = "Temperature") + workingtheme.regressions)

(model.int.wd.mean.jul.s.tempxsal <- interact_plot(model.int.wd.mean.jul.s, pred = Surface.temperature.scale, modx = Sal.max.scale, plot.points = TRUE, 
                                                   point.size = 4, point.alpha = 0.8,
                                                   y = bquote(paste("Density change (kg/", m^3, ")")),
                                                   x = "Temperature",
                                                   legend.main = "Max salinity") + workingtheme.regressions)

(model.int.wd.mean.jul.s.secchixsal <- interact_plot(model.int.wd.mean.jul.s, pred = Secchi.depth.scale, modx = Sal.max.scale, plot.points = TRUE, 
                                                     point.size = 4, point.alpha = 0.8,
                                                     y = bquote(paste("Density change (kg/", m^3, ")")),
                                                     x = "Secchi depth",
                                                     legend.main = "Max salinity") + workingtheme.regressions)

(model.int.wd.mean.jul.s.depthxsal <- interact_plot(model.int.wd.mean.jul.s, pred = Max.depth.scale, modx = Sal.max.scale, plot.points = TRUE,
                                                    point.size = 4, point.alpha = 0.8,
                                                    y = bquote(paste("Density change (kg/", m^3, ")")),
                                                    x = "Max depth",
                                                    legend.main = "Max salinity") + workingtheme.regressions)

(model.int.wd.mean.jul.s.areaxsal <- interact_plot(model.int.wd.mean.jul.s, pred = Surface.area.scale, modx = Sal.max.scale, plot.points = TRUE,
                                                   point.size = 4, point.alpha = 0.8,
                                                   y = bquote(paste("Density change (kg/", m^3, ")")),
                                                   x = "Surface area",
                                                   legend.main = "Max salinity") + workingtheme.regressions)

# Run full August regression model
model.int.wd.mean.aug.s <- glm(int.wd.mean.aug.s ~ Sal.max.scale + Surface.temperature.scale + Secchi.depth.scale + Max.depth.scale + Surface.area.scale +
                                 Sal.max.scale:Surface.temperature.scale + Sal.max.scale:Secchi.depth.scale + Sal.max.scale:Max.depth.scale + Sal.max.scale:Surface.area.scale, data = physical.aug, family = gaussian(link = "log"))

# Plot August interaction plots
(model.int.wd.mean.aug.s.salxtemp <- interact_plot(model.int.wd.mean.aug.s, pred = Sal.max.scale, modx = Surface.temperature.scale, plot.points = TRUE, 
                                                   point.size = 4, point.alpha = 0.8,
                                                   y = bquote(paste("Density change (kg/", m^3, ")")),
                                                   x = "Max salinity",
                                                   legend.main = "Temperature") + workingtheme.regressions)

(model.int.wd.mean.aug.s.tempxsal <- interact_plot(model.int.wd.mean.aug.s, pred = Surface.temperature.scale, modx = Sal.max.scale, plot.points = TRUE, 
                                                   point.size = 4, point.alpha = 0.8,
                                                   y = bquote(paste("Density change (kg/", m^3, ")")),
                                                   x = "Temperature",
                                                   legend.main = "Max salinity") + workingtheme.regressions)

(model.int.wd.mean.aug.s.secchixsal <- interact_plot(model.int.wd.mean.aug.s, pred = Secchi.depth.scale, modx = Sal.max.scale, plot.points = TRUE, 
                                                     point.size = 4, point.alpha = 0.8,
                                                     y = bquote(paste("Density change (kg/", m^3, ")")),
                                                     x = "Secchi depth",
                                                     legend.main = "Max salinity") + workingtheme.regressions)

(model.int.wd.mean.aug.s.depthxsal <- interact_plot(model.int.wd.mean.aug.s, pred = Max.depth.scale, modx = Sal.max.scale, plot.points = TRUE, 
                                                    point.size = 4, point.alpha = 0.8,
                                                    y = bquote(paste("Density change (kg/", m^3, ")")),
                                                    x = "Max depth",
                                                    legend.main = "Max salinity") + workingtheme.regressions)

(model.int.wd.mean.aug.s.areaxsal <- interact_plot(model.int.wd.mean.aug.s, pred = Surface.area.scale, modx = Sal.max.scale, plot.points = TRUE, 
                                                   point.size = 4, point.alpha = 0.8,
                                                   y = bquote(paste("Density change (kg/", m^3, ")")),
                                                   x = "Surface area",
                                                   legend.main = "Max salinity") + workingtheme.regressions)

# Arrange profile plots for vertical orientation
model.int.wd.mean.jul.s.salxtemp <- ggplotGrob(model.int.wd.mean.jul.s.salxtemp)
model.int.wd.mean.jul.s.tempxsal <- ggplotGrob(model.int.wd.mean.jul.s.tempxsal)
model.int.wd.mean.jul.s.secchixsal <- ggplotGrob(model.int.wd.mean.jul.s.secchixsal)
model.int.wd.mean.jul.s.depthxsal <- ggplotGrob(model.int.wd.mean.jul.s.depthxsal)
model.int.wd.mean.jul.s.areaxsal <- ggplotGrob(model.int.wd.mean.jul.s.areaxsal)
model.int.wd.mean.aug.s.salxtemp <- ggplotGrob(model.int.wd.mean.aug.s.salxtemp)
model.int.wd.mean.aug.s.tempxsal <- ggplotGrob(model.int.wd.mean.aug.s.tempxsal)
model.int.wd.mean.aug.s.secchixsal <- ggplotGrob(model.int.wd.mean.aug.s.secchixsal)
model.int.wd.mean.aug.s.depthxsal <- ggplotGrob(model.int.wd.mean.aug.s.depthxsal)
model.int.wd.mean.aug.s.areaxsal <- ggplotGrob(model.int.wd.mean.aug.s.areaxsal)

figS12.c1 <- rbind(model.int.wd.mean.jul.s.salxtemp,
                   model.int.wd.mean.jul.s.tempxsal,
                   model.int.wd.mean.jul.s.secchixsal,
                   model.int.wd.mean.jul.s.depthxsal,
                   model.int.wd.mean.jul.s.areaxsal)

figS12.c1$widths <- unit.pmax(model.int.wd.mean.jul.s.salxtemp$widths,
                              model.int.wd.mean.jul.s.tempxsal$widths,
                              model.int.wd.mean.jul.s.secchixsal$widths,
                              model.int.wd.mean.jul.s.depthxsal$widths,
                              model.int.wd.mean.jul.s.areaxsal$widths)

figS12.c2 <- rbind(model.int.wd.mean.aug.s.salxtemp,
                   model.int.wd.mean.aug.s.tempxsal,
                   model.int.wd.mean.aug.s.secchixsal,
                   model.int.wd.mean.aug.s.depthxsal,
                   model.int.wd.mean.aug.s.areaxsal)

figS12.c2$widths <- unit.pmax(model.int.wd.mean.aug.s.salxtemp$widths,
                              model.int.wd.mean.aug.s.tempxsal$widths,
                              model.int.wd.mean.aug.s.secchixsal$widths,
                              model.int.wd.mean.aug.s.depthxsal$widths,
                              model.int.wd.mean.aug.s.areaxsal$widths)

grid.arrange(figS12.c1, figS12.c2, ncol = 2)

###################################################### END