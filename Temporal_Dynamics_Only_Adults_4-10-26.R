library(ggplot2)
library(cowplot)
theme_set(theme_cowplot()) # Makes the ggplot look nicer (says Dave)
library(tidyverse)
library(tidyr)
library(patchwork)
library(dplyr)
library(deSolve)
library(caTools) #this package is for AUC calculations
library(sensobol)
library(data.table)
library(foreach)
library(parallel)
library(doParallel)
library(Matrix)


# compile model from C definition
setwd("~/Desktop/Rscripts/Data/Fish/")
try(dyn.unload("~/Desktop/Rscripts/Data/Fish/Fishing_Flooding_Model.so")) # unLoad dll, so for Mac
system("R CMD SHLIB ~/Desktop/Rscripts/Data/Fish/Fishing_Flooding_Model.c")
dyn.load("~/Desktop/Rscripts/Data/Fish/Fishing_Flooding_Model.so") # load dll


parameters = readRDS("simpars.RDA")

#hist(data.frame(fishchains[[4]]$samples)$d_A)

#variances = get_best_fit(fishchains)[[2]] 

parameters["latent_stages"] = 60
parameters["latent_rate"] = 4.3
parameters["lambda"] = 0.032*10
parameters["d_W"] = 0.01 #on average GW in fish live 100 days 
parameters["d_F"] = 0
parameters["convEff"] = 0 #how many fish can you build by eating one adult, temper for nauplii and juveniles (mass of n/mass over a)
parameters["ImmigrationRate"] =  0.05/50 #fish per liter per day 
parameters["FishingRate"] = 0.05 #fish fished per liter per day (fishing effort) - the average fish is caught after 1 week (1/7)
parameters["ImmigrationPeriod"] = 50
parameters["L1day"] =25 #25,150,350

#mat[, "ImmigrationPeriod"] <- qunif(mat[, "ImmigrationPeriod"], 10, 182) #units in days 
#mat[, "ImmigrationRate"] <- qunif(mat[, "ImmigrationRate"], 0.025/50, 0.5/50) #fish per liter per day
#mat[, "FishingRate"] <- qunif(mat[, "FishingRate"], 0, 0.143)



parameters = unlist(parameters)


Exposed_names = paste0("E", 1:parameters["latent_stages"])
Exposed_values = rep(0, times=parameters["latent_stages"])
names(Exposed_values) = Exposed_names
Exposed_values

#E1 <- Exposed_values[1]
#Exposed_values = Exposed_values[2:60]

Inits = c(N = 7500, J = 6000, A = 700, Exposed_values, I = 0, Preds = 0,L3F = 0, HEF = 0)/15
timespan <- seq(0, 1825, by = 1)

#####this new model runs so that introductions occur on day 150 and day 350 both within a single sim####

sim1 = lsoda(y = Inits, times=timespan, parms = parameters, func="compute_derivatives", dllname = "Fishing_Flooding_Model", initfunc="initmod", maxsteps = 1e6)


sim1Es = data.frame(rowSums(sim1[, c(5:64)])) 

colnames(sim1Es) = "Es"

sim1Esrow = sim1Es %>% mutate(time = row_number())

ggplot(data = sim1Esrow, aes(x=time,y=Es)) + geom_line() + xlim(0,1000)


sim1select = data.frame(sim1) %>% select(c(1:4,65:68))

sim1Esrowwide = sim1Esrow %>% mutate(PARAM = "Es")

colnames(sim1Esrowwide) = c("Counts","time","PARAM")

sim1summarizedlong = data.frame(sim1select) %>%
  pivot_longer(
    cols = c(2:8),
    names_to = "PARAM",
    values_to = "Counts"
  )

sim1summarizedlong = bind_rows(sim1summarizedlong,sim1Esrowwide)


# #Add together columns for uninfected stages
# sim1NonInfectedCategory = sim1summarizedlong %>% pivot_wider(names_from = PARAM,values_from = Counts) %>%
#   mutate(NI=N+J+A) %>% select(c("time","NI","Preds","L3F","HEF","I","Es")) %>%
#   pivot_longer(cols = c("NI","Preds","L3F","HEF","I","Es"), names_to = "PARAM",values_to = "Counts")
# 
#multiple HEF, I, L3F, and Preds by 100, counts per 100L 

# sim1rescaled <- sim1summarizedlong %>%
#   mutate(
#     Counts_per_100_L = case_when(
#       PARAM %in% c("HEF", "I") ~ Counts * 100,
#       PARAM %in% c("Preds", "L3F") ~ Counts * 1000,
#       TRUE ~ Counts
#     )
#   )

sim1rescaled <- sim1summarizedlong %>%
  mutate(
    Counts_per_100_L = case_when(
      PARAM %in% c("A") ~ Counts/100,
      PARAM %in% c("Es") ~ Counts/10,
      PARAM %in% c("I") ~ Counts/10,
      PARAM %in% c("Preds") ~ Counts * 100,
      TRUE ~ Counts
    )
  )

sim1rescaledAdultOnly = sim1rescaled %>% filter(PARAM != "N" & PARAM != "J")

ranges <- data.frame(
  xmin = c(1, 366, 731, 1096),
  xmax = c(50, 415, 780, 1145)
)

library(RColorBrewer)
set3_colors <- brewer.pal(12, "Set3")
my_pal <- set3_colors[c(1,4,5,7,6,10)]

#############
day350 = ggplot(data = sim1rescaledAdultOnly,
                aes(x = time, y = Counts_per_100_L,
                    group = PARAM, colour = PARAM)) + geom_rect(
                      data = ranges,
                      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
                      fill = "grey80",
                      alpha = 0.4,
                      inherit.aes = FALSE
                    ) +
  geom_line(size = 1.4, alpha = 0.8) +
  labs(y = "Counts", x = "time") +
  xlim(0,600) + 
  labs(x = "Time (days)",y = "Counts" , color = "Copepod and Fish Stages", title = "Day 350") + theme(base_size = 25) +
  scale_color_manual(
    values = my_pal,
    labels = c(
      expression('Susceptible Adult Copepods, 0.01 L'^-1),
      expression('Exposed Copepods, 0.1 L'^-1),
      expression('Cumulative Fish Exposure, L'^-1),
      expression('L3 Infected Copepods,0.1 L'^-1),
      expression('L3 Infected Fish, L'^-1),
      expression('Fish Predators, 100 L'^-1)
    )
  ) 
  

day350




day150 = ggplot(data = sim1rescaledAdultOnly,
                aes(x = time, y = Counts_per_100_L,
                    group = PARAM, colour = PARAM)) + geom_rect(
                      data = ranges,
                      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
                      fill = "grey80",
                      alpha = 0.4,
                      inherit.aes = FALSE
                    ) +
  geom_line(size = 1.4, alpha = 0.8) +
  labs(y = "Counts", x = "time") +
  xlim(0,600) + 
  labs(x = "Time (days)",y = "Counts" , color = "Copepod and Fish Stages", title = "Day 150") + theme(base_size = 25) +
  scale_color_manual(
    values = my_pal,
    labels = c(
      expression('Susceptible Adult Copepods, 0.01 L'^-1),
      expression('Exposed Copepods, 0.1 L'^-1),
      expression('Cumulative Fish Exposure, L'^-1),
      expression('L3 Infected Copepods, 0.1 L'^-1),
      expression('L3 Infected Fish, L'^-1),
      expression('Fish Predators, 100 L'^-1)
    )
  ) 

day150


day25 = ggplot(data = sim1rescaledAdultOnly,
                aes(x = time, y = Counts_per_100_L,
                    group = PARAM, colour = PARAM)) + geom_rect(
                      data = ranges,
                      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
                      fill = "grey80",
                      alpha = 0.4,
                      inherit.aes = FALSE
                    ) +
  geom_line(size = 1.4, alpha = 0.8) +
  labs(y = "Counts", x = "time") +
  xlim(0,600) +
  labs(x = "Time (days)",y = "Counts" , color = "Copepod and Fish Stages", title = "Day 25") + theme(base_size = 25) +
  scale_color_manual(
    values = my_pal,
    labels = c(
      expression('Susceptible Adult Copepods, 0.01 L'^-1),
      expression('Exposed Copepods, 0.1 L'^-1),
      expression('Cumulative Fish Exposure, L'^-1),
      expression('L3 Infected Copepods, 0.1 L'^-1),
      expression('L3 Infected Fish, L'^-1),
      expression('Fish Predators, 100 L'^-1)
    )
  ) 
ggarrange(day25, day150, day350, nrow = 1, ncol = 3, common.legend = TRUE, legend = "top")

