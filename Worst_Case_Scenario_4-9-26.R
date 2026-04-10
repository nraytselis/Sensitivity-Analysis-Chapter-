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
library(RColorBrewer)

setwd("~/Desktop/Rscripts/Data/Fish")

# compile model from C definition
setwd("~/Desktop/Rscripts/Data/Fish/")
try(dyn.unload("~/Desktop/Rscripts/Data/Fish/Fishing_Flooding_Model.so")) # unLoad dll, so for Mac
system("R CMD SHLIB ~/Desktop/Rscripts/Data/Fish/Fishing_Flooding_Model.c")
dyn.load("~/Desktop/Rscripts/Data/Fish/Fishing_Flooding_Model.so") # load dll

####Bring in parameters from RDA####
#parameters originally from
# chainsA <- readRDS("Joint_GW_full_A2_2.RDA")
# chainsB <- readRDS("Joint_GW_full_B2_2.RDA")
# chainsC <- readRDS("Joint_GW_full_C2_2.RDA")
parameters = readRDS("simpars.RDA") 
parameters["latent_stages"] = 60
parameters["latent_rate"] = 4.3
parameters["lambda"] = 0.032*10
parameters["d_W"] = 0.01 #on average GW in fish live 100 days 
parameters["d_F"] = 0
parameters["convEff"] = 0 #how many fish can you build by eating one adult, temper for nauplii and juveniles (mass of n/mass over a)
parameters["ImmigrationRate"] =  0.05/50 #fish per liter per day 
parameters["FishingRate"] = 0.05 #fish fished per liter per day (fishing effort) - the average fish is caught after 1 week (1/7)
parameters["ImmigrationPeriod"] = 50
parameters["L1day"] = 25
parameters = unlist(parameters)

####Set up Exposed States####
Exposed_names = paste0("E", 1:parameters["latent_stages"])
Exposed_values = rep(0, times=parameters["latent_stages"])
names(Exposed_values) = Exposed_names
Exposed_values

####Initial Conditions####
Inits = c(N = 7500, J = 6000, A = 700, Exposed_values, I = 0, Preds = 0,L3F = 0, HEF = 0)/15
timespan <- seq(0, 1825, by = 1)

####Run sensitivity analysis for worst case scenario - day 350####
N <- 2^13 #number of parameter combinations (recommendation from paper outlining the sensobol package)
pars <- c( "ImmigrationPeriod", "ImmigrationRate", "FishingRate","d_W")  # Choose parameters of interest
matrices <- c("A", "B", "AB", "BA") #recommendation from paper for specific analysis of interest (variation partitioning approach)
first <- total <- "azzini"
order <- "second" #what level of interactions are allowed 
R<-10^3 
type<-"percent"
conf<-0.95

mat <- sobol_matrices(matrices=matrices, N = N, params = pars, order = "third") #creates different matrices with parameter combinations (based on quantiles)

# Choose their ranges
#you have to pick the uniform because of how the sobol sensitivity is formulated 
mat[, "ImmigrationPeriod"] <- qunif(mat[, "ImmigrationPeriod"], 10, 130) #units in days 
mat[, "ImmigrationRate"] <- qunif(mat[, "ImmigrationRate"], 0.025/50, 0.5/50) #fish per liter per day
mat[, "FishingRate"] <- qunif(mat[, "FishingRate"], 0, 0.143) #fish fished per liter per day (fishing effort) - the average fish is caught between 1 week and the fish are never caught 
parameters["L1day"] = 355 #worst case under parameterization in sensitivity analysis A
mat[, "d_W"] <- qunif(mat[, "d_W"], 0, 0.1) 


sobol_fxn = function(ImmigrationPeriod,ImmigrationRate,FishingRate,d_W) { # Runs the model and calculates area under the prevalence curve (similar to a forloop)
  parameters["ImmigrationPeriod"] = ImmigrationPeriod; parameters["ImmigrationRate"] = ImmigrationRate; parameters["FishingRate"] = FishingRate;parameters["d_W"] = d_W;
  sim = lsoda(y = Inits, times=timespan, parms = parameters, func="compute_derivatives", dllname = "Fishing_Flooding_Model", initfunc="initmod", maxsteps = 500000)
  
  # Extract HEF values at the required days
  hef_731  <- sim[sim[,"time"] == 731,  "HEF"] #first day of 3rd year
  hef_1825 <- sim[sim[,"time"] == 1825, "HEF"] #last day of 5th year 
  
  HEF_diff <- hef_1825 - hef_731
  
  # Select days 731–1825
  sel <- sim[,"time"] >= 731 & sim[,"time"] <= 1825 # other way to do: sel <- 731:1825
  
  #print(sel)
  
  #print(sim[sel, "I"])
  
  
  return(list(c(
    mean_L3F = mean(sim[sel, "L3F"], na.rm = TRUE),
    mean_I   = mean(sim[sel, "I"],   na.rm = TRUE),
    HEF = as.numeric(HEF_diff),
    mean_Preds  = mean(sim[sel, "Preds"], na.rm = TRUE)
  ))) } 



sobol_fxn(ImmigrationPeriod = parameters["ImmigrationPeriod"],ImmigrationRate = parameters["ImmigrationRate"], FishingRate = parameters["FishingRate"], d_W = parameters["d_W"])

sobol_mapply_fxn = function(dt){
  return(mapply(sobol_fxn, dt[, 1], dt[, 2], dt[, 3], dt[ ,4])) } #each column in the dataframe is a parameter vector 



# mat = our inputs
# result = our model outputs
result = sobol_mapply_fxn(mat)
df = as.data.frame(do.call(rbind, result))
full.dt = data.table(mat, df) 

plot_uncertainty(Y = full.dt$mean_L3F, N = N) + labs(y = "Mean L3F", x = "y")
plot_uncertainty(Y = full.dt$mean_I, N = N) + labs(y = "Mean I", x = "y")
plot_uncertainty(Y = full.dt$HEF, N = N) + labs(y = "HEF", x = "y")
plot_uncertainty(Y = full.dt$mean_Preds, N = N) + labs(y = "Preds", x = "y")



indL3F350 <- sobol_indices(matrices = matrices, Y = df$mean_L3F, N = N, params = pars, 
                           first = first, total = total, order = "third", boot = TRUE, R = R,
                           parallel = "no", type = type, conf = conf)

indL3F350

indI350 <- sobol_indices(matrices = matrices, Y = df$mean_I, N = N, params = pars, 
                         first = first, total = total, order = "third", boot = TRUE, R = R,
                         parallel = "no", type = type, conf = conf)

indI350

indHEF350 <- sobol_indices(matrices = matrices, Y = df$HEF, N = N, params = pars, 
                           first = first, total = total, order = "third", boot = TRUE, R = R,
                           parallel = "no", type = type, conf = conf)

indHEF350


indPreds350 <- sobol_indices(matrices = matrices, Y = df$mean_Preds, N = N, params = pars, 
                             first = first, total = total, order ="third", boot = TRUE, R = R,
                             parallel = "no", type = type, conf = conf)

indPreds350


toplot.L3F350 <- subset(indL3F350$results, parameters=="ImmigrationPeriod" | parameters=="ImmigrationRate" | parameters=="FishingRate"| parameters== "d_W")
toplot.I350 <- subset(indI350$results, parameters=="ImmigrationPeriod" | parameters=="ImmigrationRate" | parameters=="FishingRate"| parameters== "d_W")
toplot.HEF350 <- subset(indHEF350$results, parameters=="ImmigrationPeriod" | parameters=="ImmigrationRate" | parameters=="FishingRate"| parameters== "d_W")
toplot.Preds350 <- subset(indPreds350$results, parameters=="ImmigrationPeriod" | parameters=="ImmigrationRate" | parameters=="FishingRate"| parameters== "d_W")

plot(indI350)
plot(indHEF350)


library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(RColorBrewer)
all_set3_colors <- brewer.pal(12, "Set3")
my_selected_colors2 <- all_set3_colors[c(4,6, 5)]


make_sobol_plot <- function(ind_obj,outcome_label = "X", title_text = NULL) {
  
  df_plot <- ind_obj$results %>%
    mutate(variable = str_split(parameters, "\\.")) %>%
    unnest_wider(variable, names_sep = "") %>%
    
    pivot_longer(starts_with("variable"),
                 values_to = "param",
                 values_drop_na = TRUE) %>%
    filter(param != "NA") %>%
    select(param, sensitivity, parameters, original) %>%
    
    mutate(original = as.numeric(original)) %>%
    
    group_by(param, sensitivity, parameters) %>%
    summarise(value = sum(original), .groups = "drop") %>%
    mutate(parameters = gsub("\\.", " x ", parameters))
  
  df_plot = df_plot %>% filter(sensitivity != "Ti")
  
  
  ggplot(df_plot, aes(param, value, fill =  fct_rev(sensitivity))) +
    geom_col() +
    theme_classic(base_size=20) +scale_fill_manual(values = my_selected_colors2, name = "Sensitivity Index") +
    ylim(0,1) + labs(x="Parameter", y = paste("Variance in", outcome_label,"explained by parameter")) +
    scale_x_discrete(labels = c("Parasite Death Rate","Fishing Rate", "Immigration Period", "Immigration Rate")) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

p_L3F350   <- make_sobol_plot(indL3F350,  "L3F")
p_I350     <- make_sobol_plot(indI350,  outcome_label = "infected copepods",  "I")
p_HEF350   <- make_sobol_plot(indHEF350, outcome_label = "exposure via fish", "HEF")
p_Preds350 <- make_sobol_plot(indPreds350,"Preds")




####FOR PUBLICATION 
library(ggpubr)
ggarrange(p_I350,
          p_HEF350,
          nrow=1,
          ncol=2,
          common.legend=TRUE, labels = c("Infected Copepods","Exposure via Fish"),label.x = 0.05, 
          label.y = 1) 



####Run simulation with worst case scenario (day 355) across parasite death rate and immigration period####

parameters["latent_stages"] = 60
parameters["latent_rate"] = 4.3
parameters["lambda"] = 0.032*10
parameters["d_W"] = 0.01#0.1 #on average GW in fish live 100 days 
parameters["d_F"] = 0
parameters["convEff"] = 0 #how many fish can you build by eating one adult, temper for nauplii and juveniles (mass of n/mass over a)
parameters["ImmigrationRate"] =  0.05/50 #fish per liter per day 
parameters["FishingRate"] = 0.05 #fish fished per liter per day (fishing effort) - the average fish is caught after 1 week (1/7)
parameters["ImmigrationPeriod"] =  50
parameters["L1day"] = 355


Inits = c(N = 7500, J = 6000, A = 700, Exposed_values, I = 0, Preds = 0,L3F = 0, HEF = 0)/15
timespan <- seq(0, 1825, by = 1)
ImmigrationPeriod = seq(0,130,by=1) 

result_vector = numeric()
mean_I = numeric()

for(i in seq_along(ImmigrationPeriod)){
  parameters["ImmigrationPeriod"] <- ImmigrationPeriod[i]
  
  sim = lsoda(y = Inits, times=timespan, parms = parameters, func="compute_derivatives", dllname = "Fishing_Flooding_Model", initfunc="initmod", maxsteps = 1e6)
  
  result_vector[i] <- sim[4*365,"HEF"] - sim[3*365,"HEF"] #dim(sim)[1] last row of object sim 
  mean_I[i] <- mean(sim[(3*365):(4*365),"I"])
}

allresults = data.frame(ImmigrationPeriod,HEF = result_vector,I = mean_I)
#allresultslong = allresults %>% pivot_longer(cols = c(2:3),names_to = "State", values_to = "Value")


####Check that it's working correctly###
# parameters["ImmigrationPeriod"] = 50
# parameters["d_W"] = 0.01
# parameters["L1day"] = 350
# sim = lsoda(y = Inits, times=timespan, parms = parameters, func="compute_derivatives", dllname = "Fishing_Flooding_Model", initfunc="initmod", maxsteps = 1e6)
# 
# ggplot(data = sim, aes(x=time,y=HEF)) + geom_line()
# ggplot(data = sim, aes(x=time,y=A)) + geom_line()
# ggplot(data = sim, aes(x=time,y=Preds)) + geom_line()

####Plotting for Publication####
set3_colors <- brewer.pal(12, "Set3")
my_pal <- set3_colors[c(4,5)]


#plot with Infectious copes, HEF, d_W = 0.01
coeff <- 0.33

p1 = ggplot(allresults, aes(x = ImmigrationPeriod)) +
  geom_line(aes(y = HEF, color = "HEF"), linewidth = 1.2) +
  geom_line(aes(y = I / coeff, color = "I"), linewidth = 1.2) +
  scale_y_continuous(
    name = expression("Average Cumulative Fish Exposure ("*L^{-1}*")"),
    sec.axis = sec_axis(~ . * coeff,
                        name = expression("Average L3 Infected Copepods ("*L^{-1}*")"))
  ) +
  scale_color_manual(
    values = c("HEF" = my_pal[1], "I" = my_pal[2]),
    name = NULL
  ) +
  labs(
    x = expression("Immigration Period")
  ) +
  theme_classic(base_size = 20) +
  theme(
    legend.position = c(0.7, 0.8),
    axis.title.y = element_text(color = my_pal[1]),
    axis.title.y.right = element_text(color = my_pal[2]),
    axis.text.y = element_text(color = my_pal[1]),
    axis.text.y.right = element_text(color = my_pal[2])
  ) + theme(legend.position = "none")

p1


#parasite death rate sim
ImmigrationPeriod = 50
d_W = seq(0,0.1,by=0.001) 

result_vectordW = numeric()
mean_IdW = numeric()

for(i in seq_along(d_W)){
  parameters["d_W"] <- d_W[i]
  
  sim = lsoda(y = Inits, times=timespan, parms = parameters, func="compute_derivatives", dllname = "Fishing_Flooding_Model", initfunc="initmod", maxsteps = 1e6)
  
  result_vectordW[i] <- sim[4*365,"HEF"] - sim[3*365,"HEF"] #dim(sim)[1] last row of object sim 
  mean_IdW[i] <- mean(sim[(3*365):(4*365),"I"])
}

allresultsdW = data.frame(d_W,HEF = result_vectordW,I = mean_IdW)


coeff <- 0.33

p2 = ggplot(allresultsdW, aes(x = d_W)) +
  geom_line(aes(y = HEF, color = "HEF"), linewidth = 1.2) +
  geom_line(aes(y = I / coeff, color = "I"), linewidth = 1.2) +
  scale_y_continuous(
    name = expression("Average Cumulative Fish Exposure ("*L^{-1}*")"),
    sec.axis = sec_axis(~ . * coeff,
                        name = expression("Average L3 Infected Copepods ("*L^{-1}*")"))
  ) +
  scale_color_manual(
    values = c("HEF" = my_pal[1], "I" = my_pal[2]),
    name = NULL
  ) +
  labs(
    x = expression("Parasite Death Rate")
  ) +
  theme_classic(base_size = 20) +
  theme(
    legend.position = c(0.7, 0.8),
    axis.title.y = element_text(color = my_pal[1]),
    axis.title.y.right = element_text(color = my_pal[2]),
    axis.text.y = element_text(color = my_pal[1]),
    axis.text.y.right = element_text(color = my_pal[2])
  ) + theme(legend.position = "none")

p2


#immigration Rate (important to infecitous copepods)
ImmigrationPeriod = 50
d_W = 0.01
ImmigrationRate = seq(0, 0.5/50, by = 0.005/50) 

result_vectorImmigrationRate = numeric()
mean_IImmigrationRate = numeric()

for(i in seq_along(ImmigrationRate)){
  parameters["ImmigrationRate"] <- ImmigrationRate[i]
  
  sim = lsoda(y = Inits, times=timespan, parms = parameters, func="compute_derivatives", dllname = "Fishing_Flooding_Model", initfunc="initmod", maxsteps = 1e6)
  
  result_vectorImmigrationRate[i] <- sim[4*365,"HEF"] - sim[3*365,"HEF"] #dim(sim)[1] last row of object sim 
  mean_IImmigrationRate[i] <- mean(sim[(3*365):(4*365),"I"])
}

allresultsImmigrationRate = data.frame(ImmigrationRate,HEF = result_vectorImmigrationRate,I = mean_IImmigrationRate)


coeff <- 0.33

p3 = ggplot(allresultsImmigrationRate, aes(x = ImmigrationRate)) +
  geom_line(aes(y = HEF, color = "HEF"), linewidth = 1.2) +
  geom_line(aes(y = I / coeff, color = "I"), linewidth = 1.2) +
  scale_y_continuous(
    name = expression("Average Cumulative Fish Exposure ("*L^{-1}*")"),
    sec.axis = sec_axis(~ . * coeff,
                        name = expression("Average L3 Infected Copepods ("*L^{-1}*")"))
  ) +
  scale_color_manual(
    values = c("HEF" = my_pal[1], "I" = my_pal[2]),
    name = NULL
  ) +
  labs(
    x = expression("ImmigrationRate")
  ) +
  theme_classic(base_size = 20) +
  theme(
    legend.position = c(0.7, 0.8),
    axis.title.y = element_text(color = my_pal[1]),
    axis.title.y.right = element_text(color = my_pal[2]),
    axis.text.y = element_text(color = my_pal[1]),
    axis.text.y.right = element_text(color = my_pal[2])
  ) + theme(legend.position = "none")

p3

ggarrange(p1,p2,p3, nrow=1,ncol=3)


#Run simulations with d_W
###Immigration Period 
allresults = data.frame(ImmigrationPeriod,HEF = result_vector)
parameters["d_W"] = 0.01 #change for each sim 0,0.1,0.01
parameters["ImmigrationRate"] =  0.05/50
Inits = c(N = 7500, J = 6000, A = 700, Exposed_values, I = 0, Preds = 0,L3F = 0, HEF = 0)/15
timespan <- seq(0, 1825, by = 1)
ImmigrationPeriod = seq(0,130,by=1) 

result_vector = numeric()
mean_I = numeric()

for(i in seq_along(ImmigrationPeriod)){
  parameters["ImmigrationPeriod"] <- ImmigrationPeriod[i]
  
  sim = lsoda(y = Inits, times=timespan, parms = parameters, func="compute_derivatives", dllname = "Fishing_Flooding_Model", initfunc="initmod", maxsteps = 1e6)
  
  result_vector[i] <- sim[4*365,"HEF"] - sim[3*365,"HEF"] #dim(sim)[1] last row of object sim 
  #mean_I[i] <- mean(sim[(3*365):(4*365),"I"])
}

allresults = data.frame(ImmigrationPeriod,HEF = result_vector)

p0 = allresults %>% pivot_longer(cols = c(2),names_to = "State", values_to = "Value")
p0.1 = allresults %>% pivot_longer(cols = c(2),names_to = "State", values_to = "Value")
p0.01 = allresults %>% pivot_longer(cols = c(2),names_to = "State", values_to = "Value")

p0 = p0 %>% mutate(dW = "0")
p0.1 = p0.1 %>% mutate(dW = "0.1")
p0.01 = p0.01 %>% mutate(dW = "0.01")

allresultsDW = rbind(p0,p0.1,p0.01)


set3_colors <- brewer.pal(12, "Set3")
my_pal2 <- set3_colors[c(1,4,3)]

p4 = ggplot(allresultsDW, aes(x=ImmigrationPeriod,y=Value, group = dW, color = dW)) + 
  geom_line(linewidth=1.5) +
  scale_color_manual(
    values = my_pal2
  ) + theme_classic(base_size=20,) +
  labs(x = expression("Immigration Period"),y = expression("Average Cumulative Fish Exposure ("*L^{-1}*")"), color = "L3 Death Rate in Fish") +
  theme(legend.position = c(0.7, 0.7))
p4

ggarrange(p1,p4,p2,p3, nrow=2,ncol=2)

library(patchwork)
(p1 + p4) / (p2 + p3)
