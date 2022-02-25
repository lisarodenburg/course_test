###############################################################################
##                  R script to analyze swelling data Exp LR-035B             ##
##                    Lisa Rodenburg, June 2021          
##                      Output: 
##        
##
###############################################################################

# load or install packages
install.packages("dplyr")
install.packages("readxl")
install.packages("gglpot2")
install.packages("docstring")

library(dplyr)
library(readxl)
library(ggplot2)
library(docstring)

###############################################################################
##                  Open files                                      ##
###############################################################################

stats_foldchange <- read_excel(
  "data/raw/LR-035B Bluewasher and Thunder 384W.xlsx",
  sheet = "stats_foldchange"   ##Location of table in Excel sheet can be changed here
  #range ="CY5:DH101"          #Location of table in Excel sheet can be changed here
)

stats_AUC <- read_excel(
  "data/raw/LR-035B Bluewasher and Thunder 384W.xlsx",
  sheet= "stats_AUC"   ##Location of table in Excel sheet can be changed here
  #range="CY5:DH101"          #Location of table in Excel sheet can be changed here
)

###############
##         DEFINE FUNCTION
##############

Plot_swelling_graph <- function(data){
  #' @tile Plot swelling graphs
  #' @description This function plots swelling of organoids over time at different conditions
  #' @param t The dataset you want to use
  
  x <- ggplot(data, aes(x = t, y = foldchange_mean2, color = Activator)) +
    geom_point(size = 1.5) +
    geom_line(size = 1) +
    geom_errorbar(aes(ymin = foldchange_mean2 - foldchange_sd2, ymax = foldchange_mean2 + foldchange_sd2), 
                  width = 0.1) +
    scale_y_continuous(limits = c(90, 150)) +
    labs(
      x = "Time",
      y = "Swelling (%)"
    ) +
    theme_classic() +
    facet_grid(cols = vars(medium)) 
  return(x)
}

docstring(Plot_swelling_graph)

?Plot_swelling_graph

###############################################################################
##                  MAKE GRAPHS: swell% in time                              ##
###############################################################################

colnames(stats_foldchange)[colnames(stats_foldchange) == "medium (uL)"] <- "medium" #remove () from colname because R does not like this
stats_foldchange$"medium" <- as.character(stats_foldchange$"medium") #change to character

stats_foldchange <-stats_foldchange %>% #remove NA values
  filter(!is.na(foldchange_mean)) %>%
  filter(well!= "D01") #exlusion of well D01 because of low amount of particles

stats_foldchange2 <- stats_foldchange %>% #Make table per condition
  group_by(t, medium, Activator) %>%
  summarise(foldchange_mean2 = mean(foldchange_mean, na.rm=T), 
            foldchange_sd2 = sd(foldchange_mean), 
            n=n())
#TODO: variable name does not five useful information, the number 2

#Plot swelling graph per donor (swell% in time)
  Plot_swelling_graph(stats_foldchange2)
    
###############################################################################
##                  MAKE GRAPHS: AUC bar graphs                              ##
###############################################################################  

colnames(stats_AUC)[colnames(stats_AUC) == "medium (uL)"] <- "medium" #remove () from colname because R does not like this
stats_AUC$"medium" <- as.character(stats_AUC$"medium") #change to character

stats_AUC <- stats_AUC %>% #remove NA values
  filter(!is.na(auc_mean))%>%
  filter(well != "D01") #exlusion of well D01 because of low amount of particles

stats_AUC2 <- stats_AUC %>% #Make table per condition
  group_by(medium, Activator) %>%
  summarise(auc_mean2 = mean(auc_mean, na.rm=T), 
            auc_sd2 = sd(auc_mean), 
            n=n())

# Make graphs (AUC)
ggplot(data = stats_AUC2, aes(x = Activator, y = auc_mean2, fill = medium))+
  geom_bar(stat = "identity", color = "black", position = "dodge") +
  geom_errorbar(aes(ymin = auc_mean2 - auc_sd2, ymax = auc_mean2 + auc_sd2), 
                width = .2, 
                color = "black", 
                position = position_dodge(width = 0.9)) +
  labs(y = "AUC (t = 60 min)", x = NULL) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))


