###############################################################################
##                  R script to analyze swelling data Exp LR-035B             ##
##                    Lisa Rodenburg, June 2021          
##                      Output: 
##        
##
###############################################################################

lint("src/R script swelling - analyse LR-035B.R")

## First specify the packages of interest
packages <- c("dplyr", "readxl", "ggplot2", "ggpubr", "ggpattern", "lintr")

## Now load or install&load all
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

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

###############################################################################
##                  MAKE GRAPHS: swell% in time                              ##
###############################################################################

rename_colname <- function (input, colname, new_colname){
  colnames(input)[colnames(input) == colname] <- new_colname
}

rename_colname(stats_foldchange, "medium (uL)", "medium")

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
  ggplot(stats_foldchange2, aes(x = t, y = foldchange_mean2, color = Activator)) +
    geom_point(size = 1.5) +
    geom_line(size = 1) +
    geom_errorbar(aes(ymin = foldchange_mean2 - foldchange_sd2, ymax = foldchange_mean2 + foldchange_sd2), 
      width = 0.1) +
    #scale_x_continuous(breaks = seq(min(Time_minutes), max(Time_minutes), by = 15)) +
    scale_y_continuous(limits = c(90, 150)) +
    labs(
      x = "Time",
      y = "Swelling (%)"
     # title = Results_foldchange_Selection$Donor
    ) +
    theme_classic() +
    #scale_color_manual(values = FDA2_color) +
    facet_grid(cols = vars(medium)) 
    
#ggsave(paste(Results_AUC$Experiment,"Swell HNEC0291,292,293 FDA2.pdf") , plot = last_plot(), device=pdf, widt = 10, height = 6, useDingbats = F)


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
  #geom_point(data = Results_combined, aes(x = Activator, y = AUC), position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = auc_mean2 - auc_sd2, ymax = auc_mean2 + auc_sd2), 
                width = .2, 
                color = "black", 
                position = position_dodge(width = 0.9)) +
  #ggtitle (Results_AUC_Donor$Donor) +
  labs(y = "AUC (t = 60 min)", x = NULL) +
  theme_classic() +
  #scale_fill_manual(values = c("black", "lightgrey", "white"), labels = medium) +
  #facet_grid(rows = vars(Condition)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

#ggsave(paste(Results_AUC_selection$Experiment,"AUC_HNEC0365 n1,n2.pdf"), device=pdf, widt = 10, height = 6, useDingbats = F)

