# Goikoetxea et al., 2022

## Code used during RIA analyses (example for E2)

install.packages("Rmisc")
install.packages("ggplot2")
install.packages("ggpubr")
install.packages("reshape")
install.packages('ggplot2', repos='http://cran.us.r-project.org')
install.packages("ggpubr")
if(!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/ggpubr")
install.packages("png")

library(Rmisc) 
library(ggplot2) 
library(ggpubr) 
library(reshape) 
library(png)
library(magrittr)

setwd("/Volumes/TITANIUM/qPCR:Spotty paper/RIA/RIA_Spotty_Decitabine_2017") 
data <- read.csv("RIA_Spotty_Cortisol_2017_by_treatment_grouped_controls.csv",sep = ",")

head(data) 
str(data) 
dim(data)

melted <- melt(data, id="Stage")
head(melted)
summary_data <- summarySE(melted, measurevar="value", Stagevars=c("Stage"))
help("summarySE")
summary_data

head(data)
data_log <- data
data_log$E2 <- log(data_log$E2)

head(data_log)

melted_log <- melt(data_log, id="Stage")
summary_data_log <- summarySE(melted_log, measurevar="value", Stagevars=c("Stage"))
help("summarySE")
summary_data_log

hist(log(data$E2)) 

kruskal.test(E2 ~ Stage, data = data_log)
?kruskal.test

install.packages("PMCMRplus")
install.packages("FSA")
install.packages("dunn.test")

library(dunn.test)
library(FSA)

require(PMCMR)
data(data_log2)
attach(data_log)


dunnTest(E2~Stage, data=data_log, method="bh")

## Code used during nanoString gene expression analyses (example for amh)

library(Rmisc) 
library(ggplot2) 
library(ggpubr) 
library(reshape) 
library(png)
library(magrittr)

data <- read.csv("amh.csv",sep = ",")

head(data) 
str(data) 
dim(data) 

melted <- melt(data, id="Treatment")
head(melted)
summary_data <- summarySE(melted, measurevar="value", groupvars=c("Treatment"))
help("summarySE")
summary_data

head(data)
data_log2 <- data
data_log2$amh <- log2(data_log2$amh)

head(data_log2)

melted_log2 <- melt(data_log2, id="Treatment")
summary_data_log2 <- summarySE(melted_log2, measurevar="value", groupvars=c("Treatment")) 

to check if by logarithmically transforming our data, the standard deviations grew closer to each other
help("summarySE")
summary_data_log2

hist(log2(data$amh)) 

p <- ggboxplot(data_log2, x = "Treatment", y = "amh", palette = c("#000000", "#000000", "#000000"),
               color = "Treatment", ylab = expression(paste("Relative normalised ", italic("amh"), " expression")), add = "jitter", order = c("CF", "TF", "TP"),
               ylim = c(0,15))+
  theme (panel.background = element_rect(fill = 'white', colour = 'white'))+
  theme (axis.line.x = element_line(colour = 'black'))+
  theme (axis.line.y = element_line(colour = 'black'))+
  theme (axis.text.x = element_text(angle=45, hjust=1))+
  theme(legend.position = "none")+
  scale_x_discrete(limits=c ("CF", "TF", "TP"), name="Treatment\n")
p

pdf (file="amh.pdf", width=5, height=4.25)
p
dev.off()

kruskal.test(amh ~ Treatment, data = data_log2)
?kruskal.test

install.packages("PMCMRplus")
install.packages("FSA")
install.packages("dunn.test")

library(dunn.test)
library(FSA)

require(PMCMR)
data(data_log2)
attach(data_log2)


dunnTest(amh~Treatment, data=data_log2, method="bh")
