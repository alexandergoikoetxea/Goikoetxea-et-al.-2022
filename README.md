# Goikoetxea-et-al.-2022

# Code used during RIA analyses

install.packages("Rmisc")
install.packages("ggplot2")
install.packages("ggpubr")
install.packages("reshape")
install.packages('ggplot2', repos='http://cran.us.r-project.org')
install.packages("ggpubr")
if(!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/ggpubr")
install.packages("png")

library(Rmisc) # Statistical test
library(ggplot2) # Graphics
library(ggpubr) # Graphics
library(reshape) # Melt
library(png)
library(magrittr)

setwd("/Volumes/TITANIUM/qPCR:Spotty paper/RIA/RIA_Spotty_Decitabine_2017") #to set your folder as the working directory
data <- read.csv("RIA_Spotty_Cortisol_2017_by_treatment_grouped_controls.csv",sep = ",")

head(data) # Check data
str(data) # Check data structure
dim(data) # Check data dimensions

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

hist(log(data$E2)) #para mirar normalidad

kruskal.test(E2 ~ Stage, data = data_log)
?kruskal.test

install.packages("PMCMRplus")
install.packages("FSA")
install.packages("dunn.test")

library(dunn.test)
library(FSA)

require(PMCMR)
#data(data_log2)
attach(data_log)


dunnTest(E2~Stage, data=data_log, method="bh")
