library(readxl)
library(plyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(tidyverse)
library(data.table)

#read in IRMS data and merge into dataframe #
#select folder where data is
#setwd("/Users/f00502n/Documents/Dartmouth/DOE Mycorrhizae/chamber/fractionation/IRMS data")
setwd("/Users/f003833/Documents/GitHub/fractionation/IRMS data")
files <- (Sys.glob("*.xlsx"))
alldata <- lapply(files, function(x) read_excel(x, sheet = 1)) 
alldata <- rbindlist(alldata, idcol = T)

#IRMS data for bulk mesh bag soil #
setwd("/Users/f003833/Documents/GitHub/fractionation")
bulk <- read.csv("chamber_corrected_IRMS_compiled.csv")[,c("Identifier.1", "PotNumber", "SampleType",
                                                        "corrected.percent.C", "C13.corrected")]
bulk <- bulk %>%
  filter(SampleType == "450 micron") %>%
  group_by(PotNumber) %>%
  dplyr::summarise(bulk.perC = mean(corrected.percent.C),
            bulk.C13 = mean(C13.corrected))

#read in chamber info and fractionation weights
chamber.info <- read.csv("Chamber Harvest Soil Bags (r).csv")
tag.id <- read.csv("Tag ID.csv")[,c(1:3,5)]
fractionation.wt <- read.csv("Chamber Fractionation.csv")[9:61, c("Sample","Wt", "FLF...OLF", "DF")]
colnames(fractionation.wt) <- c("Sample", "total.fractionation.wt", "FLF.OLF.wt", "DF.wt")
fractionation.wt$Number <- as.numeric(str_replace_all(fractionation.wt$Sample,"[A-Z, ]", "")) 

#merge chamber info, tag id info, fractionation weights, and bulk C and C13 data
chamber.info <- merge(chamber.info, tag.id, by.x = "Sample", by.y = "Number")
chamber.info <- merge(chamber.info, fractionation.wt, by.x = "Sample", by.y = "Number")
chamber.info <- merge(chamber.info, bulk, by.x = "Sample", by.y = "PotNumber")
rm(tag.id, fractionation.wt, bulk)

#remove unnecessary columns (total.wt.g = calculated total dry weight of bag)
chamber.info <- chamber.info[ , c("Sample", "mesh.bag.size", "total.wt.g", "Species", "N.treatment",
                                  "total.fractionation.wt", "FLF.OLF.wt", "DF.wt", 
                                  "bulk.perC", "bulk.C13", "Myco..Association")]


#### process IRMS data ####
# create new column that doesn't include standard numbers and replicate indicators
alldata$Identifier.2 <- ifelse(
  grepl("peach|cocoa|COCOA|spruce", alldata$Identifier.1), str_replace_all(alldata$Identifier.1, "[0-9, ]", ""),
  str_remove(alldata$Identifier.1, " REP|R1- ")
)
#remove standards
alldata <- alldata %>%
  filter(Identifier.2 != "peach" & Identifier.2 != "cocoa" & Identifier.2 != "BLANK")
#check coefficient of variation across duplicates
check.coef.var <- alldata %>%
  filter(duplicated(Identifier.2) | duplicated(Identifier.2, fromLast = TRUE))%>%
  group_by(Identifier.2) %>%
  dplyr::summarize(percent.C.Coef.Var = sd(corrected.percent.C)/ mean(corrected.percent.C) *100,
                   percent.N.Coef.Var = sd(corrected.percent.N)/ mean(corrected.percent.N) *100,
                   C13.Coef.Var = abs(sd(C13.corrected)/ mean(C13.corrected) *100),
                   N15.Coef.Var = sd(N15.corrected)/ mean(N15.corrected) *100)

#calculate average for each sample
all.avg <- alldata %>%
  group_by(Identifier.2) %>%
  dplyr::summarise(perC = mean(corrected.percent.C),
            C13 = mean(C13.corrected))
# create separate columns for ID number and fraction
all.avg$Fraction <- str_replace_all(all.avg$Identifier.2,"[0-9, ]", "") 
all.avg$Sample <- as.numeric(str_replace_all(all.avg$Identifier.2,"[A-Z,+ ]", ""))

#pivot wider
all.avg <- pivot_wider(all.avg, id_cols = Sample, names_from = Fraction, values_from = c(perC, C13))

#merge IRMS data with chamber info
alldata2 <- merge(all.avg, chamber.info, by = "Sample")

### do some math ####

# equation for calculating 
fraction <- function(C13, total.C){
  R <- (C13/1000 + 1) * 0.0112372 #PDB standard value
  F.s <- R / (1 + R) 
  g.13C <- F.s * total.C
  return(g.13C)
}

# calculate % carbon recovery 
carbon.data <- alldata2 %>%
  group_by(Sample) %>%
  mutate( 
    total.fractionated.mgC = bulk.perC/100 * total.fractionation.wt *1000,
    DF.mgC = perC_DF / 100 * DF.wt *1000,
    OLF.FLF.mgC = `perC_OLF+FLF` / 100 * FLF.OLF.wt *1000,
    recovery = (DF.mgC + OLF.FLF.mgC) / total.fractionated.mgC * 100
    )
#for some reason, C recovery is unreasonably high ??? don't know why

#calculate ratio for soil
carbon.data$mgC13.total <- fraction(carbon.data$bulk.C13, carbon.data$total.fractionated.mgC)
carbon.data$mgC13.DF <- fraction(carbon.data$C13_DF, carbon.data$DF.mgC)
carbon.data$mgC13.OLF.FLF <- fraction(carbon.data$`C13_OLF+FLF`, carbon.data$OLF.FLF.mgC)

carbon.data$recovery.C13 <- (carbon.data$mgC13.DF + carbon.data$mgC13.OLF.FLF) / carbon.data$mgC13.total *100
#this is also really high!!
carbon.data <- carbon.data %>%
  mutate( perC13.DF = mgC13.DF/ (mgC13.DF + mgC13.OLF.FLF) * 100,
          perC13.OLF.FLF = mgC13.OLF.FLF/ (mgC13.DF + mgC13.OLF.FLF) * 100
  )

# carbon.data$R.total <- (carbon.data$bulk.C13/1000 + 1) * 0.0112372 #PDB standard value
# carbon.data$F.total <- carbon.data$R.total / (1 + carbon.data$R.total) 
# carbon.data$mg.13C.v2 <- carbon.data$F.total * carbon.data$total.fractionated.mgC
# 
carbon.data$N.treatment <- factor(carbon.data$N.treatment, levels = c("low", "mid", "high "))

boxplot(carbon.data$perC13.DF~carbon.data$N.treatment)
boxplot(carbon.data$perC13.DF~carbon.data$Myco..Association+carbon.data$N.treatment)

boxplot(carbon.data$perC13.OLF.FLF~carbon.data$N.treatment)
boxplot(carbon.data$perC13.OLF.FLF~carbon.data$Myco..Association+carbon.data$N.treatment) 

ggplot() + geom_point(data = carbon.data, aes( x = carbon.data$`C13_OLF+FLF`, y = carbon.data$C13_DF)) + 
  geom_abline(slope = 1, intercept = 0)

hist(carbon.data$recovery)




