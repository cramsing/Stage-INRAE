
#Rep 2 import 

#Instead the resulting df is saved as an r file and loaded each time
# Compilation des fichiers csv All_lesions


library(dplyr)
library(readr)
list_files = list.files(path="E:/Data analysis/23.08.21 all lesions", full.names = TRUE)
a=0
b = "All_lesions.csv"

for( x in list_files) {
  if (grepl(b,x,fixed=TRUE)) {
    tmp = read_delim(file = x, delim = "\t")
    if (a == 0) {    df4 = tmp; a=1  }
    else {
      df4 = bind_rows(df4, tmp)
    }
  }
}

df4

## Renaming column names to match file names
colnames(df4)
#don't know why there is a super weird empty column. I've tried to troubleshoot it
#but am not getting anywhere so we are just going to delete it
df4$`image;leaf.number;leaf.surface;lesion.status;lesion.number;lesion.surface;lesion.perimeter;lesion.radius.mean;lesion.radius.sd;lesion.radius.min;lesion.radius.max;m.cx;m.cy;m.majoraxis;m.eccentricity;m.theta` <-NULL
df4$`<?xml version="1.0" encoding="UTF-8" standalone="yes"?>` <- NULL
names(df4)[names(df4) == "image"] <- "mutant_cultivar"

df4

## Separation of cultivar and mutant 
df4["mutant_cultivar"] = lapply(df4["mutant_cultivar"],as.character)
cs <- strsplit(df4[["mutant_cultivar"]],"_")
cs2 <- data.frame(do.call(rbind,cs),stringsAsFactors=FALSE)
names(cs2)[names(cs2) == "X1"] <- "mutant"
names(cs2)[names(cs2) == "X2"] <- "cultivar"
# 
# # 
# # 
df5 = cbind(df4, cs2)
# # 
df5$genotype =paste0(df5$mutant, ".",df5$cultivar)
# # 
#deleting lesions with 0 size
df5=df5[df5$lesion.surface>0,]
df5 <- df5[df5$lesion.status=="keep",] #only keeps lesions i've approved !
# # 
# # 
head(df5)
tail(df5)
save(df5, file="df5.RData") #saving 

load("df5.RData")

df5[ , c('genotype', 'mutant_cultivar_date')] <- list(NULL) #deleting redundant columns 
colnames(df5)

punch <- df5[, c(17, 18, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16)] 
#reordring columns and putting in a new df

# #### Dataframe clean up ####
colnames(punch)
library(dplyr)
punch <- rename(punch, sample = mutant_cultivar) #rename column
colnames(punch)



library(dplyr)
punch <- punch %>% 
  mutate(effector = recode(mutant,
                           '1' = "AVR Pia",
                           '2' = "AVR Pia",
                           '3' = 'AVR Pia',
                           '4' = 'AVR Pia',
                           '5' = 'AVR PiB',
                           '6'= 'AVR PiB',
                           '7' = 'AVR PiB',
                           '8' = 'AVR PiB',
                           '9' = 'AVR Pik',
                           '10' = 'AVR Pik',
                           '11' = 'AVR Pik',
                           '12' = 'AVR Pik',
                           '13' = 'RFP',
                           '14' = 'WT',
                           '15' = 'AVR PiB',
                           '16' = 'AVR PiB',
                           '17' = 'AVR PiB',
                           '18' = 'RFP',
                           '19' = 'AVR Pik',
                           '20' = 'AVR Pik',
                           '21' = 'AVR Pik',
                           '22' = 'AVR Pik',
                           '23' = 'RFP',
                           '24' = 'WT')) %>% #grouping effectors
  mutate(isolate = recode(mutant,
                          '1' = "FR13",
                          '2' = "FR13",
                          '3' = 'FR13',
                          '4' = 'FR13',
                          '5' = 'FR13',
                          '6'= 'FR13',
                          '7' = 'FR13',
                          '8' = 'FR13',
                          '9' = 'FR13',
                          '10' = 'FR13',
                          '11' = 'FR13',
                          '12' = 'FR13',
                          '13' = 'FR13',
                          '14' = 'FR13',
                          '15' = 'Guy11',
                          '16' = 'Guy11',
                          '17' = 'Guy11',
                          '18' = 'Guy11',
                          '19' = 'Guy11',
                          '20' = 'Guy11',
                          '21' = 'Guy11',
                          '22' = 'Guy11',
                          '23' = 'Guy11',
                          '24' = 'Guy11')) %>% # grouping isolates
  mutate(mutant = recode(mutant,
                         '1' = "mFA1",
                         '2' = "FA ectopic",
                         '3' = 'mFA2',
                         '4' = 'mFA3',
                         '5' = 'mFB1',
                         '6'= 'mFB2',
                         '7' = 'mFB3',
                         '8' = 'FB ectopic',
                         '9' = 'mFK1',
                         '10' = 'mFK2',
                         '11' = 'mFK3',
                         '12' = 'FK ectopic',
                         '13' = 'FR13 RFP',
                         '14' = 'FR13 WT',
                         '15' = 'mGB1',
                         '16' = 'mGB2',
                         '17' = 'mGB3',
                         '18' = 'Guy11 RFP_1',
                         '19' = 'mGK1',
                         '20' = 'mGK2',
                         '21' = 'mGK3',
                         '22' = 'GK ectopic',
                         '23' = 'Guy11 RFP_2',
                         '24' = 'Guy11 WT')) #renaming mutants


#rename misspelled cultivars
punch$cultivar <- as.factor(punch$cultivar)
levels(punch$cultivar) # print misspelled cultivar names
# [1] "Aichi Asahi"  "Aichi Asashi" "Aishi Asahi"  "Bl1"  "K60"  "Kasalath"   
# "Nipponbare" "Shin2" 
library(dplyr)
punch <- punch %>% 
  mutate(cultivar = recode(cultivar,
                           'Aichi Asashi' = 'Aichi Asahi',
                           'Aishi Asahi' = 'Aichi Asahi'))

levels(punch$cultivar) #verify names
# "Aichi Asahi" "Bl1"  "K60"   "Kasalath"  "Nipponbare" "Shin 2"     
#everything looks correct

punch$rep <- 'punch' #adding column for rep description
#delete superfluous columns
punch$sample <- NULL #delete sample column
punch$lesion.status <- NULL #delete lesion.status column
punch$lesion.number.1 <- NULL #delete lesion.number.1 column


#punch lesion number
punch <- rename(punch, lesion = lesion.number) #rename column because its not the number of lesions,
#it's the number OF THE lesion. confusing. 

colnames(punch)
load("rep1.RData")
colnames(rep1)#check order of column names against rep1 df

punch <- punch[, c(5, 1, 18, 17, 2, 19, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16)] 

colnames(punch)

punch <- na.omit(punch)




save(punch, file = "punch.RData") #savedf
