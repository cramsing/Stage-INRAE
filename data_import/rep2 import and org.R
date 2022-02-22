
#**Rep 2 import and organization**



#Data import ####
# Compilation des fichiers csv All_lesions


 library(dplyr)
 library(readr)
 list_files = list.files(path="E:/Data analysis/17.08.21 all lesions", full.names = TRUE)
 a=0
 b = "All_lesions.csv"
 
 for( x in list_files) {
   if (grepl(b,x,fixed=TRUE)) {
     tmp = read_delim(file = x, delim = "\t")
     if (a == 0) {    df1 = tmp; a=1  }
    else {
       df1 = bind_rows(df1, tmp)
     }
   }
 }
 
df1
  
## Renaming column names to match file names
colnames(df1)
#don't know why there is a super weird empty column. I've tried to troubleshoot it
#but am not getting anywhere so we are just going to delete it
df1$`image;leaf.number;leaf.surface;lesion.status;lesion.number;lesion.surface;lesion.perimeter;lesion.radius.mean;lesion.radius.sd;lesion.radius.min;lesion.radius.max;m.cx;m.cy;m.majoraxis;m.eccentricity;m.theta` <-NULL
df1$`<?xml version="1.0" encoding="UTF-8" standalone="yes"?>` <- NULL
names(df1)[names(df1) == "image"] <- "mutant_cultivar"

df1

## Separation of cultivar and mutant 
df1["mutant_cultivar"] = lapply(df1["mutant_cultivar"],as.character)
cs <- strsplit(df1[["mutant_cultivar"]],"_")
cs2 <- data.frame(do.call(rbind,cs),stringsAsFactors=FALSE)
names(cs2)[names(cs2) == "X1"] <- "mutant"
names(cs2)[names(cs2) == "X2"] <- "cultivar"
# 
# # 
# # 
df3 = cbind(df1, cs2)
# # 
df3$genotype =paste0(df3$mutant, ".",df3$cultivar)
# # 
#deleting lesions with 0 size
df3=df3[df3$lesion.surface>0,]
df3 <- df3[df3$lesion.status=="keep",] #only keeps lesions i've approved !
# # 
# # 
head(df3)
tail(df3)
save(df3, file="df3.RData") #saving just to be sure (and bc the above code is quite annoying to run)

load("df3.RData")

df3[ , c('genotype', 'mutant_cultivar_date')] <- list(NULL) #deleting redunant columns 
colnames(df3)

rep2 <- df3[, c(17, 18, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16)] #reordring columns and puttting in a new df
#really could've done this and previous deletions in the same step but hey 

# Dataframe clean up ####
colnames(rep2)
library(dplyr)
rep2 <- rename(rep2, sample = mutant_cultivar) #rename column
colnames(rep2)



library(dplyr)
rep2 <- rep2 %>% 
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
rep2$cultivar <- as.factor(rep2$cultivar)
levels(rep2$cultivar) # print misspelled cultivar names
# [1] "AichiAsahhi" "AichiAsahi"  "AichiAsahih" "Bl1"         "K60"      
#"Kasalath"    "Kasalth"     "Nipponbare" 
# [9] "Shin2"       "Tsuyake"     "Tsuyuake"
library(dplyr)
rep2 <- rep2 %>% 
  mutate(cultivar = recode(cultivar,
                           'AichiAsahhi' = 'Aichi Asahi',
                           'AichiAsahi' = 'Aichi Asahi',
                           'AichiAsahih' = 'Aichi Asahi',
                           'Kasalth' = 'Kasalath',
                           'Shin2' = 'Shin 2',
                           'Tsuyake' = 'Tsuyuake'))
levels(rep2$cultivar) #verify names
# "Aichi Asahi" "Bl1"  "K60"   "Kasalath"  "Nipponbare" "Shin 2"   "Tsuyuake"   
#everything looks correct

rep2$rep <- 2 #adding column for rep2
colnames(rep2)
rep2 <- rep2[, c(7, 1, 20,19, 2, 21, 3, 4, 5,6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18)] 

#delete superfluous columns
rep2$sample <- NULL #delete sample column 
rep2$lesion.number.1 <- NULL #delete lesion.number.1 column


#rep2 lesion number
rep2 <- rename(rep2, lesion = lesion.number) #rename column because its not the number of lesions,
#it's the number OF THE lesion. confusing. 


#removing lesions above 8,000 pixels (leAF Tool agglomerrated lesions so
#unfortunately we have to remove them bc they are not accurate. at the same time
#removing them is not accurate but hey)
rep2 <- subset(rep2, lesion.surface < 8000) #only keeps data with lesion surface
#values under 8,000 pixels


save (rep2, file = "rep2.RData") #saving dataframe

## rep 2 lesion count dataframe ----
#lesion count df
library(dplyr)
rep2.lc <- rep2 %>%
  group_by(mutant,leaf.number,cultivar)%>%
  tally(wt = NULL) #tallies lesions. wt=NULL makes it NOT weight the 
#lesions by number

rep2.lc <- rename(rep2.lc, lesion.count = n)

#percent lesion df, will join with same lesion count df
library(dplyr)
rep2.leaf.s <- rep2 %>%
  group_by(mutant,leaf.number,cultivar)%>%
  summarise(leaf.surface = sum(leaf.surface))

rep2.les.s <- rep2 %>%
  group_by(mutant,leaf.number,cultivar)%>%
  summarise(lesion.surface = sum(lesion.surface))

rep2.lc <- dplyr::left_join(
  rep2.lc, rep2.leaf.s, by = c("mutant","leaf.number","cultivar"))%>%
  left_join(.,rep2.les.s, by = c("mutant","leaf.number",
                                 "cultivar"))
rep2.lc <- rep2.lc %>%
  group_by(mutant,leaf.number,cultivar)%>%
  mutate(lesion.percent = lesion.surface/leaf.surface)

library(dplyr)
rep2.lc <- rep2.lc %>% 
  mutate(effector = recode(mutant,
                           'FA1' = "AVR Pia",
                           'FA ectopic' = "AVR Pia",
                           'FA2' = 'AVR Pia',
                           'FA3' = 'AVR Pia',
                           'FB1' = 'AVR PiB',
                           'FB2'= 'AVR PiB',
                           'FB3' = 'AVR PiB',
                           'FB ectopic' = 'AVR PiB',
                           'FK1' = 'AVR Pik',
                           'FK2' = 'AVR Pik',
                           'FK3' = 'AVR Pik',
                           'FK ectopic' = 'AVR Pik',
                           'FR13 RFP' = 'RFP',
                           'FR13 WT' = 'WT',
                           'GB1' = 'AVR PiB',
                           'GB2' = 'AVR PiB',
                           'GB3' = 'AVR PiB',
                           'Guy11 RFP_1' = 'RFP',
                           'GK1' = 'AVR Pik',
                           'GK2' = 'AVR Pik',
                           'GK3' = 'AVR Pik',
                           'GK ectopic' = 'AVR Pik',
                           'Guy11 RFP_2' = 'RFP',
                           'Guy11 WT' = 'WT')) %>% #grouping effectors
  mutate(isolate = recode(mutant,
                          'FA1' = "FR13",
                          'FA ectopic' = "FR13",
                          'FA2' = 'FR13',
                          'FA3' = 'FR13',
                          'FB1' = 'FR13',
                          'FB2'= 'FR13',
                          'FB3' = 'FR13',
                          'FB ectopic' = 'FR13',
                          'FK1' = 'FR13',
                          'FK2' = 'FR13',
                          'FK3' = 'FR13',
                          'FK ectopic' = 'FR13',
                          'FR13 RFP' = 'FR13',
                          'FR13 WT' = 'FR13',
                          'GB1' = 'Guy11',
                          'GB2' = 'Guy11',
                          'GB3' = 'Guy11',
                          'Guy11 RFP_1' = 'Guy11',
                          'GK1' = 'Guy11',
                          'GK2' = 'Guy11',
                          'GK3' = 'Guy11',
                          'GK ectopic' = 'Guy11',
                          'Guy11 RFP_2' = 'Guy11',
                          'Guy11 WT' = 'Guy11'))


