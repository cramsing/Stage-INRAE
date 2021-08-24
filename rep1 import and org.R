# #### Data import #### 

#Rep 1 import 
#all this will not be run because the code is too heavy to run each time. 
#Instead the resulting df is saved as an r file and loaded each time
# ## Compilation des fichiers csv All_lesions
# # list_files = list.files(path="E:/Data analysis/12.08.21 all lesions", full.names = TRUE)
# # a=0
# # b = "All_lesions.csv"
# # for( x in list_files) {
#   if (grepl(b,x,fixed=TRUE)) {
#     tmp = read_delim(file = x, delim = "\t")
#     if (a == 0) {    df = tmp; a=1  }
#     else {
#       df = bind_rows(df, tmp)
#     }
#   }
# }
# 
# df
# 
# ## Renaming column names to match file names
# colnames(df)
# names(df)[names(df) == "image"] <- "mutant_cultivar_date"
# 
# 
# df
# 
# 
# ## Separation of cultivar and mutant 
# df["mutant_cultivar"] = lapply(df["mutant_cultivar_date"],as.character)
# cs <- strsplit(df[["mutant_cultivar_date"]],"_")
# cs2 <- data.frame(do.call(rbind,cs),stringsAsFactors=FALSE)
# names(cs2)[names(cs2) == "X1"] <- "mutant"
# names(cs2)[names(cs2) == "X2"] <- "cultivar"
# names(cs2)[names(cs2) == "X3"] <- "date"
# 
# 
# df2 = cbind(df, cs2)
# 
# df2$genotype =paste0(df2$mutant, ".",df2$cultivar)
# 
# # deleting lesions with 0 size
# df2=df2[df2$lesion.surface>0,]
# 
# 
# head(df2)
# tail(df2)
# save(df2, file="df2.RData") #saving just to be sure (and bc the above code is quite annoying to run)

load("df2.RData")
df2 <- df2[d25$lesion.status=="keep",] #only keeps lesions i've approved !
df2[ , c('genotype', 'mutant_cultivar_date')] <- list(NULL) #deleting redunant columns 
colnames(df2)

rep1 <- df2[, c(16, 17, 18, 19, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)] #reordring columns and puttting in a new df
#really could've done this and previous deletions in the same step but hey 

# #### Dataframe clean up ####
colnames(rep1)
library(dplyr)
rep1 <- rename(rep1, sample = mutant_cultivar) #rename column
colnames(rep1)




# ag1 <- aggregate(df3$lesion.surface,df3["genotype_lignee"],mean) ; names(ag1)[2] <- "moyenne"
# ag2 <- aggregate(df3$lesion.surface,df3["genotype_lignee"],var) ; names(ag2)[2] <- "variance"
# ag <- merge(ag1,ag2)
# ag
# windows()
# plot(ag$moyenne,ag$variance,main=titre,xlab="Moyenne", ylab="Variance")
# windows()
# plot(ag$moyenne,sqrt(ag$variance),main=titre,xlab="Moyenne", ylab="Ecart-type")
# abline(lm(sqrt(variance)~0+moyenne,data=ag))


library(dplyr)
rep1 <- rep1 %>% 
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
                         '1' = "FA1",
                         '2' = "FA ectopic",
                         '3' = 'FA2',
                         '4' = 'FA3',
                         '5' = 'FB1',
                         '6'= 'FB2',
                         '7' = 'FB3',
                         '8' = 'FB ectopic',
                         '9' = 'FK1',
                         '10' = 'FK2',
                         '11' = 'FK3',
                         '12' = 'FK ectopic',
                         '13' = 'FR13 RFP',
                         '14' = 'FR13 WT',
                         '15' = 'GB1',
                         '16' = 'GB2',
                         '17' = 'GB3',
                         '18' = 'Guy11 RFP_1',
                         '19' = 'GK1',
                         '20' = 'GK2',
                         '21' = 'GK3',
                         '22' = 'GK ectopic',
                         '23' = 'Guy11 RFP_2',
                         '24' = 'Guy11 WT')) %>% #renaming mutants
  mutate(cultivar = recode(cultivar,
                           'A' = "Aichi Asahi", 
                           'B' = "Bl1", 
                           'N' = "Nipponbare",
                           'K' = "Kasalath",
                           'K6' = "K60", 
                           'S' = "Shin 2",
                           'T' = "Tsuyuake")) %>% #renaming cultivars
  mutate(rep = recode(date,
                      '12'= "1")) #rename date as rep (rep 1 was taped on 12.08.21 so the date was 12)

rep1$date <- NULL #delete date column
rep1$sample <- NULL #delete sample column

colnames(rep1)

rep1 <- rep1[, c(6, 1, 19, 18, 2, 20, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17)] 

colnames(rep1)

rep1$lesion.status <- NULL 




#rep1 lesion number

rep1 <- rename(rep1, lesion = lesion.number) #rename column because its not the number of lesions,
#it's the number OF THE lesion. confusing. 



#removing lesions above 8,000 pixels (leAF Tool agglomerrated lesions so
#unfortunately we have to remove them bc they are not accurate. at the same time
#removing them is not accurate but hey)
rep1 <- subset(rep1, lesion.surface < 8000) #only keeps data with lesion surface
#values under 8,000 pixels


save (rep1, file = "rep1.RData") #saving dataframe


## rep 1  lesion count dataframe---- 
#lesion count df
library(dplyr)
rep1.lc <- rep1 %>%
  group_by(mutant,leaf.number,cultivar)%>%
  tally(wt = NULL) #tallies lesions. wt=NULL makes it NOT weight the 
#lesions by number



rep1.lc <- rename(rep1.lc, lesion.count = n)

#percent lesion df, will join with same lesion count df
library(dplyr)
rep1.leaf.s <- rep1 %>%
  group_by(mutant,leaf.number,cultivar)%>%
  summarise(leaf.surface = sum(leaf.surface))

rep1.les.s <- rep1 %>%
  group_by(mutant,leaf.number,cultivar)%>%
  summarise(lesion.surface = sum(lesion.surface))

rep1.lc <- dplyr::left_join(
  rep1.lc, rep1.leaf.s, by = c("mutant","leaf.number","cultivar"))%>%
  left_join(.,rep1.les.s, by = c("mutant","leaf.number",
                                 "cultivar"))
rep1.lc <- rep1.lc %>%
  group_by(mutant,leaf.number,cultivar)%>%
  mutate(lesion.percent = lesion.surface/leaf.surface)

library(dplyr)
rep1.lc <- rep1.lc %>% 
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
