#Spray inoculation

#load replicates
load("rep1.RData")
load("rep2.RData")

spray <- rbind(rep1,rep2) #combines reps 1 and 2 

# Lesion Count ------------------------------------------------------------
## spray -----
#lesion count df
library(dplyr)
lc <- spray %>%
  group_by(mutant,leaf.number,cultivar)%>%
  tally(wt = NULL) #tallies lesions. wt=NULL makes it NOT weight the 
#lesions by number

lc <- rename(lc, lesion.count = n)

#percent lesion df, will join with same lesion count df
library(dplyr)
leaf.s <- spray %>%
  group_by(mutant,leaf.number,cultivar)%>%
  summarise(leaf.surface = sum(leaf.surface))

les.s <- spray %>%
  group_by(mutant,leaf.number,cultivar)%>%
  summarise(lesion.surface = sum(lesion.surface))

lc <- dplyr::left_join(
  lc, leaf.s, by = c("mutant","leaf.number","cultivar"))%>%
  left_join(.,les.s, by = c("mutant","leaf.number",
                            "cultivar"))
lc <- lc %>%
  group_by(mutant,leaf.number,cultivar)%>%
  mutate(lesion.percent = lesion.surface/leaf.surface)

library(dplyr)
lc <- lc %>% 
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





