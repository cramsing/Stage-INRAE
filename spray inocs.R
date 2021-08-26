#Spray inoculation

#load replicates
load("rep1.RData")
load("rep2.RData")

rep1$lesion <- NULL
rep2$lesion <- NULL

spray <- rbind(rep1,rep2) #combines reps 1 and 2 
save(spray, file = "spray.RData")
spray$lesion.status <- as.factor(spray$lesion.status)
levels(spray$lesion.status) #keep
# Lesion Count ------------------------------------------------------------
#lesion count df
library(dplyr)
lc <- spray %>%
  group_by(mutant,leaf.number,cultivar,rep)%>%
  tally(wt = NULL) #tallies lesions. wt=NULL makes it NOT weight the 
#lesions by number
lc <- subset(lc, leaf.number < 9) 
lc <- rename(lc, lesion.count = n)

# #percent lesion df, will join with same lesion count df
# library(dplyr)
# leaf.s <- spray %>%
#   group_by(mutant,leaf.number,cultivar,rep)%>%
#   summarise(leaf.surface = sum(leaf.surface))
# 
# les.s <- spray %>%
#   group_by(mutant,leaf.number,cultivar,rep)%>%
#   summarise(lesion.surface = sum(lesion.surface))
# 
# lc <- dplyr::left_join(
#   lc, leaf.s, by = c("mutant","leaf.number","cultivar"))%>%
#   left_join(.,les.s, by = c("mutant","leaf.number",
#                             "cultivar"))
# lc <- lc %>%
#   group_by(mutant,leaf.number,cultivar)%>%
#   mutate(lesion.percent = lesion.surface/leaf.surface)

library(dplyr)
lc <- lc %>% 
  mutate(effector = recode(mutant,
                           'mFA1' = "AVR Pia",
                           'FA ectopic' = "AVR Pia",
                           'mFA2' = 'AVR Pia',
                           'mFA3' = 'AVR Pia',
                           'mFB1' = 'AVR PiB',
                           'mFB2'= 'AVR PiB',
                           'mFB3' = 'AVR PiB',
                           'FB ectopic' = 'AVR PiB',
                           'mFK1' = 'AVR Pik',
                           'mFK2' = 'AVR Pik',
                           'mFK3' = 'AVR Pik',
                           'FK ectopic' = 'AVR Pik',
                           'FR13 RFP' = 'RFP',
                           'FR13 WT' = 'WT',
                           'mGB1' = 'AVR PiB',
                           'mGB2' = 'AVR PiB',
                           'mGB3' = 'AVR PiB',
                           'Guy11 RFP_1' = 'RFP',
                           'mGK1' = 'AVR Pik',
                           'mGK2' = 'AVR Pik',
                           'mGK3' = 'AVR Pik',
                           'GK ectopic' = 'AVR Pik',
                           'Guy11 RFP_2' = 'RFP',
                           'Guy11 WT' = 'WT')) %>% #grouping effectors
  mutate(isolate = recode(mutant,
                          'mFA1' = "FR13",
                          'FA ectopic' = "FR13",
                          'mFA2' = 'FR13',
                          'mFA3' = 'FR13',
                          'mFB1' = 'FR13',
                          'mFB2'= 'FR13',
                          'mFB3' = 'FR13',
                          'FB ectopic' = 'FR13',
                          'mFK1' = 'FR13',
                          'mFK2' = 'FR13',
                          'mFK3' = 'FR13',
                          'FK ectopic' = 'FR13',
                          'FR13 RFP' = 'FR13',
                          'FR13 WT' = 'FR13',
                          'mGB1' = 'Guy11',
                          'mGB2' = 'Guy11',
                          'mGB3' = 'Guy11',
                          'Guy11 RFP_1' = 'Guy11',
                          'mGK1' = 'Guy11',
                          'mGK2' = 'Guy11',
                          'mGK3' = 'Guy11',
                          'GK ectopic' = 'Guy11',
                          'Guy11 RFP_2' = 'Guy11',
                          'Guy11 WT' = 'Guy11'))
save (lc, file = "lc.RData")

# Sep by cultivar ----
load("spray.RData")
load("lc.RData")

#spray for each cultivar
A <- spray[spray$cultivar=="Aichi Asahi",] ; titre <- "Aichi Asahi"
B <- spray[spray$cultivar=="Bl1",] ; titre <- "Bl1"
N <- spray[spray$cultivar=="Nipponbare",] ; titre <- "Nipponbare"
K <- spray[spray$cultivar=="Kasalath",] ; titre <- "Kasalath"
K6 <- spray[spray$cultivar=="K60",] ; titre <- "K60"
S <- spray[spray$cultivar=="Shin 2",] ; titre <- "Shin 2"
T <- spray[spray$cultivar=="Tsuyuake",] ; titre <- "Tsuyuake"

#lesion.count for each cultivar
Alc <- lc[lc$cultivar =="Aichi Asahi",] ; titre <- "Aichi Asahi"
Blc <- lc[lc$cultivar =="Bl1",] ; titre <- "Bl1"
Nlc <- lc[lc$cultivar =="Nipponbare",] ; titre <- "Nipponbare"
Klc <- lc[lc$cultivar =="Kasalath",] ; titre <- "Kasalath"
K6lc <- lc[lc$cultivar =="K60",] ; titre <- "K60"
Slc <- lc[lc$cultivar =="Shin 2",] ; titre <- "Shin 2"
Tlc <- lc[lc$cultivar =="Tsuyuake",] ; titre <- "Tsuyuake"

# Ggplots ---- 

## Aichi Asahi ----
#ggplot lesion size
library(ggpubr)
library(ggplot2)
ap0 <- ggplot(A,aes(x=factor(mutant), y=lesion.surface, color = effector)) + 
  geom_boxplot() + geom_jitter(aes(size = lesion.surface), show.legend = TRUE) + # jitter size by lesion surface
  scale_size_continuous(range = c(0.01, 1.5))  + theme_minimal() + facet_grid(~rep)
Alsplot <- ap0 +
  labs(title = "FR13 lesion size on Aichi Asahi", x = "Mutant", y = "Lesion size") + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
# stat_compare_means(label = "p.signif", method = "wilcox.test",
#                    ref.group = ".all.", hide.ns = TRUE) 

Alsplot
#ggplot lesion number
library(ggpubr)
alc0 <- ggplot(Alc, aes(x=factor(mutant), y=lesion.count, color = effector)) + 
  geom_boxplot() + geom_jitter() + theme_minimal() + facet_grid(~rep)
Alcplot <- alc0 +
  labs(title = "FR13 lesion number on Aichi Asahi", x = "Mutant", y = "Lesion number") + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
# + stat_compare_means(label = "p.signif", method = "wilcox.test",
#                    ref.group = ".all.", hide.ns = TRUE, label.y = c(50, 300, 300,350, 50, 50, 60)) 
library(patchwork)
Aplot <- (Alsplot + Alcplot)
Aplot
ggsave(filename = "A.png", plot = Aplot, device = "png", height = 8, width = 15,
       units = "in", dpi = 500)

## Bl1 ----
