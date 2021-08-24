#Masters thesis data analysis 

# library(dplyr)
# library(readr)

# install.packages("ggforce")
# install.packages("colorspace")
# install.packages("ggtext")
# install.packages("ggsci")
#install.packages("forcats")
# 
# #### Data import #### 
#all this will not be run because the code is too heavy to run each time. 
#Instead the resulting df is saved as an r file and loaded each time
# ## Compilation des fichiers csv All_lesions
# list_files = list.files(full.names = TRUE) # making list of files
# a=0
# b = "All_lesions.csv"
# 
# #alternative code if files aren't in the same working directory
# # list_files = list.files(path="D:/Documents/folder/folder/folder", full.names = TRUE)
# # a=0
# # b = "All_lesions.csv"
# 
# for( x in list_files) {
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


# #### Dataframe clean up ####
#install.packages("patchwork")
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
load("df2.RData")

df2[ , c('genotype', 'mutant_cultivar_date')] <- list(NULL) #deleting redunant columns 
colnames(df2)

rep1 <- df2[, c(16, 17, 18, 19, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)] #reordring columns and puttting in a new df
#really could've done this and previous deletions in the same step but hey 
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


#save (rep1, file = "rep1.RData") #saving just to be sure (and bc the above code is quite annoying to run)

#load rep1
load("rep1.RData")
rep1$mutant <- as.factor(rep1$mutant)

#lesion count df
library(dplyr)
lesion.count <- rep1 %>%
  group_by(mutant,leaf.number,cultivar)%>%
  tally(wt = NULL) #tallies lesions. wt=NULL makes it NOT weight the 
#lesions by number

lesion.count <- rename(lesion.count, lesion.count = n)

#percent lesion df, will join with same lesion count df
library(dplyr)
leaf.surface <- rep1 %>%
  group_by(mutant,leaf.number,cultivar)%>%
  summarise(leaf.surface = sum(leaf.surface))

lesion.surface <- rep1 %>%
  group_by(mutant,leaf.number,cultivar)%>%
  summarise(lesion.surface = sum(lesion.surface))

lesion.count <- dplyr::left_join(
  lesion.count, leaf.surface, by = c("mutant","leaf.number","cultivar"))%>%
                                   left_join(.,lesion.surface, by = c("mutant","leaf.number",
                                                  "cultivar"))
lesion.count <- lesion.count %>%
  group_by(mutant,leaf.number,cultivar)%>%
  mutate(lesion.percent = lesion.surface/leaf.surface)

library(dplyr)
lesion.count <- lesion.count %>% 
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

#save (lesion.count, file = "lc.RData")

#rep1 for each cultivar
A <- rep1[rep1$cultivar=="Aichi Asahi",] ; titre <- "Aichi Asahi"
#save (A, file = "A.RData")
B <- rep1[rep1$cultivar=="Bl1",] ; titre <- "Bl1"
#save (B, file = "B.RData")
N <- rep1[rep1$cultivar=="Nipponbare",] ; titre <- "Nipponbare"
#save (N, file = "N.RData")
K <- rep1[rep1$cultivar=="Kasalath",] ; titre <- "Kasalath"
#save (K, file = "K.RData")
K6 <- rep1[rep1$cultivar=="K60",] ; titre <- "K60"
#save (K, file = "K6.RData")
S <- rep1[rep1$cultivar=="Shin 2",] ; titre <- "Shin 2"
#save (S, file = "S.RData")
T <- rep1[rep1$cultivar=="Tsuyuake",] ; titre <- "Tsuyuake"
#save (T, file = "T.RData")

#lesion.count for each cultivar
Alc <- lesion.count[lesion.count$cultivar =="Aichi Asahi",] ; titre <- "Aichi Asahi"
Blc <- lesion.count[lesion.count$cultivar =="Bl1",] ; titre <- "Bl1"
Nlc <- lesion.count[lesion.count$cultivar =="Nipponbare",] ; titre <- "Nipponbare"
Klc <- lesion.count[lesion.count$cultivar =="Kasalath",] ; titre <- "Kasalath"
K6lc <- lesion.count[lesion.count$cultivar =="K60",] ; titre <- "K60"
Slc <- lesion.count[lesion.count$cultivar =="Shin 2",] ; titre <- "Shin 2"
Tlc <- lesion.count[lesion.count$cultivar =="Tsuyuake",] ; titre <- "Tsuyuake"

#### Aichi Asahi ####

##Lesion Number

#Anova lesion surface
a0 <- aov(lesion.surface ~ mutant, data = A) 
anova(a0)
#p-value = < 2.2e-16 ***, significant
#normality of residuals check
#distrubution/normality
par(mfrow=c(1,2))
#Histogramme
hist(a0$residuals) #histogram to check normality
# Q-Q Plot
qqnorm(a0$residuals) # qqplot to check normality 
qqline(a0$residuals)
#very anormal
shapiro.test(a0$residuals)
#too many entries to run shapiro but clearly anormal
kruskal.test(lesion.surface ~ mutant, data = A)
#p-value <2e-16, significant !

#Post-hoc 
#post-hoc test Dunn
library(FSA)
dunnTest(lesion.surface ~ mutant, data = A, method ="bh")
# #post-hoc test Wilcoxon rank sum test
#will allow us to easily add significance levels to our ggplots
library(ggpubr)
atest<- compare_means(lesion.surface ~ mutant, A, method="wilcox.test",
              ref.group = ".all.")

#Anova lesion number

a1 <- aov(lesion.count ~ mutant, data = Alc) 
anova(a1)
#p-value = <2e-16 ***, significant
#normality of residuals check
#distrubution/normality
par(mfrow=c(1,2))
#Histogramme
hist(a1$residuals)
# Q-Q Plot
qqnorm(a1$residuals)
qqline(a1$residuals)
#vvvv not normal 
shapiro.test(a1$residuals)
#1e-07, not normal
kruskal.test(lesion.count ~ mutant, data = Alc)
#p = 5e-11, significant

#Post-hoc 

#post-hoc test Dunn
library(FSA)
dunnTest(lesion.count ~ mutant, data = Alc, method ="bh")
# #post-hoc test Wilcoxon rank sum test
#will allow us to easily add significance levels to our ggplots
library(ggpubr)
atest <- compare_means(lesion.count ~ mutant, Alc, method="wilcox.test",
              ref.group = ".all.")


#Anova lesion percent
a2 <- aov(lesion.percent ~ mutant, data = Alc) 
anova(a2)
#p-value = 0.64, not significant
#normality of residuals check
#distrubution/normality
par(mfrow=c(1,2))
#Histogramme
hist(a2$residuals)
# Q-Q Plot
qqnorm(a2$residuals)
qqline(a2$residuals)
#vvvv not normal 
shapiro.test(a2$residuals)
#"Error in shapiro.test(a0$residuals) : 
#la taille de l'échantillon doit être comprise entre 3 et 5000"
#error in shapiro test but histogram does not look good so we will use kruskal wallis
kruskal.test(lesion.percent ~ mutant, data = Alc)
#p  = 2e-07, significant

#Post-hoc 

#post-hoc test Dunn
library(FSA)
dunnTest(lesion.percent ~ mutant, data = Alc, method ="bh")
# #post-hoc test Wilcoxon rank sum test
#will allow us to easily add significance levels to our ggplots
library(ggpubr)
compare_means(lesion.percent ~ mutant, Alc, method="wilcox.test",
              ref.group = ".all.")

### ####Aichi ggplots####
#ggplot lesion size
library(ggpubr)
library(ggplot2)
ap0 <- ggplot(A,aes(x=factor(mutant), y=lesion.surface, color = effector)) + 
  geom_boxplot() + geom_jitter(aes(size = lesion.surface), show.legend = TRUE) + # jitter size by lesion surface
  scale_size_continuous(range = c(0.01, 1.5))  + theme_minimal()
Alsplot <- ap0 +
  labs(title = "FR13 lesion size on Aichi Asahi", x = "Mutant", y = "Lesion size") + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
# stat_compare_means(label = "p.signif", method = "wilcox.test",
#                    ref.group = ".all.", hide.ns = TRUE) 

Alsplot
#ggplot lesion number
library(ggpubr)
alc0 <- ggplot(Alc, aes(x=factor(mutant), y=lesion.count, color = effector)) + 
  geom_boxplot() + geom_jitter() + theme_minimal()
Alcplot <- alc0 +
  labs(title = "FR13 lesion number on Aichi Asahi", x = "Mutant", y = "Lesion number") + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) +
  stat_compare_means(label = "p.signif", method = "wilcox.test",
                    ref.group = ".all.", hide.ns = TRUE, label.y = c(50, 300, 300,350, 50, 50, 60)) 
# #ggplot lesion percent
# library(ggpubr)
# alp0 <- ggplot(Alc, aes(x=factor(mutant), y=lesion.percent, color = effector)) + 
#   geom_boxplot() + geom_jitter() + theme_minimal()
# Alpplot <- alp0 +
#   labs(title = "Lesion percent on Aichi Asahi", x = "Mutant", y = "Lesion percent") + 
#   theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
# # stat_compare_means(label = "p.signif", method = "wilcox.test",
# #                    ref.group = ".all.", hide.ns = TRUE) 
# Alpplot

library(patchwork)
Aplot <- (Alsplot + Alcplot)

ggsave(filename = "A.png", plot = Aplot, device = "png", height = 8, width = 15,
                      units = "in", dpi = 500)

### #### Bl1 ####

#ggplot lesion size
BF <- dplyr::filter(B, isolate == "FR13")
BG <- dplyr::filter(B, isolate == "Guy11")
library(ggpubr)
bf.ls.p <- ggplot(BF,aes(x=factor(mutant), y=lesion.surface, color = effector)) + geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface), show.legend = TRUE) + # jitter size by lesion surface
  scale_size_continuous(range = c(0.01, 2))  + theme_minimal() + scale_x_discrete(limits=rev)
BF.ls.plot <- bf.ls.p +
  labs(title = "FR13 lesion surface", x = "Mutant", y = "Lesion surface (in pixels)") +
  # theme(axis.title.y = element_blank(), axis.title.x = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
# stat_compare_means(label = "p.signif", method = "wilcox.test",
#                    ref.group = ".all.", hide.ns = TRUE) 
bg.ls.p <- ggplot(BG,aes(x=factor(mutant), y=lesion.surface, color = effector)) + geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface), show.legend = TRUE) + # jitter size by lesion surface
  scale_size_continuous(range = c(0.01, 2))  + theme_minimal() + scale_x_discrete(limits=rev)
BG.ls.plot <- bg.ls.p +
  labs(title = "Guy11 lesion surface", x = "Mutant", y = "Lesion surface (in pixels)") +   
  #theme(axis.title.y = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0))   
library(patchwork)
Blsplot <- (BF.ls.plot + BG.ls.plot) + plot_annotation(
  title = 'Lesion surface and number on Bl1')
Blsplot
ggsave(filename = "Bls.png", plot = Blsplot, device = "png", height = 10, width = 15,
       units = "in", dpi = 500)
#ggplot lesion number
BFlc <- dplyr::filter(Blc, isolate == "FR13")
BGlc <- dplyr::filter(Blc, isolate == "Guy11")
library(ggpubr)
bf.lc0 <- ggplot(BFlc, aes(x=factor(mutant), y=lesion.count, color = effector)) + 
  geom_boxplot() + geom_jitter() + theme_minimal()
BFlc <- bf.lc0 +
  labs(title = "FR13 lesion number", x = "Mutant", y = "Lesion number") + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) + 
 stat_compare_means(label = "p.signif", method = "wilcox.test",
                    ref.group = ".all.", hide.ns = TRUE, label.y = c(40, 60, 200, 30,30)) 
bg.lc0 <- ggplot(BGlc, aes(x=factor(mutant), y=lesion.count, color = effector)) + 
  geom_boxplot() + geom_jitter() + theme_minimal()
Bglc <- bg.lc0 +
  labs(title = "Guy11 lesion number", x = "Mutant", y = "Lesion number") + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) +
 stat_compare_means(label = "p.signif", method = "wilcox.test",
                    ref.group = ".all.", hide.ns = TRUE, label.y = c(125, 170, 170, 30,30)) 
Bglc
library(patchwork)
Bplot <- ((BF.ls.plot + BG.ls.plot)/(BFlc+Bglc)) + plot_annotation(
  title = 'Lesion surface and number on Bl1')
Bplot

ggsave(filename = "B.png", plot = Bplot, device = "png", height = 10, width = 10,
       units = "in", dpi = 500)

### #### Nipponbare####
#ggplot lesion size

NF <- dplyr::filter(N, isolate == "FR13")
NG <- dplyr::filter(N, isolate == "Guy11")
library(ggpubr)
Nf.ls.p <- ggplot(NF,aes(x=factor(mutant), y=lesion.surface, color = effector)) + geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface), show.legend = TRUE) + # jitter size by lesion surface
  scale_size_continuous(range = c(0.01, 2.5))  + theme_minimal()  + scale_x_discrete(limits=rev)
NF.ls.plot <- Nf.ls.p +
  labs(title = "FR13 lesion surface", x = "Mutant", y = "Lesion surface (in pixels)") +
  # theme(axis.title.y = element_blank(), axis.title.x = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
# stat_compare_means(label = "p.signif", method = "wilcox.test",
#                    ref.group = ".all.", hide.ns = TRUE) 
Ng.ls.p <- ggplot(NG,aes(x=factor(mutant), y=lesion.surface, color = effector)) + geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface), show.legend = TRUE) + # jitter size by lesion surface
  scale_size_continuous(range = c(0.01, 2.5))  + theme_minimal() + scale_x_discrete(limits=rev)
NG.ls.plot <- Ng.ls.p +
  labs(title = "Guy11 lesion surface", x = "Mutant", y = "Lesion surface (in pixels)") +   
  #theme(axis.title.y = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0))   
library(patchwork)
Nlsplot <- (NF.ls.plot + NG.ls.plot) + plot_annotation(
  title = 'Lesion surface and number on Nipponbare')
Nlsplot
ggsave(filename = "Nls.png", plot = Nlsplot, device = "png", height = 10, width = 15,
       units = "in", dpi = 500)


#ggplot lesion number
NFlc <- dplyr::filter(Nlc, isolate == "FR13")
NGlc <- dplyr::filter(Nlc, isolate == "Guy11")
library(ggpubr)
Nf.lc0 <- ggplot(NFlc, aes(x=factor(mutant), y=lesion.count, color = effector)) + 
  geom_boxplot() + geom_jitter() + theme_minimal()
NFlc <- Nf.lc0 +
  labs(title = "FR13 lesion number", x = "Mutant", y = "Lesion number") + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
# stat_compare_means(label = "p.signif", method = "wilcox.test",
#                    ref.group = ".all.", hide.ns = TRUE) 
Ng.lc0 <- ggplot(NGlc, aes(x=factor(mutant), y=lesion.count, color = effector)) + 
  geom_boxplot() + geom_jitter() + theme_minimal()
Nglc <- Ng.lc0 +
  labs(title = "Guy11 lesion number", x = "Mutant", y = "Lesion number") + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
# stat_compare_means(label = "p.signif", method = "wilcox.test",
#                    ref.group = ".all.", hide.ns = TRUE) 

library(patchwork)
Nplot <- ((NF.ls.plot + NG.ls.plot)/(NFlc+Nglc)) + plot_annotation(
  title = 'Lesion surface and number on Nipponbare')
Nplot

ggsave(filename = "N.png", plot = Nplot, device = "png", height = 10, width = 10,
       units = "in", dpi = 500)

### #### Kasalath####
#ggplot lesion size
KF <- dplyr::filter(K, isolate == "FR13")
KG <- dplyr::filter(K, isolate == "Guy11")
library(ggpubr)
Kf.ls.p <- ggplot(KF,aes(x=factor(mutant), y=lesion.surface, color = effector)) + geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface), show.legend = TRUE) + # jitter size by lesion surface
  scale_size_continuous(range = c(0.01, 1.5))  + theme_minimal() 
KF.ls.plot <- Kf.ls.p +
  labs(title = "FR13 lesion surface", x = "Mutant", y = "Lesion surface (in pixels") +
  # theme(axis.title.y = element_blank(), axis.title.x = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
# stat_compare_means(label = "p.signif", method = "wilcox.test",
#                    ref.group = ".all.", hide.ns = TRUE) 
Kg.ls.p <- ggplot(KG,aes(x=factor(mutant), y=lesion.surface, color = effector)) + geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface), show.legend = TRUE) + # jitter size by lesion surface
  scale_size_continuous(range = c(0.01, 1.5))  + theme_minimal()
KG.ls.plot <- Kg.ls.p +
  labs(title = "Guy11 lesion surface", x = "Mutant", y = "Lesion surface (in pixels") +   
  #theme(axis.title.y = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0))   

#ggplot lesion number
KFlc <- dplyr::filter(Klc, isolate == "FR13")
KGlc <- dplyr::filter(Klc, isolate == "Guy11")
library(ggpubr)
Kf.lc0 <- ggplot(KFlc, aes(x=factor(mutant), y=lesion.count, color = effector)) + 
  geom_boxplot() + geom_jitter() + theme_minimal()
KFlc <- Kf.lc0 +
  labs(title = "FR13 lesion number", x = "Mutant", y = "Lesion number") + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
# stat_compare_means(label = "p.signif", method = "wilcox.test",
#                    ref.group = ".all.", hide.ns = TRUE) 
Kg.lc0 <- ggplot(KGlc, aes(x=factor(mutant), y=lesion.count, color = effector)) + 
  geom_boxplot() + geom_jitter() + theme_minimal()
Kglc <- Kg.lc0 +
  labs(title = "Guy11 lesion number", x = "Mutant", y = "Lesion number") + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
# stat_compare_means(label = "p.signif", method = "wilcox.test",
#                    ref.group = ".all.", hide.ns = TRUE) 

library(patchwork)
Kplot <- ((KF.ls.plot + KG.ls.plot)/(KFlc+Kglc)) + plot_annotation(
  title = 'Lesion surface and number on Kasalath')
Kplot

ggsave(filename = "K.png", plot = Kplot, device = "png", height = 10, width = 10,
       units = "in", dpi = 500)

### #### K6####

#ggplot lesion size
K6F <- dplyr::filter(K6, isolate == "FR13")
K6G <- dplyr::filter(K6, isolate == "Guy11")
library(ggpubr)
K6f.ls.p <- ggplot(K6F,aes(x=factor(mutant), y=lesion.surface, color = effector)) + geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface), show.legend = TRUE) + # jitter size by lesion surface
  scale_size_continuous(range = c(0.01, 1.5))  + theme_minimal() 
K6F.ls.plot <- K6f.ls.p +
  labs(title = "FR13 lesion surface", x = "Mutant", y = "Lesion surface (in pixels") +
  # theme(axis.title.y = element_blank(), axis.title.x = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
# stat_compare_means(label = "p.signif", method = "wilcox.test",
#                    ref.group = ".all.", hide.ns = TRUE) 
K6g.ls.p <- ggplot(K6G,aes(x=factor(mutant), y=lesion.surface, color = effector)) + geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface), show.legend = TRUE) + # jitter size by lesion surface
  scale_size_continuous(range = c(0.01, 1.5))  + theme_minimal()
K6G.ls.plot <- K6g.ls.p +
  labs(title = "Guy11 lesion surface", x = "Mutant", y = "Lesion surface (in pixels") +   
  #theme(axis.title.y = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0))   

#ggplot lesion number
K6Flc <- dplyr::filter(K6lc, isolate == "FR13")
K6Glc <- dplyr::filter(K6lc, isolate == "Guy11")
library(ggpubr)
K6f.lc0 <- ggplot(K6Flc, aes(x=factor(mutant), y=lesion.count, color = effector)) + 
  geom_boxplot() + geom_jitter() + theme_minimal()
K6Flc <- K6f.lc0 +
  labs(title = "FR13 lesion number", x = "Mutant", y = "Lesion number") + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
# stat_compare_means(label = "p.signif", method = "wilcox.test",
#                    ref.group = ".all.", hide.ns = TRUE) 
K6g.lc0 <- ggplot(K6Glc, aes(x=factor(mutant), y=lesion.count, color = effector)) + 
  geom_boxplot() + geom_jitter() + theme_minimal()
K6glc <- K6g.lc0 +
  labs(title = "Guy11 lesion number", x = "Mutant", y = "Lesion number") + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
# stat_compare_means(label = "p.signif", method = "wilcox.test",
#                    ref.group = ".all.", hide.ns = TRUE) 

library(patchwork)
K6plot <- ((K6F.ls.plot + K6G.ls.plot)/(K6Flc+K6glc)) + plot_annotation(
  title = 'Lesion surface and number on K60')
K6plot

ggsave(filename = "K6.png", plot = K6plot, device = "png", height = 10, width = 10,
       units = "in", dpi = 500)

### #### Shin 2####
#ggplot lesion size
SF <- dplyr::filter(S, isolate == "FR13")
SG <- dplyr::filter(S, isolate == "Guy11")
library(ggpubr)
Sf.ls.p <- ggplot(SF,aes(x=factor(mutant), y=lesion.surface, color = effector)) + geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface), show.legend = TRUE) + # jitter size by lesion surface
  scale_size_continuous(range = c(0.01, 1.5))  + theme_minimal() 
SF.ls.plot <- Sf.ls.p +
  labs(title = "FR13 lesion surface", x = "Mutant", y = "Lesion surface (in pixels") +
  # theme(axis.title.y = element_blank(), axis.title.x = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
# stat_compare_means(label = "p.signif", method = "wilcox.test",
#                    ref.group = ".all.", hide.ns = TRUE) 
Sg.ls.p <- ggplot(SG,aes(x=factor(mutant), y=lesion.surface, color = effector)) + geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface), show.legend = TRUE) + # jitter size by lesion surface
  scale_size_continuous(range = c(0.01, 1.5))  + theme_minimal()
SG.ls.plot <- Sg.ls.p +
  labs(title = "Guy11 lesion surface", x = "Mutant", y = "Lesion surface (in pixels") +   
  #theme(axis.title.y = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0))   

#ggplot lesion number
SFlc <- dplyr::filter(Slc, isolate == "FR13")
SGlc <- dplyr::filter(Slc, isolate == "Guy11")
library(ggpubr)
Sf.lc0 <- ggplot(SFlc, aes(x=factor(mutant), y=lesion.count, color = effector)) + 
  geom_boxplot() + geom_jitter() + theme_minimal()
SFlc <- Sf.lc0 +
  labs(title = "FR13 lesion number", x = "Mutant", y = "Lesion number") + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
# stat_compare_means(label = "p.signif", method = "wilcox.test",
#                    ref.group = ".all.", hide.ns = TRUE) 
Sg.lc0 <- ggplot(SGlc, aes(x=factor(mutant), y=lesion.count, color = effector)) + 
  geom_boxplot() + geom_jitter() + theme_minimal()
Sglc <- Sg.lc0 +
  labs(title = "Guy11 lesion number", x = "Mutant", y = "Lesion number") + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
# stat_compare_means(label = "p.signif", method = "wilcox.test",
#                    ref.group = ".all.", hide.ns = TRUE) 

library(patchwork)
Splot <- ((SF.ls.plot + SG.ls.plot)/(SFlc+Sglc)) + plot_annotation(
  title = 'Lesion surface and number on Shin 2')
Splot

ggsave(filename = "S.png", plot = Splot, device = "png", height = 10, width = 10,
       units = "in", dpi = 500)



### #### Tsuyuake####
t1 <- aov(lesion.surface ~ mutant, data = T) 
anova(t1)
#p-value = 4.8e-12 ***, significant
#normality of residuals check
#distrubution/normality
par(mfrow=c(1,2))
#Histogramme
hist(t1$residuals)
# Q-Q Plot
qqnorm(t1$residuals)
qqline(t1$residuals)
#vvvv not normal 
shapiro.test(t1$residuals)
#"Error in shapiro.test(a0$residuals) : 
#la taille de l'échantillon doit être comprise entre 3 et 5000"
#error in shapiro test but histogram does not look good so we will use kruskal wallis
kruskal.test(lesion.surface ~ mutant, data = T)
#p < 2.2e-16, significant


#ggplot lesion size
library(ggpubr)
library(ggplot2)
tp0 <- ggplot(T,aes(x=factor(mutant), y=lesion.surface, color = effector)) + 
  geom_jitter(aes(size = lesion.surface), show.legend = TRUE) + # jitter size by lesion surface
  scale_size_continuous(range = c(0.01, 1.5))  + theme_minimal()
Tlsplot <- tp0 +
  labs(title = "Lesion surface on Tsuyuake", x = "Mutant", y = "Lesion surface") + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
# stat_compare_means(label = "p.signif", method = "wilcox.test",
#                    ref.group = ".all.", hide.ns = TRUE) 
#ggplot lesion number
library(ggpubr)
tlc0 <- ggplot(Tlc, aes(x=factor(mutant), y=lesion.count, color = effector)) + 
  geom_boxplot() + geom_jitter() + theme_minimal()
Tlcplot <- tlc0 +
  labs(title = "Lesion number on Tsuyuake", x = "Mutant", y = "Lesion number") + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
# stat_compare_means(label = "p.signif", method = "wilcox.test",
#                    ref.group = ".all.", hide.ns = TRUE) 

#ggplot lesion percent
library(ggpubr)
tlp0 <- ggplot(Tlc, aes(x=factor(mutant), y=lesion.percent, color = effector)) + 
  geom_boxplot() + geom_jitter() + theme_minimal()
Tlpplot <- tlp0 +
  labs(title = "Lesion percent on Tsuyuake", x = "Mutant", y = "Lesion percent") + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
# stat_compare_means(label = "p.signif", method = "wilcox.test",
#                    ref.group = ".all.", hide.ns = TRUE) 
Tlpplot

library(patchwork)
Tplot <- (Tlsplot + Tlcplot/Tlpplot)

ggsave("T.png", Tplot)



#### jfkdla####


library(ggplot2)
library(ggpubr)
AAplot <- ggplot(AA) + geom_boxplot(aes(x=factor(mutant), y=lesion.number)) +
  labs(title = "Aichi Asahi", x = "mutant", y = "lesion number") 
AAplot



ggplot(A, aes(x=factor(mutant),y=lesion.surface,fill=effector)) +
  theme_classic() + coord_flip() + theme(legend.position="none",axis.text = element_text(face="bold",size=10)) +
  scale_x_discrete(limits=levels(A$mutant)) +
  labs(x="", y="lesion surface (pixels)") +
  geom_boxplot() + geom_point(cex=1.5,pch=1) +
  stat_summary(fun.y=mean, geom="point", shape=20, size=3, color="red", fill="Black")


blep <- ggplot(A, aes(x=mutant, y=lesion.surface, fill=effector)) +
  geom_boxplot(alpha=0.4,width=0.5,lwd=1,outlier.shape=NA)+
  geom_jitter(position=position_jitter(0.1), cex=1, pch=1)+
  stat_summary(fun.y=mean, geom="point", shape=20, size=3, color="Red", fill="Black")
blep + theme_classic() +
  theme(legend.position="none")+ 
  ggtitle("Aichi Asahi")+
  ylab("lesion surface (pixels)") +
  theme(plot.title = element_text(color = "black", size = 17, face = ,hjust =0.5),
        #plot.background = element_rect(colour = 'Black', size = 3),
        axis.text.x = element_text(color = "black", size = 15, face = "bold",angle = 90, hjust = 1 ),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color = "black", size = 15, face = "bold",angle = 90, hjust = ))





