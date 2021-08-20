#Masters thesis data analysis 

# library(dplyr)
# library(readr)
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
                         '1' = "Pia",
                         '2' = "Pia",
                         '3' = 'Pia',
                         '4' = 'Pia',
                         '5' = 'PiB',
                         '6'= 'PiB',
                         '7' = 'PiB',
                         '8' = 'PiB',
                         '9' = 'Pik',
                         '10' = 'Pik',
                         '11' = 'Pik',
                         '12' = 'Pik',
                         '13' = 'pCR17',
                         '14' = 'WT',
                         '15' = 'PiB',
                         '16' = 'PiB',
                         '17' = 'PiB',
                         '18' = 'pCR17',
                         '19' = 'Pik',
                         '20' = 'Pik',
                         '21' = 'Pik',
                         '22' = 'Pik',
                         '23' = 'pCR17',
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
                         '13' = 'FR13 pCR17',
                         '14' = 'FR13 WT',
                         '15' = 'GB1',
                         '16' = 'GB2',
                         '17' = 'GB3',
                         '18' = 'Guy11 pCR17_1',
                         '19' = 'GK1',
                         '20' = 'GK2',
                         '21' = 'GK3',
                         '22' = 'GK ectopic',
                         '23' = 'Guy11 pCR17_2',
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

colnames(rep1)

rep1 <- rep1[, c(1, 19, 20, 2, 3, 21, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16)] 

colnames(rep1)

A <- rep1[rep1$cultivar=="Aichi Asahi",] ; titre <- "Aichi Asahi"
B <- rep1[rep1$cultivar=="Bl1",] ; titre <- "Bl1"
N <- rep1[rep1$cultivar=="Nipponbare",] ; titre <- "Nipponbare"
K <- rep1[rep1$cultivar=="Kasalath",] ; titre <- "Kasalath"
K6 <- rep1[rep1$cultivar=="K60",] ; titre <- "K60"
S <- rep1[rep1$cultivar=="Shin 2",] ; titre <- "Shin 2"
T <- rep1[rep1$cultivar=="Tsuyuake",] ; titre <- "Tsuyuake"

#### Aichi Asahi ####

##Lesion Number

#Anova lesion number
a0 <- aov(lesion.number ~ mutant, data = AA) 
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
#fairly normal! 
shapiro.test(a0$residuals)
#"Error in shapiro.test(a0$residuals) : 
#la taille de l'échantillon doit être comprise entre 3 et 5000"
#error in shapiro test but histogram looks good so we will stick with anova

#Post-hoc lesion number
# #post-hoc test Dunn
# library(FSA)
# dunnTest(percent.lesions ~ mutant, data = AA, method ="bh")
# # #post-hoc test Wilcoxon rank sum test
# #will allow us to easily add significance levels to our ggplots
# library(ggpubr)
# a.w.test <- compare_means(percent.lesions ~ mutant, rep1, method="wilcox.test",
#                           ref.group = ".all.")
# a.w.test <- dplyr::rename(a.w.test, 
#                           mutant = "group2",
#                           Ap = "p.signif") #rename columns
# library(dplyr)
# A.p.sig <- a.w.test %>%
#   select(mutant, Ap)
# a.w.test <- subset(a.w.test, p.signif!="ns") #remove non sig rows
# a.w.test <- left_join(AA, a.w.test, by = "mutant") #joins with original data
# colnames(a.w.test) #check column names
# a.w.test <- a.w.test[, -c(2:5,8:11)] #remove unecessary columns (aka "leaf.number"     "leaf.surface"    "lesion.nb"      
# # [5] "lesion.surface")
# a.w.test <-na.omit(a.w.test) #removes na 
# a.w.test #view table
#will change all this bc its normal 

#ggplot lesion number

library(ggpubr)
ap0 <- ggplot(A,aes(x=factor(mutant), y=lesion.number, color = effector)) + geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface), show.legend = TRUE) + # jitter size by lesion surface
  scale_size_continuous(range = c(0.01, 1.5))  + theme_minimal()
Aplot <- ap0 +
  labs(title = "Lesion number on Aichi Asahi", x = "Mutant", y = "Lesion number") + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
  # stat_compare_means(label = "p.signif", method = "wilcox.test",
  #                    ref.group = ".all.", hide.ns = TRUE) 
Aplot
ggsave("A.png",Aplot)

### #### Bl1 ####

BF <- dplyr::filter(B, isolate == "FR13")
BG <- dplyr::filter(B, isolate == "Guy11")


library(ggpubr)
bfp0 <- ggplot(BF,aes(x=factor(mutant), y=lesion.number, color = effector)) + geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface), show.legend = TRUE) + # jitter size by lesion surface
  scale_size_continuous(range = c(0.01, 1.5))  + theme_minimal() 
BFplot <- bfp0 +
  labs(x = "Mutant", y = "Lesion number") +
  # theme(axis.title.y = element_blank(), axis.title.x = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
# stat_compare_means(label = "p.signif", method = "wilcox.test",
#                    ref.group = ".all.", hide.ns = TRUE) 
bgp0 <- ggplot(BG,aes(x=factor(mutant), y=lesion.number, color = effector)) + geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface), show.legend = TRUE) + # jitter size by lesion surface
  scale_size_continuous(range = c(0.01, 1.5))  + theme_minimal()
BGplot <- bgp0 +
  labs(x = "Mutant", y = "Lesion number") +   
  # theme(axis.title.y = element_blank(), axis.title.x = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0))   


library(patchwork)
bfgplot <- (BFplot + BGplot)
Bplot <- nfgplot + plot_annotation(
  title = 'Lesion number on Bl1',
)

ggsave(filename = "B.png", plot = Bplot, device = "png", height = 5, width = 10,
       units = "in", dpi = 500)

### #### Nipponbare####

NF <- dplyr::filter(N, isolate == "FR13")
NG <- dplyr::filter(N, isolate == "Guy11")


library(ggpubr)
nfp0 <- ggplot(NF,aes(x=factor(mutant), y=lesion.number, color = effector)) + geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface), show.legend = TRUE) + # jitter size by lesion surface
  scale_size_continuous(range = c(0.01, 1.5))  + theme_minimal() 
NFplot <- nfp0 +
  labs(x = "Mutant", y = "Lesion number") +
  # theme(axis.title.y = element_blank(), axis.title.x = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
# stat_compare_means(label = "p.signif", method = "wilcox.test",
#                    ref.group = ".all.", hide.ns = TRUE) 
ngp0 <- ggplot(NG,aes(x=factor(mutant), y=lesion.number, color = effector)) + geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface), show.legend = TRUE) + # jitter size by lesion surface
  scale_size_continuous(range = c(0.01, 1.5))  + theme_minimal()
NGplot <- ngp0 +
  labs(x = "Mutant", y = "Lesion number") +   
  # theme(axis.title.y = element_blank(), axis.title.x = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0))   


library(patchwork)
nfgplot <- (NFplot + NGplot)
Nplot <- nfgplot + plot_annotation(
  title = 'Lesion number on Nipponbare',
)

ggsave(filename = "N.png", plot = Nplot, device = "png", height = 5, width = 10,
       units = "in", dpi = 500)


### #### Nipponbare####

NF <- dplyr::filter(N, isolate == "FR13")
NG <- dplyr::filter(N, isolate == "Guy11")


library(ggpubr)
nfp0 <- ggplot(NF,aes(x=factor(mutant), y=lesion.number, color = effector)) + geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface), show.legend = TRUE) + # jitter size by lesion surface
  scale_size_continuous(range = c(0.01, 1.5))  + theme_minimal() 
NFplot <- nfp0 +
  labs(x = "Mutant", y = "Lesion number") +
  # theme(axis.title.y = element_blank(), axis.title.x = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
# stat_compare_means(label = "p.signif", method = "wilcox.test",
#                    ref.group = ".all.", hide.ns = TRUE) 
ngp0 <- ggplot(NG,aes(x=factor(mutant), y=lesion.number, color = effector)) + geom_boxplot() + 
    geom_jitter(aes(size = lesion.surface), show.legend = TRUE) + # jitter size by lesion surface
    scale_size_continuous(range = c(0.01, 1.5))  + theme_minimal()
NGplot <- ngp0 +
    labs(x = "Mutant", y = "Lesion number") +   
  # theme(axis.title.y = element_blank(), axis.title.x = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0))   


library(patchwork)
nfgplot <- (NFplot + NGplot)
Nplot <- nfgplot + plot_annotation(
  title = 'Lesion number on Nipponbare',
)

ggsave(filename = "N.png", plot = Nplot, device = "png", height = 5, width = 10,
       units = "in", dpi = 500)

### #### Kasalath####

KF <- dplyr::filter(K, isolate == "FR13")
KG <- dplyr::filter(K, isolate == "Guy11")


library(ggpubr)
Kfp0 <- ggplot(KF,aes(x=factor(mutant), y=lesion.number, color = effector)) + geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface), show.legend = TRUE) + # jitter size by lesion surface
  scale_size_continuous(range = c(0.01, 1.5))  + theme_minimal() 
KFplot <- Kfp0 +
  labs(x = "Mutant", y = "Lesion number") +
  # theme(axis.title.y = element_blank(), axis.title.x = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
# stat_compare_means(label = "p.signif", method = "wilcox.test",
#                    ref.group = ".all.", hide.ns = TRUE) 
Kgp0 <- ggplot(KG,aes(x=factor(mutant), y=lesion.number, color = effector)) + geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface), show.legend = TRUE) + # jitter size by lesion surface
  scale_size_continuous(range = c(0.01, 1.5))  + theme_minimal()
KGplot <- Kgp0 +
  labs(x = "Mutant", y = "Lesion number") +   
  # theme(axis.title.y = element_blank(), axis.title.x = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0))   


library(patchwork)
Kfgplot <- (KFplot + KGplot)
Kplot <- Kfgplot + plot_annotation(
  title = 'Lesion number on Kasalath',
)

ggsave(filename = "K.png", plot = Kplot, device = "png", height = 5, width = 10,
       units = "in", dpi = 500)


### #### K6####

K6F <- dplyr::filter(K6, isolate == "FR13")
K6G <- dplyr::filter(K6, isolate == "Guy11")


library(ggpubr)
K6fp0 <- ggplot(K6F,aes(x=factor(mutant), y=lesion.number, color = effector)) + geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface), show.legend = TRUE) + # jitter size by lesion surface
  scale_size_continuous(range = c(0.01, 1.5))  + theme_minimal() 
K6Fplot <- K6fp0 +
  labs(x = "Mutant", y = "Lesion number") +
  # theme(axis.title.y = element_blank(), axis.title.x = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
# stat_compare_means(label = "p.signif", method = "wilcox.test",
#                    ref.group = ".all.", hide.ns = TRUE) 
K6gp0 <- ggplot(K6G,aes(x=factor(mutant), y=lesion.number, color = effector)) + geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface), show.legend = TRUE) + # jitter size by lesion surface
  scale_size_continuous(range = c(0.01, 1.5))  + theme_minimal()
K6Gplot <- K6gp0 +
  labs(x = "Mutant", y = "Lesion number") +   
  # theme(axis.title.y = element_blank(), axis.title.x = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0))   


library(patchwork)
K6fgplot <- (K6Fplot + K6Gplot)
K6plot <- K6fgplot + plot_annotation(
  title = 'Lesion number on K60',
)

ggsave(filename = "K6.png", plot = K6plot, device = "png", height = 5, width = 10,
       units = "in", dpi = 500)


### #### Shin 2####

SF <- dplyr::filter(S, isolate == "FR13")
SG <- dplyr::filter(S, isolate == "Guy11")


library(ggpubr)
Sfp0 <- ggplot(SF,aes(x=factor(mutant), y=lesion.number, color = effector)) + geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface), show.legend = TRUE) + # jitter size by lesion surface
  scale_size_continuous(range = c(0.01, 1.5))  + theme_minimal() 
SFplot <- Sfp0 +
  labs(x = "Mutant", y = "Lesion number") +
  # theme(axis.title.y = element_blank(), axis.title.x = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
# stat_compare_means(label = "p.signif", method = "wilcox.test",
#                    ref.group = ".all.", hide.ns = TRUE) 
Sgp0 <- ggplot(SG,aes(x=factor(mutant), y=lesion.number, color = effector)) + geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface), show.legend = TRUE) + # jitter size by lesion surface
  scale_size_continuous(range = c(0.01, 1.5))  + theme_minimal()
SGplot <- Sgp0 +
  labs(x = "Mutant", y = "Lesion number") +   
  # theme(axis.title.y = element_blank(), axis.title.x = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0))   


library(patchwork)
Sfgplot <- (SFplot + SGplot)
Splot <- Sfgplot + plot_annotation(
  title = 'Lesion number on Shin 2',
)

ggsave(filename = "S.png", plot = Splot, device = "png", height = 5, width = 10,
       units = "in", dpi = 500)

### #### K6####

K6F <- dplyr::filter(K6, isolate == "FR13")
K6G <- dplyr::filter(K6, isolate == "Guy11")


library(ggpubr)
K6fp0 <- ggplot(K6F,aes(x=factor(mutant), y=lesion.number, color = effector)) + geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface), show.legend = TRUE) + # jitter size by lesion surface
  scale_size_continuous(range = c(0.01, 1.5))  + theme_minimal() 
K6Fplot <- K6fp0 +
  labs(x = "Mutant", y = "Lesion number") +
  # theme(axis.title.y = element_blank(), axis.title.x = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
# stat_compare_means(label = "p.signif", method = "wilcox.test",
#                    ref.group = ".all.", hide.ns = TRUE) 
K6gp0 <- ggplot(K6G,aes(x=factor(mutant), y=lesion.number, color = effector)) + geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface), show.legend = TRUE) + # jitter size by lesion surface
  scale_size_continuous(range = c(0.01, 1.5))  + theme_minimal()
K6Gplot <- K6gp0 +
  labs(x = "Mutant", y = "Lesion number") +   
  # theme(axis.title.y = element_blank(), axis.title.x = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0))   


library(patchwork)
K6fgplot <- (K6Fplot + K6Gplot)
K6plot <- K6fgplot + plot_annotation(
  title = 'Lesion number on K60',
)

ggsave(filename = "K6.png", plot = K6plot, device = "png", height = 5, width = 10,
       units = "in", dpi = 500)

### #### Tsuyuake####


library(ggpubr)
Tp1 <- ggplot(T,aes(x=factor(mutant), y=lesion.number, color = effector)) + geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface), show.legend = TRUE) + # jitter size by lesion surface
  scale_size_continuous(range = c(0.01, 1.5))  + theme_minimal()
Tplot <- Tp1 +
  labs(title = "Lesion number on Tsuyuake", x = "Mutant", y = "Lesion number") + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
# stat_compare_means(label = "p.signif", method = "wilcox.test",
#                    ref.group = ".all.", hide.ns = TRUE) 
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



#Anova lesion surface

a1 <- aov(lesion.surface ~ mutant, data = AA) 
anova(a1)
#p-value = 8.58e-05 ***, significant
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
#"Error in shapiro.test(a0$residuals) : 
#la taille de l'échantillon doit être comprise entre 3 et 5000"
#error in shapiro test but histogram does not look good so we will use kruskal wallis
kruskal.test(lesion.surface ~ mutant, data = AA)
#p < 2.2e-16, significant

#Post-hoc 

#post-hoc test Dunn
library(FSA)
dunnTest(lesion.surface ~ mutant, data = AA, method ="bh")
# #post-hoc test Wilcoxon rank sum test
#will allow us to easily add significance levels to our ggplots
library(ggpubr)
compare_means(lesion.surface ~ mutant, rep1, method="wilcox.test",
                          ref.group = ".all.")

