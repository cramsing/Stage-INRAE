#Master Inoculation Sheet
# install.packages("xlsx")
# install.packages("agricolae")

#Spray compiling ----
#load replicates
load("rep1.RData")
load("rep2.RData")

spray <- rbind(rep1,rep2) #combines reps 1 and 2 
save(spray, file = "spray.RData")

# Lesion Count ------------------------------------------------------------
#lesion count df
library(dplyr)
lc <- spray %>%
  group_by(mutant,leaf.number,cultivar,rep)%>%
  tally(wt = NULL) #tallies lesions. wt=NULL makes it NOT weight the 
#lesions by number

lc <- rename(lc, lesion.count = n)

#percent lesion df, will join with same lesion count df
library(dplyr)
leaf.s <- spray %>%
  group_by(mutant,leaf.number,cultivar,rep)%>%
  summarise(leaf.surface = sum(leaf.surface))

les.s <- spray %>%
  group_by(mutant,leaf.number,cultivar,rep)%>%
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
save (lc, file = "lc.RData")

# Sep by cultivar ----
load("spray.RData")
load("lc.RData")
load("punch.RData")

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

#punch inoculation for each cultivar
A.p <- punch[punch$cultivar=="Aichi Asahi",] ; titre <- "Aichi Asahi"
B.p <- punch[punch$cultivar=="Bl1",] ; titre <- "Bl1"
K.p <- punch[punch$cultivar=="Kasalath",] ; titre <- "Kasalath"
K6.p <- punch[punch$cultivar=="K60",] ; titre <- "K60"
N.p <- punch[punch$cultivar=="Nipponbare",] ; titre <- "Nipponbare"
S.p <- punch[punch$cultivar=="Shin2",] ; titre <- "Shin 2"

# Stats ----

## Aichi Asahi ---- 

### Lesion surface ----
#Initial stats lesion surface
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
#p-value < 2.2e-16, significant !

#Post-hoc 
#post-hoc test Dunn
library(FSA)
dunnTest(lesion.surface ~ mutant, data = A, method ="bh")
library(rcompanion)
A.d = dunnTest(lesion.surface ~ mutant, data = A, method ="bh")
A.d=A.d$res
cldList(comparison = A.d$Comparison,
                p.value    = A.d$P.adj,
                threshold  = 0.05)
# #post-hoc test Wilcoxon rank sum test
#will allow us to easily add significance levels to our ggplots
library(ggpubr)
A.ls.w <- compare_means(lesion.surface ~ mutant, A, method="wilcox.test",
                        ref.group = ".all.")
A.ls.w <- dplyr::rename(A.ls.w, 
                        mutant = "group2",
                        A.pvalue = "p.format",
                        A.sig = "p.signif") #rename columns
colnames(A.ls.w) #check column names
A.ls.w <- A.ls.w[, c(3,6,7)] # keep only necessary columns
A.ls.w
### Lesion number ---- 
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
#kind of normal looking
shapiro.test(a1$residuals)
#5e-09, not normal
kruskal.test(lesion.count ~ mutant, data = Alc)
#p =  < 2.2e-16, significant

#Post-hoc 

#post-hoc test Dunn
library(FSA)
dunnTest(lesion.count ~ mutant, data = Alc, method ="bh")
library(rcompanion)
Alc.d = dunnTest(lesion.count ~ mutant, data = Alc, method ="bh")
Alc.d=Alc.d$res
cldList(comparison = Alc.d$Comparison,
        p.value    = Alc.d$P.adj,
        threshold  = 0.05)
# #post-hoc test Wilcoxon rank sum test
#will allow us to easily add significance levels to our ggplots
library(ggpubr)
A.lc.w <- compare_means(lesion.count ~ mutant, Alc, method="wilcox.test",
                        ref.group = ".all.") #issue with this code, will come back to it

### Punch inoculation ----
a3 <- aov(lesion.surface ~ mutant, data = A.p) 
anova(a3)
#p-value = 9.7e-15 ***, significant
#normality of residuals check
#distrubution/normality
par(mfrow=c(1,2))
#Histogramme
hist(a3$residuals)
# Q-Q Plot
qqnorm(a3$residuals)
qqline(a3$residuals)
#kind of normal looking
shapiro.test(a3$residuals)
#5e-09, not normal
kruskal.test(lesion.surface ~ mutant, data = A.p)
#p  <2e-16, significant


#Post-hoc 

#post-hoc test Dunn
library(FSA)
Apd <- dunnTest(lesion.surface ~ mutant, data = A.p, method ="bh")
Apd=Apd$res
cldList(comparison = Apd$Comparison,
        p.value    = Apd$P.adj,
        threshold  = 0.05)
# #post-hoc test Wilcoxon rank sum test
#will allow us to easily add significance levels to our ggplots
library(ggpubr)
A.p.w <- compare_means(lesion.surface ~ mutant, data = A.p, method="wilcox.test",
                       ref.group = ".all.") #issue with this code, will come back to it
A.p.w

## Bl1 ----
### Lesion surface ----
#Initial stats lesion surface
b0 <- aov(lesion.surface ~ mutant, data = B) 
anova(b0)
#p-value = < 2.2e-16 ***, significant
#normality of residuals check
#distrubution/normality
par(mfrow=c(1,2))
#Histogramme
hist(b0$residuals) #histogram to check normality
# Q-Q Plot
qqnorm(b0$residuals) # qqplot to check normality 
qqline(b0$residuals)
#very anormal
shapiro.test(b0$residuals)
#too many entries to run shapiro but clearly anormal
kruskal.test(lesion.surface ~ mutant, data = B)
#p-value <2e-16, significant !

BF <- dplyr::filter(B, isolate == "FR13")
bf0 <- aov(lesion.surface ~ mutant, data = BF) 
anova(bf0)
#p-value = < 2.2e-16 ***, significant
#normality of residuals check
#distrubution/normality
par(mfrow=c(1,2))
#Histogramme
hist(bf0$residuals) #histogram to check normality
# Q-Q Plot
qqnorm(bf0$residuals) # qqplot to check normality 
qqline(bf0$residuals)
#very anormal
shapiro.test(bf0$residuals)
#too many entries to run shapiro but clearly anormal
kruskal.test(lesion.surface ~ mutant, data = BF)

#guy11
BG <- dplyr::filter(B, isolate == "Guy11")
bg0 <- aov(lesion.surface ~ mutant, data = BG) 
anova(bg0)
#p-value = < 2.2e-16 ***, significant
#normality of residuals check
#distrubution/normality
par(mfrow=c(1,2))
#Histogramme
hist(bg0$residuals) #histogram to check normality
# Q-Q Plot
qqnorm(bg0$residuals) # qqplot to check normality 
qqline(bg0$residuals)
#very anormal
shapiro.test(bg0$residuals)
#too many entries to run shapiro but clearly anormal
kruskal.test(lesion.surface ~ mutant, data = BG)

#Post-hoc 
#post-hoc test Dunn
library(FSA)
library(rcompanion)
dunnTest(lesion.surface ~ mutant, data = B, method ="bh")
bfd <- dunnTest(lesion.surface ~ mutant, data = BF, method ="bh")
bfd=bfd$res
cldList(comparison = bfd$Comparison,
        p.value    = bfd$P.adj,
        threshold  = 0.05)
bgd <- dunnTest(lesion.surface ~ mutant, data = BG, method ="bh")
bgd=bgd$res
cldList(comparison = bgd$Comparison,
        p.value    = bgd$P.adj,
        threshold  = 0.05)
# #post-hoc test Wilcoxon rank sum test
#will allow us to easily add significance levels to our ggplots
library(ggpubr)
B.ls.w <- compare_means(lesion.surface ~ mutant, B, method="wilcox.test",
                        ref.group = ".all.")
B.ls.w <- dplyr::rename(B.ls.w, 
                        mutant = "group2",
                        B.pvalue = "p.format",
                        B.sig = "p.signif") #rename columns
colnames(B.ls.w) #check column names
B.ls.w <- B.ls.w[, c(3,6,7)] # keep only necessary columns
B.ls.w

### Lesion number ---- 
b1 <- aov(lesion.count ~ mutant, data = Blc) 
anova(b1)
#p-value = <2e-16 ***, significant
#normality of residuals check
#distrubution/normality
par(mfrow=c(1,2))
#Histogramme
hist(b1$residuals)
# Q-Q Plot
qqnorm(b1$residuals)
qqline(b1$residuals)
#kind of normal looking
shapiro.test(b1$residuals)
#5e-09, not normal
kruskal.test(lesion.count ~ mutant, data = Blc)
#p <2e-16, significant

BFlc <- dplyr::filter(Blc, isolate == "FR13")
bf1 <- aov(lesion.count ~ mutant, data = BFlc) 
anova(bf1)
#p-value = <2e-16 ***, significant
#normality of residuals check
shapiro.test(bf1$residuals)
#6.016e-13, not normal
kruskal.test(lesion.count ~ mutant, data = BFlc)
BGlc <- dplyr::filter(Blc, isolate == "Guy11")
bg1 <- aov(lesion.count ~ mutant, data = BGlc) 
anova(bg1)
#p-value = 2.385e-12 ***, significant
#normality of residuals check
shapiro.test(bg1$residuals)
#6.922e-08, not normal
kruskal.test(lesion.count ~ mutant, data = BGlc)
#p-value = 1.717e-12

#Post-hoc 
#post-hoc test Dunn
library(FSA)
dunnTest(lesion.count ~ mutant, data = Blc, method ="bh")
bfcd <- dunnTest(lesion.count ~ mutant, data = BFlc, method ="bh")
bfcd=bfcd$res
cldList(comparison = bfcd$Comparison,
        p.value    = bfcd$P.adj,
        threshold  = 0.05)
bgcd <- dunnTest(lesion.count ~ mutant, data = BGlc, method ="bh")
bgcd=bgcd$res
cldList(comparison = bgcd$Comparison,
        p.value    = bgcd$P.adj,
        threshold  = 0.05)
# #post-hoc test Wilcoxon rank sum test
#will allow us to easily add significance levels to our ggplots
library(ggpubr)
A.lc.w <- compare_means(lesion.count ~ mutant, Blc, method="wilcox.test",
                        ref.group = ".all.") #issue with this code, will come back to it

#Wilcoxon test with only rep 1
compare_means(lesion.count ~ mutant, Blc.1, method="wilcox.test",
              ref.group = ".all.") 
#Wilcoxon test with only rep 2
compare_means(lesion.count ~ mutant, Blc.2, method="wilcox.test",
              ref.group = ".all.") 

### Punch inoculation ----
b3 <- aov(lesion.surface ~ mutant, data = B.p) 
anova(b3)
#p-value = 7.6e-05 ***, significant
#normality of residuals check
#distrubution/normality
par(mfrow=c(1,2))
#Histogramme
hist(b3$residuals)
# Q-Q Plot
qqnorm(b3$residuals)
qqline(b3$residuals)
#kind of normal looking
shapiro.test(b3$residuals)
#3e-04, not normal
kruskal.test(lesion.surface ~ mutant, data = B.p)
#p = 7e-04, significant

BFp <- dplyr::filter(B.p, isolate == "FR13")
bf3 <- aov(lesion.surface ~ mutant, data = BFp) 
anova(bf3)
#p-value = 2.385e-12 ***, significant
#normality of residuals check
shapiro.test(bf3$residuals)
#0.06990 normal, anova stands!


BGp <- dplyr::filter(B.p, isolate == "Guy11")
bg3 <- aov(lesion.surface ~ mutant, data = BGp) 
anova(bg3)
#p-value = 0.00651 **, significant
#normality of residuals check
shapiro.test(bg3$residuals)
#0.02258,  normal anova stands

#post-hoc test Dunn
library(FSA)
TukeyHSD(bf3, "mutant", ordered = FALSE, conf.level = 0.95)
TukeyHSD(bg3, "mutant", ordered = FALSE, conf.level = 0.95)
# #post-hoc test Wilcoxon rank sum test
#will allow us to easily add significance levels to our ggplots
library(ggpubr)
B.p.w <- compare_means(lesion.surface ~ mutant, data = B.p, method="wilcox.test",
                       ref.group = ".all.") #issue with this code, will come back to it
B.p.w


## Kasalath ----

### Lesion surface ----
#Initial stats lesion surface
k0 <- aov(lesion.surface ~ mutant, data = K) 
anova(k0)
#p-value = < 2.2e-16 ***, significant
#normality of residuals check
#distrubution/normality
par(mfrow=c(1,2))
#Histogramme
hist(k0$residuals) #histogram to check normality
# Q-Q Plot
qqnorm(k0$residuals) # qqplot to check normality 
qqline(k0$residuals)
#very anormal
shapiro.test(k0$residuals)
#too many entries to run shapiro but clearly anormal
kruskal.test(lesion.surface ~ mutant, data = K)
#p-value <2e-16, significant !
KF <- dplyr::filter(K, isolate == "FR13")
kf0 <- aov(lesion.surface ~ mutant, data = KF) 
anova(kf0)
#p-value = < 2.2e-16 ***, significant
shapiro.test(kf0$residuals)
#too many entries to run shapiro but clearly anormal
kruskal.test(lesion.surface ~ mutant, data = KF)

KG <- dplyr::filter(K, isolate == "Guy11")
kg0 <- aov(lesion.surface ~ mutant, data = KG) 
anova(kg0)
#p-value = < 2.2e-16 ***, significant
shapiro.test(kg0$residuals)
#too many entries to run shapiro but clearly anormal
kruskal.test(lesion.surface ~ mutant, data = KG)
#Post-hoc 
#post-hoc test Dunn
library(FSA)
dunnTest(lesion.surface ~ mutant, data = K, method ="bh")
kfd <- dunnTest(lesion.surface ~ mutant, data = KF, method ="bh")
kfd=kfd$res
cldList(comparison = kfd$Comparison,
        p.value    = kfd$P.adj,
        threshold  = 0.05)
kgd <- dunnTest(lesion.surface ~ mutant, data = KG, method ="bh")
kgd=kgd$res
cldList(comparison = kgd$Comparison,
        p.value    = kgd$P.adj,
        threshold  = 0.05)
# #post-hoc test Wilcoxon rank sum test
#will allow us to easily add significance levels to our ggplots
library(ggpubr)
K.ls.w <- compare_means(lesion.surface ~ mutant, K, method="wilcox.test",
                        ref.group = ".all.")
K.ls.w <- dplyr::rename(K.ls.w, 
                         mutant = "group2",
                        K.pvalue = "p.format",
                         K.sig = "p.signif") #rename columns
colnames(K.ls.w) #check column names
K.ls.w <- K.ls.w[, c(3,6,7)] # keep only necessary columns
K.ls.w
### Lesion number ---- 
k1 <- aov(lesion.count ~ mutant, data = Klc) 
anova(k1)
#p-value =< 2.2e-16 ***, significant
#normality of residuals check
#distrubution/normality
par(mfrow=c(1,2))
#Histogramme
hist(k1$residuals)
# Q-Q Plot
qqnorm(k1$residuals)
qqline(k1$residuals)
#kind of normal looking
shapiro.test(k1$residuals)
#<2e-16, not normal
kruskal.test(lesion.count ~ mutant, data = Klc)
#p <2e-16, significant
KFlc <- dplyr::filter(Klc, isolate == "FR13")
kf1 <- aov(lesion.count ~ mutant, data = Klc) 
anova(kf1)
#p-value =< 2.2e-16 ***, significant
#normality of residuals check
shapiro.test(kf1$residuals)
# 9.973e-08 not normal
kruskal.test(lesion.count ~ mutant, data = KFlc)
#1.721e-09
KGlc <- dplyr::filter(Klc, isolate == "Guy11")
kg1 <- aov(lesion.count ~ mutant, data = Klc) 
anova(kg1)
#p-value =< 2.2e-16 ***, significant
#normality of residuals check
shapiro.test(kg1$residuals)
# 9.973e-08 not normal
kruskal.test(lesion.count ~ mutant, data = KGlc)
# 5.976e-08


#Post-hoc 
#post-hoc test Dunn
library(FSA)
dunnTest(lesion.count ~ mutant, data = Klc, method ="bh")
kfcd <- dunnTest(lesion.count ~ mutant, data = KFlc, method ="bh")
kfcd=kfcd$res
cldList(comparison = kfcd$Comparison,
        p.value    = kfcd$P.adj,
        threshold  = 0.05)
kgcd <- dunnTest(lesion.count ~ mutant, data = KGlc, method ="bh")
kgcd=kgcd$res
cldList(comparison = kgcd$Comparison,
        p.value    = kgcd$P.adj,
        threshold  = 0.05)
# #post-hoc test Wilcoxon rank sum test
#will allow us to easily add significance levels to our ggplots
library(ggpubr)
K.lc.w <- compare_means(lesion.count ~ mutant, Klc, method="wilcox.test",
                        ref.group = ".all.") #issue with this code, will come back to it

### Punch inoculation ----
k2 <- aov(lesion.surface ~ mutant, data = K.p) 
anova(k2)
#p-value = 2.9e-09 ***, significant
#normality of residuals check
#distrubution/normality
par(mfrow=c(1,2))
#Histogramme
hist(k2$residuals)
# Q-Q Plot
qqnorm(k2$residuals)
qqline(k2$residuals)
#kind of normal looking
shapiro.test(k2$residuals)
#0.3, normal, anova stands!
#fr13
KFp <- dplyr::filter(K.p, isolate == "FR13")
kf2 <- aov(lesion.surface ~ mutant, data = KFp) 
anova(kf2)
#p-value = 0.0009117 ***, significant
#normality of residuals check
shapiro.test(kf2$residuals)
#0.5678, normal 
#guy11
KGp <- dplyr::filter(K.p, isolate == "Guy11")
kg2 <- aov(lesion.surface ~ mutant, data = KGp) 
anova(kg2)
#p-value = 1.895e-05 ***, significant
#normality of residuals check
shapiro.test(kg2$residuals)
#0.06516, normal anova stands!

#Post-hoc 
#post-hoc test Dunn
library(FSA)
TukeyHSD(kf2, "mutant", ordered = FALSE, conf.level = 0.95)
TukeyHSD(kg2, "mutant", ordered = FALSE, conf.level = 0.95)
# #post-hoc test Wilcoxon rank sum test
#will allow us to easily add significance levels to our ggplots
library(ggpubr)
K.p.w <- compare_means(lesion.surface ~ mutant, data = K.p, method="wilcox.test",
                       ref.group = ".all.") #issue with this code, will come back to it
K.p.w

##K60 ----
### Lesion surface ----
#Initial stats lesion surface
K60 <- aov(lesion.surface ~ mutant, data = K6) 
anova(K60)
#p-value = < 2.2e-16 ***, significant
#normality of residuals check
#distrubution/normality
par(mfrow=c(1,2))
#Histogramme
hist(K60$residuals) #histogram to check normality
# Q-Q Plot
qqnorm(K60$residuals) # qqplot to check normality 
qqline(K60$residuals)
#very anormal
shapiro.test(K60$residuals)
#too many entries to run shapiro but clearly anormal
kruskal.test(lesion.surface ~ mutant, data = K6)
#p-value <2e-16, significant !
#fr13
K6F <- dplyr::filter(K6, isolate == "FR13")
K6f0 <- aov(lesion.surface ~ mutant, data = K6F) 
anova(K6f0)
#p-value = < 2.2e-16 ***, significant
#normality of residuals check
shapiro.test(K6f0$residuals)
#too many entries to run shapiro but clearly anormal
kruskal.test(lesion.surface ~ mutant, data = K6F)
#p-value < 2.2e-16
#guy11
K6G <- dplyr::filter(K6, isolate == "Guy11")
K6g0 <- aov(lesion.surface ~ mutant, data = K6G) 
anova(K6g0)
#p-value = < 2.2e-16 ***, significant
#normality of residuals check
shapiro.test(K6g0$residuals)
#too many entries to run shapiro but clearly anormal
kruskal.test(lesion.surface ~ mutant, data = K6G)
#p-value < 2.2e-16

#Post-hoc 
#post-hoc test Dunn
library(FSA)
library(rcompanion)
dunnTest(lesion.surface ~ mutant, data = K6, method ="bh")
k6fd <- dunnTest(lesion.surface ~ mutant, data = K6F, method ="bh")
k6fd=k6fd$res
cldList(comparison = k6fd$Comparison,
        p.value    = k6fd$P.adj,
        threshold  = 0.05)
k6gd <- dunnTest(lesion.surface ~ mutant, data = K6G, method ="bh")
k6gd=k6gd$res
cldList(comparison = k6gd$Comparison,
        p.value    = k6gd$P.adj,
        threshold  = 0.05)
# #post-hoc test Wilcoxon rank sum test
#will allow us to easily add significance levels to our ggplots
library(ggpubr)
K6.ls.w <- compare_means(lesion.surface ~ mutant, K6, method="wilcox.test",
                         ref.group = ".all.")
K6.ls.w <- dplyr::rename(K6.ls.w, 
                        mutant = "group2",
                        k6.pvalue = "p.format",
                        k6.sig = "p.signif") #rename columns
colnames(K6.ls.w) #check column names
K6.ls.w <- K6.ls.w[, c(3,6,7)] # keep only necessary columns
K6.ls.w

### Lesion number ---- 
k61 <- aov(lesion.count ~ mutant, data = K6lc) 
anova(k61)
#p-value = <2e-16 ***, significant
#normality of residuals check
#distrubution/normality
par(mfrow=c(1,2))
#Histogramme
hist(k61$residuals)
# Q-Q Plot
qqnorm(k61$residuals)
qqline(k61$residuals)
#kind of normal looking
shapiro.test(k61$residuals)
#<2e-16, not normal
kruskal.test(lesion.count ~ mutant, data = K6lc)
#p <2e-16, significant
#fr13
K6Flc <- dplyr::filter(K6lc, isolate == "FR13")
k6f1 <- aov(lesion.count ~ mutant, data = K6Flc) 
anova(k6f1)
#p-value = 1.087e-09 ***, significant
#normality of residuals check
shapiro.test(k6f1$residuals)
#1.035e-12, not normal
kruskal.test(lesion.count ~ mutant, data = K6Flc)
# 8.365e-13
#guy11
K6Glc <- dplyr::filter(K6lc, isolate == "Guy11")
k6g1 <- aov(lesion.count ~ mutant, data = K6Glc) 
anova(k6g1)
#p-value = 1.087e-09 ***, significant
#normality of residuals check
shapiro.test(k6g1$residuals)
#1.035e-12, not normal
kruskal.test(lesion.count ~ mutant, data = K6Glc)

#Post-hoc 
#post-hoc test Dunn
library(FSA)
dunnTest(lesion.count ~ mutant, data = K6lc, method ="bh")
k6fcd <- dunnTest(lesion.count ~ mutant, data = K6Flc, method ="bh")
k6fcd=k6fcd$res
cldList(comparison = k6fcd$Comparison,
        p.value    = k6fcd$P.adj,
        threshold  = 0.05)
k6gcd <- dunnTest(lesion.count ~ mutant, data = K6Glc, method ="bh")
k6gcd=k6gcd$res
cldList(comparison = k6gcd$Comparison,
        p.value    = k6gcd$P.adj,
        threshold  = 0.05)
# #post-hoc test Wilcoxon rank sum test
#will allow us to easily add significance levels to our ggplots
library(ggpubr)
K6.lc.w <- compare_means(lesion.count ~ mutant, K6lc, method="wilcox.test",
                         ref.group = ".all.") #issue with this code, will come back to it


### Punch inoculation ----
k62 <- aov(lesion.surface ~ mutant, data = K6.p) 
anova(k62)
#p-value =  0.00032 ***, significant
#normality of residuals check
#distrubution/normality
par(mfrow=c(1,2))
#Histogramme
hist(k62$residuals)
# Q-Q Plot
qqnorm(k62$residuals)
qqline(k62$residuals)
#kind of normal looking
shapiro.test(k62$residuals)
#1e-07, not normal
kruskal.test(lesion.surface ~ mutant, data = K6.p)
#p = 0.003, barely significant
#fr13
K6Fp <- dplyr::filter(K6.p, isolate == "FR13")
k6f2 <- aov(lesion.surface ~ mutant, data = K6Fp) 
anova(k6f2)
#p-value =  0.01846 * ***, significant
#normality of residuals check
shapiro.test(k6f2$residuals)
#0.004125, not normal
kruskal.test(lesion.surface ~ mutant, data = K6Fp)
#p-value = 0.02815
K6Gp <- dplyr::filter(K6.p, isolate == "Guy11")
k6g2 <- aov(lesion.surface ~ mutant, data = K6Gp) 
anova(k6g2)
#p-value =  0.01564 *, significant
#normality of residuals check
shapiro.test(k6g2$residuals)
#0.1654,  normal



#Post-hoc 
#post-hoc test Dunn
library(FSA)
#fr13
dunnTest(lesion.surface ~ mutant, data = K6Fp, method ="bh")
#guy11
TukeyHSD(k6g2, "mutant", ordered = FALSE, conf.level = 0.95)
library(agricolae)
HSD.test(k6g2, "mutant", group = TRUE, console = TRUE) 
#no sig diffs
# #post-hoc test Wilcoxon rank sum test
#will allow us to easily add significance levels to our ggplots
library(ggpubr)
K6.p.w <- compare_means(lesion.surface ~ mutant, data = K6.p, method="wilcox.test",
                        ref.group = ".all.") #issue with this code, will come back to it
K6.p.w

##Nipponbare ----
### Lesion surface ----
#Initial stats lesion surface
N0 <- aov(lesion.surface ~ mutant, data = N) 
anova(N0)
#p-value = < 2.2e-16 ***, significant
#normality of residuals check
#distrubution/normality
par(mfrow=c(1,2))
#Histogramme
hist(N0$residuals) #histogram to check normality
# Q-Q Plot
qqnorm(N0$residuals) # qqplot to check normality 
qqline(N0$residuals)
#very anormal
shapiro.test(N0$residuals)
#too many entries to run shapiro but clearly anormal
kruskal.test(lesion.surface ~ mutant, data = N)
#p-value <2e-16, significant !
#fr13
NF <- dplyr::filter(N, isolate == "FR13")
Nf0 <- aov(lesion.surface ~ mutant, data = NF) 
anova(Nf0)
#p-value = < 2.2e-16 ***, significant
#normality of residuals check
#very anormal
shapiro.test(Nf0$residuals)
#too many entries to run shapiro but clearly anormal
kruskal.test(lesion.surface ~ mutant, data = NF)
#guy11
NG <- dplyr::filter(N, isolate == "Guy11")
Ng0 <- aov(lesion.surface ~ mutant, data = NG) 
anova(Ng0)
#p-value = < 2.2e-16 ***, significant
#normality of residuals check
#very anormal
shapiro.test(Ng0$residuals)
#too many entries to run shapiro but clearly anormal
kruskal.test(lesion.surface ~ mutant, data = NG)

#Post-hoc 
#post-hoc test Dunn
library(FSA)
dunnTest(lesion.surface ~ mutant, data = N, method ="bh")
nfd <- dunnTest(lesion.surface ~ mutant, data = NF, method ="bh")
nfd=nfd$res
cldList(comparison = nfd$Comparison,
        p.value    = nfd$P.adj,
        threshold  = 0.05)
ngd <- dunnTest(lesion.surface ~ mutant, data = NG, method ="bh")
ngd=ngd$res
cldList(comparison = ngd$Comparison,
        p.value    = ngd$P.adj,
        threshold  = 0.05)
# #post-hoc test Wilcoxon rank sum test
#will allow us to easily add significance levels to our ggplots
library(ggpubr)
N.ls.w <- compare_means(lesion.surface ~ mutant, N, method="wilcox.test",
                        ref.group = ".all.")
N.ls.w <- dplyr::rename(N.ls.w, 
                        mutant = "group2",
                        N.pvalue = "p.format",
                        N.sig = "p.signif") #rename columns
colnames(N.ls.w) #check column names
N.ls.w <- N.ls.w[, c(3,6,7)] # keep only necessary columns

N.ls.w

### Lesion number ---- 
N1 <- aov(lesion.count ~ mutant, data = Nlc) 
anova(N1)
#p-value = <2e-16 ***, significant
#normality of residuals check
#distrubution/normality
par(mfrow=c(1,2))
#Histogramme
hist(N1$residuals)
# Q-Q Plot
qqnorm(N1$residuals)
qqline(N1$residuals)
#kind of normal looking
shapiro.test(N1$residuals)
#<2e-16, not normal
kruskal.test(lesion.count ~ mutant, data = Nlc)
#p <2e-16, significant
#fr13
NFlc <- dplyr::filter(Nlc, isolate == "FR13")
Nf1 <- aov(lesion.count ~ mutant, data = NFlc) 
anova(Nf1)
#p-value = < 2.2e-16 ***, significant
#normality of residuals check
#very anormal
shapiro.test(Nf1$residuals)
#too many entries to run shapiro but clearly anormal
kruskal.test(lesion.count ~ mutant, data = NFlc)
#guy11
NGlc <- dplyr::filter(Nlc, isolate == "Guy11")
Ng1 <- aov(lesion.count ~ mutant, data = NGlc) 
anova(Ng1)
#p-value = < 2.2e-16 ***, significant
#normality of residuals check
#very anormal
shapiro.test(Ng1$residuals)
#too many entries to run shapiro but clearly anormal
kruskal.test(lesion.count ~ mutant, data = NGlc)

#Post-hoc 
#post-hoc test Dunn
library(FSA)
dunnTest(lesion.count ~ mutant, data = Nlc, method ="bh")
nfcd <- dunnTest(lesion.count ~ mutant, data = NFlc, method ="bh")
nfcd=nfcd$res
cldList(comparison = nfcd$Comparison,
        p.value    = nfcd$P.adj,
        threshold  = 0.05)
dunnTest(lesion.count ~ mutant, data = NGlc, method ="bh")
ngcd=ngcd$res
cldList(comparison = ngcd$Comparison,
        p.value    = ngcd$P.adj,
        threshold  = 0.05)
# #post-hoc test Wilcoxon rank sum test
#will allow us to easily add significance levels to our ggplots
library(ggpubr)
N.lc.w <- compare_means(lesion.count ~ mutant, Nlc, method="wilcox.test",
                        ref.group = ".all.") #issue with this code, will come back to it
N.lc.w
### Punch inoculation ----
N2 <- aov(lesion.surface ~ mutant, data = N.p) 
anova(N2)
#p-value =   <2e-16 ***, significant
#normality of residuals check
#distrubution/normality
par(mfrow=c(1,2))
#Histogramme
hist(N2$residuals)
# Q-Q Plot
qqnorm(N2$residuals)
qqline(N2$residuals)
#kind of normal looking
shapiro.test(N2$residuals)
#0.7,  normal!!! anova stands!
#fr13
NFp <- dplyr::filter(N.p, isolate == "FR13")
Nf2 <- aov(lesion.surface ~ mutant, data = NFp) 
anova(Nf2)
#p-value =   9.465e-08 ***, significant
#normality of residuals check
shapiro.test(Nf2$residuals)
#0.9649 normal, annova stands
#guy11
NGp <- dplyr::filter(N.p, isolate == "Guy11")
Ng2 <- aov(lesion.surface ~ mutant, data = NGp) 
anova(Ng2)
#p-value =   9.465e-08 ***, significant
#normality of residuals check
shapiro.test(Ng2$residuals)
#0.0006698, not normal 
kruskal.test(lesion.surface ~ mutant, data = NGp)

#Post-hoc 
TukeyHSD(Nf2, "mutant", ordered = FALSE, conf.level = 0.95)
library(FSA)
ngpd <- dunnTest(lesion.surface ~ mutant, data = NGp, method ="bh")
ngpd=ngpd$res
cldList(comparison = ngpd$Comparison,
        p.value    = ngpd$P.adj,
        threshold  = 0.05)
# #post-hoc test Wilcoxon rank sum test
#will allow us to easily add significance levels to our ggplots
library(ggpubr)
N.p.w <- compare_means(lesion.surface ~ mutant, data = N.p, method="wilcox.test",
                       ref.group = ".all.") #issue with this code, will come back to it
N.p.w


##Shin 2 ----
### Lesion surface ----
#Initial stats lesion surface
S0 <- aov(lesion.surface ~ mutant, data = S) 
anova(S0)
#p-value = < 2.2e-16 ***, significant
#normality of residuals check
#distrubution/normality
par(mfrow=c(1,2))
#Histogramme
hist(S0$residuals) #histogram to check normality
# Q-Q Plot
qqnorm(S0$residuals) # qqplot to check normality 
qqline(S0$residuals)
#very anormal
shapiro.test(S0$residuals)
#too many entries to run shapiro but clearly anormal
kruskal.test(lesion.surface ~ mutant, data = S)
#p-value <2e-16, significant !
#fr13
SF <- dplyr::filter(S, isolate == "FR13")
Sf0 <- aov(lesion.surface ~ mutant, data = SF) 
anova(Sf0)
#p-value = < 2.2e-16 ***, significant
#normality of residuals check
shapiro.test(Sf0$residuals)
#too many entries to run shapiro but clearly anormal
kruskal.test(lesion.surface ~ mutant, data = SF)
#p-value <2e-16, significant !
#guy11
SG <- dplyr::filter(S, isolate == "Guy11")
Sg0 <- aov(lesion.surface ~ mutant, data = SG) 
anova(Sg0)
#p-value = < 2.2e-16 ***, significant
#normality of residuals check
shapiro.test(Sg0$residuals)
#too many entries to run shapiro but clearly anormal
kruskal.test(lesion.surface ~ mutant, data = SG)
#p-value <2e-16, significant !


#Post-hoc 
#post-hoc test Dunn
library(FSA)
dunnTest(lesion.surface ~ mutant, data = S, method ="bh")
sfd <- dunnTest(lesion.surface ~ mutant, data = SF, method ="bh")
sfd=sfd$res
cldList(comparison = sfd$Comparison,
        p.value    = sfd$P.adj,
        threshold  = 0.05)
sgd <- dunnTest(lesion.surface ~ mutant, data = SG, method ="bh")
sgd=sgd$res
cldList(comparison = sgd$Comparison,
        p.value    = sgd$P.adj,
        threshold  = 0.05)
# #post-hoc test Wilcoxon rank sum test
#will allow us to easily add significance levels to our ggplots
library(ggpubr)
S.ls.w <- compare_means(lesion.surface ~ mutant, S, method="wilcox.test",
                        ref.group = ".all.")
S.ls.w <- dplyr::rename(S.ls.w, 
                        mutant = "group2",
                        S.pvalue = "p.format",
                        S.sig = "p.signif")
colnames(S.ls.w) #check column names
S.ls.w <- S.ls.w[, c(3,6,7)] 
S.ls.w
### Lesion number ---- 
S1 <- aov(lesion.count ~ mutant, data = Slc) 
anova(S1)
#p-value = <2e-16 ***, significant
#normality of residuals check
#distrubution/normality
par(mfrow=c(1,2))
#Histogramme
hist(S1$residuals)
# Q-Q Plot
qqnorm(S1$residuals)
qqline(S1$residuals)
#kind of normal looking
shapiro.test(S1$residuals)
#<2e-16, not normal
kruskal.test(lesion.count ~ mutant, data = Slc)
#p <2e-16, significant
#fr13
SFlc <- dplyr::filter(Slc, isolate == "FR13")
Sf1 <- aov(lesion.count ~ mutant, data = SFlc) 
anova(Sf1)
#p-value = <2e-16 ***, significant
#normality of residuals check
shapiro.test(Sf1$residuals)
#<2e-16, not normal
kruskal.test(lesion.count ~ mutant, data = SFlc)
#guy11
SGlc <- dplyr::filter(Slc, isolate == "Guy11")
Sg1 <- aov(lesion.count ~ mutant, data = SGlc) 
anova(Sg1)
#p-value = <2e-16 ***, significant
#normality of residuals check
shapiro.test(Sg1$residuals)
#<2e-16, not normal
kruskal.test(lesion.count ~ mutant, data = SGlc)

#Post-hoc 
#post-hoc test Dunn
library(FSA)
sfcd <- dunnTest(lesion.count ~ mutant, data = SFlc, method ="bh")
sfcd=sfcd$res
cldList(comparison = sfcd$Comparison,
        p.value    = sfcd$P.adj,
        threshold  = 0.05)
sgcd <- dunnTest(lesion.count ~ mutant, data = SGlc, method ="bh")
sgcd=sgcd$res
cldList(comparison = sgcd$Comparison,
        p.value    = sgcd$P.adj,
        threshold  = 0.05)
# #post-hoc test Wilcoxon rank sum test
#will allow us to easily add significance levels to our ggplots
library(ggpubr)
S.lc.w <- compare_means(lesion.count ~ mutant, Slc, method="wilcox.test",
                        ref.group = ".all.") #issue with this code, will come back to it

### Punch inoculation ----
S2 <- aov(lesion.surface ~ mutant, data = S.p) 
anova(S2)
#p-value =  1.3e-07 ***, significant
#normality of residuals check
#distrubution/normality
par(mfrow=c(1,2))
#Histogramme
hist(S2$residuals)
# Q-Q Plot
qqnorm(S2$residuals)
qqline(S2$residuals)
#kind of normal looking
shapiro.test(S2$residuals)
# 6e-04, not normal
kruskal.test(lesion.surface ~ mutant, data = S.p)
#p = 1e-05, significant
#fr13
SFp <- dplyr::filter(S.p, isolate == "FR13")
Sf2 <- aov(lesion.surface ~ mutant, data = SFp) 
anova(Sf2)
#p-value =  1.3e-07 ***, significant
#normality of residuals check
shapiro.test(Sf2$residuals)
# 0.06389, normal

#gguy11
SGp <- dplyr::filter(S.p, isolate == "Guy11")
Sg2 <- aov(lesion.surface ~ mutant, data = SGp) 
anova(Sg2)
#p-value =  1.3e-07 ***, significant
#normality of residuals check
shapiro.test(Sg2$residuals)
# 0.004701, not normal
kruskal.test(lesion.surface ~ mutant, data = SGp)

#Post-hoc 
#post-hoc tests
TukeyHSD(Sf2, "mutant", ordered = FALSE, conf.level = 0.95)
library(agricolae)
HSD.test(Sf2, "mutant", group = TRUE, console = TRUE) 
library(FSA)
sgpd <- dunnTest(lesion.surface ~ mutant, data = SGp, method ="bh")
sgpd=sgpd$res
cldList(comparison = sgpd$Comparison,
        p.value    = sgpd$P.adj,
        threshold  = 0.05)
# #post-hoc test Wilcoxon rank sum test
#will allow us to easily add significance levels to our ggplots
library(ggpubr)
S.p.w <- compare_means(lesion.surface ~ mutant, data = S.p, method="wilcox.test",
                       ref.group = ".all.") #issue with this code, will come back to it
S.p.w


##Tsuyuake ----
### Lesion surface ----
#Initial stats lesion surface
T0 <- aov(lesion.surface ~ mutant, data = T) 
anova(T0)
#p-value = 2.5e-07 ***, significant
#normality of residuals check
#distrubution/normality
par(mfrow=c(1,2))
#Histogramme
hist(T0$residuals) #histogram to check normality
# Q-Q Plot
qqnorm(T0$residuals) # qqplot to check normality 
qqline(T0$residuals)
#very anormal
shapiro.test(T0$residuals)
#too many entries to run shapiro but clearly anormal
kruskal.test(lesion.surface ~ mutant, data = T)
#p-value <2e-16, significant !

#Post-hoc 
#post-hoc test Dunn
library(FSA)
td <- dunnTest(lesion.surface ~ mutant, data = T, method ="bh")
td=td$res
cldList(comparison = td$Comparison,
        p.value    = td$P.adj,
        threshold  = 0.05)
# #post-hoc test Wilcoxon rank sum test
#will allow us to easily add significance levels to our ggplots
library(ggpubr)
T.ls.w <- compare_means(lesion.surface ~ mutant, T, method="wilcox.test",
                        ref.group = ".all.")
T.ls.w <- dplyr::rename(T.ls.w, 
                          mutant = "group2",
                        T.pvalue = "p.format",
                          T.sig = "p.signif") #rename columns
colnames(T.ls.w) #check column names
T.ls.w <- T.ls.w[, c(3,6,7)] # keep only necessary columns
### Lesion number ---- 
T1 <- aov(lesion.count ~ mutant, data = Tlc) 
anova(T1)
#p-value = 0.00049 ***, significant
#normality of residuals check
#distrubution/normality
par(mfrow=c(1,2))
#Histogramme
hist(T1$residuals)
# Q-Q Plot
qqnorm(T1$residuals)
qqline(T1$residuals)
#kind of normal looking
shapiro.test(T1$residuals)
#4e-12, not normal
kruskal.test(lesion.count ~ mutant, data = Tlc)
#p 0.1083, not significant

#Post-hoc 

#post-hoc test Dunn
library(FSA)
dunnTest(lesion.count ~ mutant, data = Tlc, method ="bh")
#no sig diffs
# #post-hoc test Wilcoxon rank sum test
#will allow us to easily add significance levels to our ggplots
library(ggpubr)
T.lc.w <- compare_means(lesion.count ~ mutant, Tlc, method="wilcox.test",
                        ref.group = ".all.") #issue with this code, will come back to it


##Wilcoxon master ----
library(dplyr)
p.sig.ls <- dplyr::left_join(K.ls.w, B.ls.w, by = "mutant") %>%
  left_join(., A.ls.w, by = "mutant" )%>%
  left_join(., K6.ls.w, by = "mutant" )%>%
  left_join(., N.ls.w, by = "mutant" )%>%
  left_join(., S.ls.w, by = "mutant" )%>%
  left_join(., T.ls.w, by = "mutant" )
p.sig.ls
write.csv2(p.sig.ls, file = "p.sig.ls.csv")
# library(xlsx)
# write.xlsx(p.sig.ls, file = "wilcoxon.xlsx",
#            sheetName = "lesion.surface", append = FALSE)






# Ggplots ---- 
#All ggplots
blue.palette <- c("#f0f9e8","#bae4bc", "#7bccc4","#43a2ca","#0868ac")
f.level_order <- c("FR13 WT","FR13 RFP" ,"mFA1","mFA2","mFA3","FA ectopic", "mFB1",
                   "mFB2", "mFB3","FB ectopic", "mFK1","mFK2", "mFK3","FK ectopic")
g.level_order <- c("Guy11 WT", "Guy11 RFP_1", "Guy11 RFP_2", "mGB1", "mGB2","mGB3",
                   "mGK1","mGK2","mGK3", "GK ectopic")
layout <- "
AAAA#
BBBCC"
## Aichi Asahi ----
A$mutant <- as.factor(A$mutant)
levels(A$mutant)
#re-order for plots
f.level_order <- c("FR13 WT","FR13 RFP" ,"mFA1","mFA2","mFA3","FA ectopic", "mFB1",
                   "mFB2", "mFB3","FB ectopic", "mFK1","mFK2", "mFK3","FK ectopic")

blue.palette <- c("#f0f9e8","#bae4bc", "#7bccc4","#43a2ca","#0868ac")
layout <- "
AAAA#
BBBCC"
#ggplot lesion size
library(ggpubr)
library(ggplot2)
#color by boxplot
als0 <- ggplot(A, aes(x=factor(mutant, level = f.level_order), y=lesion.surface, fill = effector)) + 
  geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface,color = rep, shape = rep), show.legend = TRUE) + # jitter size by lesion surface
  scale_size_continuous(range = c(0.1, 2.0))  + theme_minimal() + scale_x_discrete(limits=rev) +
  scale_y_continuous(trans = "log10", breaks = c(10, 100, 500, 1000, 5000, 10000)) 
Alsplot <- als0 + scale_colour_grey() + scale_fill_manual(values = blue.palette) +
  labs(title = "Spray lesion size", x = "Isolate", y = "Lesion size (in log10 pixels)") +
  # theme(axis.title.y = element_blank(), axis.title.x = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) +  coord_flip() +
  theme(text=element_text(size=18)) + theme(legend.position = "none")


ggsave(filename = "Als.png", plot = Alsplot, device = "png", height = 20, width = 30,
        units = "cm", dpi = 500)

#ggplot lesion number
library(ggpubr)
alc0 <- ggplot(Alc, aes(x=factor(mutant, level = f.level_order), y=lesion.count, fill = effector)) + 
  geom_boxplot(color = "darkslategray") + geom_jitter(size=2.0, aes(shape = rep, color = rep)) + 
  theme_minimal() 
Alcplot <- alc0 +  scale_colour_grey() + scale_fill_manual(values = blue.palette) +
  labs(title = "Spray lesion number", x = "Isolate", y = "Lesion number") + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) + theme(legend.position = "none")+  
  theme(text=element_text(size=18))


ggsave(filename = "Alc.png", plot = Alcplot, device = "png", height = 20, width = 30,
        units = "cm", dpi = 500)

# library(patchwork)
# Aplot <- ((Alcplot + plot_spacer()) /Alsplot ) + plot_annotation(tag_levels = 'A') 
# ggsave(filename = "A.png", plot = Aplot, device = "png", height = 35, width = 40,
#        units = "cm", dpi = 500)
#punch inoculation 
library(ggpubr)
library(ggplot2)
ap0 <- ggplot(A.p,aes(x=factor(mutant, level = f.level_order), y=lesion.surface, fill = effector)) + 
  geom_boxplot() + geom_jitter(aes(size = lesion.surface), show.legend = TRUE) + # jitter size by lesion surface
  scale_size_continuous(range = c(0.01, 3))  + theme_minimal() 
A.punch <- ap0 + scale_fill_manual(values = blue.palette) +
  labs(title = "Punch inoculation lesion size", x = "Isolate", y = "Lesion size (in pixels)") + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) + theme(text=element_text(size=18))

ggsave(filename = "A.punch.png", plot = A.punch, device = "png", height = 15, width = 35,
units = "cm", dpi = 500)

library(patchwork)
# Aplot <- (Alsplot/(Alcplot + A.punch)) + plot_annotation(tag_levels = 'A') + 
#   plot_layout(guides = 'collect')
Aplot <- (Alsplot + Alcplot + A.punch) +  
  plot_layout(design = layout) + plot_annotation(tag_levels = 'A')  
ggsave(filename = "A.png", plot = Aplot, device = "png", height = 40, width = 45,
       units = "cm", dpi = 500)
## Bl1 ----
### ggplot lesion size ----
BF <- dplyr::filter(B, isolate == "FR13")
BG <- dplyr::filter(B, isolate == "Guy11")
library(ggpubr)
bf.ls.p <- ggplot(BF,aes(x=factor(mutant, level = f.level_order), y=lesion.surface, fill = effector)) +   
  geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface, shape = rep, color = rep), show.legend = TRUE) + 
  scale_size_continuous(range = c(0.1, 2))+ theme_minimal() + scale_x_discrete(limits=rev) + 
  scale_y_continuous(trans = "log10", breaks = c(10, 100, 500, 1000, 5000, 10000)) 
leg.bfls <- get_legend(bf.ls.p)#get legend so we can run without the legend 

BF.ls.plot <- bf.ls.p + scale_colour_grey() + scale_fill_manual(values = blue.palette) +
  labs(title = "Spray lesion size", x = "Isolate", y = "Lesion size (in log10 pixels)") +
  # theme(axis.title.y = element_blank(), axis.title.x = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) +  coord_flip() + 
  theme(legend.position = "none")+ theme(text=element_text(size=18))

bg.ls.p <- ggplot(BG,aes(x=factor(mutant, level = g.level_order), y=lesion.surface, fill = effector)) + 
  geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface, shape = rep, color = rep), show.legend = TRUE) + 
  scale_size_continuous(range = c(0.1, 2))  + theme_minimal() + scale_x_discrete(limits=rev) +
  scale_y_continuous(trans = "log10", breaks = c(10, 100, 500, 1000, 5000, 10000)) 
leg.bgls <- get_legend(bg.ls.p)#get legend so we can run without the legend
BG.ls.plot <- bg.ls.p + scale_colour_grey() + scale_fill_manual(values = blue.palette) +
  labs(title = "Spray lesion size", x = "Isolate", y = "Lesion size (in log10 pixels)") +   
  #theme(axis.title.y = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) +  coord_flip() + 
  theme(legend.position = "none") + theme(text=element_text(size=18))
# library(patchwork)
# Blsplot <- (BF.ls.plot/BG.ls.plot) + plot_annotation(title = 'Lesion surface on Bl1')
# Blsplot
ggsave(filename = "Bfls.png", plot = BF.ls.plot, device = "png", height = 20, width = 30,
       units = "cm", dpi = 500)
ggsave(filename = "Bgls.png", plot = BG.ls.plot, device = "png", height = 20, width = 30,
       units = "cm", dpi = 500)

#ggplot lesion number
BFlc <- dplyr::filter(Blc, isolate == "FR13")
BGlc <- dplyr::filter(Blc, isolate == "Guy11")
library(ggpubr)
bf.lc0 <- ggplot(BFlc, aes(x=factor(mutant, level = f.level_order), y=lesion.count, fill = effector)) + 
  geom_boxplot() + geom_jitter(size=2.0, aes(shape = rep, color = rep)) + 
  theme_minimal()
BFlc <- bf.lc0 +  scale_colour_grey() + scale_fill_manual(values = blue.palette) +
  labs(title = "Spray lesion number", x = "Isolate", y = "Lesion number") + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) + 
  theme(legend.position = "none") + theme(text=element_text(size=18))
bg.lc0 <- ggplot(BGlc, aes(x=factor(mutant, level = g.level_order), y=lesion.count, fill = effector)) + 
  geom_boxplot() + geom_jitter(size=2.0, aes(shape = rep, color = rep)) + 
  theme_minimal()
Bglc <- bg.lc0 +  scale_colour_grey() + scale_fill_manual(values = blue.palette) +
  labs(title = "Spray lesion number", x = "Isolate", y = "Lesion number") + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) + 
  theme(legend.position = "none") + theme(text=element_text(size=18))

ggsave(filename = "Bflc.png", plot = BFlc, device = "png", height = 20, width = 30,
       units = "cm", dpi = 500)
ggsave(filename = "Bglc.png", plot = Bglc, device = "png", height = 20, width = 30,
       units = "cm", dpi = 500)
# library(patchwork)
# Bplot <- ((BF.ls.plot + BFlc) / (BG.ls.plot + Bglc)) + plot_annotation(
#   title = 'Lesion size and count on Bl1', tag_levels = 'A') 

#ggplot punch
BFp <- dplyr::filter(B.p, isolate == "FR13")
BGp <- dplyr::filter(B.p, isolate == "Guy11")
library(ggpubr)
library(ggplot2)
bf.p.p <- ggplot(BFp,aes(x=factor(mutant, level = f.level_order), y=lesion.surface, fill = effector)) + geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface), show.legend = TRUE) + # jitter size by lesion surface
  scale_size_continuous(range = c(0.01, 2))  + theme_minimal()
BF.p.plot <- bf.p.p + scale_fill_manual(values = blue.palette) +
  labs(title = "Punch inoculation lesion size", x = "Isolate", y = "Lesion size (in pixels)") +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) + theme(text=element_text(size=18))
bg.p.p <- ggplot(BGp, aes(x=factor(mutant, level = g.level_order), y=lesion.surface, fill = effector)) + geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface), show.legend = TRUE) + # jitter size by lesion surface
  scale_size_continuous(range = c(0.01, 2))  + theme_minimal()
BG.p.plot <- bg.p.p + scale_fill_manual(values = blue.palette) +
  labs(title = "Punch inoculation lesion size", x = "Isolate", y = "Lesion size (in pixels)") +   
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) + theme(text=element_text(size=18))

ggsave(filename = "Bfp.png", plot = BF.p.plot, device = "png", height = 15, width = 20,
       units = "cm", dpi = 500)
ggsave(filename = "Bgp.png", plot = BG.p.plot, device = "png", height = 15, width = 20,
       units = "cm", dpi = 500)

library(patchwork)
Bgplot <- (BG.ls.plot + Bglc + BG.p.plot) +  
  plot_layout(design = layout) + plot_annotation(tag_levels = 'A') 
ggsave(filename = "Bg.png", plot = Bgplot, device = "png", height = 40, width = 45,
       units = "cm", dpi = 500)
library(patchwork)
Bfplot <- (BF.ls.plot + BFlc + BF.p.plot) +  
  plot_layout(design = layout) + plot_annotation(tag_levels = 'A')  
ggsave(filename = "Bf.png", plot = Bfplot, device = "png", height = 40, width = 45,
       units = "cm", dpi = 500)


## K60 ----
### ggplot lesion size ----
K6F <- dplyr::filter(K6, isolate == "FR13")
library(ggpubr)
K6f.ls.p <- ggplot(K6F,aes(x=factor(mutant, level = f.level_order), y=lesion.surface, fill = effector)) +   
  geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface, shape = rep, color = rep), show.legend = TRUE) + 
  scale_size_continuous(range = c(0.1, 2))+ theme_minimal() + scale_x_discrete(limits=rev) + 
  scale_y_continuous(trans = "log10", breaks = c(10, 100, 500, 1000, 5000, 10000)) 
leg.k6fls <- get_legend(K6f.ls.p)#get legend so we can run without the legend 
K6F.ls.plot <- K6f.ls.p + scale_colour_grey() + scale_fill_manual(values = blue.palette) +
  labs(title = "Spray lesion size", x = "Isolate", y = "Lesion size (in log10 pixels)") +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) +  coord_flip() +
  theme(legend.position = "none") + theme(text=element_text(size=18))


ggsave(filename = "K6fls.png", plot = K6F.ls.plot, device = "png", height = 20, width = 30,
       units = "cm", dpi = 500)


### ggplot lesion number ----
K6Flc <- dplyr::filter(K6lc, isolate == "FR13")
library(ggpubr)
K6f.lc0 <- ggplot(K6Flc, aes(x=factor(mutant, level = f.level_order), y=lesion.count, fill = effector)) + 
  geom_boxplot() + geom_jitter(size=2.0, aes(shape = rep, color = rep)) + 
  theme_minimal() 
K6Flc <- K6f.lc0 +  scale_colour_grey() + scale_fill_manual(values = blue.palette) +
  labs(title = "Spray lesion number", x = "Isolate", y = "Lesion number") + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) +
  theme(legend.position = "none")+ theme(text=element_text(size=18))


ggsave(filename = "K6flc.png", plot = K6Flc, device = "png", height = 20, width = 30,
       units = "cm", dpi = 500)



### ggplot punch ----
K6Fp <- dplyr::filter(K6.p, isolate == "FR13")
library(ggpubr)
library(ggplot2)
K6f.p.p <- ggplot(K6Fp,aes(x=factor(mutant, level = f.level_order), y=lesion.surface, fill = effector)) + geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface), show.legend = TRUE) + # jitter size by lesion surface
  scale_size_continuous(range = c(0.01, 2))  + theme_minimal() 
K6F.p.plot <- K6f.p.p + scale_fill_manual(values = blue.palette) +
  labs(title = "Punch inoculation lesion size", x = "Isolate", y = "Lesion size (in pixels)") +
  # theme(axis.title.y = element_blank(), axis.title.x = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) + theme(text=element_text(size=18))


library(patchwork)
K6fplot <- (K6F.ls.plot + K6Flc +K6F.p.plot)+ plot_layout(design = layout) + 
  plot_annotation(tag_levels = 'A') 
ggsave(filename = "K6f.png", plot = K6fplot, device = "png",  height = 40, width = 45,
       units = "cm", dpi = 500)

## Tsuyuake ----
#ggplot lesion size
library(ggpubr)
t.ls.p <- ggplot(T,aes(x=factor(mutant, level = f.level_order), y=lesion.surface, fill = effector)) +   
  geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface, shape = rep, color = rep), show.legend = TRUE) + 
  scale_size_continuous(range = c(0.1, 2))+ theme_minimal() + scale_x_discrete(limits=rev) + 
  scale_y_continuous(trans = "log10", breaks = c(10, 100, 500, 1000, 5000, 10000))  
t.ls.p <- get_legend(t.ls.p)#get legend so we can run without the legend 
T.ls.plot <- t.ls.p + scale_color_grey() + scale_fill_manual(values = blue.palette) +
  labs(title = "Spray lesion size", x = "Isolate", y = "Lesion surface (in log10 pixels)") +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) +  coord_flip() + 
  theme(text=element_text(size=18))

ggsave(filename = "Tls.png", plot = T.ls.plot, device = "png", height = 20, width = 30,
       units = "cm", dpi = 500)
ggsave(filename = "Tls.leg.png", plot = t.ls.p)

#ggplot lesion number
library(ggpubr)
t.lc0 <- ggplot(Tlc, aes(x=factor(mutant, level = f.level_order), y=lesion.count, fill = effector)) + 
  geom_boxplot() + geom_jitter(size=2.0, aes(shape = rep, color = rep)) + 
  theme_minimal() 
Tlc <- t.lc0 +  scale_colour_grey() + scale_fill_manual(values = blue.palette) +
  labs(title = "Spray lesion number", x = "Isolate", y = "Lesion number") + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) + theme(legend.position = "none")
+ theme(text=element_text(size=18))

ggsave(filename = "Tlc.png", plot = Tlc, device = "png", height = 20, width = 30,
       units = "cm", dpi = 500)


library(patchwork)
Tplot <- (T.ls.plot/(Tlc + plot_spacer())) + plot_annotation(tag_levels = 'A') 
ggsave(filename = "T.png", plot = Tplot, device = "png",  height = 40, width = 45,
       units = "cm", dpi = 500)


## Susceptible ----
### Kasalath ----
#### ggplot lesion size ----
KF <- dplyr::filter(K, isolate == "FR13")
KG <- dplyr::filter(K, isolate == "Guy11")
library(ggpubr)
kf.ls.p <- ggplot(KF,aes(x=factor(mutant, level = f.level_order), y=lesion.surface, fill = effector)) +   
  geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface, shape = rep, color = rep), show.legend = TRUE) + 
  scale_size_continuous(range = c(0.1, 1.5))+ theme_minimal() + 
  scale_y_continuous(trans = "log10", breaks = c(10, 100, 500, 1000, 5000, 10000))  
kF.ls.plot <- kf.ls.p + scale_colour_grey() + scale_fill_manual(values = blue.palette) +
  labs(title = "Spray lesion size", x = "Isolate", y = "Lesion size (in log10 pixels)") +
  # theme(axis.title.y = element_blank(), axis.title.x = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) + theme(text=element_text(size=15))
# stat_compare_means(label = "p.signif", method = "wilcox.test",
#                    ref.group = ".all.", hide.ns = TRUE) 
kg.ls.p <- ggplot(KG,aes(x=factor(mutant, level = g.level_order), y=lesion.surface, fill = effector)) + 
  geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface, shape = rep, color = rep), show.legend = TRUE) + 
  scale_size_continuous(range = c(0.1, 2))  + theme_minimal() + 
  scale_y_continuous(trans = "log10", breaks = c(10, 100, 500, 1000, 5000, 10000)) 
kG.ls.plot <- kg.ls.p + scale_colour_grey() + scale_fill_manual(values = blue.palette) +
  labs(title = "Spray lesion size", x = "Isolate", y = "Lesion size (in log10 pixels)") +   
  #theme(axis.title.y = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0))+ theme(text=element_text(size=15))
# library(patchwork)
# klsplot <- (kF.ls.plot/kG.ls.plot) + plot_annotation(title = 'Lesion surface on Kasalath')
# klsplot
ggsave(filename = "kfls.png", plot = kF.ls.plot, device = "png", height = 20, width = 30,
       units = "cm", dpi = 500)
ggsave(filename = "kgls.png", plot = kG.ls.plot, device = "png", height = 20, width = 30,
       units = "cm", dpi = 500)

#### ggplot lesion number ----
KFlc <- dplyr::filter(Klc, isolate == "FR13")
KGlc <- dplyr::filter(Klc, isolate == "Guy11")
library(ggpubr)
kf.lc0 <- ggplot(KFlc, aes(x=factor(mutant, level = f.level_order), y=lesion.count, fill = effector)) + 
  geom_boxplot() + geom_jitter(size=2.0, aes(shape = rep, color = rep)) + 
  theme_minimal() 
kFlc <- kf.lc0 +  scale_colour_grey() + scale_fill_manual(values = blue.palette) +
  labs(title = "Spray lesion number", x = "Isolate", y = "Lesion number") + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) + theme(text=element_text(size=15))
# + stat_compare_means(label = "p.signif", method = "wilcox.test",
#                      ref.group = ".all.", hide.ns = TRUE, label.y = c(40, 60, 200, 30,30)) 
kg.lc0 <- ggplot(KGlc, aes(x=factor(mutant, level = g.level_order), y=lesion.count, fill = effector)) + 
  geom_boxplot() + geom_jitter(size=2.0, aes(shape = rep, color = rep)) + 
  theme_minimal() 
kglc <- kg.lc0 +  scale_colour_grey() + scale_fill_manual(values = blue.palette) +
  labs(title = "Spray lesion number", x = "Isolate", y = "Lesion number") + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) + theme(text=element_text(size=15))
# + stat_compare_means(label = "p.signif", method = "wilcox.test",
#                      ref.group = ".all.", hide.ns = TRUE, label.y = c(125, 170, 170, 30,30)) 
# library(patchwork)
# klcplot <- (kFlc / kglc) + plot_annotation(
#   title = 'Lesion number on Kasalath')
# klcplot
ggsave(filename = "kflc.png", plot = kFlc, device = "png", height = 20, width = 30,
       units = "cm", dpi = 500)
ggsave(filename = "kglc.png", plot = kglc, device = "png", height = 20, width = 30,
       units = "cm", dpi = 500)
# library(patchwork)
# kplot <- ((kF.ls.plot + kFlc) / (kG.ls.plot + kglc)) + plot_annotation(
#   title = 'Lesion size and count on Kasalath', tag_levels = 'A') 
library(patchwork)
kgplot <- ((kglc + plot_spacer())/kG.ls.plot ) + plot_annotation(title = 'Guy11 lesion size and count on Kasalath', tag_levels = 'A') 
ggsave(filename = "kg.png", plot = kgplot, device = "png", height = 25, width = 30,
       units = "cm", dpi = 500)
#### ggplot punch ----
KFp <- dplyr::filter(K.p, isolate == "FR13")
KGp <- dplyr::filter(K.p, isolate == "Guy11")
library(ggpubr)
library(ggplot2)
kf.p.p <- ggplot(KFp,aes(x=factor(mutant, level = f.level_order), y=lesion.surface, fill = effector)) + geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface), show.legend = TRUE) + # jitter size by lesion surface
  scale_size_continuous(range = c(0.01, 2))  + theme_minimal()
kF.p.plot <- kf.p.p + scale_fill_manual(values = blue.palette) +
  labs(title = "Punch inoculation lesion size", x = "Isolate", y = "Lesion size (in pixels)") +
  # theme(axis.title.y = element_blank(), axis.title.x = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) + theme(text=element_text(size=15))
kg.p.p <- ggplot(KGp, aes(x=factor(mutant, level = g.level_order), y=lesion.surface, fill = effector)) + geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface), show.legend = TRUE) + # jitter size by lesion surface
  scale_size_continuous(range = c(0.01, 2))  + theme_minimal() 
kG.p.plot <- kg.p.p + scale_fill_manual(values = blue.palette) +
  labs(title = "Punch inoculation lesion size", x = "Isolate", y = "Lesion size (in pixels)") +   
  #theme(axis.title.y = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) + theme(text=element_text(size=15))
library(patchwork)
kpplot <- (kF.p.plot / kG.p.plot) + plot_annotation(
  title = 'Punch inoculation on Kasalath')
kpplot
ggsave(filename = "kp.png", plot = kpplot, device = "png", height = 20, width = 30,
       units = "cm", dpi = 500)

### K60 Guy11 ----
#### Lesion size ----
K6G <- dplyr::filter(K6, isolate == "Guy11")
K6g.ls.p <- ggplot(K6G,aes(x=factor(mutant, level = g.level_order), y=lesion.surface, fill = effector)) + 
  geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface, shape = rep, color = rep), show.legend = TRUE) + 
  scale_size_continuous(range = c(0.1, 2))  + theme_minimal() + 
  scale_y_continuous(trans = "log10", breaks = c(10, 100, 500, 1000, 5000, 10000)) 
leg.k6gls <- get_legend(K6g.ls.p)#get legend so we can run without the legend 
K6G.ls.plot <- K6g.ls.p + scale_colour_grey() + scale_fill_manual(values = blue.palette) +
  labs(title = "Spray lesion size", x = "Isolate", y = "Lesion surface (in log10 pixels)") +   
  theme(axis.text.x=element_text(angle = -90, hjust = 0))  +
theme(text=element_text(size=15))
#### Lesion number ----
K6Glc <- dplyr::filter(K6lc, isolate == "Guy11")
K6g.lc0 <- ggplot(K6Glc, aes(x=factor(mutant, level = g.level_order), y=lesion.count, fill = effector)) + 
  geom_boxplot() + geom_jitter(size=2.0, aes(shape = rep, color = rep)) + 
  theme_minimal() 
K6glc <- K6g.lc0 +  scale_colour_grey() + scale_fill_manual(values = blue.palette) +
  labs(title = "Spray lesion number", x = "Isolate", y = "Lesion number") + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) +
  theme(legend.position = "none")+ theme(text=element_text(size=15))
#### Punch ----
K6Gp <- dplyr::filter(K6.p, isolate == "Guy11")
K6g.p.p <- ggplot(K6Gp, aes(x=factor(mutant, level = g.level_order), y=lesion.surface, fill = effector)) + geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface), show.legend = TRUE) + # jitter size by lesion surface
  scale_size_continuous(range = c(0.01, 2))  + theme_minimal() 
K6G.p.plot <- K6g.p.p + scale_fill_manual(values = blue.palette) +
  labs(title = "Punch inoculation lesion size", x = "Isolate", y = "Lesion surface (in pixels)") +   
  #theme(axis.title.y = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0))+ theme(text=element_text(size=15))   

ggsave(filename = "K6gls.png", plot = K6G.ls.plot, device = "png",  height = 20, width = 40,
        units = "cm", dpi = 500)

### Nipponbare ----
#### ggplot lesion size ----
NF <- dplyr::filter(N, isolate == "FR13")
NG <- dplyr::filter(N, isolate == "Guy11")
library(ggpubr)
Nf.ls.p <- ggplot(NF,aes(x=factor(mutant, level = f.level_order), y=lesion.surface, fill = effector)) +   
  geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface, shape = rep, color = rep), show.legend = TRUE) + 
  scale_size_continuous(range = c(0.1, 1.5))+ theme_minimal() +  
  scale_y_continuous(trans = "log10", breaks = c(10, 100, 500, 1000, 5000, 10000))  
NF.ls.plot <- Nf.ls.p + scale_colour_grey() + scale_fill_manual(values = blue.palette) +
  labs(title = "Spray lesion size", x = "Isolate", y = "Lesion size (in log10 pixels)") +
  # theme(axis.title.y = element_blank(), axis.title.x = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) + theme(text=element_text(size=15))
# stat_compare_means(label = "p.signif", method = "wilcox.test",
#                    ref.group = ".all.", hide.ns = TRUE) 
Ng.ls.p <- ggplot(NG,aes(x=factor(mutant, level = g.level_order), y=lesion.surface, fill = effector)) + 
  geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface, shape = rep, color = rep), show.legend = TRUE) + 
  scale_size_continuous(range = c(0.1, 2))  + theme_minimal() + 
  scale_y_continuous(trans = "log10", breaks = c(10, 100, 500, 1000, 5000, 10000)) 
NG.ls.plot <- Ng.ls.p + scale_colour_grey() + scale_fill_manual(values = blue.palette) +
  labs(title = "Spray lesion size", x = "Isolate", y = "Lesion size (in log10 pixels)") +   
  #theme(axis.title.y = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0))+ theme(text=element_text(size=15))
# library(patchwork)
# Nlsplot <- (NF.ls.plot/NG.ls.plot) + plot_annotation(title = 'Lesion surface on Nipponbare')
# Nlsplot
ggsave(filename = "Nfls.png", plot = NF.ls.plot, device = "png", height = 20, width = 30,
       units = "cm", dpi = 500)
ggsave(filename = "Ngls.png", plot = NG.ls.plot, device = "png", height = 20, width = 30,
       units = "cm", dpi = 500)

#### ggplot lesion number ----
NFlc <- dplyr::filter(Nlc, isolate == "FR13")
NGlc <- dplyr::filter(Nlc, isolate == "Guy11")
library(ggpubr)
Nf.lc0 <- ggplot(NFlc, aes(x=factor(mutant, level = f.level_order), y=lesion.count, fill = effector)) + 
  geom_boxplot() + geom_jitter(size=2.0, aes(shape = rep, color = rep)) + 
  theme_minimal() 
NFlc <- Nf.lc0 +  scale_colour_grey() + scale_fill_manual(values = blue.palette) +
  labs(title = "Spray lesion number", x = "Isolate", y = "Lesion number") + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) + theme(text=element_text(size=15))
# + stat_compare_means(label = "p.signif", method = "wilcox.test",
#                      ref.group = ".all.", hide.ns = TRUE, label.y = c(40, 60, 200, 30,30)) 
Ng.lc0 <- ggplot(NGlc, aes(x=factor(mutant, level = g.level_order), y=lesion.count, fill = effector)) + 
  geom_boxplot() + geom_jitter(size=2.0, aes(shape = rep, color = rep)) + 
  theme_minimal() 
Nglc <- Ng.lc0 +  scale_colour_grey() + scale_fill_manual(values = blue.palette) +
  labs(title = "Spray lesion number", x = "Isolate", y = "Lesion number") + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) + theme(text=element_text(size=15))
# + stat_compare_means(label = "p.signif", method = "wilcox.test",
#                      ref.group = ".all.", hide.ns = TRUE, label.y = c(125, 170, 170, 30,30)) 
# library(patchwork)
# Nlcplot <- (NFlc / Nglc) + plot_annotation(
#   title = 'Lesion number on Nipponbare')
# Nlcplot
ggsave(filename = "Nflc.png", plot = NFlc, device = "png", height = 20, width = 30,
       units = "cm", dpi = 500)
ggsave(filename = "Nglc.png", plot = Nglc, device = "png", height = 20, width = 30,
       units = "cm", dpi = 500)
# library(patchwork)
# Nplot <- ((NF.ls.plot + NFlc) / (NG.ls.plot + Nglc)) + plot_annotation(
#   title = 'Lesion size and count on Nipponbare', tag_levels = 'A') 
library(patchwork)
Ngplot <- ((Nglc + plot_spacer())/NG.ls.plot ) + plot_annotation(title = 'Guy11 lesion size and count on Nipponbare', tag_levels = 'A') 
ggsave(filename = "Ng.png", plot = Ngplot, device = "png", height = 25, width = 30,
       units = "cm", dpi = 500)
#### ggplot punch ----
NFp <- dplyr::filter(N.p, isolate == "FR13")
NGp <- dplyr::filter(N.p, isolate == "Guy11")
library(ggpubr)
library(ggplot2)
Nf.p.p <- ggplot(NFp,aes(x=factor(mutant, level = f.level_order), y=lesion.surface, fill = effector)) + geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface), show.legend = TRUE) + # jitter size by lesion surface
  scale_size_continuous(range = c(0.01, 2))  + theme_minimal() 
NF.p.plot <- Nf.p.p + scale_fill_manual(values = blue.palette) +
  labs(title = "Punch inoculation lesion size", x = "Isolate", y = "Lesion size (in pixels)") +
  # theme(axis.title.y = element_blank(), axis.title.x = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) + theme(text=element_text(size=15))
Ng.p.p <- ggplot(NGp, aes(x=factor(mutant, level = g.level_order), y=lesion.surface, fill = effector)) + geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface), show.legend = TRUE) + # jitter size by lesion surface
  scale_size_continuous(range = c(0.01, 2))  + theme_minimal() 
NG.p.plot <- Ng.p.p + scale_fill_manual(values = blue.palette) +
  labs(title = "Punch inoculation lesion size", x = "Isolate", y = "Lesion size (in pixels)") +   
  #theme(axis.title.y = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0))+ theme(text=element_text(size=15)) 
library(patchwork)
Npplot <- (NF.p.plot / NG.p.plot) + plot_annotation(
  title = 'Punch inoculation on Nipponbare')
Npplot
ggsave(filename = "Np.png", plot = Npplot, device = "png", height = 20, width = 30,
       units = "cm", dpi = 500)
### Shin 2 ----
##### ggplot lesion size ----
SF <- dplyr::filter(S, isolate == "FR13")
SG <- dplyr::filter(S, isolate == "Guy11")
library(ggpubr)
Sf.ls.p <- ggplot(SF,aes(x=factor(mutant, level = f.level_order), y=lesion.surface, fill = effector)) +   
  geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface, shape = rep, color = rep), show.legend = TRUE) + 
  scale_size_continuous(range = c(0.1, 2))+ theme_minimal() + 
  scale_y_continuous(trans = "log10", breaks = c(10, 100, 500, 1000, 5000, 10000))  
SF.ls.plot <- Sf.ls.p + scale_colour_grey() + scale_fill_manual(values = blue.palette) +
  labs(title = "Spray lesion size", x = "Isolate", y = "Lesion surface (in log10 pixels)") +
  # theme(axis.title.y = element_blank(), axis.title.x = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) + theme(text=element_text(size=15))
# stat_compare_means(label = "p.signif", method = "wilcox.test",
#                    ref.group = ".all.", hide.ns = TRUE) 
Sg.ls.p <- ggplot(SG,aes(x=factor(mutant, level = g.level_order), y=lesion.surface, fill = effector)) + 
  geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface, shape = rep, color = rep), show.legend = TRUE) + 
  scale_size_continuous(range = c(0.1, 2))  + theme_minimal() + 
  scale_y_continuous(trans = "log10", breaks = c(10, 100, 500, 1000, 5000, 10000)) 
SG.ls.plot <- Sg.ls.p + scale_colour_grey() + scale_fill_manual(values = blue.palette) +
  labs(title = "Spray lesion size", x = "Isolate", y = "Lesion surface (in log10 pixels)") +   
  #theme(axis.title.y = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) + theme(text=element_text(size=15))
# library(patchwork)
# Slsplot <- (SF.ls.plot/SG.ls.plot) + plot_annotation(title = 'Lesion surface on Shin 2')
# Slsplot
ggsave(filename = "Sfls.png", plot = SF.ls.plot, device = "png", height = 20, width = 30,
       units = "cm", dpi = 500)
ggsave(filename = "Sgls.png", plot = SG.ls.plot, device = "png", height = 20, width = 30,
       units = "cm", dpi = 500)

##### ggplot lesion number -----
SFlc <- dplyr::filter(Slc, isolate == "FR13")
SGlc <- dplyr::filter(Slc, isolate == "Guy11")
library(ggpubr)
Sf.lc0 <- ggplot(SFlc, aes(x=factor(mutant, level = f.level_order), y=lesion.count, fill = effector)) + 
  geom_boxplot() + geom_jitter(size=2.0, aes(shape = rep, color = rep)) + 
  theme_minimal() 
SFlc <- Sf.lc0 +  scale_colour_grey() + scale_fill_manual(values = blue.palette) +
  labs(title = "Spray lesion number", x = "Isolate", y = "Lesion number") + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) + theme(text=element_text(size=15))
# + stat_compare_means(label = "p.signif", method = "wilcox.test",
#                      ref.group = ".all.", hide.ns = TRUE, label.y = c(40, 60, 200, 30,30)) 
Sg.lc0 <- ggplot(SGlc, aes(x=factor(mutant, level = g.level_order), y=lesion.count, fill = effector)) + 
  geom_boxplot() + geom_jitter(size=2.0, aes(shape = rep, color = rep)) + 
  theme_minimal() 
Sglc <- Sg.lc0 +  scale_colour_grey() + scale_fill_manual(values = blue.palette) +
  labs(title = "Spray lesion number", x = "Isolate", y = "Lesion number") + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) + theme(text=element_text(size=15))
# + stat_compare_means(label = "p.signif", method = "wilcox.test",
#                      ref.group = ".all.", hide.ns = TRUE, label.y = c(125, 170, 170, 30,30)) 
# library(patchwork)
# Slcplot <- (SFlc / Sglc) + plot_annotation(
#   title = 'Lesion number on Shin 2')
# Slcplot
ggsave(filename = "Sflc.png", plot = SFlc, device = "png", height = 20, width = 30,
       units = "cm", dpi = 500)
ggsave(filename = "Sglc.png", plot = Sglc, device = "png", height = 20, width = 30,
       units = "cm", dpi = 500)
# library(patchwork)
# Splot <- ((SF.ls.plot + SFlc) / (SG.ls.plot + Sglc)) + plot_annotation(
#   title = 'Lesion size and count on Shin 2', tag_levels = 'A') 
library(patchwork)
Sgplot <- ((Sglc + plot_spacer())/SG.ls.plot ) + plot_annotation(title = 'Guy11 lesion size and count on Shin 2', tag_levels = 'A') 
ggsave(filename = "Sg.png", plot = Sgplot, device = "png", height = 25, width = 30,
       units = "cm", dpi = 500)
##### ggplot punch ----
SFp <- dplyr::filter(S.p, isolate == "FR13")
SGp <- dplyr::filter(S.p, isolate == "Guy11")
library(ggpubr)
library(ggplot2)
Sf.p.p <- ggplot(SFp,aes(x=factor(mutant, level = f.level_order), y=lesion.surface, fill = effector)) + geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface), show.legend = TRUE) + # jitter size by lesion surface
  scale_size_continuous(range = c(0.01, 2))  + theme_minimal() 
SF.p.plot <- Sf.p.p + scale_fill_manual(values = blue.palette) +
  labs(title = "Punch inoculation lesion size", x = "Isolate", y = "Lesion size (in pixels)") +
  # theme(axis.title.y = element_blank(), axis.title.x = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) + theme(text=element_text(size=15))
Sg.p.p <- ggplot(SGp, aes(x=factor(mutant, level = g.level_order), y=lesion.surface, fill = effector)) + geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface), show.legend = TRUE) + # jitter size by lesion surface
  scale_size_continuous(range = c(0.01, 2))  + theme_minimal() 
SG.p.plot <- Sg.p.p + scale_fill_manual(values = blue.palette) +
  labs(title = "Punch inoculation lesion size", x = "Isolate", y = "Lesion size (in pixels)") +   
  #theme(axis.title.y = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0))+ theme(text=element_text(size=15))
library(patchwork)
Spplot <- (SF.p.plot / SG.p.plot) + plot_annotation(
  title = 'Punch inoculation on Shin 2')
Spplot
ggsave(filename = "Sp.png", plot = Spplot, device = "png", height = 20, width = 30,
       units = "cm", dpi = 500)


### compiled  virulence ggplots ----

#FR13
library(patchwork)
#Kasalath
KFplot <- (kF.ls.plot | kFlc | kF.p.plot) + plot_layout(guides = 'collect') + plot_annotation(tag_levels = 'A') 
ggsave(filename = "KF.png", plot = KFplot, device = "png", height =  20, width = 40,
       units = "cm", dpi = 500)
#Nipponbare
NFplot <- (NF.ls.plot | NFlc | NF.p.plot) + plot_layout(guides = 'collect') + 
  plot_annotation(tag_levels = 'A') 
ggsave(filename = "NF.png", plot = NFplot, device = "png", height =  20, width = 40,
       units = "cm", dpi = 500)
#Shin2
SFplot <- (SF.ls.plot | SFlc | SF.p.plot) + plot_layout(guides = 'collect') + plot_annotation(tag_levels = 'A') 
ggsave(filename = "SF.png", plot = SFplot, device = "png", height =  20, width = 40,
       units = "cm", dpi = 500)

#Guy11
library(patchwork)
#Kasalath
KGplot <- (kG.ls.plot | kglc | kG.p.plot) + plot_layout(guides = 'collect') + plot_annotation(tag_levels = 'A') 
ggsave(filename = "KG.png", plot = KGplot, device = "png", height =  20, width = 40,
       units = "cm", dpi = 500)
#K60
K6gplot <- (K6G.ls.plot + K6glc + K6G.p.plot)  + plot_layout(guides = 'collect') +
  plot_annotation(tag_levels = 'A') 
ggsave(filename = "K6G.png", plot = K6gplot, device = "png", height =  20, width = 40,
       units = "cm", dpi = 500)
#Nipponbare
library(patchwork)
NGplot <- (NG.ls.plot + Nglc + NG.p.plot) + plot_layout(guides = 'collect') + 
  plot_annotation(tag_levels = 'A') 
ggsave(filename = "NG.png", plot = NGplot, device = "png", height =  20, width = 40,
       units = "cm", dpi = 500)
#Shin2
SGplot <- (SG.ls.plot | SGlc | SG.p.plot) + plot_layout(guides = 'collect') + plot_annotation(tag_levels = 'A') 
ggsave(filename = "SG.png", plot = SGplot, device = "png", height =  20, width = 40,
       units = "cm", dpi = 500)


#FR13
# #Lesion size
# library(patchwork)
# F.ls <- (kF.ls.plot | NF.ls.plot | SF.ls.plot ) + plot_layout(guides = 'collect') + plot_annotation(tag_levels = 'A') 
# ggsave(filename = "F.ls.png", plot = F.ls, device = "png", height =  20, width = 40,
#        units = "cm", dpi = 500)
# #lesion count
# F.lc <- (kFlc | NFlc | SFlc ) + plot_layout(guides = 'collect')+ plot_annotation(tag_levels = 'A') 
# ggsave(filename = "F.lc.png", plot = F.lc, device = "png", height = 20, width = 40,
#        units = "cm", dpi = 500)
# #punch
# F.p <- (kF.p.plot | NF.p.plot | SF.p.plot ) + plot_layout(guides = 'collect') + plot_annotation(tag_levels = 'A') 
# ggsave(filename = "F.p.png", plot = F.p, device = "png",height = 20, width = 40,
#        units = "cm", dpi = 500)
# Guy11
# #Lesion size
# library(patchwork)
# G.ls <- (kG.ls.plot | NG.ls.plot | SG.ls.plot ) + plot_layout(guides = 'collect')+ plot_annotation(tag_levels = 'A') 
# ggsave(filename = "G.ls.png", plot = G.ls, device = "png", height =  20, width = 40,
#        units = "cm", dpi = 500)
# #lesion count
# G.lc <- (kglc | Nglc | Sglc ) + plot_layout(guides = 'collect') 
# ggsave(filename = "G.lc.png", plot = G.lc, device = "png", height = 20, width = 40,
#        units = "cm", dpi = 500)
# #punch
# G.p <- (kG.p.plot | NG.p.plot | SG.p.plot ) + plot_layout(guides = 'collect')
# ggsave(filename = "G.p.png", plot = G.p, device = "png", height = 20, width = 40,
#        units = "cm", dpi = 500)
# 
