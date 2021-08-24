#Master Inoculation Sheet
install.packages("xlsx")

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
#p-value <2e-16, significant !

#Post-hoc 
#post-hoc test Dunn
library(FSA)
dunnTest(lesion.surface ~ mutant, data = A, method ="bh")
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
#p =  2e-12, significant

#Post-hoc 

#post-hoc test Dunn
library(FSA)
dunnTest(lesion.count ~ mutant, data = Alc, method ="bh")
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
dunnTest(lesion.surface ~ mutant, data = A.p, method ="bh")
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

#Post-hoc 
#post-hoc test Dunn
library(FSA)
dunnTest(lesion.surface ~ mutant, data = B, method ="bh")
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

#Kruskal test with only rep 1
kruskal.test(lesion.count ~ mutant, data = Blc.1) #p-value <2e-16
#Kruskal test with only rep 2
kruskal.test(lesion.count ~ mutant, data = Blc.2) #p-value <2e-16
#same values for both

#Post-hoc 

#post-hoc test Dunn
library(FSA)
dunnTest(lesion.count ~ mutant, data = Blc, method ="bh")
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


#Post-hoc 

#post-hoc test Dunn
library(FSA)
dunnTest(lesion.surface ~ mutant, data = B.p, method ="bh")
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

#Post-hoc 
#post-hoc test Dunn
library(FSA)
dunnTest(lesion.surface ~ mutant, data = K, method ="bh")
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
#p-value = <2e-16 ***, significant
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


#Post-hoc 
#post-hoc test Dunn
library(FSA)
dunnTest(lesion.count ~ mutant, data = Klc, method ="bh")
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

#Post-hoc 
#post-hoc test Dunn
library(FSA)
dunnTest(lesion.surface ~ mutant, data = K.p, method ="bh")
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

#Post-hoc 
#post-hoc test Dunn
library(FSA)
dunnTest(lesion.surface ~ mutant, data = K6, method ="bh")
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

#Post-hoc 
#post-hoc test Dunn
library(FSA)
dunnTest(lesion.count ~ mutant, data = K6lc, method ="bh")
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


#Post-hoc 

#post-hoc test Dunn
library(FSA)
dunnTest(lesion.surface ~ mutant, data = K6.p, method ="bh")
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


#Post-hoc 
#post-hoc test Dunn
library(FSA)
dunnTest(lesion.surface ~ mutant, data = N, method ="bh")
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

#Post-hoc 
#post-hoc test Dunn
library(FSA)
dunnTest(lesion.count ~ mutant, data = Nlc, method ="bh")
# #post-hoc test Wilcoxon rank sum test
#will allow us to easily add significance levels to our ggplots
library(ggpubr)
N.lc.w <- compare_means(lesion.count ~ mutant, Nlc, method="wilcox.test",
                        ref.group = ".all.") #issue with this code, will come back to it
N.lc.w
## Punch inoculation ----
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

#Post-hoc 
#post-hoc test Dunn
library(FSA)
dunnTest(lesion.surface ~ mutant, data = N.p, method ="bh")
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

#Kruskal test with only rep 1
kruskal.test(lesion.surface ~ mutant, data = S.1) #p-value <2e-16
#Kruskal test with only rep 2
kruskal.test(lesion.surface ~ mutant, data = S.2) #p-value <2e-16
#same values for both

#Post-hoc 
#post-hoc test Dunn
library(FSA)
dunnTest(lesion.surface ~ mutant, data = S, method ="bh")
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

#Kruskal test with only rep 1
kruskal.test(lesion.count ~ mutant, data = Slc.1) #p-value <2e-16
#Kruskal test with only rep 2
kruskal.test(lesion.count ~ mutant, data = Slc.2) #p-value <2e-16
#same values for both

#Post-hoc 
#post-hoc test Dunn
library(FSA)
dunnTest(lesion.count ~ mutant, data = Slc, method ="bh")
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

#Post-hoc 
#post-hoc test Dunn
library(FSA)
dunnTest(lesion.surface ~ mutant, data = S.p, method ="bh")
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
dunnTest(lesion.surface ~ mutant, data = T, method ="bh")
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
#p 1e-05, significant

#Post-hoc 

#post-hoc test Dunn
library(FSA)
dunnTest(lesion.count ~ mutant, data = Tlc, method ="bh")
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
