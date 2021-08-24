#Aichi Asahi
load("spray.RData")
load("lc.RData")
load("punch.RData")

A <- spray[spray$cultivar=="Aichi Asahi",] ; titre <- "Aichi Asahi"
A.1 <- dplyr::filter(A, rep == "1") #A for rep 1
A.2 <- dplyr::filter(A, rep == "2") #A for rep 2
Alc <- lc[lc$cultivar =="Aichi Asahi",] ; titre <- "Aichi Asahi"
Alc.1 <-dplyr::filter(Alc, rep == "1") # Alc for rep 1
Alc.2 <-dplyr::filter(Alc, rep == "2") # Alc for rep 2
A.p <- punch[punch$cultivar=="Aichi Asahi",] ; titre <- "Aichi Asahi"

# Stats ----

## Lesion surface ----
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

#Kruskal test with only rep 1
kruskal.test(lesion.surface ~ mutant, data = A.1) #p-value <2e-16
#Kruskal test with only rep 2
kruskal.test(lesion.surface ~ mutant, data = A.2) #p-value <2e-16
#same values for both

#Post-hoc 
#post-hoc test Dunn
library(FSA)
dunnTest(lesion.surface ~ mutant, data = A, method ="bh")
# #post-hoc test Wilcoxon rank sum test
#will allow us to easily add significance levels to our ggplots
library(ggpubr)
A.ls.w <- compare_means(lesion.surface ~ mutant, A, method="wilcox.test",
                      ref.group = ".all.")

#Wilcoxon test with only rep 1
compare_means(lesion.surface ~ mutant, A.1, method="wilcox.test",
              ref.group = ".all.") 
#Wilcoxon test with only rep 2
compare_means(lesion.surface ~ mutant, A.2, method="wilcox.test",
              ref.group = ".all.") 
#actually fairly different results... 

# Lesion number ---- 
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

#Kruskal test with only rep 1
kruskal.test(lesion.surface ~ mutant, data = Alc.1) #p-value <2e-16
#Kruskal test with only rep 2
kruskal.test(lesion.surface ~ mutant, data = Alc.2) #p-value <2e-16
#same values for both

#Post-hoc 

#post-hoc test Dunn
library(FSA)
dunnTest(lesion.count ~ mutant, data = Alc, method ="bh")
# #post-hoc test Wilcoxon rank sum test
#will allow us to easily add significance levels to our ggplots
library(ggpubr)
A.lc.w <- compare_means(lesion.count ~ mutant, Alc, method="wilcox.test",
                       ref.group = ".all.") #issue with this code, will come back to it

#Wilcoxon test with only rep 1
compare_means(lesion.count ~ mutant, Alc.1, method="wilcox.test",
              ref.group = ".all.") 
#Wilcoxon test with only rep 2
compare_means(lesion.count ~ mutant, Alc.2, method="wilcox.test",
              ref.group = ".all.") 

## Punch inoculation ----
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

### ####Aichi ggplots####
#ggplot lesion size
library(ggpubr)
library(ggplot2)
als0 <- ggplot(A,aes(x=factor(mutant), y=lesion.surface, color = effector)) + 
  geom_boxplot() + geom_jitter(aes(size = lesion.surface), show.legend = TRUE) + # jitter size by lesion surface
  scale_size_continuous(range = c(0.01, 1.5))  + theme_minimal() + scale_x_discrete(limits=rev) + facet_grid(~rep)
Alsplot <- als0 +
  labs(title = "FR13 lesion size on Aichi Asahi", x = "Mutant", y = "Lesion size") + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
# stat_compare_means(label = "p.signif", method = "wilcox.test",
#                    ref.group = ".all.", hide.ns = TRUE) 

Alsplot
#ggplot lesion number
library(ggpubr)
alc0 <- ggplot(Alc, aes(x=factor(mutant), y=lesion.count, color = effector)) + 
  geom_boxplot() + geom_jitter() + theme_minimal() + scale_x_discrete(limits=rev) + facet_grid(~rep)
Alcplot <- alc0 +
  labs(title = "FR13 lesion number on Aichi Asahi", x = "Mutant", y = "Lesion number") + 
   theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
  # + stat_compare_means(label = "p.signif", method = "wilcox.test",
  #                    ref.group = ".all.", hide.ns = TRUE, label.y = c(50, 300, 300,350, 50, 50, 60)) 
library(patchwork)
Aplot <- (Alsplot + Alcplot)
Aplot
ggsave(filename = "A.png", plot = Aplot, device = "png", height = 15, width = 35,
       units = "cm", dpi = 500)


library(ggpubr)
library(ggplot2)
ap0 <- ggplot(A.p,aes(x=factor(mutant), y=lesion.surface, color = effector)) + 
  geom_boxplot() + geom_jitter(aes(size = lesion.surface), show.legend = TRUE) + # jitter size by lesion surface
  scale_size_continuous(range = c(0.01, 3))  + theme_minimal() + scale_x_discrete(limits=rev)
A.punch <- ap0 +
  labs(title = "Punch inoculation lesion surface on Aichi Asahi", x = "Mutant", y = "Lesion size") + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
A.punch
ggsave(filename = "A.punch.png", plot = A.punch, device = "png", height = 15, width = 35,
       units = "cm", dpi = 500)