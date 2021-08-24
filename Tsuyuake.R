#Tsuyuake
load("spray.RData")
load("lc.RData")
#no Tsuyuake used in the punch inoculations

T <- spray[spray$cultivar=="Tsuyuake",] ; titre <- "Tsuyuake"
T.1 <- dplyr::filter(T, rep == "1") #B for rep 1
T.2 <- dplyr::filter(T, rep == "2") #B for rep 2
Tlc <- lc[lc$cultivar =="Tsuyuake",] ; titre <- "Tsuyuake"
Tlc.1 <-dplyr::filter(Tlc, rep == "1") # Blc for rep 1
Tlc.2 <-dplyr::filter(Tlc, rep == "2") # Blc for rep 2

# Stats ----

## Lesion surface ----
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

#Kruskal test with only rep 1
kruskal.test(lesion.surface ~ mutant, data = T.1) #p-value <2e-16
#Kruskal test with only rep 2
kruskal.test(lesion.surface ~ mutant, data = T.2) #p-value 1e-11
#not the same values for both

#Post-hoc 
#post-hoc test Dunn
library(FSA)
dunnTest(lesion.surface ~ mutant, data = T, method ="bh")
# #post-hoc test Wilcoxon rank sum test
#will allow us to easily add significance levels to our ggplots
library(ggpubr)
T.ls.w <- compare_means(lesion.surface ~ mutant, T, method="wilcox.test",
                         ref.group = ".all.")
T.ls.w
#Wilcoxon test with only rep 1
compare_means(lesion.surface ~ mutant, T.1, method="wilcox.test",
              ref.group = ".all.") 
#Wilcoxon test with only rep 2
compare_means(lesion.surface ~ mutant, T.2, method="wilcox.test",
              ref.group = ".all.") 
#actually fairly different results... 

# Lesion number ---- 
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

#Kruskal test with only rep 1
kruskal.test(lesion.count ~ mutant, data = Tlc.1) #p-value <2e-16
#Kruskal test with only rep 2
kruskal.test(lesion.count ~ mutant, data = Tlc.2) #p-value <2e-16
#same values for both

#Post-hoc 

#post-hoc test Dunn
library(FSA)
dunnTest(lesion.count ~ mutant, data = Tlc, method ="bh")
# #post-hoc test Wilcoxon rank sum test
#will allow us to easily add significance levels to our ggplots
library(ggpubr)
T.lc.w <- compare_means(lesion.count ~ mutant, Tlc, method="wilcox.test",
                         ref.group = ".all.") #issue with this code, will come back to it

#Wilcoxon test with only rep 1
compare_means(lesion.count ~ mutant, Tlc.1, method="wilcox.test",
              ref.group = ".all.") 
#Wilcoxon test with only rep 2
compare_means(lesion.count ~ mutant, Tlc.2, method="wilcox.test",
              ref.group = ".all.") 



# Ggplots ----

#ggplot lesion size
library(ggpubr)
library(ggplot2)
tls0 <- ggplot(T,aes(x=factor(mutant), y=lesion.surface, color = effector)) + 
  geom_boxplot() + geom_jitter(aes(size = lesion.surface), show.legend = TRUE) + # jitter size by lesion surface
  scale_size_continuous(range = c(0.01, 1.5))  + theme_minimal() + scale_x_discrete(limits=rev) + facet_grid(~rep)
Tlsplot <- tls0 +
  labs(title = "FR13 lesion size on Tsuyuake", x = "Mutant", y = "Lesion size") + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
# stat_compare_means(label = "p.signif", method = "wilcox.test",
#                    ref.group = ".all.", hide.ns = TRUE) 

Tlsplot
#ggplot lesion number
library(ggpubr)
tlc0 <- ggplot(Tlc, aes(x=factor(mutant), y=lesion.count, color = effector)) + 
  geom_boxplot() + geom_jitter() + theme_minimal() + scale_x_discrete(limits=rev) + facet_grid(~rep)
Tlcplot <- tlc0 +
  labs(title = "FR13 lesion number on Tsuyuake", x = "Mutant", y = "Lesion number") + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
# + stat_compare_means(label = "p.signif", method = "wilcox.test",
#                    ref.group = ".all.", hide.ns = TRUE, label.y = c(50, 300, 300,350, 50, 50, 60)) 
library(patchwork)
Tplot <- (Tlsplot + Tlcplot)
Tplot
ggsave(filename = "T.png", plot = Tplot, device = "png", height = 15, width = 35,
       units = "cm", dpi = 500)
