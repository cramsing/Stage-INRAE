#Shin 2
load("spray.RData")
load("lc.RData")
load("punch.RData")

S <- spray[spray$cultivar=="Shin 2",] ; titre <- "Shin 2"
S.1 <- dplyr::filter(S, rep == "1") #B for rep 1
S.2 <- dplyr::filter(S, rep == "2") #B for rep 2
Slc <- lc[lc$cultivar =="Shin 2",] ; titre <- "Shin 2"
Slc.1 <-dplyr::filter(Slc, rep == "1") # Blc for rep 1
Slc.2 <-dplyr::filter(Slc, rep == "2") # Blc for rep 2
S.p <- punch[punch$cultivar=="Shin2",] ; titre <- "Shin 2"
# Stats ----

## Lesion surface ----
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
S.ls.w
#Wilcoxon test with only rep 1
compare_means(lesion.surface ~ mutant, S.1, method="wilcox.test",
              ref.group = ".all.") 
#Wilcoxon test with only rep 2
compare_means(lesion.surface ~ mutant, S.2, method="wilcox.test",
              ref.group = ".all.") 
#actually fairly different results... 

# Lesion number ---- 
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

#Wilcoxon test with only rep 1
compare_means(lesion.count ~ mutant, Slc.1, method="wilcox.test",
              ref.group = ".all.") 
#Wilcoxon test with only rep 2
compare_means(lesion.count ~ mutant, Slc.2, method="wilcox.test",
              ref.group = ".all.") 

## Punch inoculation ----
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

# Ggplots ----

#ggplot lesion size
SF <- dplyr::filter(S, isolate == "FR13")
SG <- dplyr::filter(S, isolate == "Guy11")
library(ggpubr)
Sf.ls.p <- ggplot(SF,aes(x=factor(mutant), y=lesion.surface, color = effector)) + geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface), show.legend = TRUE) + # jitter size by lesion surface
  scale_size_continuous(range = c(0.01, 2))  + theme_minimal() + scale_x_discrete(limits=rev) + 
  facet_grid(~rep)
SF.ls.plot <- Sf.ls.p +
  labs(title = "FR13 lesion surface", x = "Mutant", y = "Lesion surface (in pixels)") +
  # theme(axis.title.y = element_blank(), axis.title.x = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
# stat_compare_means(label = "p.signif", method = "wilcox.test",
#                    ref.group = ".all.", hide.ns = TRUE) 
Sg.ls.p <- ggplot(SG,aes(x=factor(mutant), y=lesion.surface, color = effector)) + geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface), show.legend = TRUE) + # jitter size by lesion surface
  scale_size_continuous(range = c(0.01, 2))  + theme_minimal() + scale_x_discrete(limits=rev) +
  facet_grid(~rep)
SG.ls.plot <- Sg.ls.p +
  labs(title = "Guy11 lesion surface", x = "Mutant", y = "Lesion surface (in pixels)") +   
  #theme(axis.title.y = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0))   
library(patchwork)
Slsplot <- (SF.ls.plot / SG.ls.plot) + plot_annotation(
  title = 'Lesion surface on Shin 2')
Slsplot
ggsave(filename = "Sls.png", plot = Slsplot, device = "png", height = 20, width = 30,
       units = "cm", dpi = 500)

#ggplot lesion number
SFlc <- dplyr::filter(Slc, isolate == "FR13")
SGlc <- dplyr::filter(Slc, isolate == "Guy11")
library(ggpubr)
Sf.lc0 <- ggplot(SFlc, aes(x=factor(mutant), y=lesion.count, color = effector)) + 
  geom_boxplot() + geom_jitter() + theme_minimal() + scale_x_discrete(limits=rev) +
  facet_grid(~rep)
SFlc <- Sf.lc0 +
  labs(title = "FR13 lesion number", x = "Mutant", y = "Lesion number") + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
# + stat_compare_means(label = "p.signif", method = "wilcox.test",
#                      ref.group = ".all.", hide.ns = TRUE, label.y = c(40, 60, 200, 30,30)) 
Sg.lc0 <- ggplot(SGlc, aes(x=factor(mutant), y=lesion.count, color = effector)) + 
  geom_boxplot() + geom_jitter() + theme_minimal() + scale_x_discrete(limits=rev) +
  facet_grid(~rep)
Sglc <- Sg.lc0 +
  labs(title = "Guy11 lesion number", x = "Mutant", y = "Lesion number") + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
# + stat_compare_means(label = "p.signif", method = "wilcox.test",
#                      ref.group = ".all.", hide.ns = TRUE, label.y = c(125, 170, 170, 30,30)) 
library(patchwork)
Slcplot <- (SFlc / Sglc) + plot_annotation(
  title = 'Lesion number on Shin 2')
Slcplot
ggsave(filename = "Slc.png", plot = Slcplot, device = "png", height = 20, width = 30,
       units = "cm", dpi = 500)


#ggplot punch
SFp <- dplyr::filter(S.p, isolate == "FR13")
SGp <- dplyr::filter(S.p, isolate == "Guy11")
library(ggpubr)
library(ggplot2)
Sf.p.p <- ggplot(SFp,aes(x=factor(mutant), y=lesion.surface, color = effector)) + geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface), show.legend = TRUE) + # jitter size by lesion surface
  scale_size_continuous(range = c(0.01, 2))  + theme_minimal() + scale_x_discrete(limits=rev)
SF.p.plot <- Sf.p.p +
  labs(title = "FR13 punch inoculation", x = "Mutant", y = "Lesion surface (in pixels)") +
  # theme(axis.title.y = element_blank(), axis.title.x = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
Sg.p.p <- ggplot(SGp,aes(x=factor(mutant), y=lesion.surface, color = effector)) + geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface), show.legend = TRUE) + # jitter size by lesion surface
  scale_size_continuous(range = c(0.01, 2))  + theme_minimal() + scale_x_discrete(limits=rev)
SG.p.plot <- Sg.p.p +
  labs(title = "Guy11 punch inoculation", x = "Mutant", y = "Lesion surface (in pixels)") +   
  #theme(axis.title.y = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0))   
library(patchwork)
Spplot <- (SF.p.plot / SG.p.plot) + plot_annotation(
  title = 'Punch inoculation on Shin 2')
Spplot
ggsave(filename = "Sp.png", plot = Spplot, device = "png", height = 20, width = 30,
       units = "cm", dpi = 500)