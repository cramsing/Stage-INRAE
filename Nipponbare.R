#Nipponbare
load("spray.RData")
load("lc.RData")
load("punch.RData")

N <- spray[spray$cultivar=="Nipponbare",] ; titre <- "Nipponbare"
N.1 <- dplyr::filter(N, rep == "1") #B for rep 1
N.2 <- dplyr::filter(N, rep == "2") #B for rep 2
Nlc <- lc[lc$cultivar =="Nipponbare",] ; titre <- "Nipponbare"
Nlc.1 <-dplyr::filter(Nlc, rep == "1") # Blc for rep 1
Nlc.2 <-dplyr::filter(Nlc, rep == "2") # Blc for rep 2
N.p <- punch[punch$cultivar=="Nipponbare",] ; titre <- "Nipponbare"
# Stats ----

## Lesion surface ----
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

#Kruskal test with only rep 1
kruskal.test(lesion.surface ~ mutant, data = N.1) #p-value <2e-16
#Kruskal test with only rep 2
kruskal.test(lesion.surface ~ mutant, data = N.2) #p-value <2e-16
#same values for both

#Post-hoc 
#post-hoc test Dunn
library(FSA)
dunnTest(lesion.surface ~ mutant, data = N, method ="bh")
# #post-hoc test Wilcoxon rank sum test
#will allow us to easily add significance levels to our ggplots
library(ggpubr)
N.ls.w <- compare_means(lesion.surface ~ mutant, N, method="wilcox.test",
                         ref.group = ".all.")
N.ls.w
#Wilcoxon test with only rep 1
compare_means(lesion.surface ~ mutant, N.1, method="wilcox.test",
              ref.group = ".all.") 
#Wilcoxon test with only rep 2
compare_means(lesion.surface ~ mutant, N.2, method="wilcox.test",
              ref.group = ".all.") 
#actually fairly different results... 

# Lesion number ---- 
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

#Kruskal test with only rep 1
kruskal.test(lesion.count ~ mutant, data = Nlc.1) #p-value <2e-16
#Kruskal test with only rep 2
kruskal.test(lesion.count ~ mutant, data = Nlc.2) #p-value <2e-16
#same values for both

#Post-hoc 

#post-hoc test Dunn
library(FSA)
dunnTest(lesion.count ~ mutant, data = Nlc, method ="bh")
# #post-hoc test Wilcoxon rank sum test
#will allow us to easily add significance levels to our ggplots
library(ggpubr)
N.lc.w <- compare_means(lesion.count ~ mutant, Nlc, method="wilcox.test",
                         ref.group = ".all.") #issue with this code, will come back to it

#Wilcoxon test with only rep 1
compare_means(lesion.count ~ mutant, Nlc.1, method="wilcox.test",
              ref.group = ".all.") 
#Wilcoxon test with only rep 2
compare_means(lesion.count ~ mutant, Nlc.2, method="wilcox.test",
              ref.group = ".all.") 

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

# Ggplots ----

#ggplot lesion size
NF <- dplyr::filter(N, isolate == "FR13")
NG <- dplyr::filter(N, isolate == "Guy11")
library(ggpubr)
Nf.ls.p <- ggplot(NF,aes(x=factor(mutant), y=lesion.surface, color = effector)) + geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface), show.legend = TRUE) + # jitter size by lesion surface
  scale_size_continuous(range = c(0.01, 2))  + theme_minimal() + scale_x_discrete(limits=rev) + 
  facet_grid(~rep)
NF.ls.plot <- Nf.ls.p +
  labs(title = "FR13 lesion surface", x = "Mutant", y = "Lesion surface (in pixels)") +
  # theme(axis.title.y = element_blank(), axis.title.x = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
# stat_compare_means(label = "p.signif", method = "wilcox.test",
#                    ref.group = ".all.", hide.ns = TRUE) 
Ng.ls.p <- ggplot(NG,aes(x=factor(mutant), y=lesion.surface, color = effector)) + geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface), show.legend = TRUE) + # jitter size by lesion surface
  scale_size_continuous(range = c(0.01, 2))  + theme_minimal() + scale_x_discrete(limits=rev) +
  facet_grid(~rep)
NG.ls.plot <- Ng.ls.p +
  labs(title = "Guy11 lesion surface", x = "Mutant", y = "Lesion surface (in pixels)") +   
  #theme(axis.title.y = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0))   
library(patchwork)
Nlsplot <- (NF.ls.plot / NG.ls.plot) + plot_annotation(
  title = 'Lesion surface on Nipponbare')
Nlsplot
ggsave(filename = "Nls.png", plot = Nlsplot, device = "png", height = 20, width = 30,
       units = "cm", dpi = 500)

#ggplot lesion number
NFlc <- dplyr::filter(Nlc, isolate == "FR13")
NGlc <- dplyr::filter(Nlc, isolate == "Guy11")
library(ggpubr)
Nf.lc0 <- ggplot(NFlc, aes(x=factor(mutant), y=lesion.count, color = effector)) + 
  geom_boxplot() + geom_jitter() + theme_minimal() + scale_x_discrete(limits=rev) +
  facet_grid(~rep)
NFlc <- Nf.lc0 +
  labs(title = "FR13 lesion number", x = "Mutant", y = "Lesion number") + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
# + stat_compare_means(label = "p.signif", method = "wilcox.test",
#                      ref.group = ".all.", hide.ns = TRUE, label.y = c(40, 60, 200, 30,30)) 
Ng.lc0 <- ggplot(NGlc, aes(x=factor(mutant), y=lesion.count, color = effector)) + 
  geom_boxplot() + geom_jitter() + theme_minimal() + scale_x_discrete(limits=rev) +
  facet_grid(~rep)
Nglc <- Ng.lc0 +
  labs(title = "Guy11 lesion number", x = "Mutant", y = "Lesion number") + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
# + stat_compare_means(label = "p.signif", method = "wilcox.test",
#                      ref.group = ".all.", hide.ns = TRUE, label.y = c(125, 170, 170, 30,30)) 
library(patchwork)
Nlcplot <- (NFlc / Nglc) + plot_annotation(
  title = 'Lesion number on Nipponbare')
Nlcplot
ggsave(filename = "Nlc.png", plot = Nlcplot, device = "png", height = 20, width = 30,
       units = "cm", dpi = 500)


#ggplot punch
NFp <- dplyr::filter(N.p, isolate == "FR13")
NGp <- dplyr::filter(N.p, isolate == "Guy11")
library(ggpubr)
library(ggplot2)
Nf.p.p <- ggplot(NFp,aes(x=factor(mutant), y=lesion.surface, color = effector)) + geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface), show.legend = TRUE) + # jitter size by lesion surface
  scale_size_continuous(range = c(0.01, 2))  + theme_minimal() + scale_x_discrete(limits=rev)
NF.p.plot <- Nf.p.p +
  labs(title = "FR13 punch inoculation", x = "Mutant", y = "Lesion surface (in pixels)") +
  # theme(axis.title.y = element_blank(), axis.title.x = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
Ng.p.p <- ggplot(NGp,aes(x=factor(mutant), y=lesion.surface, color = effector)) + geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface), show.legend = TRUE) + # jitter size by lesion surface
  scale_size_continuous(range = c(0.01, 2))  + theme_minimal() + scale_x_discrete(limits=rev)
NG.p.plot <- Ng.p.p +
  labs(title = "Guy11 punch inoculation", x = "Mutant", y = "Lesion surface (in pixels)") +   
  #theme(axis.title.y = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0))   
library(patchwork)
Npplot <- (NF.p.plot / NG.p.plot) + plot_annotation(
  title = 'Punch inoculation on Nipponbare')
Npplot
ggsave(filename = "Np.png", plot = Npplot, device = "png", height = 20, width = 30,
       units = "cm", dpi = 500)