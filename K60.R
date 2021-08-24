#K60
load("spray.RData")
load("lc.RData")
load("punch.RData")

K6 <- spray[spray$cultivar=="K60",] ; titre <- "K60"
K6.1 <- dplyr::filter(K6, rep == "1") #B for rep 1
K6.2 <- dplyr::filter(K6, rep == "2") #B for rep 2
K6lc <- lc[lc$cultivar =="K60",] ; titre <- "K60"
K6lc.1 <-dplyr::filter(K6lc, rep == "1") # Blc for rep 1
K6lc.2 <-dplyr::filter(K6lc, rep == "2") # Blc for rep 2
K6.p <- punch[punch$cultivar=="K60",] ; titre <- "K60"
# Stats ----

## Lesion surface ----
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

#Kruskal test with only rep 1
kruskal.test(lesion.surface ~ mutant, data = K6.1) #p-value <2e-16
#Kruskal test with only rep 2
kruskal.test(lesion.surface ~ mutant, data = K6.2) #p-value <2e-16
#same values for both

#Post-hoc 
#post-hoc test Dunn
library(FSA)
dunnTest(lesion.surface ~ mutant, data = K6, method ="bh")
# #post-hoc test Wilcoxon rank sum test
#will allow us to easily add significance levels to our ggplots
library(ggpubr)
K6.ls.w <- compare_means(lesion.surface ~ mutant, K6, method="wilcox.test",
                        ref.group = ".all.")
K6.ls.w
#Wilcoxon test with only rep 1
compare_means(lesion.surface ~ mutant, K6.1, method="wilcox.test",
              ref.group = ".all.") 
#Wilcoxon test with only rep 2
compare_means(lesion.surface ~ mutant, K6.2, method="wilcox.test",
              ref.group = ".all.") 
#actually fairly different results... 

# Lesion number ---- 
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

#Kruskal test with only rep 1
kruskal.test(lesion.count ~ mutant, data = K6lc.1) #p-value <2e-16
#Kruskal test with only rep 2
kruskal.test(lesion.count ~ mutant, data = K6lc.2) #p-value <2e-16
#same values for both

#Post-hoc 

#post-hoc test Dunn
library(FSA)
dunnTest(lesion.count ~ mutant, data = K6lc, method ="bh")
# #post-hoc test Wilcoxon rank sum test
#will allow us to easily add significance levels to our ggplots
library(ggpubr)
K6.lc.w <- compare_means(lesion.count ~ mutant, K6lc, method="wilcox.test",
                        ref.group = ".all.") #issue with this code, will come back to it

#Wilcoxon test with only rep 1
compare_means(lesion.count ~ mutant, K6lc.1, method="wilcox.test",
              ref.group = ".all.") 
#Wilcoxon test with only rep 2
compare_means(lesion.count ~ mutant, K6lc.2, method="wilcox.test",
              ref.group = ".all.") 

## Punch inoculation ----
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

# Ggplots ----

#ggplot lesion size
K6F <- dplyr::filter(K6, isolate == "FR13")
K6G <- dplyr::filter(K6, isolate == "Guy11")
library(ggpubr)
K6f.ls.p <- ggplot(K6F,aes(x=factor(mutant), y=lesion.surface, color = effector)) + geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface), show.legend = TRUE) + # jitter size by lesion surface
  scale_size_continuous(range = c(0.01, 2))  + theme_minimal() + scale_x_discrete(limits=rev) + 
  facet_grid(~rep)
K6F.ls.plot <- K6f.ls.p +
  labs(title = "FR13 lesion surface", x = "Mutant", y = "Lesion surface (in pixels)") +
  # theme(axis.title.y = element_blank(), axis.title.x = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
# stat_compare_means(label = "p.signif", method = "wilcox.test",
#                    ref.group = ".all.", hide.ns = TRUE) 
K6g.ls.p <- ggplot(K6G,aes(x=factor(mutant), y=lesion.surface, color = effector)) + geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface), show.legend = TRUE) + # jitter size by lesion surface
  scale_size_continuous(range = c(0.01, 2))  + theme_minimal() + scale_x_discrete(limits=rev) +
  facet_grid(~rep)
K6G.ls.plot <- K6g.ls.p +
  labs(title = "Guy11 lesion surface", x = "Mutant", y = "Lesion surface (in pixels)") +   
  #theme(axis.title.y = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0))   
library(patchwork)
K6lsplot <- (K6F.ls.plot / K6G.ls.plot) + plot_annotation(
  title = 'Lesion surface on K60')
K6lsplot
ggsave(filename = "K6ls.png", plot = K6lsplot, device = "png", height = 20, width = 30,
       units = "cm", dpi = 500)

#ggplot lesion number
K6Flc <- dplyr::filter(K6lc, isolate == "FR13")
K6Glc <- dplyr::filter(K6lc, isolate == "Guy11")
library(ggpubr)
k6f.lc0 <- ggplot(K6Flc, aes(x=factor(mutant), y=lesion.count, color = effector)) + 
  geom_boxplot() + geom_jitter() + theme_minimal() + scale_x_discrete(limits=rev) +
  facet_grid(~rep)
K6Flc <- k6f.lc0 +
  labs(title = "FR13 lesion number", x = "Mutant", y = "Lesion number") + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
# + stat_compare_means(label = "p.signif", method = "wilcox.test",
#                      ref.group = ".all.", hide.ns = TRUE, label.y = c(40, 60, 200, 30,30)) 
k6g.lc0 <- ggplot(K6Glc, aes(x=factor(mutant), y=lesion.count, color = effector)) + 
  geom_boxplot() + geom_jitter() + theme_minimal() + scale_x_discrete(limits=rev) +
  facet_grid(~rep)
K6glc <- k6g.lc0 +
  labs(title = "Guy11 lesion number", x = "Mutant", y = "Lesion number") + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
# + stat_compare_means(label = "p.signif", method = "wilcox.test",
#                      ref.group = ".all.", hide.ns = TRUE, label.y = c(125, 170, 170, 30,30)) 
library(patchwork)
K6lcplot <- (K6Flc / K6glc) + plot_annotation(
  title = 'Lesion number on K60')
K6lcplot
ggsave(filename = "K6lc.png", plot = K6lcplot, device = "png", height = 20, width = 30,
       units = "cm", dpi = 500)


#ggplot punch
K6Fp <- dplyr::filter(K6.p, isolate == "FR13")
K6Gp <- dplyr::filter(K6.p, isolate == "Guy11")
library(ggpubr)
library(ggplot2)
k6f.p.p <- ggplot(K6Fp,aes(x=factor(mutant), y=lesion.surface, color = effector)) + geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface), show.legend = TRUE) + # jitter size by lesion surface
  scale_size_continuous(range = c(0.01, 2))  + theme_minimal() + scale_x_discrete(limits=rev)
K6F.p.plot <- k6f.p.p +
  labs(title = "FR13 punch inoculation", x = "Mutant", y = "Lesion surface (in pixels)") +
  # theme(axis.title.y = element_blank(), axis.title.x = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
k6g.p.p <- ggplot(K6Gp,aes(x=factor(mutant), y=lesion.surface, color = effector)) + geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface), show.legend = TRUE) + # jitter size by lesion surface
  scale_size_continuous(range = c(0.01, 2))  + theme_minimal() + scale_x_discrete(limits=rev)
K6G.p.plot <- k6g.p.p +
  labs(title = "Guy11 punch inoculation", x = "Mutant", y = "Lesion surface (in pixels)") +   
  #theme(axis.title.y = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0))   
library(patchwork)
K6pplot <- (K6F.p.plot / K6G.p.plot) + plot_annotation(
  title = 'Punch inoculation on K60')
K6pplot
ggsave(filename = "K6p.png", plot = K6pplot, device = "png", height = 20, width = 30,
       units = "cm", dpi = 500)