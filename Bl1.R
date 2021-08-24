#Bl1
load("spray.RData")
load("lc.RData")
load("punch.RData")

B <- spray[spray$cultivar=="Bl1",] ; titre <- "BL1"
B.1 <- dplyr::filter(B, rep == "1") #B for rep 1
B.2 <- dplyr::filter(B, rep == "2") #B for rep 2
Blc <- lc[lc$cultivar =="Bl1",] ; titre <- "BL1"
Blc.1 <-dplyr::filter(Blc, rep == "1") # Blc for rep 1
Blc.2 <-dplyr::filter(Blc, rep == "2") # Blc for rep 2
B.p <- punch[punch$cultivar=="Bl1",] ; titre <- "Bl1"
# Stats ----

## Lesion surface ----
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

#Kruskal test with only rep 1
kruskal.test(lesion.surface ~ mutant, data = B.1) #p-value <2e-16
#Kruskal test with only rep 2
kruskal.test(lesion.surface ~ mutant, data = B.2) #p-value <2e-16
#same values for both

#Post-hoc 
#post-hoc test Dunn
library(FSA)
dunnTest(lesion.surface ~ mutant, data = B, method ="bh")
# #post-hoc test Wilcoxon rank sum test
#will allow us to easily add significance levels to our ggplots
library(ggpubr)
B.ls.w <- compare_means(lesion.surface ~ mutant, B, method="wilcox.test",
                        ref.group = ".all.")
B.ls.w
#Wilcoxon test with only rep 1
compare_means(lesion.surface ~ mutant, B.1, method="wilcox.test",
              ref.group = ".all.") 
#Wilcoxon test with only rep 2
compare_means(lesion.surface ~ mutant, B.2, method="wilcox.test",
              ref.group = ".all.") 
#actually fairly different results... 

# Lesion number ---- 
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

## Punch inoculation ----
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

# Ggplots ----

#ggplot lesion size
BF <- dplyr::filter(B, isolate == "FR13")
BG <- dplyr::filter(B, isolate == "Guy11")
library(ggpubr)
bf.ls.p <- ggplot(BF,aes(x=factor(mutant), y=lesion.surface, color = effector)) + geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface), show.legend = TRUE) + # jitter size by lesion surface
  scale_size_continuous(range = c(0.01, 2))  + theme_minimal() + scale_x_discrete(limits=rev) + 
  facet_grid(~rep)
BF.ls.plot <- bf.ls.p +
  labs(title = "FR13 lesion surface", x = "Mutant", y = "Lesion surface (in pixels)") +
  # theme(axis.title.y = element_blank(), axis.title.x = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
# stat_compare_means(label = "p.signif", method = "wilcox.test",
#                    ref.group = ".all.", hide.ns = TRUE) 
bg.ls.p <- ggplot(BG,aes(x=factor(mutant), y=lesion.surface, color = effector)) + geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface), show.legend = TRUE) + # jitter size by lesion surface
  scale_size_continuous(range = c(0.01, 2))  + theme_minimal() + scale_x_discrete(limits=rev) +
  facet_grid(~rep)
BG.ls.plot <- bg.ls.p +
  labs(title = "Guy11 lesion surface", x = "Mutant", y = "Lesion surface (in pixels)") +   
  #theme(axis.title.y = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0))   
library(patchwork)
Blsplot <- (BF.ls.plot / BG.ls.plot) + plot_annotation(
  title = 'Lesion surface on Bl1')
Blsplot
ggsave(filename = "Bls.png", plot = Blsplot, device = "png", height = 20, width = 30,
       units = "cm", dpi = 500)

#ggplot lesion number
BFlc <- dplyr::filter(Blc, isolate == "FR13")
BGlc <- dplyr::filter(Blc, isolate == "Guy11")
library(ggpubr)
bf.lc0 <- ggplot(BFlc, aes(x=factor(mutant), y=lesion.count, color = effector)) + 
  geom_boxplot() + geom_jitter() + theme_minimal() + scale_x_discrete(limits=rev) +
  facet_grid(~rep)
BFlc <- bf.lc0 +
  labs(title = "FR13 lesion number", x = "Mutant", y = "Lesion number") + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
# + stat_compare_means(label = "p.signif", method = "wilcox.test",
#                      ref.group = ".all.", hide.ns = TRUE, label.y = c(40, 60, 200, 30,30)) 
bg.lc0 <- ggplot(BGlc, aes(x=factor(mutant), y=lesion.count, color = effector)) + 
  geom_boxplot() + geom_jitter() + theme_minimal() + scale_x_discrete(limits=rev) +
  facet_grid(~rep)
Bglc <- bg.lc0 +
  labs(title = "Guy11 lesion number", x = "Mutant", y = "Lesion number") + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
# + stat_compare_means(label = "p.signif", method = "wilcox.test",
#                      ref.group = ".all.", hide.ns = TRUE, label.y = c(125, 170, 170, 30,30)) 
library(patchwork)
Blcplot <- (BFlc / Bglc) + plot_annotation(
  title = 'Lesion number on Bl1')
Blcplot
ggsave(filename = "Blc.png", plot = Blcplot, device = "png", height = 20, width = 30,
       units = "cm", dpi = 500)


#ggplot punch
BFp <- dplyr::filter(B.p, isolate == "FR13")
BGp <- dplyr::filter(B.p, isolate == "Guy11")
library(ggpubr)
library(ggplot2)
bf.p.p <- ggplot(BFp,aes(x=factor(mutant), y=lesion.surface, color = effector)) + geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface), show.legend = TRUE) + # jitter size by lesion surface
  scale_size_continuous(range = c(0.01, 2))  + theme_minimal() + scale_x_discrete(limits=rev)
BF.p.plot <- bf.p.p +
  labs(title = "FR13 punch inoculation", x = "Mutant", y = "Lesion surface (in pixels)") +
  # theme(axis.title.y = element_blank(), axis.title.x = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
bg.p.p <- ggplot(BGp,aes(x=factor(mutant), y=lesion.surface, color = effector)) + geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface), show.legend = TRUE) + # jitter size by lesion surface
  scale_size_continuous(range = c(0.01, 2))  + theme_minimal() + scale_x_discrete(limits=rev)
BG.p.plot <- bg.p.p +
  labs(title = "Guy11 punch inoculation", x = "Mutant", y = "Lesion surface (in pixels)") +   
  #theme(axis.title.y = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0))   
library(patchwork)
Bpplot <- (BF.p.plot / BG.p.plot) + plot_annotation(
  title = 'Punch inoculation on Bl1')
Bpplot
ggsave(filename = "Bp.png", plot = Bpplot, device = "png", height = 20, width = 30,
       units = "cm", dpi = 500)