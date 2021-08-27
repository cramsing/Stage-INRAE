#Kasalath
load("spray.RData")
load("lc.RData")
load("punch.RData")

K <- spray[spray$cultivar=="Kasalath",] ; titre <- "Kasalath"
K.1 <- dplyr::filter(K, rep == "1") #B for rep 1
K.2 <- dplyr::filter(K, rep == "2") #B for rep 2
Klc <- lc[lc$cultivar =="Kasalath",] ; titre <- "Kasalath"
Klc.1 <-dplyr::filter(Klc, rep == "1") # Blc for rep 1
Klc.2 <-dplyr::filter(Klc, rep == "2") # Blc for rep 2
K.p <- punch[punch$cultivar=="Kasalath",] ; titre <- "Kasalath"
# Stats ----

## Lesion surface ----
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
#p-value <2.2e-16, significant !

#Kruskal test with only rep 1
kruskal.test(lesion.surface ~ mutant, data = K.1) #p-value <2e-16
#Kruskal test with only rep 2
kruskal.test(lesion.surface ~ mutant, data = K.2) #p-value <2e-16
#same values for both

#Post-hoc 
#post-hoc test Dunn
library(FSA)
dunnTest(lesion.surface ~ mutant, data = K, method ="bh")
# #post-hoc test Wilcoxon rank sum test
#will allow us to easily add significance levels to our ggplots
library(ggpubr)
K.ls.w <- compare_means(lesion.surface ~ mutant, K, method="wilcox.test",
                        ref.group = ".all.")
K.ls.w
#Wilcoxon test with only rep 1
compare_means(lesion.surface ~ mutant, K.1, method="wilcox.test",
              ref.group = ".all.") 
#Wilcoxon test with only rep 2
compare_means(lesion.surface ~ mutant, K.2, method="wilcox.test",
              ref.group = ".all.") 
#actually fairly similar.. ish... results 

# Lesion number ---- 
k1 <- aov(lesion.count ~ mutant, data = Klc) 
anova(k1)
#p-value = <2.2e-16 ***, significant
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
#9.973e-08, not normal
kruskal.test(lesion.count ~ mutant, data = Klc)
#p = 2.296e-16, significant

#Kruskal test with only rep 1
kruskal.test(lesion.count ~ mutant, data = Klc.1) #p-value <2e-16
#Kruskal test with only rep 2
kruskal.test(lesion.count ~ mutant, data = Klc.2) #p-value <2e-16
#same values for both

#Post-hoc 

#post-hoc test Dunn
library(FSA)
dunnTest(lesion.count ~ mutant, data = Klc, method ="bh")
# #post-hoc test Wilcoxon rank sum test
#will allow us to easily add significance levels to our ggplots
library(ggpubr)
K.lc.w <- compare_means(lesion.count ~ mutant, Klc, method="wilcox.test",
                        ref.group = ".all.") #issue with this code, will come back to it

#Wilcoxon test with only rep 1
compare_means(lesion.count ~ mutant, Klc.1, method="wilcox.test",
              ref.group = ".all.") 
#Wilcoxon test with only rep 2
compare_means(lesion.count ~ mutant, Klc.2, method="wilcox.test",
              ref.group = ".all.") 

## Punch inoculation ----
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

#check for each isolate
kFp <- dplyr::filter(K.p, isolate == "FR13")
kGp <- dplyr::filter(K.p, isolate == "Guy11")

k3 <- aov(lesion.surface ~ mutant, data = kFp) 
anova(k3)
#p-value = 0.0009117 ***, significant
#normality of residuals check
#distrubution/normality
par(mfrow=c(1,2))
#Histogramme
hist(k3$residuals)
# Q-Q Plot
qqnorm(k3$residuals)
qqline(k3$residuals)
#kind of normal looking
shapiro.test(k3$residuals)
# 0.5678, normal, anova stands!

k4 <- aov(lesion.surface ~ mutant, data = kGp) 
anova(k4)
#p-value =  1.895e-05 ***, significant
#normality of residuals check
#distrubution/normality
par(mfrow=c(1,2))
#Histogramme
hist(k4$residuals)
# Q-Q Plot
qqnorm(k4$residuals)
qqline(k4$residuals)
#kind of normal looking
shapiro.test(k4$residuals)
#0.06516 normal!


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


## Compiled kruskal groups ----
KF <- dplyr::filter(K, isolate == "FR13")
KG <- dplyr::filter(K, isolate == "Guy11")
#lesion surface kruskal groups
library(agricolae)
library(dplyr)
K.f.ls.groups <- kruskal(KF$lesion.surface,
                          KF$mutant, group = TRUE) %>%
  .$groups %>%
  as_tibble(rownames = "mutant") %>%
  rename("K.f.ls.groups" = "groups")
K.f.ls.groups$`KF$lesion.surface` <- NULL
K.g.ls.groups <- kruskal(KG$lesion.surface,
                          KG$mutant, group = TRUE) %>%
  .$groups %>%
  as_tibble(rownames = "mutant") %>%
  rename("K.g.ls.groups" = "groups")
K.g.ls.groups$`KG$lesion.surface` <- NULL


#lesion count kruskal groups
KFlc <- dplyr::filter(Klc, isolate == "FR13")
KGlc <- dplyr::filter(Klc, isolate == "Guy11")
library(agricolae)
library(dplyr)
K.f.lc.groups <- kruskal(KFlc$lesion.count,
                          KFlc$mutant, group = TRUE) %>%
  .$groups %>%
  as_tibble(rownames = "mutant") %>%
  rename("K.f.lc.groups" = "groups")
K.f.lc.groups$`KFlc$lesion.count`<- NULL
K.g.lc.groups <- kruskal(KGlc$lesion.count,
                          KGlc$mutant, group = TRUE) %>%
  .$groups %>%
  as_tibble(rownames = "mutant") %>%
  rename("K.g.lc.groups" = "groups")
K.g.lc.groups$`KGlc$lesion.count` <- NULL

#punch inoculation groups
# should use t-test bc the data is normal!
kFp <- dplyr::filter(K.p, isolate == "FR13")
kGp <- dplyr::filter(K.p, isolate == "Guy11")

k3 <- aov(lesion.surface ~ mutant, data = kFp) 
anova(k3)
library(agricolae)
K.f.p.groups <- HSD.test(k3, "mutant", group = TRUE) %>%
  .$groups %>%
  as_tibble(rownames = "mutant") %>%
  rename("K.f.p.groups" = "groups")
K.f.p.groups$lesion.surface<- NULL

k4 <- aov(lesion.surface ~ mutant, data = kGp) 
anova(k4)
K.g.p.groups <- HSD.test(k4, "mutant", group = TRUE) %>%
  .$groups %>%
  as_tibble(rownames = "mutant") %>%
  rename("K.g.p.groups" = "groups")
K.g.p.groups$lesion.surface<- NULL

#compilation
library(dplyr)
K.f.groups <- dplyr::left_join(K.f.ls.groups, K.f.lc.groups, by = "mutant")%>%
  left_join(., K.f.p.groups, by = "mutant" )
K.f.groups <- K.f.groups[order(K.f.groups$mutant),]
write.csv2(K.f.groups, file = "K.f.groups.csv")
K.g.groups <- dplyr::left_join(K.g.ls.groups, K.g.lc.groups, by = "mutant")%>%
  left_join(., K.f.p.groups, by = "mutant" )
K.g.groups <- K.g.p.groups[order(K.g.groups$mutant),] 
write.csv2(K.g.groups, file = "K.g.groups.csv")

# Ggplots ----
#ggplot lesion size
kF <- dplyr::filter(k, isolate == "FR13")
kG <- dplyr::filter(k, isolate == "Guy11")
library(ggpubr)
kf.ls.p <- ggplot(kF,aes(x=factor(mutant, level = f.level_order), y=lesion.surface, fill = effector)) +   
  geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface, shape = rep, color = rep), show.legend = TRUE) + 
  scale_size_continuous(range = c(0.1, 2))+ theme_minimal() + scale_x_discrete(limits=rev) + 
  scale_y_continuous(trans = "log10", breaks = c(10, 100, 500, 1000, 5000, 10000))  
kF.ls.plot <- kf.ls.p + scale_colour_grey() + scale_fill_manual(values = blue.palette) +
  labs(title = "FR13 lesion surface", x = "Isolate", y = "Lesion surface (in log10 pixels)") +
  # theme(axis.title.y = element_blank(), axis.title.x = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) +  coord_flip()
# stat_compare_means(label = "p.signif", method = "wilcox.test",
#                    ref.group = ".all.", hide.ns = TRUE) 
kg.ls.p <- ggplot(kG,aes(x=factor(mutant, level = g.level_order), y=lesion.surface, fill = effector)) + 
  geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface, shape = rep, color = rep), show.legend = TRUE) + 
  scale_size_continuous(range = c(0.1, 2))  + theme_minimal() + scale_x_discrete(limits=rev) +
  scale_y_continuous(trans = "log10", breaks = c(10, 100, 500, 1000, 5000, 10000)) 
kG.ls.plot <- kg.ls.p + scale_colour_grey() + scale_fill_manual(values = blue.palette) +
  labs(title = "Guy11 lesion surface", x = "Isolate", y = "Lesion surface (in log10 pixels)") +   
  #theme(axis.title.y = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) +  coord_flip() 
# library(patchwork)
# klsplot <- (kF.ls.plot/kG.ls.plot) + plot_annotation(title = 'Lesion surface on Kasalath')
# klsplot
ggsave(filename = "kfls.png", plot = kF.ls.plot, device = "png", height = 20, width = 30,
       units = "cm", dpi = 500)
ggsave(filename = "kgls.png", plot = kG.ls.plot, device = "png", height = 20, width = 30,
       units = "cm", dpi = 500)

#ggplot lesion number
kFlc <- dplyr::filter(klc, isolate == "FR13")
kGlc <- dplyr::filter(klc, isolate == "Guy11")
library(ggpubr)
kf.lc0 <- ggplot(kFlc, aes(x=factor(mutant, level = f.level_order), y=lesion.count, fill = effector)) + 
  geom_boxplot() + geom_jitter(size=2.0, aes(shape = rep, color = rep)) + 
  theme_minimal() + scale_x_discrete(limits=rev) + coord_flip()
kFlc <- kf.lc0 +  scale_colour_grey() + scale_fill_manual(values = blue.palette) +
  labs(title = "FR13 lesion number", x = "Isolate", y = "Lesion number") + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
# + stat_compare_means(label = "p.signif", method = "wilcox.test",
#                      ref.group = ".all.", hide.ns = TRUE, label.y = c(40, 60, 200, 30,30)) 
kg.lc0 <- ggplot(kGlc, aes(x=factor(mutant, level = g.level_order), y=lesion.count, fill = effector)) + 
  geom_boxplot() + geom_jitter(size=2.0, aes(shape = rep, color = rep)) + 
  theme_minimal() + scale_x_discrete(limits=rev) + coord_flip()
kglc <- kg.lc0 +  scale_colour_grey() + scale_fill_manual(values = blue.palette) +
  labs(title = "Guy11 lesion number", x = "Isolate", y = "Lesion number") + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
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
#ggplot punch
kFp <- dplyr::filter(k.p, isolate == "FR13")
kGp <- dplyr::filter(k.p, isolate == "Guy11")
library(ggpubr)
library(ggplot2)
kf.p.p <- ggplot(kFp,aes(x=factor(mutant, level = f.level_order), y=lesion.surface, fill = effector)) + geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface), show.legend = TRUE) + # jitter size by lesion surface
  scale_size_continuous(range = c(0.01, 2))  + theme_minimal() + scale_x_discrete(limits=rev)
kF.p.plot <- kf.p.p + scale_fill_manual(values = blue.palette) +
  labs(title = "FR13 punch inoculation", x = "Isolate", y = "Lesion surface (in pixels)") +
  # theme(axis.title.y = element_blank(), axis.title.x = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
kg.p.p <- ggplot(kGp, aes(x=factor(mutant, level = g.level_order), y=lesion.surface, fill = effector)) + geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface), show.legend = TRUE) + # jitter size by lesion surface
  scale_size_continuous(range = c(0.01, 2))  + theme_minimal() + scale_x_discrete(limits=rev)
kG.p.plot <- kg.p.p + scale_fill_manual(values = blue.palette) +
  labs(title = "Guy11 punch inoculation", x = "Isolate", y = "Lesion surface (in pixels)") +   
  #theme(axis.title.y = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0))   
library(patchwork)
kpplot <- (kF.p.plot / kG.p.plot) + plot_annotation(
  title = 'Punch inoculation on Kasalath')
kpplot
ggsave(filename = "kp.png", plot = kpplot, device = "png", height = 20, width = 30,
       units = "cm", dpi = 500)