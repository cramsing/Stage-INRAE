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
#p-value <2e-16, significant !

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

#punch inoculation kruskal groups

KFp <- dplyr::filter(K.p, isolate == "FR13")
KGp <- dplyr::filter(K.p, isolate == "Guy11")
library(agricolae)
library(dplyr)
K.f.p.groups <- kruskal(KFp$lesion.surface,
                         KFp$mutant, group = TRUE) %>%
  .$groups %>%
  as_tibble(rownames = "mutant") %>%
  rename("K.f.p.groups" = "groups")
K.f.p.groups$`KFp$lesion.surface`<- NULL
K.g.p.groups <- kruskal(KGp$lesion.surface,
                         KGp$mutant, group = TRUE) %>%
  .$groups %>%
  as_tibble(rownames = "mutant") %>%
  rename("K.g.p.groups" = "groups")
K.g.p.groups$`KGp$lesion.surface` <- NULL

#compilation
library(dplyr)
K.f.groups <- dplyr::left_join(K.f.ls.groups, K.f.lc.groups, by = "mutant") %>%
  left_join(., K.f.p.groups, by = "mutant" )
K.f.groups <- K.f.groups[order(K.f.groups$mutant),]
write.csv2(K.f.groups, file = "K.f.groups.csv")
K.g.groups <- dplyr::left_join(K.g.ls.groups, K.g.lc.groups, by = "mutant") %>%
  left_join(., K.g.p.groups, by = "mutant" )
K.g.groups <- K.g.groups[order(K.g.groups$mutant),] 
write.csv2(K.g.groups, file = "K.g.groups.csv")

# Ggplots ----

#ggplot lesion size
KF <- dplyr::filter(K, isolate == "FR13")
KG <- dplyr::filter(K, isolate == "Guy11")
library(ggpubr)
Kf.ls.p <- ggplot(KF,aes(x=factor(mutant), y=lesion.surface, color = effector)) + 
  geom_boxplot(color = "darkslategray") + 
  geom_jitter(aes(size = lesion.surface, shape = rep), show.legend = TRUE) + # jitter size by lesion surface
  scale_size_continuous(range = c(0.01, 1.5))  + theme_minimal() + scale_x_discrete(limits=rev) + 
  scale_y_continuous(trans = "log10") 
KF.ls.plot <- Kf.ls.p +
  labs(title = "FR13 lesion surface", x = "Mutant", y = "Lesion size (in log10 pixels)") +
  # theme(axis.title.y = element_blank(), axis.title.x = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
# stat_compare_means(label = "p.signif", method = "wilcox.test",
#                    ref.group = ".all.", hide.ns = TRUE) 
Kg.ls.p <- ggplot(KG,aes(x=factor(mutant), y=lesion.surface, color = effector)) + 
  geom_boxplot(color = "darkslategray") + 
  geom_jitter(aes(size = lesion.surface, shape = rep), show.legend = TRUE) + # jitter size by lesion surface
  scale_size_continuous(range = c(0.01, 1.5))  + theme_minimal() + scale_x_discrete(limits=rev) + 
  scale_y_continuous(trans = "log10") 
KG.ls.plot <- Kg.ls.p +
  labs(title = "Guy11 lesion surface", x = "Mutant", y = "Lesion size (in log10 pixels)") +   
  #theme(axis.title.y = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0))   
# library(patchwork)
# Klsplot <- (KF.ls.plot / KG.ls.plot) + plot_annotation(
#   title = 'Lesion surface on Kasalath')
# Klsplot
# ggsave(filename = "Kls.png", plot = Klsplot, device = "png", height = 20, width = 30,
#        units = "cm", dpi = 500)

#ggplot lesion number
KFlc <- dplyr::filter(Klc, isolate == "FR13")
KGlc <- dplyr::filter(Klc, isolate == "Guy11")
library(ggpubr)
kf.lc0 <- ggplot(KFlc, aes(x=factor(mutant), y=lesion.count, color = effector)) + 
  geom_boxplot(color = "darkslategray") + geom_jitter(size=1.0, aes(shape = rep.x)) + 
  theme_minimal() + scale_x_discrete(limits=rev)
KFlc <- kf.lc0 +
  labs(title = "FR13 lesion number", x = "Mutant", y = "Lesion number") + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
# + stat_compare_means(label = "p.signif", method = "wilcox.test",
#                      ref.group = ".all.", hide.ns = TRUE, label.y = c(40, 60, 200, 30,30)) 
kg.lc0 <- ggplot(KGlc, aes(x=factor(mutant), y=lesion.count, color = effector)) + 
  geom_boxplot(color = "darkslategray") + geom_jitter(size=1.0, aes(shape = rep.x)) + 
  theme_minimal() + scale_x_discrete(limits=rev)
Kglc <- kg.lc0 +
  labs(title = "Guy11 lesion number", x = "Mutant", y = "Lesion number") + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
# + stat_compare_means(label = "p.signif", method = "wilcox.test",
#                      ref.group = ".all.", hide.ns = TRUE, label.y = c(125, 170, 170, 30,30)) 
# library(patchwork)
# Klcplot <- (KFlc / Kglc) + plot_annotation(
#   title = 'Lesion number on Kasalath')
# Klcplot
# ggsave(filename = "Klc.png", plot = Klcplot, device = "png", height = 20, width = 30,
#        units = "cm", dpi = 500)

library(patchwork)
Kplot <- ((KF.ls.plot + KFlc) / (KG.ls.plot + Kglc)) + plot_annotation(
  title = 'Lesion size and count on Kasalath', tag_levels = 'A') 
Kplot
ggsave(filename = "K.png", plot = Kplot, device = "png", height = 20, width = 40,
       units = "cm", dpi = 500)


#ggplot punch
KFp <- dplyr::filter(K.p, isolate == "FR13")
KGp <- dplyr::filter(K.p, isolate == "Guy11")
library(ggpubr)
library(ggplot2)
kf.p.p <- ggplot(KFp,aes(x=factor(mutant), y=lesion.surface, color = effector)) + geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface), show.legend = TRUE) + # jitter size by lesion surface
  scale_size_continuous(range = c(0.01, 2))  + theme_minimal() + scale_x_discrete(limits=rev)
KF.p.plot <- kf.p.p +
  labs(title = "FR13 punch inoculation", x = "Mutant", y = "Lesion surface (in pixels)") +
  # theme(axis.title.y = element_blank(), axis.title.x = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
kg.p.p <- ggplot(KGp,aes(x=factor(mutant), y=lesion.surface, color = effector)) + geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface), show.legend = TRUE) + # jitter size by lesion surface
  scale_size_continuous(range = c(0.01, 2))  + theme_minimal() + scale_x_discrete(limits=rev)
KG.p.plot <- kg.p.p +
  labs(title = "Guy11 punch inoculation", x = "Mutant", y = "Lesion surface (in pixels)") +   
  #theme(axis.title.y = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0))   
library(patchwork)
Kpplot <- (KF.p.plot / KG.p.plot) + plot_annotation(
  title = 'Punch inoculation on Kasalath')
Kpplot
ggsave(filename = "Kp.png", plot = Kpplot, device = "png", height = 20, width = 30,
       units = "cm", dpi = 500)