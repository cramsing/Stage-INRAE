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


K6F <- dplyr::filter(K6, isolate == "FR13")
K6G <- dplyr::filter(K6, isolate == "Guy11")
#kruskal groups
library(agricolae)
library(dplyr)
k6.f.ls.groups <- kruskal(K6F$lesion.surface,
                           K6F$mutant, group = TRUE) %>%
  .$groups %>%
  as_tibble(rownames = "mutant") %>%
  rename("k6.f.ls.groups" = "groups")
k6.f.ls.groups$`K6F$lesion.surface` <- NULL
k6.g.ls.groups <- kruskal(K6G$lesion.surface,
                          K6G$mutant, group = TRUE) %>%
  .$groups %>%
  as_tibble(rownames = "mutant") %>%
  rename("k6.g.ls.groups" = "groups")
k6.g.ls.groups$`K6G$lesion.surface` <- NULL


#standard error
Pch.D.Error <- PchConc %>%
  group_by(diam.class) %>%
  summarise(Pch.D.Error = standard_error(Conc_pg.ul))
Pch.D.Error <- dplyr::rename(Pch.D.Error, "diameter class" = "diam.class")
Pch.D.Error$`diameter class`<- as.character(Pch.D.Error$`diameter class`)
pch.diam.groups <- left_join(pch.diam.groups, Pch.D.Error, by = "diameter class")
print(pch.diam.groups)

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


#kruskal groups
#ggplot lesion number
K6Flc <- dplyr::filter(K6lc, isolate == "FR13")
K6Glc <- dplyr::filter(K6lc, isolate == "Guy11")
library(agricolae)
library(dplyr)
k6.f.lc.groups <- kruskal(K6Flc$lesion.count,
                          K6Flc$mutant, group = TRUE) %>%
  .$groups %>%
  as_tibble(rownames = "mutant") %>%
  rename("k6.f.lc.groups" = "groups")
k6.f.lc.groups$`K6Flc$lesion.count`<- NULL
k6.g.lc.groups <- kruskal(K6Glc$lesion.count,
                          K6Glc$mutant, group = TRUE) %>%
  .$groups %>%
  as_tibble(rownames = "mutant") %>%
  rename("k6.g.lc.groups" = "groups")
k6.g.lc.groups$`K6Glc$lesion.count` <- NULL
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

K6Fp <- dplyr::filter(K6.p, isolate == "FR13")
K6Gp <- dplyr::filter(K6.p, isolate == "Guy11")
library(agricolae)
library(dplyr)
k6.f.p.groups <- kruskal(K6Fp$lesion.surface,
                         K6Fp$mutant, group = TRUE) %>%
  .$groups %>%
  as_tibble(rownames = "mutant") %>%
  rename("k6.f.p.groups" = "groups")
k6.f.p.groups$`K6Fp$lesion.surface`<- NULL
k6.g.p.groups <- kruskal(K6Gp$lesion.surface,
                         K6Gp$mutant, group = TRUE) %>%
  .$groups %>%
  as_tibble(rownames = "mutant") %>%
  rename("k6.g.p.groups" = "groups")
k6.g.p.groups$`K6Gp$lesion.surface` <- NULL


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

## Compiled kruskal groups ----
library(dplyr)
k6.f.groups <- dplyr::left_join(k6.f.ls.groups, k6.f.lc.groups, by = "mutant") %>%
  left_join(., k6.f.p.groups, by = "mutant" )
k6.f.groups <- k6.f.groups[order(k6.f.groups$mutant),]
write.csv2(k6.f.groups, file = "k6.f.groups.csv")
k6.g.groups <- dplyr::left_join(k6.g.ls.groups, k6.g.lc.groups, by = "mutant") %>%
  left_join(., k6.g.p.groups, by = "mutant" )
k6.g.groups <- k6.g.groups[order(k6.g.groups$mutant),] 
write.csv2(k6.g.groups, file = "k6.g.groups.csv")


# Ggplots ----

#ggplot lesion size
K6F <- dplyr::filter(K6, isolate == "FR13")
K6G <- dplyr::filter(K6, isolate == "Guy11")
library(ggpubr)
library(ggplot2)

K6f.ls.p <- ggplot(K6F, aes(x=factor(mutant), y=lesion.surface, color = effector)) + 
  geom_boxplot(color = "darkslategray") + 
  geom_jitter(aes(size = lesion.surface, shape = rep), show.legend = TRUE) + # jitter size by lesion surface
  scale_size_continuous(range = c(0.01, 1.5))  + theme_minimal() + scale_x_discrete(limits=rev) + 
  scale_y_continuous(trans = "log10") 
K6F.ls.plot <- K6f.ls.p +
  labs(title = "FR13 lesion size", x = "Mutant", y = "Lesion size (in log10 pixels)") +
  # theme(axis.title.y = element_blank(), axis.title.x = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
# stat_compare_means(label = "p.signif", method = "wilcox.test",
#                    ref.group = ".all.", hide.ns = TRUE) 
K6g.ls.p <- ggplot(K6G,aes(x=factor(mutant), y=lesion.surface, color = effector)) + 
  geom_boxplot(color = "darkslategray") + 
  geom_jitter(aes(size = lesion.surface, shape = rep), show.legend = TRUE) + # jitter size by lesion surface
  scale_size_continuous(range = c(0.01, 1.5))  + theme_minimal() + scale_x_discrete(limits=rev) + 
  scale_y_continuous(trans = "log10") 
K6G.ls.plot <- K6g.ls.p +
  labs(title = "Guy11 lesion size", x = "Mutant", y = "Lesion size (in log10 pixels)") +   
  #theme(axis.title.y = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0))   
# library(patchwork)
# K6lsplot <- (K6F.ls.plot / K6G.ls.plot) + plot_annotation(
#   title = 'Lesion surface on K60')
# K6lsplot
# ggsave(filename = "K6ls.png", plot = K6lsplot, device = "png", height = 20, width = 30,
#        units = "cm", dpi = 500)

#ggplot lesion number
K6Flc <- dplyr::filter(K6lc, isolate == "FR13")
K6Glc <- dplyr::filter(K6lc, isolate == "Guy11")
library(ggpubr)
k6f.lc0 <- ggplot(K6Flc, aes(x=factor(mutant), y=lesion.count, color = effector)) + 
  geom_boxplot(color = "darkslategray") + geom_jitter(size=1.0, aes(shape = rep.x)) + 
  theme_minimal() + scale_x_discrete(limits=rev) 
K6Flc <- k6f.lc0 +
  labs(title = "FR13 lesion number", x = "Mutant", y = "Lesion number") + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
# + stat_compare_means(label = "p.signif", method = "wilcox.test",
#                      ref.group = ".all.", hide.ns = TRUE, label.y = c(40, 60, 200, 30,30)) 
k6g.lc0 <- ggplot(K6Glc, aes(x=factor(mutant), y=lesion.count, color = effector)) + 
  geom_boxplot(color = "darkslategray") + geom_jitter(size=1.0, aes(shape = rep.x)) + 
  theme_minimal() + scale_x_discrete(limits=rev) 
K6glc <- k6g.lc0 +
  labs(title = "Guy11 lesion number", x = "Mutant", y = "Lesion number") + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
# + stat_compare_means(label = "p.signif", method = "wilcox.test",
#                      ref.group = ".all.", hide.ns = TRUE, label.y = c(125, 170, 170, 30,30)) 
# library(patchwork)
# K6lcplot <- (K6Flc / K6glc) + plot_annotation(
#   title = 'Lesion number on K60')
# K6lcplot
# ggsave(filename = "K6lc.png", plot = K6lcplot, device = "png", height = 20, width = 30,
#        units = "cm", dpi = 500)



library(patchwork)
K6plot <- ((K6F.ls.plot + K6Flc) / (K6G.ls.plot + K6glc)) + plot_annotation(
  title = 'Lesion size and count on K60', tag_levels = 'A') 
K6plot
ggsave(filename = "K6.png", plot = K6plot, device = "png", height = 20, width = 40,
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