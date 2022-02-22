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
#p-value <2.2e-16, significant !

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
#p-value = 4.604e-13 ***, significant
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
#1.001e-12, not normal
kruskal.test(lesion.count ~ mutant, data = Slc)
#p 1.481e-08, significant

#Kruskal test with only rep 1
kruskal.test(lesion.count ~ mutant, data = Slc.1) #p-value 8.194e-08
#Kruskal test with only rep 2
kruskal.test(lesion.count ~ mutant, data = Slc.2) #p-value 7.443e-13
#very different values

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


## Compiled kruskal groups ----
SF <- dplyr::filter(S, isolate == "FR13")
SG <- dplyr::filter(S, isolate == "Guy11")
#lesion surface kruskal groups
library(agricolae)
library(dplyr)
S.f.ls.groups <- kruskal(SF$lesion.surface,
                          SF$mutant, group = TRUE) %>%
  .$groups %>%
  as_tibble(rownames = "mutant") %>%
  rename("S.f.ls.groups" = "groups")
S.f.ls.groups$`SF$lesion.surface` <- NULL
S.g.ls.groups <- kruskal(SG$lesion.surface,
                          SG$mutant, group = TRUE) %>%
  .$groups %>%
  as_tibble(rownames = "mutant") %>%
  rename("S.g.ls.groups" = "groups")
S.g.ls.groups$`SG$lesion.surface` <- NULL


#lesion count kruskal groups
SFlc <- dplyr::filter(Slc, isolate == "FR13")
SGlc <- dplyr::filter(Slc, isolate == "Guy11")
library(agricolae)
library(dplyr)
S.f.lc.groups <- kruskal(SFlc$lesion.count,
                          SFlc$mutant, group = TRUE) %>%
  .$groups %>%
  as_tibble(rownames = "mutant") %>%
  rename("S.f.lc.groups" = "groups")
S.f.lc.groups$`SFlc$lesion.count`<- NULL
S.g.lc.groups <- kruskal(SGlc$lesion.count,
                          SGlc$mutant, group = TRUE) %>%
  .$groups %>%
  as_tibble(rownames = "mutant") %>%
  rename("S.g.lc.groups" = "groups")
S.g.lc.groups$`SGlc$lesion.count` <- NULL

#punch inoculation kruskal groups

SFp <- dplyr::filter(S.p, isolate == "FR13")
SGp <- dplyr::filter(S.p, isolate == "Guy11")
library(agricolae)
library(dplyr)
S.f.p.groups <- kruskal(SFp$lesion.surface,
                         SFp$mutant, group = TRUE) %>%
  .$groups %>%
  as_tibble(rownames = "mutant") %>%
  rename("S.f.p.groups" = "groups")
S.f.p.groups$`SFp$lesion.surface`<- NULL
S.g.p.groups <- kruskal(SGp$lesion.surface,
                         SGp$mutant, group = TRUE) %>%
  .$groups %>%
  as_tibble(rownames = "mutant") %>%
  rename("S.g.p.groups" = "groups")
S.g.p.groups$`SGp$lesion.surface` <- NULL

#compilation
library(dplyr)
S.f.groups <- dplyr::left_join(S.f.ls.groups, S.f.lc.groups, by = "mutant") %>%
  left_join(., S.f.p.groups, by = "mutant" )
S.f.groups <- S.f.groups[order(S.f.groups$mutant),]
write.csv2(S.f.groups, file = "S.f.groups.csv")
S.g.groups <- dplyr::left_join(S.g.ls.groups, S.g.lc.groups, by = "mutant") %>%
  left_join(., S.g.p.groups, by = "mutant" )
S.g.groups <- S.g.groups[order(S.g.groups$mutant),] 
write.csv2(S.g.groups, file = "S.g.groups.csv")


# Ggplots ----

#ggplot lesion size
#ggplot lesion size
SF <- dplyr::filter(S, isolate == "FR13")
SG <- dplyr::filter(S, isolate == "Guy11")
library(ggpubr)
Sf.ls.p <- ggplot(SF,aes(x=factor(mutant, level = f.level_order), y=lesion.surface, fill = effector)) +   
  geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface, shape = rep, color = rep), show.legend = TRUE) + 
  scale_size_continuous(range = c(0.1, 2))+ theme_minimal() + scale_x_discrete(limits=rev) + 
  scale_y_continuous(trans = "log10", breaks = c(10, 100, 500, 1000, 5000, 10000))  
SF.ls.plot <- Sf.ls.p + scale_colour_grey() + scale_fill_manual(values = blue.palette) +
  labs(title = "FR13 lesion surface", x = "Isolate", y = "Lesion surface (in log10 pixels)") +
  # theme(axis.title.y = element_blank(), axis.title.x = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) +  coord_flip()
# stat_compare_means(label = "p.signif", method = "wilcox.test",
#                    ref.group = ".all.", hide.ns = TRUE) 
Sg.ls.p <- ggplot(SG,aes(x=factor(mutant, level = g.level_order), y=lesion.surface, fill = effector)) + 
  geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface, shape = rep, color = rep), show.legend = TRUE) + 
  scale_size_continuous(range = c(0.1, 2))  + theme_minimal() + scale_x_discrete(limits=rev) +
  scale_y_continuous(trans = "log10", breaks = c(10, 100, 500, 1000, 5000, 10000)) 
SG.ls.plot <- Sg.ls.p + scale_colour_grey() + scale_fill_manual(values = blue.palette) +
  labs(title = "Guy11 lesion surface", x = "Isolate", y = "Lesion surface (in log10 pixels)") +   
  #theme(axis.title.y = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) +  coord_flip() 
# library(patchwork)
# Slsplot <- (SF.ls.plot/SG.ls.plot) + plot_annotation(title = 'Lesion surface on Shin 2')
# Slsplot
ggsave(filename = "Sfls.png", plot = SF.ls.plot, device = "png", height = 20, width = 30,
       units = "cm", dpi = 500)
ggsave(filename = "Sgls.png", plot = SG.ls.plot, device = "png", height = 20, width = 30,
       units = "cm", dpi = 500)

#ggplot lesion number
SFlc <- dplyr::filter(Slc, isolate == "FR13")
SGlc <- dplyr::filter(Slc, isolate == "Guy11")
library(ggpubr)
Sf.lc0 <- ggplot(SFlc, aes(x=factor(mutant, level = f.level_order), y=lesion.count, fill = effector)) + 
  geom_boxplot() + geom_jitter(size=2.0, aes(shape = rep, color = rep)) + 
  theme_minimal() + scale_x_discrete(limits=rev) + coord_flip()
SFlc <- Sf.lc0 +  scale_colour_grey() + scale_fill_manual(values = blue.palette) +
  labs(title = "FR13 lesion number", x = "Isolate", y = "Lesion number") + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
# + stat_compare_means(label = "p.signif", method = "wilcox.test",
#                      ref.group = ".all.", hide.ns = TRUE, label.y = c(40, 60, 200, 30,30)) 
Sg.lc0 <- ggplot(SGlc, aes(x=factor(mutant, level = g.level_order), y=lesion.count, fill = effector)) + 
  geom_boxplot() + geom_jitter(size=2.0, aes(shape = rep, color = rep)) + 
  theme_minimal() + scale_x_discrete(limits=rev) + coord_flip()
Sglc <- Sg.lc0 +  scale_colour_grey() + scale_fill_manual(values = blue.palette) +
  labs(title = "Guy11 lesion number", x = "Isolate", y = "Lesion number") + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
# + stat_compare_means(label = "p.signif", method = "wilcox.test",
#                      ref.group = ".all.", hide.ns = TRUE, label.y = c(125, 170, 170, 30,30)) 
# library(patchwork)
# Slcplot <- (SFlc / Sglc) + plot_annotation(
#   title = 'Lesion number on Shin 2')
# Slcplot
ggsave(filename = "Sflc.png", plot = SFlc, device = "png", height = 20, width = 30,
       units = "cm", dpi = 500)
ggsave(filename = "Sglc.png", plot = Sglc, device = "png", height = 20, width = 30,
       units = "cm", dpi = 500)
# library(patchwork)
# Splot <- ((SF.ls.plot + SFlc) / (SG.ls.plot + Sglc)) + plot_annotation(
#   title = 'Lesion size and count on Shin 2', tag_levels = 'A') 
library(patchwork)
Sgplot <- ((Sglc + plot_spacer())/SG.ls.plot ) + plot_annotation(title = 'Guy11 lesion size and count on Shin 2', tag_levels = 'A') 
ggsave(filename = "Sg.png", plot = Sgplot, device = "png", height = 25, width = 30,
       units = "cm", dpi = 500)
#ggplot punch
SFp <- dplyr::filter(S.p, isolate == "FR13")
SGp <- dplyr::filter(S.p, isolate == "Guy11")
library(ggpubr)
library(ggplot2)
Sf.p.p <- ggplot(SFp,aes(x=factor(mutant, level = f.level_order), y=lesion.surface, fill = effector)) + geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface), show.legend = TRUE) + # jitter size by lesion surface
  scale_size_continuous(range = c(0.01, 2))  + theme_minimal() + scale_x_discrete(limits=rev)
SF.p.plot <- Sf.p.p + scale_fill_manual(values = blue.palette) +
  labs(title = "FR13 punch inoculation", x = "Isolate", y = "Lesion surface (in pixels)") +
  # theme(axis.title.y = element_blank(), axis.title.x = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
Sg.p.p <- ggplot(SGp, aes(x=factor(mutant, level = g.level_order), y=lesion.surface, fill = effector)) + geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface), show.legend = TRUE) + # jitter size by lesion surface
  scale_size_continuous(range = c(0.01, 2))  + theme_minimal() + scale_x_discrete(limits=rev)
SG.p.plot <- Sg.p.p + scale_fill_manual(values = blue.palette) +
  labs(title = "Guy11 punch inoculation", x = "Isolate", y = "Lesion surface (in pixels)") +   
  #theme(axis.title.y = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0))   
library(patchwork)
Spplot <- (SF.p.plot / SG.p.plot) + plot_annotation(
  title = 'Punch inoculation on Shin 2')
Spplot
ggsave(filename = "Sp.png", plot = Spplot, device = "png", height = 20, width = 30,
       units = "cm", dpi = 500)
