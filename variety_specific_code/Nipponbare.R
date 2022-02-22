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
#p-value <2.2e-16, significant !

#Kruskal test with only rep 1
kruskal.test(lesion.surface ~ mutant, data = N.1) #p-value <2.2e-16
#Kruskal test with only rep 2
kruskal.test(lesion.surface ~ mutant, data = N.2) #p-value <2.2e-16
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
#p-value = <2.2e-16 ***, significant
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
#2.169e-11, not normal
kruskal.test(lesion.count ~ mutant, data = Nlc)
#p 1.767e-13, significant

#Kruskal test with only rep 1
kruskal.test(lesion.count ~ mutant, data = Nlc.1) #p-value 3.072e-10
#Kruskal test with only rep 2
kruskal.test(lesion.count ~ mutant, data = Nlc.2) #p-value 1.119e-11
#very different values

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

#check for each isolate
NFp <- dplyr::filter(N.p, isolate == "FR13")
NGp <- dplyr::filter(N.p, isolate == "Guy11")

N3 <- aov(lesion.surface ~ mutant, data = NFp) 
anova(N3)
#p-value =   <2e-16 ***, significant
#normality of residuals check
#distrubution/normality
par(mfrow=c(1,2))
#Histogramme
hist(N2$residuals)
# Q-Q Plot
qqnorm(N3$residuals)
qqline(N3$residuals)
#kind of normal looking
shapiro.test(N3$residuals)

N4 <- aov(lesion.surface ~ mutant, data = N.p) 
anova(N4)
#p-value =   <2e-16 ***, significant
#normality of residuals check
#distrubution/normality
par(mfrow=c(1,2))
#Histogramme
hist(N4$residuals)
# Q-Q Plot
qqnorm(N4$residuals)
qqline(N4$residuals)
#kind of normal looking
shapiro.test(N4$residuals)

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

## Compiled dunn groups ----
NF <- dplyr::filter(N, isolate == "FR13")
NG <- dplyr::filter(N, isolate == "Guy11")
#lesion surface kruskal groups
library(rcompanion)
nf.d = dunnTest(lesion.surface ~ mutant, data = NF, method ="bh")
nf.d=nf.d$res
nf.d <- cldList(comparison = nf.d$Comparison,
        p.value    = nf.d$P.adj,
        threshold  = 0.05)
nf.d$MonoLetter <- NULL
nf.d <- dplyr::rename(nf.d,"mutant" = "Group")
ng.d = dunnTest(lesion.surface ~ mutant, data = NG, method ="bh")
ng.d=nf.d$res
ng.d <- cldList(comparison = ng.d$Comparison,
                p.value    = ng.d$P.adj,
                threshold  = 0.05)
ng.d$MonoLetter <- NULL
ng.d <- dplyr::rename(ng.d,"mutant" = "Group")

#lesion count kruskal groups
NFlc <- dplyr::filter(Nlc, isolate == "FR13")
NGlc <- dplyr::filter(Nlc, isolate == "Guy11")
library(rcompanion)
nflc.d = dunnTest(lesion.surface ~ mutant, data = NFlc, method ="bh")
nflc.d=nflc.d$res
nflc.d <- cldList(comparison = nflc.d$Comparison,
                p.value    = nflc.d$P.adj,
                threshold  = 0.05)
nflc.d$MonoLetter <- NULL
nflc.d <- dplyr::rename(nflc.d,"mutant" = "Group")
nglc.d = dunnTest(lesion.surface ~ mutant, data = NGlc, method ="bh")
nglc.d=nglc.d$res
nglc.d <- cldList(comparison = nglc.d$Comparison,
                p.value    = nglc.d$P.adj,
                threshold  = 0.05)
nglc.d$MonoLetter <- NULL
nglc.d <- dplyr::rename(nglc.d,"mutant" = "Group")

#punch inoculation kruskal groups

NFp <- dplyr::filter(N.p, isolate == "FR13")
NGp <- dplyr::filter(N.p, isolate == "Guy11")
library(agricolae)
library(dplyr)

k3 <- aov(lesion.surface ~ mutant, data = NFp) 
anova(k3)
N.f.p.groups <- HSD.test(k3, "mutant", group = TRUE) %>%
  .$groups %>%
  as_tibble(rownames = "mutant") %>%
  rename("N.f.p.groups" = "groups")
N.f.p.groups$lesion.surface<- NULL
N.g.p.groups <-HSD.test(k4, "mutant", group = TRUE) %>%
  .$groups %>%
  as_tibble(rownames = "mutant") %>%
  rename("N.g.p.groups" = "groups")
N.g.p.groups$lesion.surface<- NULL

#compilation
library(dplyr)
N.f.groups <- dplyr::left_join(nf.d, nflc.d, by = "mutant") %>%
  left_join(., N.f.p.groups, by = "mutant" )
write.csv2(N.f.groups, file = "N.f.groups.csv")
N.g.groups <- dplyr::left_join(ng.d, nglc.d,by = "mutant") %>%
  left_join(., N.g.p.groups, by = "mutant" )
write.csv2(N.g.groups, file = "N.g.groups.csv")

# Ggplots ----

#ggplot lesion size
NF <- dplyr::filter(N, isolate == "FR13")
NG <- dplyr::filter(N, isolate == "Guy11")
library(ggpubr)
Nf.ls.p <- ggplot(NF,aes(x=factor(mutant, level = f.level_order), y=lesion.surface, fill = effector)) +   
  geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface, shape = rep, color = rep), show.legend = TRUE) + 
  scale_size_continuous(range = c(0.1, 2))+ theme_minimal() + scale_x_discrete(limits=rev) + 
  scale_y_continuous(trans = "log10", breaks = c(10, 100, 500, 1000, 5000, 10000))  
NF.ls.plot <- Nf.ls.p + scale_colour_grey() + scale_fill_manual(values = blue.palette) +
  labs(title = "FR13 lesion surface", x = "Isolate", y = "Lesion surface (in log10 pixels)") +
  # theme(axis.title.y = element_blank(), axis.title.x = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) +  coord_flip()
# stat_compare_means(label = "p.signif", method = "wilcox.test",
#                    ref.group = ".all.", hide.ns = TRUE) 
Ng.ls.p <- ggplot(NG,aes(x=factor(mutant, level = g.level_order), y=lesion.surface, fill = effector)) + 
  geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface, shape = rep, color = rep), show.legend = TRUE) + 
  scale_size_continuous(range = c(0.1, 2))  + theme_minimal() + scale_x_discrete(limits=rev) +
  scale_y_continuous(trans = "log10", breaks = c(10, 100, 500, 1000, 5000, 10000)) 
NG.ls.plot <- Ng.ls.p + scale_colour_grey() + scale_fill_manual(values = blue.palette) +
  labs(title = "Guy11 lesion surface", x = "Isolate", y = "Lesion surface (in log10 pixels)") +   
  #theme(axis.title.y = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) +  coord_flip() 
# library(patchwork)
# Nlsplot <- (NF.ls.plot/NG.ls.plot) + plot_annotation(title = 'Lesion surface on Nipponbare')
# Nlsplot
ggsave(filename = "Nfls.png", plot = NF.ls.plot, device = "png", height = 20, width = 30,
       units = "cm", dpi = 500)
ggsave(filename = "Ngls.png", plot = NG.ls.plot, device = "png", height = 20, width = 30,
       units = "cm", dpi = 500)

#ggplot lesion number
NFlc <- dplyr::filter(Nlc, isolate == "FR13")
NGlc <- dplyr::filter(Nlc, isolate == "Guy11")
library(ggpubr)
Nf.lc0 <- ggplot(NFlc, aes(x=factor(mutant, level = f.level_order), y=lesion.count, fill = effector)) + 
  geom_boxplot() + geom_jitter(size=2.0, aes(shape = rep, color = rep)) + 
  theme_minimal() + scale_x_discrete(limits=rev) + coord_flip()
NFlc <- Nf.lc0 +  scale_colour_grey() + scale_fill_manual(values = blue.palette) +
  labs(title = "FR13 lesion number", x = "Isolate", y = "Lesion number") + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
# + stat_compare_means(label = "p.signif", method = "wilcox.test",
#                      ref.group = ".all.", hide.ns = TRUE, label.y = c(40, 60, 200, 30,30)) 
Ng.lc0 <- ggplot(NGlc, aes(x=factor(mutant, level = g.level_order), y=lesion.count, fill = effector)) + 
  geom_boxplot() + geom_jitter(size=2.0, aes(shape = rep, color = rep)) + 
  theme_minimal() + scale_x_discrete(limits=rev) + coord_flip()
Nglc <- Ng.lc0 +  scale_colour_grey() + scale_fill_manual(values = blue.palette) +
  labs(title = "Guy11 lesion number", x = "Isolate", y = "Lesion number") + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
# + stat_compare_means(label = "p.signif", method = "wilcox.test",
#                      ref.group = ".all.", hide.ns = TRUE, label.y = c(125, 170, 170, 30,30)) 
# library(patchwork)
# Nlcplot <- (NFlc / Nglc) + plot_annotation(
#   title = 'Lesion number on Nipponbare')
# Nlcplot
ggsave(filename = "Nflc.png", plot = NFlc, device = "png", height = 20, width = 30,
       units = "cm", dpi = 500)
ggsave(filename = "Nglc.png", plot = Nglc, device = "png", height = 20, width = 30,
       units = "cm", dpi = 500)
# library(patchwork)
# Nplot <- ((NF.ls.plot + NFlc) / (NG.ls.plot + Nglc)) + plot_annotation(
#   title = 'Lesion size and count on Nipponbare', tag_levels = 'A') 
library(patchwork)
Ngplot <- ((Nglc + plot_spacer())/NG.ls.plot ) + plot_annotation(title = 'Guy11 lesion size and count on Nipponbare', tag_levels = 'A') 
ggsave(filename = "Ng.png", plot = Ngplot, device = "png", height = 25, width = 30,
       units = "cm", dpi = 500)
#ggplot punch
NFp <- dplyr::filter(N.p, isolate == "FR13")
NGp <- dplyr::filter(N.p, isolate == "Guy11")
library(ggpubr)
library(ggplot2)
Nf.p.p <- ggplot(NFp,aes(x=factor(mutant, level = f.level_order), y=lesion.surface, fill = effector)) + geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface), show.legend = TRUE) + # jitter size by lesion surface
  scale_size_continuous(range = c(0.01, 2))  + theme_minimal() + scale_x_discrete(limits=rev)
NF.p.plot <- Nf.p.p + scale_fill_manual(values = blue.palette) +
  labs(title = "FR13 punch inoculation", x = "Isolate", y = "Lesion surface (in pixels)") +
  # theme(axis.title.y = element_blank(), axis.title.x = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
Ng.p.p <- ggplot(NGp, aes(x=factor(mutant, level = g.level_order), y=lesion.surface, fill = effector)) + geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface), show.legend = TRUE) + # jitter size by lesion surface
  scale_size_continuous(range = c(0.01, 2))  + theme_minimal() + scale_x_discrete(limits=rev)
NG.p.plot <- Ng.p.p + scale_fill_manual(values = blue.palette) +
  labs(title = "Guy11 punch inoculation", x = "Isolate", y = "Lesion surface (in pixels)") +   
  #theme(axis.title.y = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0))   
library(patchwork)
Npplot <- (NF.p.plot / NG.p.plot) + plot_annotation(
  title = 'Punch inoculation on Nipponbare')
Npplot
ggsave(filename = "Np.png", plot = Npplot, device = "png", height = 20, width = 30,
       units = "cm", dpi = 500)
