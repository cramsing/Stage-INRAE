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


## Compiled kruskal groups ----
BF <- dplyr::filter(B, isolate == "FR13")
BG <- dplyr::filter(B, isolate == "Guy11")
#lesion surface kruskal groups
library(agricolae)
library(dplyr)
B.f.ls.groups <- kruskal(BF$lesion.surface,
                          BF$mutant, group = TRUE) %>%
  .$groups %>%
  as_tibble(rownames = "mutant") %>%
  rename("B.f.ls.groups" = "groups")
B.f.ls.groups$`BF$lesion.surface` <- NULL
B.g.ls.groups <- kruskal(BG$lesion.surface,
                          BG$mutant, group = TRUE) %>%
  .$groups %>%
  as_tibble(rownames = "mutant") %>%
  rename("B.g.ls.groups" = "groups")
B.g.ls.groups$`BG$lesion.surface` <- NULL


#lesion count kruskal groups
BFlc <- dplyr::filter(Blc, isolate == "FR13")
BGlc <- dplyr::filter(Blc, isolate == "Guy11")
library(agricolae)
library(dplyr)
B.f.lc.groups <- kruskal(BFlc$lesion.count,
                          BFlc$mutant, group = TRUE) %>%
  .$groups %>%
  as_tibble(rownames = "mutant") %>%
  rename("B.f.lc.groups" = "groups")
B.f.lc.groups$`BFlc$lesion.count`<- NULL
B.g.lc.groups <- kruskal(BGlc$lesion.count,
                          BGlc$mutant, group = TRUE) %>%
  .$groups %>%
  as_tibble(rownames = "mutant") %>%
  rename("B.g.lc.groups" = "groups")
B.g.lc.groups$`BGlc$lesion.count` <- NULL

#punch inoculation kruskal groups

BFp <- dplyr::filter(B.p, isolate == "FR13")
BGp <- dplyr::filter(B.p, isolate == "Guy11")
library(agricolae)
library(dplyr)
B.f.p.groups <- kruskal(BFp$lesion.surface,
                         BFp$mutant, group = TRUE) %>%
  .$groups %>%
  as_tibble(rownames = "mutant") %>%
  rename("B.f.p.groups" = "groups")
B.f.p.groups$`BFp$lesion.surface`<- NULL
B.g.p.groups <- kruskal(BGp$lesion.surface,
                         BGp$mutant, group = TRUE) %>%
  .$groups %>%
  as_tibble(rownames = "mutant") %>%
  rename("B.g.p.groups" = "groups")
B.g.p.groups$`BGp$lesion.surface` <- NULL

#compilation
library(dplyr)
B.f.groups <- dplyr::left_join(B.f.ls.groups, B.f.lc.groups, by = "mutant") %>%
  left_join(., B.f.p.groups, by = "mutant" )
B.f.groups <- B.f.groups[order(B.f.groups$mutant),]
write.csv2(B.f.groups, file = "B.f.groups.csv")
B.g.groups <- dplyr::left_join(B.g.ls.groups, B.g.lc.groups, by = "mutant") %>%
  left_join(., B.g.p.groups, by = "mutant" )
B.g.groups <- B.g.groups[order(B.g.groups$mutant),] 
write.csv2(B.g.groups, file = "B.g.groups.csv")


# Ggplots ----
blue.palette <- c("#f0f9e8","#bae4bc", "#7bccc4","#43a2ca","#0868ac")
f.level_order <- c("FR13 WT","FR13 RFP" ,"mFA1","mFA2","mFA3","FA ectopic", "mFB1",
                   "mFB2", "mFB3","FB ectopic", "mFK1","mFK2", "mFK3","FK ectopic")
g.level_order <- c("Guy11 WT", "Guy11 RFP_1", "Guy11 RFP_2", "mGB1", "mGB2","mGB3",
                   "mGK1","mGK2","mGK3", "GK ectopic")

#ggplot lesion size
BF <- dplyr::filter(B, isolate == "FR13")
BG <- dplyr::filter(B, isolate == "Guy11")
library(ggpubr)
bf.ls.p <- ggplot(BF,aes(x=factor(mutant, level = f.level_order), y=lesion.surface, fill = effector)) +   
  geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface, shape = rep, color = rep), show.legend = TRUE) + 
  scale_size_continuous(range = c(0.1, 2))+ theme_minimal() + scale_x_discrete(limits=rev) + 
  scale_y_continuous(trans = "log10", breaks = c(10, 100, 500, 1000, 5000, 10000))  
BF.ls.plot <- bf.ls.p + scale_colour_grey() + scale_fill_manual(values = blue.palette) +
  labs(x = "Isolate", y = "Lesion surface (in log10 pixels)") +
  # theme(axis.title.y = element_blank(), axis.title.x = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) +  coord_flip()
# stat_compare_means(label = "p.signif", method = "wilcox.test",
#                    ref.group = ".all.", hide.ns = TRUE) 
bg.ls.p <- ggplot(BG,aes(x=factor(mutant, level = g.level_order), y=lesion.surface, fill = effector)) + 
  geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface, shape = rep, color = rep), show.legend = TRUE) + 
  scale_size_continuous(range = c(0.1, 2))  + theme_minimal() + scale_x_discrete(limits=rev) +
  scale_y_continuous(trans = "log10", breaks = c(10, 100, 500, 1000, 5000, 10000)) 
BG.ls.plot <- bg.ls.p + scale_colour_grey() + scale_fill_manual(values = blue.palette) +
  labs(title = "Guy11 lesion surface", x = "Isolate", y = "Lesion surface (in log10 pixels)") +   
  #theme(axis.title.y = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) +  coord_flip() 
# library(patchwork)
# Blsplot <- (BF.ls.plot/BG.ls.plot) + plot_annotation(title = 'Lesion surface on Bl1')
# Blsplot
ggsave(filename = "Bfls.png", plot = BF.ls.plot, device = "png", height = 20, width = 30,
        units = "cm", dpi = 500)
ggsave(filename = "Bgls.png", plot = BG.ls.plot, device = "png", height = 20, width = 30,
       units = "cm", dpi = 500)

#ggplot lesion number
BFlc <- dplyr::filter(Blc, isolate == "FR13")
BGlc <- dplyr::filter(Blc, isolate == "Guy11")
library(ggpubr)
bf.lc0 <- ggplot(BFlc, aes(x=factor(mutant, level = f.level_order), y=lesion.count, fill = effector)) + 
  geom_boxplot() + geom_jitter(size=2.0, aes(shape = rep, color = rep)) + 
  theme_minimal() + scale_x_discrete(limits=rev) + coord_flip()
BFlc <- bf.lc0 +  scale_colour_grey() + scale_fill_manual(values = blue.palette) +
  labs(x = "Isolate", y = "Lesion number") + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) + theme(legend.position = "none")
# + stat_compare_means(label = "p.signif", method = "wilcox.test",
#                      ref.group = ".all.", hide.ns = TRUE, label.y = c(40, 60, 200, 30,30)) 
bg.lc0 <- ggplot(BGlc, aes(x=factor(mutant, level = g.level_order), y=lesion.count, fill = effector)) + 
  geom_boxplot() + geom_jitter(size=2.0, aes(shape = rep, color = rep)) + 
  theme_minimal() + scale_x_discrete(limits=rev) + coord_flip()
Bglc <- bg.lc0 +  scale_colour_grey() + scale_fill_manual(values = blue.palette) +
  labs(title = "Guy11 lesion surface", x = "Isolate", y = "Lesion number") + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) + theme(legend.position = "none")
# + stat_compare_means(label = "p.signif", method = "wilcox.test",
#                      ref.group = ".all.", hide.ns = TRUE, label.y = c(125, 170, 170, 30,30)) 
# library(patchwork)
# Blcplot <- (BFlc / Bglc) + plot_annotation(
#   title = 'Lesion number on Bl1')
# Blcplot
ggsave(filename = "Bflc.png", plot = BFlc, device = "png", height = 20, width = 30,
        units = "cm", dpi = 500)
ggsave(filename = "Bglc.png", plot = Bglc, device = "png", height = 20, width = 30,
       units = "cm", dpi = 500)
# library(patchwork)
# Bplot <- ((BF.ls.plot + BFlc) / (BG.ls.plot + Bglc)) + plot_annotation(
#   title = 'Lesion size and count on Bl1', tag_levels = 'A') 
library(patchwork)
Bgplot <- ((Bglc + plot_spacer())/BG.ls.plot ) + #plot_layout(guides = 'collect') + 
  plot_annotation( tag_levels = 'A') 
ggsave(filename = "Bg.png", plot = Bgplot, device = "png", height = 25, width = 30,
       units = "cm", dpi = 500)
Bfplot <- ((BFlc + plot_spacer())/ BF.ls.plot ) + #plot_layout(guides = 'collect') +
  plot_annotation( tag_levels = 'A') 
ggsave(filename = "Bf.png", plot = Bfplot, device = "png", height = 25, width = 30,
       units = "cm", dpi = 500)
#ggplot punch
BFp <- dplyr::filter(B.p, isolate == "FR13")
BGp <- dplyr::filter(B.p, isolate == "Guy11")
library(ggpubr)
library(ggplot2)
bf.p.p <- ggplot(BFp,aes(x=factor(mutant, level = f.level_order), y=lesion.surface, fill = effector)) + geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface), show.legend = TRUE) + # jitter size by lesion surface
  scale_size_continuous(range = c(0.01, 2))  + theme_minimal() + scale_x_discrete(limits=rev)
BF.p.plot <- bf.p.p + scale_fill_manual(values = blue.palette) +
  labs(title = "FR13 punch inoculation", x = "Isolate", y = "Lesion surface (in pixels)") +
  # theme(axis.title.y = element_blank(), axis.title.x = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
bg.p.p <- ggplot(BGp, aes(x=factor(mutant, level = g.level_order), y=lesion.surface, fill = effector)) + geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface), show.legend = TRUE) + # jitter size by lesion surface
  scale_size_continuous(range = c(0.01, 2))  + theme_minimal() + scale_x_discrete(limits=rev)
BG.p.plot <- bg.p.p + scale_fill_manual(values = blue.palette) +
  labs(title = "Guy11 punch inoculation", x = "Isolate", y = "Lesion surface (in pixels)") +   
  #theme(axis.title.y = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0))   
library(patchwork)
Bpplot <- (BF.p.plot / BG.p.plot) + plot_annotation(
  title = 'Punch inoculation on Bl1')
Bpplot
ggsave(filename = "Bp.png", plot = Bpplot, device = "png", height = 20, width = 30,
       units = "cm", dpi = 500)
