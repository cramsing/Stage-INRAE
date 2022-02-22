#Tsuyuake
load("spray.RData")
load("lc.RData")
#no Tsuyuake used in the punch inoculations

T <- spray[spray$cultivar=="Tsuyuake",] ; titre <- "Tsuyuake"
T.1 <- dplyr::filter(T, rep == "1") #B for rep 1
T.2 <- dplyr::filter(T, rep == "2") #B for rep 2
Tlc <- lc[lc$cultivar =="Tsuyuake",] ; titre <- "Tsuyuake"
Tlc.1 <-dplyr::filter(Tlc, rep == "1") # Blc for rep 1
Tlc.2 <-dplyr::filter(Tlc, rep == "2") # Blc for rep 2

# Stats ----

## Lesion surface ----
#Initial stats lesion surface
T0 <- aov(lesion.surface ~ mutant, data = T) 
anova(T0)
#p-value = 3.577e-07 ***, significant
#normality of residuals check
#distrubution/normality
par(mfrow=c(1,2))
#Histogramme
hist(T0$residuals) #histogram to check normality
# Q-Q Plot
qqnorm(T0$residuals) # qqplot to check normality 
qqline(T0$residuals)
#very anormal
shapiro.test(T0$residuals)
#too many entries to run shapiro but clearly anormal
kruskal.test(lesion.surface ~ mutant, data = T)
#p-value <2.2e-16, significant !

#Kruskal test with only rep 1
kruskal.test(lesion.surface ~ mutant, data = T.1) #p-value <2.2e-16
#Kruskal test with only rep 2
kruskal.test(lesion.surface ~ mutant, data = T.2) #p-value 1.044e-11
#not the same values for both

#Post-hoc 
#post-hoc test Dunn
library(FSA)
dunnTest(lesion.surface ~ mutant, data = T, method ="bh")
# #post-hoc test Wilcoxon rank sum test
#will allow us to easily add significance levels to our ggplots
library(ggpubr)
T.ls.w <- compare_means(lesion.surface ~ mutant, T, method="wilcox.test",
                         ref.group = ".all.")
T.ls.w
#Wilcoxon test with only rep 1
compare_means(lesion.surface ~ mutant, T.1, method="wilcox.test",
              ref.group = ".all.") 
#Wilcoxon test with only rep 2
compare_means(lesion.surface ~ mutant, T.2, method="wilcox.test",
              ref.group = ".all.") 
#actually fairly different results... 

# Lesion number ---- 
T1 <- aov(lesion.count ~ mutant, data = Tlc) 
anova(T1)
#p-value = 0.2113, not significant
#normality of residuals check
#distrubution/normality
par(mfrow=c(1,2))
#Histogramme
hist(T1$residuals)
# Q-Q Plot
qqnorm(T1$residuals)
qqline(T1$residuals)
#kind of normal looking
shapiro.test(T1$residuals)
#0.0001963, slightly not normal
kruskal.test(lesion.count ~ mutant, data = Tlc)
#p = 0.1083, not significant

#Kruskal test with only rep 1
kruskal.test(lesion.count ~ mutant, data = Tlc.1) #p-value  0.01895, not significant
#Kruskal test with only rep 2
kruskal.test(lesion.count ~ mutant, data = Tlc.2) #p-value 0.01751, not significant
#same values for both

#Post-hoc 

#post-hoc test Dunn
library(FSA)
dunnTest(lesion.count ~ mutant, data = Tlc, method ="bh")
# #post-hoc test Wilcoxon rank sum test
#will allow us to easily add significance levels to our ggplots
library(ggpubr)
T.lc.w <- compare_means(lesion.count ~ mutant, Tlc, method="wilcox.test",
                         ref.group = ".all.") #issue with this code, will come back to it

#Wilcoxon test with only rep 1
compare_means(lesion.count ~ mutant, Tlc.1, method="wilcox.test",
              ref.group = ".all.") 
#Wilcoxon test with only rep 2
compare_means(lesion.count ~ mutant, Tlc.2, method="wilcox.test",
              ref.group = ".all.") 


## Compiled kruskal groups ----
#lesion surface kruskal groups
library(agricolae)
library(dplyr)
T.ls.groups <- kruskal(T$lesion.surface,
                          T$mutant, group = TRUE) %>%
  .$groups %>%
  as_tibble(rownames = "mutant") %>%
  rename("T.ls.groups" = "groups")
T.ls.groups$`T$lesion.surface` <- NULL
#lesion count kruskal groups
T.lc.groups <- kruskal(Tlc$lesion.count,
                          Tlc$mutant, group = TRUE) %>%
  .$groups %>%
  as_tibble(rownames = "mutant") %>%
  rename("T..lc.groups" = "groups")
#compilation
library(dplyr)
T.groups <- dplyr::left_join(T.ls.groups, T.lc.groups, by = "mutant") 
T.groups <- T.groups[order(T.groups$mutant),]
write.csv2(T.groups, file = "T.groups.csv")


# Ggplots ----

#ggplot lesion size
library(ggpubr)
t.ls.p <- ggplot(T,aes(x=factor(mutant, level = f.level_order), y=lesion.surface, fill = effector)) +   
  geom_boxplot() + 
  geom_jitter(aes(size = lesion.surface, shape = rep, color = rep), show.legend = TRUE) + 
  scale_size_continuous(range = c(0.1, 2))+ theme_minimal() + scale_x_discrete(limits=rev) + 
  scale_y_continuous(trans = "log10", breaks = c(10, 100, 500, 1000, 5000, 10000))  
T.ls.plot <- t.ls.p + scale_colour_grey() + scale_fill_manual(values = blue.palette) +
  labs(title = "FR13 lesion surface", x = "Isolate", y = "Lesion surface (in log10 pixels)") +
  # theme(axis.title.y = element_blank(), axis.title.x = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) +  coord_flip()

ggsave(filename = "Tls.png", plot = T.ls.plot, device = "png", height = 20, width = 30,
       units = "cm", dpi = 500)


#ggplot lesion number
library(ggpubr)
t.lc0 <- ggplot(Tlc, aes(x=factor(mutant, level = f.level_order), y=lesion.count, fill = effector)) + 
  geom_boxplot() + geom_jitter(size=2.0, aes(shape = rep, color = rep)) + 
  theme_minimal() + scale_x_discrete(limits=rev) + coord_flip()
Tlc <- t.lc0 +  scale_colour_grey() + scale_fill_manual(values = blue.palette) +
  labs(title = "FR13 lesion number", x = "Isolate", y = "Lesion number") + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 

ggsave(filename = "Tlc.png", plot = Tlc, device = "png", height = 20, width = 30,
       units = "cm", dpi = 500)


library(patchwork)
Tplot <- (Tlc /T.ls.plot ) + plot_annotation(title = 'FR13 sion size and count on Tsuyuake', tag_levels = 'A') 
ggsave(filename = "T.png", plot = Tplot, device = "png", height = 25, width = 30,
       units = "cm", dpi = 500)
