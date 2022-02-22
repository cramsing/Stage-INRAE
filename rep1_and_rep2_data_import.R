#Masters thesis data analysis 

# library(dplyr)
# library(readr)

# install.packages("ggforce")
# install.packages("colorspace")
# install.packages("ggtext")
# install.packages("ggsci")
#install.packages("forcats")
#install.packages("patchwork")
# 
# #### Data import #### 

#Rep 1 import 
#all this will not be run because the code is too heavy to run each time. 
#Instead the resulting df is saved as an r file and loaded each time
# ## Compilation des fichiers csv All_lesions
# list_files = list.files(full.names = TRUE) # making list of files
# a=0
# b = "All_lesions.csv"
# 
# #alternative code if files aren't in the same working directory
# # list_files = list.files(path="D:/Documents/folder/folder/folder", full.names = TRUE)
# # a=0
# # b = "All_lesions.csv"
# 
# for( x in list_files) {
#   if (grepl(b,x,fixed=TRUE)) {
#     tmp = read_delim(file = x, delim = "\t")
#     if (a == 0) {    df = tmp; a=1  }
#     else {
#       df = bind_rows(df, tmp)
#     }
#   }
# }
# 
# df
# 
# ## Renaming column names to match file names
# colnames(df)
# names(df)[names(df) == "image"] <- "mutant_cultivar_date"
# 
# 
# df
# 
# 
# ## Separation of cultivar and mutant 
# df["mutant_cultivar"] = lapply(df["mutant_cultivar_date"],as.character)
# cs <- strsplit(df[["mutant_cultivar_date"]],"_")
# cs2 <- data.frame(do.call(rbind,cs),stringsAsFactors=FALSE)
# names(cs2)[names(cs2) == "X1"] <- "mutant"
# names(cs2)[names(cs2) == "X2"] <- "cultivar"
# names(cs2)[names(cs2) == "X3"] <- "date"
# 
# 
# df2 = cbind(df, cs2)
# 
# df2$genotype =paste0(df2$mutant, ".",df2$cultivar)
# 
# # deleting lesions with 0 size
# df2=df2[df2$lesion.surface>0,]
# 
# 
# head(df2)
# tail(df2)
# save(df2, file="df2.RData") #saving just to be sure (and bc the above code is quite annoying to run)

load("df2.RData")

df2[ , c('genotype', 'mutant_cultivar_date')] <- list(NULL) #deleting redunant columns 
colnames(df2)

rep1 <- df2[, c(16, 17, 18, 19, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)] #reordring columns and puttting in a new df
#really could've done this and previous deletions in the same step but hey 



# #Rep 2 import 
# #all this will not be run because the code is too heavy to run each time. 
# #Instead the resulting df is saved as an r file and loaded each time
# # ## Compilation des fichiers csv All_lesions
#  # list_files = list.files(full.names = TRUE) # making list of files
#  # a=0
#  # b = "All_lesions.csv"
#  
# # #alternative code if files aren't in the same working directory
# library(dplyr)
# library(readr)
# list_files = list.files(path="E:/results 23.08", full.names = TRUE)
# a=0
# b = "All_lesions.csv"
# 
# for( x in list_files) {
#   if (grepl(b,x,fixed=TRUE)) {
#     tmp = read_delim(file = x, delim = "\t")
#     if (a == 0) {    df1 = tmp; a=1  }
#     else {
#       df1 = bind_rows(df1, tmp)
#     }
#   }
# }
# 
# df1
#  
# # ## Renaming column names to match file names
# colnames(df1)
# names(df1)[names(df1) == "image"] <- "mutant_cultivar"
# # 
# # 
# df1
# # 
# # 
# # ## Separation of cultivar and mutant 
# df1["mutant_cultivar"] = lapply(df1["mutant_cultivar"],as.character)
# cs <- strsplit(df1[["mutant_cultivar"]],"_")
# cs2 <- data.frame(do.call(rbind,cs),stringsAsFactors=FALSE)
# names(cs2)[names(cs2) == "X1"] <- "mutant"
# names(cs2)[names(cs2) == "X2"] <- "cultivar"
# 
# # 
# # 
# df3 = cbind(df1, cs2)
# # 
# df3$genotype =paste0(df3$mutant, ".",df3$cultivar)
# # 
# #deleting lesions with 0 size
# df3=df3[df3$lesion.surface>0,]
# # 
# # 
# head(df3)
# tail(df3)
# save(df3, file="df3.RData") #saving just to be sure (and bc the above code is quite annoying to run)
