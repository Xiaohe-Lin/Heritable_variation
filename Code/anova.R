# Load necessary R packages for data analysis
library(tidyverse)
library(lme4)
library(car)
library(lmerTest)

# Function to read data and prepare factors
read_and_prepare_data <- function(file_path){
  matr <- read.csv(file_path, header = TRUE)
  matr$P.layer <- factor(matr$P.layer)
  matr$layer <- factor(matr$layer)
  matr<-matr[which(matr$rm==FALSE),]
  return(matr)
}

# Function to perform for specified generations
perform_test <- function(matr, generations,  file_prefix){
  aa <- data.frame()
  bb <- data.frame()
  
  for (Generation in generations){
    matr.sub.gen <- matr[matr$Generation == Generation,]
    print(Generation)
    for (trait in c("flt", "height", "diam", "biomass", "fruit")){
      if (Generation == "Anc") {
        formu <- as.formula(paste(trait, "~ecotype+treat.1+ecotype:treat.1+(1|P.layer)+(1|P.row)"))
      } else {
        formu<-as.formula(paste(trait, "~layer+ecotype+treat.1+ecotype:treat.1+(1|P.layer)+(1|P.row)"))
      }
      print(formu)
      mod <- lmer(formu, na.action = na.omit, data = matr.sub.gen)
      a <- Anova(mod, test = "F", type = "II")
      b <- ranova(mod)
      a <- cbind(trait, Generation, a)
      b <- cbind(trait, Generation, b)
      aa <- rbind(aa, a)
      bb <- rbind(bb, b)
    }
  }
  
  write.csv(aa, paste(file_prefix, "_Anova_each_generations.csv", sep = ""))
  write.csv(bb, paste(file_prefix, "_ranova_each_generations.csv", sep = ""))
}

# Function to perform for all generations
perform_test_all_generations <- function(matr, generations, file_prefix){
  matr_unite <- unite(matr, "Layer", c("Generation", "layer"), sep = "_", remove = FALSE)
  matr_sub <- subset(matr_unite, Generation %in% generations)
  aa <- data.frame()
  bb <- data.frame()
  for (trait in c("flt", "height", "diam", "biomass", "fruit")){
    if (file_prefix == 'TestII') {
      formu <- as.formula(paste(trait, "~layer+ecotype+treat.1+Generation+ecotype*treat.1+ecotype*Generation+treat.1*Generation+ecotype*treat.1*Generation+(1|id)+(1|P.layer)+(1|P.row)"))
    } else{
      formu <- as.formula(paste(trait, "~Layer+ecotype+treat.1+ecotype*treat.1+ecotype*Generation+treat.1*Generation+ecotype*treat.1*Generation+(1|id)+(1|P.layer)+(1|P.row)"))
    }
    print(file_prefix)
    print(formu)
    mod <- lmer(formu, na.action = na.omit, data = matr_sub)
    a <- Anova(mod, test = "F", type = "II")
    b <- ranova(mod)
    a <- cbind(trait, a)
    aa <- rbind(aa, a)
    b <- cbind(trait, b)
    bb <- rbind(bb, b)
  }
  
  write.csv(aa, paste(file_prefix, "_Anova_all_generations.csv", sep = ""))
  write.csv(bb, paste(file_prefix, "_ranova_all_generations.csv", sep = ""))
}

# Read and prepare data
matr1 <- read_and_prepare_data("TestI_raw.csv")
matr2 <- read_and_prepare_data("TestII_raw.csv")
matr_Anc <- read_and_prepare_data("Anc_raw.csv")

perform_test(matr1, c("F1", "F2", "F3", "F4"), "TestI")
perform_test(matr2, c("F1", "F2", "F3", "F4"), "TestII")
perform_test(matr_Anc, c("Anc"), "Anc")


perform_test_all_generations(matr1, c("F1", "F2", "F3", "F4"), "TestI")
perform_test_all_generations(matr1, c("F1", "F2", "F3", "F4"), "TestII")
