Calculate_and_export_ef_size <- function(input_file, output_file, treat_list, gen_list, trait_list) {
  # Load necessary R packages for data analysis
  library(lme4)
  library(lmerTest)
  library(stringr)
  
  # Function to read data and prepare factors
  read_and_prepare_data <- function(file_path){
    matr <- read.csv(file_path, header = TRUE)
    matr<-matr[which(matr$rm==FALSE),]
    matr$P.layer <- factor(matr$P.layer)
    matr$layer <- factor(matr$layer)
    return(matr)
  }
  
  # Function to calculate effect size
  Calculate_ef_size <- function(matr, treat_list, gen_list, trait_list) {
    result_list <- lapply(trait_list, function(trait) {
      trait_results <- lapply(treat_list, function(treat) {
        if (treat == 'JA') {
          treat_level = 'JActr'
        } else {
          treat_level = 'actr'
        }
        matr.sub <- matr[matr$treat.1 %in% treat | matr$treat.1 == treat_level, ]
        matr.sub$treat.1 <- relevel(factor(matr.sub$treat.1), ref = treat_level)
        result <- lapply(gen_list, function(gen) {
          if (gen == 'Anc') {
            formu <- as.formula(paste(trait, "~ ecotype + ecotype:treat.1 + (1|P.layer) + (1|P.row)"))
          } else{
            formu <- as.formula(paste(trait, "~ layer + ecotype + ecotype:treat.1 + (1|P.layer) + (1|P.row)")) 
          }
          print(paste0('Processing trait ',trait, ' of treatment ', treat,' in ', gen, ' generation'))
          print(formu)
          mod <- lmer(formu, na.action = na.omit, data = matr.sub[matr.sub$Generation == gen, ])
          coef_summary <- coef(summary(mod))
          data.frame(treat = treat, gen = gen, trait = trait, coef_summary)
        })
        do.call(rbind, result)
      })
      do.call(rbind, trait_results)
    })
    aa <- do.call(rbind, result_list)
    return(aa)
  }
  
  # Read and prepare data
  matr1 <- read_and_prepare_data(input_file)
  
  # Calculate effect size
  result <- Calculate_ef_size(matr1, treat_list, gen_list, trait_list)
  
  # Rename columns
  colnames(result) <- c("treat", "gen", "trait", "Estimate", "Std.Error", "df", "tvalue", "Pvalue")
  
  # Filter and process data
  result <- as.data.frame(result)
  result <- result[grepl('ecotype', rownames(result)) & grepl('treat', rownames(result)), ]
  result$ecotype <- gsub(".*ecotype(.+):treat.*", "\\1", rownames(result))
  result <- result[, c(9, 1:8)]
  
  # Export result to CSV
  write.csv(result, output_file, row.names = FALSE, quote = FALSE)
}

# Define your custom treatment, generation, and trait lists
custom_treat_list <- c("Drought", "Salt1", "Salt2", "Cd1", "Cd2", "Nutrient1", "Nutrient2", "Leaf", "JA")
custom_gen_list <- c("F1", "F2", "F3", "F4")
custom_trait_list <-  c("flt", "height", "diam", "biomass", "fruit")

# Call the function with custom lists
Calculate_and_export_ef_size('TestI_raw.1.csv', 'TestI_effect.size.csv', custom_treat_list, custom_gen_list, custom_trait_list)

custom_treat_list <- c("Drought", "Salt1", "Salt2", "Cd1", "Cd2", "Nutrient1", "Nutrient2", "Leaf", "JA")
custom_gen_list <- c("Anc")
custom_trait_list <-  c("flt", "height", "diam", "biomass", "fruit")
Calculate_and_export_ef_size('Anc_raw.csv', 'Anc_effect.size.csv', custom_treat_list, custom_gen_list, custom_trait_list)

custom_treat_list <- c("Drought", "Salt1", "Salt2", "Cd1", "Cd2")
custom_gen_list <- c("F1", "F2", "F3", "F4")
custom_trait_list <-  c("flt", "height", "diam", "biomass", "fruit")
Calculate_and_export_ef_size('TestII_raw.csv', 'TestII_effect.size.csv', custom_treat_list, custom_gen_list, custom_trait_list)