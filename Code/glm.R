# Load necessary R packages for data analysis
library(lme4)
library(lmtest)
library(car)
library(multcomp)
library(reshape2)

# Read the effect size data from a CSV file
ef.size <- read.csv('TestI_effect.size.csv', header = TRUE)

# Check the frequency distribution of categorical variables
table(ef.size$ecotype)
table(ef.size$treat)
table(ef.size$gen)
table(ef.size$trait)

# Replace some ecotype names for better readability
ef.size <- ef.size %>%
  mutate(ecotype = str_replace_all(ecotype, c('aaCol-0'='Col-0', 'Ang-1'='Ang-0', 'Olympia-0'='Oly-2')))

# Create a binary variable indicating whether Pvalue < 0.05
ef.size$value_sig <- ifelse(ef.size$Pvalue < 0.05, 1, 0)
ef.size <- ef.size[ef.size$treat!='actr1',]
# Create a new dataset with selected columns
d1 <- ef.size[, c(1:4, 10)]
d1$value_sig <- as.numeric(d1$value_sig)
# Reshape the data to wide format using 'dcast' and then back to long format using 'melt'
d2 <- d1 |> 
  dcast(ecotype + treat + gen ~ trait, value.var = 'value_sig')
d2$group <- 1:487
d3 <- d2 |> 
  melt(id.vars = c('ecotype', 'treat', 'gen', 'group'))
colnames(d3) <- c('ecotype', 'treat', 'gen', 'group', 'trait', 'value_sig')
d3$group <- factor(d3$group)

# Create a generalized linear mixed-effects model (glmer) with a binary outcome variable
mod <- glmer(value_sig ~ ecotype + treat + gen + trait + (1|group),
             family = binomial, data = d3)

# Perform Anova on the mixed-effects model
Anova(mod)

# Fit a generalized linear model (glm) to investigate the effect of 'trait'
d1$trait <- factor(d1$trait, levels = c('flt', 'fruit', 'biomass', 'diam', 'height'))
modle_trait <- glm(value_sig ~ trait, family = binomial, data = d1)
m1_trait <- summary(glht(modle_trait, linfct = mcp("trait" = "Tukey")))
m1_trait

# Fit a generalized linear model (glm) to investigate the effect of 'gen'
d1$gen <- as.factor(d1$gen)
modle_gen <- glm(value_sig ~ gen, family = binomial, data = d1)
m1_gen <- summary(glht(modle_gen, linfct = mcp("gen" = "Tukey")))
m1_gen

# Fit a generalized linear model (glm) to investigate the effect of 'treat'
d1$treat <- as.factor(d1$treat)
levels(d1$treat)
modle_treat <- glm(value_sig ~ treat, family = binomial, data = d1)
m1_treat <- summary(glht(modle_treat, linfct = mcp(treat = c(0, 0, 0, 0, 0, 0, 0, 1, -1))))
m3_treat <- summary(glht(modle_treat, linfct = mcp(treat = c(1, -1, 0, 0, 0, 0, 0, 0, 0))))
m1_treat
m3_treat
