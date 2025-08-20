# Mediation effects of Circadian Rest-activity Rhythm on the association between light durations and dementia
# Using Uk Biobank data
# Author: Wei Wang, Dr.rer.hum.biol
# Email: wei.wang@gzhmu.edu.cn
# © 2024 Wei Wang. All rights reserved.

library(broom)
library(dplyr)
library(readr)
library(survival)
library(mediation)
library(car)

# Load data
data.lightCRAR <- read_csv("your/data/path.csv", show_col_types = FALSE)

# Converting covariates to factors
data.lightCRAR$assessment.centre54 <- factor(data.lightCRAR$assessment.centre54)
data.lightCRAR$start.time.of.wear90010.season <- factor(data.lightCRAR$start.time.of.wear90010.season)
data.lightCRAR$highest.hour <- factor(data.lightCRAR$highest.hour)
data.lightCRAR$highest.twohours <- factor(data.lightCRAR$highest.twohours)

covariates <- c("age.acti", "sex31", "ethnic.white21000.2gp", "photoperiod", "VitaminDsupp.2gp","chronotype.4gp","townsend.deprivation.index189", "assessment.centre54",
                "Final_Education_FULL", "start.time.of.wear90010.season", "pm2.5absorbance.24007", "Final_Healthy_diet_score_FULL",
                "Final_Smoking_status_FULL", "Final_Alcohol_frequency_categories_FULL", "MVPA1005M400.Weekly.0.23",
                "Social.Network.Index3gp", "obesity", "Diabetes.comorbidities.FINAL", 
                "Hypertension.comorbidities.FINAL", "TraumaticBrainInjury.acti", "HearingLoss.acti")

# Loop through independent variables
for (iv_index in 97:100) {
  # Extract the independent variable name
  iv_name <- names(data.lightCRAR)[iv_index]
  
  # Construct the file path for the results
  directory_path <- "your/path/"
  full_file_path <- paste0(directory_path, "regression/results_", iv_name, ".csv")
  
  # Read the CSV file containing regression results
  regression_results <- read.csv(full_file_path)
  
  # Filter significant dependent variables
  significant_dependent_vars <- regression_results %>%
    filter(p.value < 0.05) %>%
    dplyr::select(dependent_var)
  
  # Initialize a data frame to store mediation analysis results for this independent variable
  mediation_results_iv <- data.frame()
  
  # Loop through significant dependent variables
  for (mediator_var in significant_dependent_vars$dependent_var) {
    # Check if there is valid data
    if (!all(is.na(data.lightCRAR[[mediator_var]])) && length(data.lightCRAR[[mediator_var]]) == nrow(data.lightCRAR)) {
      
      # Ensure the variable is numeric
      if (!is.numeric(data.lightCRAR[[mediator_var]])) {
        data.lightCRAR[[mediator_var]] <- as.numeric(as.factor(data.lightCRAR[[mediator_var]]))
      }
      
      # Fit the mediator model
      mediator_formula <- paste("scale(", mediator_var, ") ~ ", iv_name, " + ", paste(covariates, collapse = " + "))
      model_mediator <- lm(as.formula(mediator_formula), data = data.lightCRAR)
      
      # Fit the dependent variable model
      model_dv_formula <- paste("Surv(FUduration.acti.AllDementia.HDC.ICD10.y, Incidental.acti.AllDementia.HDC.ICD10.y) ~ ", iv_name, " + ", "scale(", mediator_var, ")", " + ", paste(covariates, collapse = " + "))
      model_dv <- survreg(as.formula(model_dv_formula), data = data.lightCRAR)
      
      # Check if the model coefficients contain NA
      if (any(is.na(coef(model_dv)))) {
        warning(paste("NA coefficients detected in the model for mediator variable", mediator_var, "and independent variable", iv_name))
        next
      }
      
      # Perform mediation analysis
      mediator_name <- paste("scale(", mediator_var, ")", sep="")
      mediation_results <- mediate(model_mediator, model_dv, treat = iv_name, mediator = mediator_name, sims = 5000)
      
      # Summarize mediation analysis results
      mediation_results_iv <- rbind(mediation_results_iv, data.frame(
        Mediator = mediator_name,
        ACME = mediation_results$d.avg,
        p_ACME = mediation_results$d.avg.p,
        ADE = mediation_results$z.avg,
        p_ADE = mediation_results$z.avg.p,
        Mediated = mediation_results$n.avg,
        p_Mediated = mediation_results$n.avg.p,
        lower.ACME = mediation_results$d.avg.ci[1],
        upp.ACME = mediation_results$d.avg.ci[2],
        lower.ADE = mediation_results$z.avg.ci[1],
        upp.ADE = mediation_results$z.avg.ci[2],
        lower.pro = mediation_results$n.avg.ci[1],
        upp.pro = mediation_results$n.avg.ci[2]
      ))
      med_sum <- summary(mediation_results)
      med_sum
    }
  }
  # Store the results in the list, with the key as iv_name
  all_results[[iv_name]] <- mediation_results_iv
}

directory_path2 <- "your/saving/path/"

# Write the results to a CSV file
for (iv_name in names(all_results)) {
  filename <- paste0(directory_path2, "MediationCRAR_", iv_name, ".csv")
  write.csv(all_results[[iv_name]], filename, row.names = FALSE)
}


