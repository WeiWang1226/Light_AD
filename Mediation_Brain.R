# Mediation effects of Brain structures on the association between light durations and dementia
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
data_with_residuals <- read_csv("data/with/DPIs.csv", show_col_types = FALSE)

# Convert covariates that need to be factors into factors
data_with_residuals$assessment.centre54 <- factor(data_with_residuals$assessment.centre54)
data_with_residuals$start.time.of.wear90010.season <- factor(data_with_residuals$start.time.of.wear90010.season)
data_with_residuals$highest.hour <- factor(data_with_residuals$highest.hour)
data_with_residuals$highest.twohours <- factor(data_with_residuals$highest.twohours)

# Define covariates
covariates <- c("days_int", "age.acti", "sex31", "ethnic.white21000.2gp", "photoperiod", "VitaminDsupp.2gp","chronotype.4gp","townsend.deprivation.index189", "assessment.centre54",
                "Final_Education_FULL", "start.time.of.wear90010.season", "pm2.5absorbance.24007", "Final_Healthy_diet_score_FULL",
                "Final_Smoking_status_FULL", "Final_Alcohol_frequency_categories_FULL", "MVPA1005M400.Weekly.0.23",
                "Social.Network.Index3gp", "obesity", "Diabetes.comorbidities.FINAL", 
                "Hypertension.comorbidities.FINAL", "TraumaticBrainInjury.acti", "HearingLoss.acti")

# Loop through the indices of the independent variables of interest
for (iv_index in 685:688) {
  
  # Extract the independent variable name
  iv_name <- names(data_with_residuals)[iv_index]
  
  # Construct the file path for the results
  directory_path <- "your/path/"
  full_file_path <- paste0(directory_path, "regression/result_", iv_name, ".csv")
  
  # Read the CSV file containing regression results
  regression_results <- read.csv(full_file_path)
  
  # Filter significant dependent variables
  significant_dependent_vars <- regression_results %>%
    filter(p.value < 0.05) %>%
    dplyr::select(dependent_var)
  
  # Initialize a data frame to store mediation analysis results for this iv_name
  mediation_results_iv <- data.frame()
  
  # Loop through significant dependent variables
  for (mediator_var in significant_dependent_vars$dependent_var) {
    # Check if valid data exists
    if (!all(is.na(data_with_residuals[[mediator_var]])) && length(data_with_residuals[[mediator_var]]) == nrow(data_with_residuals)) {
      
      # Ensure the variable is numeric
      if (!is.numeric(data_with_residuals[[mediator_var]])) {
        data_with_residuals[[mediator_var]] <- as.numeric(as.factor(data_with_residuals[[mediator_var]]))
      }
      
      # Fit the mediator model
      mediator_formula <- paste("scale(", mediator_var, ") ~ ", iv_name, " + ", paste(covariates, collapse = " + "))
      model_mediator <- lm(as.formula(mediator_formula), data = data_with_residuals)
      
      # Fit the dependent variable model
      model_dv_formula <- paste("Surv(FUduration.acti.AllDementia.HDC.ICD10.y, Incidental.acti.AllDementia.HDC.ICD10.y) ~ ", iv_name, " + ", "scale(", mediator_var, ")", " + ", paste(covariates, collapse = " + "))
      model_dv <- survreg(as.formula(model_dv_formula), data = data_with_residuals)
      
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
  
  # Store the results in the list with the key as iv_name
  all_results[[iv_name]] <- mediation_results_iv
}

# Write the results to a CSV file
for (iv_name in names(all_results)) {
  filename <- paste0(directory_path, "mediation_", iv_name, ".csv")
  write.csv(all_results[[iv_name]], filename, row.names = FALSE)
}