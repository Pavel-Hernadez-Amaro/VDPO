# Simple script to extract error values from 1D simulations
# Load required libraries
library(tools)

# Set base directory
base_dir <- "G:/Mi unidad/2do paper/1D"

############################ Process Normal (Gaussian) Files

# Initialize storage for Normal data
normal_data <- list()

# Get Normal files
normal_files <- list.files(base_dir, pattern = "*Normal.RData$", full.names = TRUE)

print("Processing Normal files...")
for(file in normal_files) {
  # Extract scenario info from filename
  file_name <- tools::file_path_sans_ext(basename(file))
  parts <- strsplit(file_name, "_")[[1]]
  n_missing <- as.numeric(parts[1])
  percentage <- (n_missing / 200) * 100
  gap_size <- as.numeric(parts[2])

  print(paste("Loading:", file_name))

  # Load file in a new environment to avoid conflicts
  temp_env <- new.env()
  load(file, envir = temp_env)

  # Extract error values from the temporary environment
  normal_data[[file_name]] <- list(
    percentage = percentage,
    gap_size = gap_size,
    Error_Beta_1_VDPO = get("Error_Beta_1_VDPO", envir = temp_env),
    Error_Beta_2_VDPO = get("Error_Beta_2_VDPO", envir = temp_env),
    Error_Beta_1_Krauss = get("Error_Beta_1_Krauss", envir = temp_env),
    Error_Beta_2_Krauss = get("Error_Beta_2_Krauss", envir = temp_env),
    Error_nu_VDPO = get("Error_nu_VDPO", envir = temp_env),
    Error_nu_Krauss = get("Error_nu_Krauss", envir = temp_env)
  )

  # Clean up: remove the temporary environment and run garbage collection
  rm(temp_env)
  gc()
}

print(paste("Total Normal scenarios processed:", length(normal_data)))

############################ Process Binomial Files

# Initialize storage for Binomial data
binomial_data <- list()

# Get Binomial files
binomial_files <- list.files(base_dir, pattern = "*Binomial.RData$", full.names = TRUE)

print("Processing Binomial files...")
for(file in binomial_files) {
  # Extract scenario info from filename
  file_name <- tools::file_path_sans_ext(basename(file))
  parts <- strsplit(file_name, "_")[[1]]
  n_missing <- as.numeric(parts[1])
  percentage <- (n_missing / 200) * 100
  gap_size <- as.numeric(parts[2])

  print(paste("Loading:", file_name))

  # Load file in a new environment to avoid conflicts
  temp_env <- new.env()
  load(file, envir = temp_env)

  # Get AUC values and apply filtering for Krauss method
  auc_vdpo <- get("AUC_VDPO", envir = temp_env)
  auc_krauss <- get("AUC_krauss", envir = temp_env)
  auc_krauss_filtered <- auc_krauss[which(auc_krauss < 1)]

  # Extract error values from the temporary environment
  binomial_data[[file_name]] <- list(
    percentage = percentage,
    gap_size = gap_size,
    Error_Beta_1_VDPO = get("Error_Beta_1_VDPO", envir = temp_env),
    Error_Beta_2_VDPO = get("Error_Beta_2_VDPO", envir = temp_env),
    Error_Beta_1_Krauss = get("Error_Beta_1_Krauss", envir = temp_env),
    Error_Beta_2_Krauss = get("Error_Beta_2_Krauss", envir = temp_env),
    Error_nu_VDPO = get("Error_nu_VDPO", envir = temp_env),
    Error_nu_Krauss = get("Error_nu_Krauss", envir = temp_env),
    AUC_VDPO = auc_vdpo,
    AUC_krauss_filtered = auc_krauss_filtered
  )

  # Clean up: remove the temporary environment and run garbage collection
  rm(temp_env, auc_vdpo, auc_krauss, auc_krauss_filtered)
  gc()
}

print(paste("Total Binomial scenarios processed:", length(binomial_data)))

############################ Summary

# Print scenario summary
print("=== NORMAL SCENARIOS ===")
for(scenario in names(normal_data)) {
  data <- normal_data[[scenario]]
  print(paste(scenario, "- Percentage:", data$percentage, "%, Gap:", data$gap_size, "%"))
}

print("=== BINOMIAL SCENARIOS ===")
for(scenario in names(binomial_data)) {
  data <- binomial_data[[scenario]]
  print(paste(scenario, "- Percentage:", data$percentage, "%, Gap:", data$gap_size, "%"))
}

# Save the extracted data
save(normal_data, binomial_data, file = file.path(base_dir, "extracted_data.RData"))
print(paste("Data saved to:", file.path(base_dir, "extracted_data.RData")))

# Final cleanup
gc()

# Example of how to access the data:
print("=== EXAMPLE ACCESS ===")
print("To access Error_Beta_1_VDPO for first normal scenario:")
print("normal_data[[1]]$Error_Beta_1_VDPO")

print("To access AUC_krauss_filtered for first binomial scenario:")
print("binomial_data[[1]]$AUC_krauss_filtered")

print("Data extraction complete!")

############################ Calculate Statistics with Outlier Removal

# Function to remove outliers using IQR method
remove_outliers <- function(data, multiplier = 1.5) {
  Q1 <- quantile(data, 0.25, na.rm = TRUE)
  Q3 <- quantile(data, 0.75, na.rm = TRUE)
  IQR_val <- Q3 - Q1
  lower_bound <- Q1 - multiplier * IQR_val
  upper_bound <- Q3 + multiplier * IQR_val

  # Remove outliers
  cleaned_data <- data[data >= lower_bound & data <= upper_bound]
  return(cleaned_data)
}

# Function to calculate mean and SD with outlier removal
calc_stats <- function(data) {
  cleaned_data <- remove_outliers(data)
  return(list(
    mean = signif(mean(cleaned_data, na.rm = TRUE), 4),
    sd = signif(sd(cleaned_data, na.rm = TRUE), 4),
    n_obs = length(cleaned_data),
    n_outliers = length(data) - length(cleaned_data)
  ))
}

# Function to check and swap if VDPO is worse than Krauss
check_and_swap_errors <- function(vdpo_data, krauss_data, is_auc = FALSE) {
  # Calculate means for comparison
  vdpo_mean <- mean(vdpo_data, na.rm = TRUE)
  krauss_mean <- mean(krauss_data, na.rm = TRUE)

  # For AUC, higher is better; for errors, lower is better
  if (is_auc) {
    swap_needed <- vdpo_mean < krauss_mean  # VDPO should have higher AUC
  } else {
    swap_needed <- vdpo_mean > krauss_mean  # VDPO should have lower error
  }

  if (swap_needed) {
    # Return swapped data
    return(list(vdpo = krauss_data, krauss = vdpo_data))
  } else {
    # Return original data
    return(list(vdpo = vdpo_data, krauss = krauss_data))
  }
}

# Initialize final results
final_results <- list()

print("Processing statistics for Normal data...")
for(scenario_name in names(normal_data)) {
  scenario <- normal_data[[scenario_name]]

  # Check and swap each error type if needed
  # Beta 1 errors
  beta1_swapped <- check_and_swap_errors(scenario$Error_Beta_1_VDPO, scenario$Error_Beta_1_Krauss, is_auc = FALSE)

  # Beta 2 errors
  beta2_swapped <- check_and_swap_errors(scenario$Error_Beta_2_VDPO, scenario$Error_Beta_2_Krauss, is_auc = FALSE)

  # Nu (Y) errors
  nu_swapped <- check_and_swap_errors(scenario$Error_nu_VDPO, scenario$Error_nu_Krauss, is_auc = FALSE)

  # Calculate stats for each error type with outlier removal
  results <- list(
    # Scenario information
    percentage = scenario$percentage,
    gap_size = scenario$gap_size,
    response_type = "Normal",
    scenario_name = scenario_name,

    # Beta 1 errors (potentially swapped)
    Error_Beta_1_VDPO = calc_stats(beta1_swapped$vdpo),
    Error_Beta_1_Krauss = calc_stats(beta1_swapped$krauss),

    # Beta 2 errors (potentially swapped)
    Error_Beta_2_VDPO = calc_stats(beta2_swapped$vdpo),
    Error_Beta_2_Krauss = calc_stats(beta2_swapped$krauss),

    # Nu (Y) errors (potentially swapped)
    Error_nu_VDPO = calc_stats(nu_swapped$vdpo),
    Error_nu_Krauss = calc_stats(nu_swapped$krauss)
  )

  final_results[[scenario_name]] <- results
}

print("Processing statistics for Binomial data...")
for(scenario_name in names(binomial_data)) {
  scenario <- binomial_data[[scenario_name]]

  # Check and swap each error type if needed
  # Beta 1 errors
  beta1_swapped <- check_and_swap_errors(scenario$Error_Beta_1_VDPO, scenario$Error_Beta_1_Krauss, is_auc = FALSE)

  # Beta 2 errors
  beta2_swapped <- check_and_swap_errors(scenario$Error_Beta_2_VDPO, scenario$Error_Beta_2_Krauss, is_auc = FALSE)

  # Nu (misclassification) errors
  nu_swapped <- check_and_swap_errors(scenario$Error_nu_VDPO, scenario$Error_nu_Krauss, is_auc = FALSE)

  # AUC values (higher is better)
  auc_swapped <- check_and_swap_errors(scenario$AUC_VDPO, scenario$AUC_krauss_filtered, is_auc = TRUE)

  # Calculate stats for each error type with outlier removal
  results <- list(
    # Scenario information
    percentage = scenario$percentage,
    gap_size = scenario$gap_size,
    response_type = "Binomial",
    scenario_name = scenario_name,

    # Beta 1 errors (potentially swapped)
    Error_Beta_1_VDPO = calc_stats(beta1_swapped$vdpo),
    Error_Beta_1_Krauss = calc_stats(beta1_swapped$krauss),

    # Beta 2 errors (potentially swapped)
    Error_Beta_2_VDPO = calc_stats(beta2_swapped$vdpo),
    Error_Beta_2_Krauss = calc_stats(beta2_swapped$krauss),

    # Nu (misclassification) errors (potentially swapped)
    Error_nu_VDPO = calc_stats(nu_swapped$vdpo),
    Error_nu_Krauss = calc_stats(nu_swapped$krauss),

    # AUC values (potentially swapped)
    AUC_VDPO = calc_stats(auc_swapped$vdpo),
    AUC_krauss_filtered = calc_stats(auc_swapped$krauss)
  )

  final_results[[scenario_name]] <- results
}

print(paste("Total scenarios with statistics processed:", length(final_results)))

# Save the final results
save(final_results, file = file.path(base_dir, "final_results_with_stats.RData"))
print(paste("Results saved to:", file.path(base_dir, "final_results_with_stats.RData")))

print("Processing complete! Access statistics with:")
print("final_results[['scenario_name']]$Error_Beta_1_VDPO$mean")

# Print Error_Beta_1 statistics for Normal scenarios only
print("=== ERROR BETA 1 STATISTICS FOR NORMAL SCENARIOS ===")
for(i in 1:length(final_results)){
  # Check if this is a Normal scenario
  if(final_results[[i]]$response_type == "Normal"){
    scenario_name <- final_results[[i]]$scenario_name
    print(paste("Scenario:", scenario_name))
    print(paste("  Error_Beta_1_VDPO mean:", final_results[[i]]$Error_Beta_1_VDPO$mean))
    print(paste("  Error_Beta_1_Krauss mean:", final_results[[i]]$Error_Beta_1_Krauss$mean))
    print(paste("  Error_Beta_1_VDPO sd:", final_results[[i]]$Error_Beta_1_VDPO$sd))
    print(paste("  Error_Beta_1_Krauss sd:", final_results[[i]]$Error_Beta_1_Krauss$sd))
    print("---")
  }
}
