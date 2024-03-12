#Install required packages----

#List of packages
pkg_list <- c("ggplot2", "dplyr", "readxl", "smatr", "cowplot")

#Function to check and install missing packages
check_and_install <- function(pkg){
  if (!require(pkg, character.only = TRUE)){
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

#Apply the function to the list of packages
sapply(pkg_list, check_and_install)

#Load the required packages

library(ggplot2)
library(dplyr)
library(readxl)
library(smatr)
library(cowplot)

#Import dataset----

path_to_dataset<- "insert path file for your machine" 

df_full <- read_excel(path_to_dataset)

#Subset aboveground data for full dataset

df_full_aboveground <- df_full[which(df_full$Organ=='Leaf'|df_full$Organ=='Twig'| df_full$Organ=='Branch'| df_full$Organ=='Stem'), ]
df_full_aboveground$Organ<- factor(df_full_aboveground$Organ, levels=(c("Leaf", "Twig", "Branch", "Stem"))) 

##Estimating minimum population size needed to estimate population mean for different organs----

#Summarize by organ

z_crit_95percCI<- 1.960
margin_of_error<- 0.05

df_full_summary <- df_full %>%
  group_by(Organ) %>%
  summarize(
    pop_std_dev = sd(DAVG)*sqrt((length(DAVG)-1)/length(DAVG)),
    min_sample_size = ((z_crit_95percCI^2)*(pop_std_dev^2))/margin_of_error^2
  )

#Simulating the effects of data preprocessing decisions----

##Bootstrapping data at different sample sizes, averaging the data, and fitting a model----

## Function to perform random sampling, fit SMA model, and extract coefficients

run_iteration <- function(df_full_aboveground, sample_size) {
  # Random sampling of D across each and every L
  bootstrapped_data <- df_full_aboveground %>%
    group_by(L) %>%
    slice_sample(n = sample_size, replace = TRUE) %>%
    ungroup()
  
  # Calculate the average DAVG for each L
  averaged_data <- bootstrapped_data %>%
    group_by(L) %>%
    summarize(DAVG = mean(DAVG))
  
  # Fit a SMA model to averaged_data
  bootstrapped_scaling <- sma(data = averaged_data, DAVG ~ L, log = "XY", method = "SMA")
  
  # Extract coefficients
  bootstrapped_scaling_summary <- bootstrapped_scaling$groupsummary
  bootstrapped_results <- c(
    slope = bootstrapped_scaling_summary$Slope,
    y_intercept = bootstrapped_scaling_summary$Int,
    slope_lowCI = bootstrapped_scaling_summary$Slope_lowCI,
    slope_highCI = bootstrapped_scaling_summary$Slope_highCI
  )
  
  # Return a list of coefficients
  return(bootstrapped_results)
}

## Specify the sample sizes and number of iterations
sample_sizes <- c(1, 5, 10, 25, 50, 75, 100,1000)
num_iterations <- 1000  

## Create an empty dataframe to store the results
bootstrapped_results_df <- data.frame(
  iteration = integer(),
  mean_slope = double(),
  mean_y_intercept = double(),
  mean_slope_lowCI = double(),
  mean_slope_highCI = double(),
  se_slope = double(),
  se_y_intercept = double(),
  se_slope_lowCI = double(),
  se_slope_highCI = double(),
  cv_slope= double(),
  cv_slope_lowCI = double(),
  cv_slope_highCI = double(),
  sample_size = integer()
)

## Loop through different sample sizes
for (sample_size in sample_sizes) {
  for (num_iteration in num_iterations) {
    slopes <- numeric(num_iteration)
    y_intercepts <- numeric(num_iteration)
    slope_lowCIs <- numeric(num_iteration)
    slope_highCIs <- numeric(num_iteration)
    
    # Run multiple iterations
    for (iteration in 1:num_iteration) {
      # Set seed for this iteration based on the sample size and iteration
      set.seed(iteration + sample_size)
      
      # Run the iteration and extract the coefficients
      coefficients_iteration <- run_iteration(df_full_aboveground, sample_size)
      
      # Store the results
      slopes[iteration] <- coefficients_iteration["slope"]
      y_intercepts[iteration] <- coefficients_iteration["y_intercept"]
      slope_lowCIs[iteration] <- coefficients_iteration["slope_lowCI"]
      slope_highCIs[iteration] <- coefficients_iteration["slope_highCI"]
    }
    
    # Calculate mean and standard error for each coefficient
    mean_slope <- mean(slopes)
    se_slope <- sd(slopes) / sqrt(num_iteration)
    
    mean_y_intercept <- mean(y_intercepts)
    se_y_intercept <- sd(y_intercepts) / sqrt(num_iteration)
    
    mean_slope_lowCI <- mean(slope_lowCIs)
    se_slope_lowCI <- sd(slope_lowCIs) / sqrt(num_iteration)
    
    mean_slope_highCI <- mean(slope_highCIs)
    se_slope_highCI <- sd(slope_highCIs) / sqrt(num_iteration)
    
    cv_slope<- sd(slopes) / mean(slopes)
    cv_slope_lowCI<- sd(slope_lowCIs) / mean(slope_lowCIs)
    cv_slope_highCI<- sd(slope_highCIs) / mean(slope_highCIs)
    
    
    
    # Store the results in the dataframe
    bootstrapped_results_df <- rbind(bootstrapped_results_df, c(iteration, mean_slope, mean_y_intercept, mean_slope_lowCI, mean_slope_highCI, se_slope, se_y_intercept, se_slope_lowCI, se_slope_highCI, cv_slope, cv_slope_lowCI, cv_slope_highCI, sample_size))
  }
}

## Rename the columns in bootstrapped_results_df
colnames(bootstrapped_results_df) <- c("iteration", "mean_slope", "mean_y_intercept", "mean_slope_lowCI", "mean_slope_highCI", "se_slope", "se_y_intercept", "se_slope_lowCI", "se_slope_highCI","cv_slope", "cv_slope_lowCI", "cv_slope_highCI",  "sample_size")

##Store the results

bootstrap_simulation_arith_mean_summary<- bootstrapped_results_df

##Bootstrapping data at different sample sizes, averaging the data and weighing by conductance (Dh), and fitting a model----

## Function to perform random sampling, fit SMA model, and extract coefficients

run_iteration <- function(df_full_aboveground, sample_size) {
  # Random sampling of D across each and every L
  bootstrapped_data <- df_full_aboveground %>%
    group_by(L) %>%
    slice_sample(n = sample_size, replace = TRUE) %>%
    ungroup()
  
  # Calculate the average DAVG for each L
  averaged_data <- bootstrapped_data %>%
    group_by(L) %>%
    summarize(DhAVG = sum(DAVG^5)/sum(DAVG^4))
  
  # Fit a SMA model to averaged_data
  bootstrapped_scaling <- sma(data = averaged_data, DhAVG ~ L, log = "XY", method = "SMA")
  
  # Extract coefficients
  bootstrapped_scaling_summary <- bootstrapped_scaling$groupsummary
  bootstrapped_results <- c(
    slope = bootstrapped_scaling_summary$Slope,
    y_intercept = bootstrapped_scaling_summary$Int,
    slope_lowCI = bootstrapped_scaling_summary$Slope_lowCI,
    slope_highCI = bootstrapped_scaling_summary$Slope_highCI
  )
  
  # Return a list of coefficients
  return(bootstrapped_results)
}

## Specify the sample sizes and number of iterations
sample_sizes <- c(1, 5, 10, 25, 50, 75, 100,1000)
num_iterations <- 1000  

## Create an empty dataframe to store the results
bootstrapped_results_df <- data.frame(
  iteration = integer(),
  mean_slope = double(),
  mean_y_intercept = double(),
  mean_slope_lowCI = double(),
  mean_slope_highCI = double(),
  se_slope = double(),
  se_y_intercept = double(),
  se_slope_lowCI = double(),
  se_slope_highCI = double(),
  cv_slope= double(),
  cv_slope_lowCI = double(),
  cv_slope_highCI = double(),
  sample_size = integer()
)

## Loop through different sample sizes
for (sample_size in sample_sizes) {
  for (num_iteration in num_iterations) {
    slopes <- numeric(num_iteration)
    y_intercepts <- numeric(num_iteration)
    slope_lowCIs <- numeric(num_iteration)
    slope_highCIs <- numeric(num_iteration)
    
    # Run multiple iterations
    for (iteration in 1:num_iteration) {
      # Set seed for this iteration based on the sample size and iteration
      set.seed(iteration + sample_size)
      
      # Run the iteration and extract the coefficients
      coefficients_iteration <- run_iteration(df_full_aboveground, sample_size)
      
      # Store the results
      slopes[iteration] <- coefficients_iteration["slope"]
      y_intercepts[iteration] <- coefficients_iteration["y_intercept"]
      slope_lowCIs[iteration] <- coefficients_iteration["slope_lowCI"]
      slope_highCIs[iteration] <- coefficients_iteration["slope_highCI"]
    }
    
    # Calculate mean and standard error for each coefficient
    mean_slope <- mean(slopes)
    se_slope <- sd(slopes) / sqrt(num_iteration)
    
    mean_y_intercept <- mean(y_intercepts)
    se_y_intercept <- sd(y_intercepts) / sqrt(num_iteration)
    
    mean_slope_lowCI <- mean(slope_lowCIs)
    se_slope_lowCI <- sd(slope_lowCIs) / sqrt(num_iteration)
    
    mean_slope_highCI <- mean(slope_highCIs)
    se_slope_highCI <- sd(slope_highCIs) / sqrt(num_iteration)
    
    cv_slope<- sd(slopes) / mean(slopes)
    cv_slope_lowCI<- sd(slope_lowCIs) / mean(slope_lowCIs)
    cv_slope_highCI<- sd(slope_highCIs) / mean(slope_highCIs)
    
    
    
    # Store the results in the dataframe
    bootstrapped_results_df <- rbind(bootstrapped_results_df, c(iteration, mean_slope, mean_y_intercept, mean_slope_lowCI, mean_slope_highCI, se_slope, se_y_intercept, se_slope_lowCI, se_slope_highCI, cv_slope, cv_slope_lowCI, cv_slope_highCI, sample_size))
  }
}

## Rename the columns in bootstrapped_results_df
colnames(bootstrapped_results_df) <- c("iteration", "mean_slope", "mean_y_intercept", "mean_slope_lowCI", "mean_slope_highCI", "se_slope", "se_y_intercept", "se_slope_lowCI", "se_slope_highCI","cv_slope", "cv_slope_lowCI", "cv_slope_highCI",  "sample_size")

##Store the results

bootstrap_simulation_hydr_mean_summary<- bootstrapped_results_df

##Bootstrapping data at different sample sizes and fitting a model without averaging----

## Function to perform random sampling, fit SMA model, and extract coefficients
run_iteration <- function(df_full_aboveground, sample_size) {
  # Random sampling of D across each and every L
  bootstrapped_data <- df_full_aboveground %>%
    group_by(L) %>%
    slice_sample(n = sample_size, replace = TRUE) %>%
    ungroup()
  
  # Fit a SMA model to bootstrapped_data
  bootstrapped_scaling <- sma(data = bootstrapped_data, DAVG ~ L, log = "XY", method = "SMA")
  
  # Extract coefficients
  bootstrapped_scaling_summary <- bootstrapped_scaling$groupsummary
  bootstrapped_results <- c(
    slope = bootstrapped_scaling_summary$Slope,
    y_intercept = bootstrapped_scaling_summary$Int,
    slope_lowCI = bootstrapped_scaling_summary$Slope_lowCI,
    slope_highCI = bootstrapped_scaling_summary$Slope_highCI
  )
  
  # Return a list of coefficients
  return(bootstrapped_results)
}

## Specify the sample sizes and number of iterations
sample_sizes <- c(1, 5, 10, 25, 50, 75, 100,1000)
num_iterations <- 1000  # You can adjust the number of iterations

## Create an empty dataframe to store the results
bootstrapped_results_df <- data.frame(
  iteration = integer(),
  mean_slope = double(),
  mean_y_intercept = double(),
  mean_slope_lowCI = double(),
  mean_slope_highCI = double(),
  se_slope = double(),
  se_y_intercept = double(),
  se_slope_lowCI = double(),
  se_slope_highCI = double(),
  sample_size = integer()
)

## Loop through different sample sizes
for (sample_size in sample_sizes) {
  for (num_iteration in num_iterations) {
    slopes <- numeric(num_iteration)
    y_intercepts <- numeric(num_iteration)
    slope_lowCIs <- numeric(num_iteration)
    slope_highCIs <- numeric(num_iteration)
    
    # Run multiple iterations
    for (iteration in 1:num_iteration) {
      # Set seed for this iteration based on the sample size and iteration
      set.seed(iteration + sample_size)
      
      # Run the iteration and extract the coefficients
      coefficients_iteration <- run_iteration(df_full_aboveground, sample_size)
      
      # Store the results
      slopes[iteration] <- coefficients_iteration["slope"]
      y_intercepts[iteration] <- coefficients_iteration["y_intercept"]
      slope_lowCIs[iteration] <- coefficients_iteration["slope_lowCI"]
      slope_highCIs[iteration] <- coefficients_iteration["slope_highCI"]
    }
    
    # Calculate mean and standard error for each coefficient
    mean_slope <- mean(slopes)
    se_slope <- sd(slopes) / sqrt(num_iteration)
    
    mean_y_intercept <- mean(y_intercepts)
    se_y_intercept <- sd(y_intercepts) / sqrt(num_iteration)
    
    mean_slope_lowCI <- mean(slope_lowCIs)
    se_slope_lowCI <- sd(slope_lowCIs) / sqrt(num_iteration)
    
    mean_slope_highCI <- mean(slope_highCIs)
    se_slope_highCI <- sd(slope_highCIs) / sqrt(num_iteration)
    
    cv_slope<- sd(slopes) / mean(slopes)
    cv_slope_lowCI<- sd(slope_lowCIs) / mean(slope_lowCIs)
    cv_slope_highCI<- sd(slope_highCIs) / mean(slope_highCIs)
    
    # Store the results in the dataframe
    bootstrapped_results_df <- rbind(bootstrapped_results_df, c(iteration, mean_slope, mean_y_intercept, mean_slope_lowCI, mean_slope_highCI, se_slope, se_y_intercept, se_slope_lowCI, se_slope_highCI, cv_slope, cv_slope_lowCI, cv_slope_highCI, sample_size))
  }
}

## Rename the columns in bootstrapped_results_df
colnames(bootstrapped_results_df) <- c("iteration", "mean_slope", "mean_y_intercept", "mean_slope_lowCI", "mean_slope_highCI", "se_slope", "se_y_intercept", "se_slope_lowCI", "se_slope_highCI", "cv_slope", "cv_slope_lowCI", "cv_slope_highCI", "sample_size")

## Store the summary in separate dataframe

bootstrap_simulation_raw_summary<- bootstrapped_results_df

###Plot the results----

##Make fake plot to extract legend from

Model_data<- c("Mean conduit diameter (arithmetic)", "Mean conduit diameter (hydraulic)", "Raw data")
Value<- c(1,2,3)

df_data<- data.frame(Model_data, Value)

df_data$Model_data<- factor(df_data$Model_data, levels = c("Mean conduit diameter (arithmetic)", "Mean conduit diameter (hydraulic)", "Raw data"))

data_grp<- levels(df_data$Model_data)
data_color<- c("red", "yellow2", "blue3")
names(data_color) <- data_grp

data_fake_plot<- ggplot(data=df_data, aes(x=Model_data, y=Value))+
  geom_bar(aes(fill= Model_data), stat = "identity")+
  scale_fill_manual(values = data_color)+
  labs(fill = "Data preprocessing")+
  theme(legend.position = "right",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

data_fake_plot

data_fake_plot_legend<- cowplot::get_legend(data_fake_plot)

##Build actual plot

bootstrap_simulation_V2_plot<- ggplot()+
  geom_path(data=bootstrap_simulation_arith_mean_summary, aes(y=mean_slope, x=sample_size), color="red", linewidth=1.5, linetype = "solid")+
  geom_ribbon(data=bootstrap_simulation_arith_mean_summary, aes(x=sample_size, ymin=mean_slope_lowCI, ymax=mean_slope), fill="red", alpha=0.15)+
  geom_ribbon(data=bootstrap_simulation_arith_mean_summary, aes(x=sample_size, ymin=mean_slope, ymax=mean_slope_highCI), fill="red", alpha=0.15)+
  geom_path(data=bootstrap_simulation_hydr_mean_summary, aes(y=mean_slope, x=sample_size), color="yellow2", linewidth=1.5, linetype = "solid")+
  geom_ribbon(data=bootstrap_simulation_hydr_mean_summary, aes(x=sample_size, ymin=mean_slope_lowCI, ymax=mean_slope), fill="yellow2", alpha=0.15)+
  geom_ribbon(data=bootstrap_simulation_hydr_mean_summary, aes(x=sample_size, ymin=mean_slope, ymax=mean_slope_highCI), fill="yellow2", alpha=0.15)+
  geom_path(data=bootstrap_simulation_raw_summary, aes(y=mean_slope, x=sample_size), color="blue3", linewidth=1.5, linetype = "solid")+
  geom_ribbon(data=bootstrap_simulation_raw_summary, aes(x=sample_size, ymin=mean_slope_lowCI, ymax=mean_slope), fill="blue3", alpha=0.15)+
  geom_ribbon(data=bootstrap_simulation_raw_summary, aes(x=sample_size, ymin=mean_slope, ymax=mean_slope_highCI), fill="blue3", alpha=0.15)+
  xlab("Number of conduits measured per sample (dimensionless)")+
  ylab(expression(paste("Fitted scaling exponent, ", italic("b"), " (dimensionless)")))+
  scale_x_log10(breaks=c(1, 5, 10, 25, 50, 100,1000))+
  theme_classic()+
  theme(aspect.ratio = 1,
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))

bootstrap_simulation_V2_plot

bootstrap_simulation_V2_plot_combined<- ggdraw() +
  draw_plot(bootstrap_simulation_V2_plot, x = 0, y = 0, width = 1, height = 1)+
  draw_grob(data_fake_plot_legend, x = 0.6, y = 0.85, width = 0.1, height = 0.1)

bootstrap_simulation_V2_plot_combined
    
#Save plot!

ggsave("bootstrap_simulation_plot_1.JPG", plot = bootstrap_simulation_V2_plot_combined, width=6, height = 6, dpi = 600)
ggsave("bootstrap_simulation_plot_1.svg", plot = bootstrap_simulation_V2_plot_combined, width=6, height = 6, dpi = 600)

#Simulating the effects of regression model choices----

##Bootstrapping data at sample size of 1000 at each L, averaging the data, and fitting an OLS model----

## Function to perform random sampling, fit SMA model, and extract coefficients

run_iteration <- function(df_full_aboveground, sample_size) {
  # Random sampling of D across each and every L
  bootstrapped_data <- df_full_aboveground %>%
    group_by(L) %>%
    slice_sample(n = sample_size, replace = TRUE) %>%
    ungroup()
  
  # Calculate the average DAVG for each z
  averaged_data <- bootstrapped_data %>%
    group_by(L) %>%
    summarize(DAVG = mean(DAVG))
  
  # Fit a SMA model to averaged_data
  bootstrapped_scaling <- sma(data = averaged_data, DAVG ~ L, log = "XY", method = "OLS")
  
  # Extract coefficients
  bootstrapped_scaling_summary <- bootstrapped_scaling$groupsummary
  bootstrapped_results <- c(
    slope = bootstrapped_scaling_summary$Slope,
    y_intercept = bootstrapped_scaling_summary$Int,
    slope_lowCI = bootstrapped_scaling_summary$Slope_lowCI,
    slope_highCI = bootstrapped_scaling_summary$Slope_highCI
  )
  
  # Return a list of coefficients
  return(bootstrapped_results)
}

## Specify the sample sizes and number of iterations
sample_sizes <- 1000
num_iterations <- 1000  # You can adjust the number of iterations

## Create an empty dataframe to store the results
bootstrapped_results_df <- data.frame(
  iteration = integer(),
  mean_slope = double(),
  mean_y_intercept = double(),
  mean_slope_lowCI = double(),
  mean_slope_highCI = double(),
  se_slope = double(),
  se_y_intercept = double(),
  se_slope_lowCI = double(),
  se_slope_highCI = double(),
  cv_slope= double(),
  cv_slope_lowCI = double(),
  cv_slope_highCI = double(),
  sample_size = integer()
)

## Loop through different sample sizes
for (sample_size in sample_sizes) {
  for (num_iteration in num_iterations) {
    slopes <- numeric(num_iteration)
    y_intercepts <- numeric(num_iteration)
    slope_lowCIs <- numeric(num_iteration)
    slope_highCIs <- numeric(num_iteration)
    
    # Run multiple iterations
    for (iteration in 1:num_iteration) {
      # Set seed for this iteration based on the sample size and iteration
      set.seed(iteration + sample_size)
      
      # Run the iteration and extract the coefficients
      coefficients_iteration <- run_iteration(df_full_aboveground, sample_size)
      
      # Store the results
      slopes[iteration] <- coefficients_iteration["slope"]
      y_intercepts[iteration] <- coefficients_iteration["y_intercept"]
      slope_lowCIs[iteration] <- coefficients_iteration["slope_lowCI"]
      slope_highCIs[iteration] <- coefficients_iteration["slope_highCI"]
    }
    
    # Calculate mean and standard error for each coefficient
    mean_slope <- mean(slopes)
    se_slope <- sd(slopes) / sqrt(num_iteration)
    
    mean_y_intercept <- mean(y_intercepts)
    se_y_intercept <- sd(y_intercepts) / sqrt(num_iteration)
    
    mean_slope_lowCI <- mean(slope_lowCIs)
    se_slope_lowCI <- sd(slope_lowCIs) / sqrt(num_iteration)
    
    mean_slope_highCI <- mean(slope_highCIs)
    se_slope_highCI <- sd(slope_highCIs) / sqrt(num_iteration)
    
    cv_slope<- sd(slopes) / mean(slopes)
    cv_slope_lowCI<- sd(slope_lowCIs) / mean(slope_lowCIs)
    cv_slope_highCI<- sd(slope_highCIs) / mean(slope_highCIs)
    
    
    
    # Store the results in the dataframe
    bootstrapped_results_df <- rbind(bootstrapped_results_df, c(iteration, mean_slope, mean_y_intercept, mean_slope_lowCI, mean_slope_highCI, se_slope, se_y_intercept, se_slope_lowCI, se_slope_highCI, cv_slope, cv_slope_lowCI, cv_slope_highCI, sample_size))
  }
}

## Rename the columns in bootstrapped_results_df
colnames(bootstrapped_results_df) <- c("iteration", "mean_slope", "mean_y_intercept", "mean_slope_lowCI", "mean_slope_highCI", "se_slope", "se_y_intercept", "se_slope_lowCI", "se_slope_highCI","cv_slope", "cv_slope_lowCI", "cv_slope_highCI",  "sample_size")

##Store the results

bootstrap_simulation_arith_mean_summary_OLS<- bootstrapped_results_df

##Bootstrapping data at a sample size of 1000 at each L, averaging the data and weighing by conductance (Dh), and fitting a model----

## Function to perform random sampling, fit SMA model, and extract coefficients

run_iteration <- function(df_full_aboveground, sample_size) {
  # Random sampling of D across each and every L
  bootstrapped_data <- df_full_aboveground %>%
    group_by(L) %>%
    slice_sample(n = sample_size, replace = TRUE) %>%
    ungroup()
  
  # Calculate the average DAVG for each z
  averaged_data <- bootstrapped_data %>%
    group_by(L) %>%
    summarize(DhAVG = sum(DAVG^5)/sum(DAVG^4))
  
  # Fit a SMA model to averaged_data
  bootstrapped_scaling <- sma(data = averaged_data, DhAVG ~ L, log = "XY", method = "OLS")
  
  # Extract coefficients
  bootstrapped_scaling_summary <- bootstrapped_scaling$groupsummary
  bootstrapped_results <- c(
    slope = bootstrapped_scaling_summary$Slope,
    y_intercept = bootstrapped_scaling_summary$Int,
    slope_lowCI = bootstrapped_scaling_summary$Slope_lowCI,
    slope_highCI = bootstrapped_scaling_summary$Slope_highCI
  )
  
  # Return a list of coefficients
  return(bootstrapped_results)
}

## Specify the sample sizes and number of iterations
sample_sizes <- 1000
num_iterations <- 1000  # You can adjust the number of iterations

## Create an empty dataframe to store the results
bootstrapped_results_df <- data.frame(
  iteration = integer(),
  mean_slope = double(),
  mean_y_intercept = double(),
  mean_slope_lowCI = double(),
  mean_slope_highCI = double(),
  se_slope = double(),
  se_y_intercept = double(),
  se_slope_lowCI = double(),
  se_slope_highCI = double(),
  cv_slope= double(),
  cv_slope_lowCI = double(),
  cv_slope_highCI = double(),
  sample_size = integer()
)

## Loop through different sample sizes
for (sample_size in sample_sizes) {
  for (num_iteration in num_iterations) {
    slopes <- numeric(num_iteration)
    y_intercepts <- numeric(num_iteration)
    slope_lowCIs <- numeric(num_iteration)
    slope_highCIs <- numeric(num_iteration)
    
    # Run multiple iterations
    for (iteration in 1:num_iteration) {
      # Set seed for this iteration based on the sample size and iteration
      set.seed(iteration + sample_size)
      
      # Run the iteration and extract the coefficients
      coefficients_iteration <- run_iteration(df_full_aboveground, sample_size)
      
      # Store the results
      slopes[iteration] <- coefficients_iteration["slope"]
      y_intercepts[iteration] <- coefficients_iteration["y_intercept"]
      slope_lowCIs[iteration] <- coefficients_iteration["slope_lowCI"]
      slope_highCIs[iteration] <- coefficients_iteration["slope_highCI"]
    }
    
    # Calculate mean and standard error for each coefficient
    mean_slope <- mean(slopes)
    se_slope <- sd(slopes) / sqrt(num_iteration)
    
    mean_y_intercept <- mean(y_intercepts)
    se_y_intercept <- sd(y_intercepts) / sqrt(num_iteration)
    
    mean_slope_lowCI <- mean(slope_lowCIs)
    se_slope_lowCI <- sd(slope_lowCIs) / sqrt(num_iteration)
    
    mean_slope_highCI <- mean(slope_highCIs)
    se_slope_highCI <- sd(slope_highCIs) / sqrt(num_iteration)
    
    cv_slope<- sd(slopes) / mean(slopes)
    cv_slope_lowCI<- sd(slope_lowCIs) / mean(slope_lowCIs)
    cv_slope_highCI<- sd(slope_highCIs) / mean(slope_highCIs)
    
    
    
    # Store the results in the dataframe
    bootstrapped_results_df <- rbind(bootstrapped_results_df, c(iteration, mean_slope, mean_y_intercept, mean_slope_lowCI, mean_slope_highCI, se_slope, se_y_intercept, se_slope_lowCI, se_slope_highCI, cv_slope, cv_slope_lowCI, cv_slope_highCI, sample_size))
  }
}

## Rename the columns in bootstrapped_results_df
colnames(bootstrapped_results_df) <- c("iteration", "mean_slope", "mean_y_intercept", "mean_slope_lowCI", "mean_slope_highCI", "se_slope", "se_y_intercept", "se_slope_lowCI", "se_slope_highCI","cv_slope", "cv_slope_lowCI", "cv_slope_highCI",  "sample_size")

##Store the results

bootstrap_simulation_hydr_mean_summary_OLS<- bootstrapped_results_df

##Bootstrapping data at a sample size of 1000 and fitting OLS model----

## Function to perform random sampling, fit SMA model, and extract coefficients
run_iteration <- function(df_full_aboveground_clean, sample_size) {
  # Random sampling of D across each and every L
  bootstrapped_data <- df_full_aboveground_clean %>%
    group_by(z) %>%
    slice_sample(n = sample_size, replace = TRUE) %>%
    ungroup()
  
  # Fit a SMA model to bootstrapped_data
  bootstrapped_scaling <- sma(data = bootstrapped_data, DAVG ~ z, log = "XY", method = "OLS")
  
  # Extract coefficients
  bootstrapped_scaling_summary <- bootstrapped_scaling$groupsummary
  bootstrapped_results <- c(
    slope = bootstrapped_scaling_summary$Slope,
    y_intercept = bootstrapped_scaling_summary$Int,
    slope_lowCI = bootstrapped_scaling_summary$Slope_lowCI,
    slope_highCI = bootstrapped_scaling_summary$Slope_highCI
  )
  
  # Return a list of coefficients
  return(bootstrapped_results)
}

## Specify the sample sizes and number of iterations
sample_sizes <- 1000
num_iterations <- 1000  # You can adjust the number of iterations

## Create an empty dataframe to store the results
bootstrapped_results_df <- data.frame(
  iteration = integer(),
  mean_slope = double(),
  mean_y_intercept = double(),
  mean_slope_lowCI = double(),
  mean_slope_highCI = double(),
  se_slope = double(),
  se_y_intercept = double(),
  se_slope_lowCI = double(),
  se_slope_highCI = double(),
  sample_size = integer()
)

## Loop through different sample sizes
for (sample_size in sample_sizes) {
  for (num_iteration in num_iterations) {
    slopes <- numeric(num_iteration)
    y_intercepts <- numeric(num_iteration)
    slope_lowCIs <- numeric(num_iteration)
    slope_highCIs <- numeric(num_iteration)
    
    # Run multiple iterations
    for (iteration in 1:num_iteration) {
      # Set seed for this iteration based on the sample size and iteration
      set.seed(iteration + sample_size)
      
      # Run the iteration and extract the coefficients
      coefficients_iteration <- run_iteration(df_full_aboveground_clean, sample_size)
      
      # Store the results
      slopes[iteration] <- coefficients_iteration["slope"]
      y_intercepts[iteration] <- coefficients_iteration["y_intercept"]
      slope_lowCIs[iteration] <- coefficients_iteration["slope_lowCI"]
      slope_highCIs[iteration] <- coefficients_iteration["slope_highCI"]
    }
    
    # Calculate mean and standard error for each coefficient
    mean_slope <- mean(slopes)
    se_slope <- sd(slopes) / sqrt(num_iteration)
    
    mean_y_intercept <- mean(y_intercepts)
    se_y_intercept <- sd(y_intercepts) / sqrt(num_iteration)
    
    mean_slope_lowCI <- mean(slope_lowCIs)
    se_slope_lowCI <- sd(slope_lowCIs) / sqrt(num_iteration)
    
    mean_slope_highCI <- mean(slope_highCIs)
    se_slope_highCI <- sd(slope_highCIs) / sqrt(num_iteration)
    
    cv_slope<- sd(slopes) / mean(slopes)
    cv_slope_lowCI<- sd(slope_lowCIs) / mean(slope_lowCIs)
    cv_slope_highCI<- sd(slope_highCIs) / mean(slope_highCIs)
    
    # Store the results in the dataframe
    bootstrapped_results_df <- rbind(bootstrapped_results_df, c(iteration, mean_slope, mean_y_intercept, mean_slope_lowCI, mean_slope_highCI, se_slope, se_y_intercept, se_slope_lowCI, se_slope_highCI, cv_slope, cv_slope_lowCI, cv_slope_highCI, sample_size))
  }
}

## Rename the columns in bootstrapped_results_df
colnames(bootstrapped_results_df) <- c("iteration", "mean_slope", "mean_y_intercept", "mean_slope_lowCI", "mean_slope_highCI", "se_slope", "se_y_intercept", "se_slope_lowCI", "se_slope_highCI", "cv_slope", "cv_slope_lowCI", "cv_slope_highCI", "sample_size")

## Store the summary in separate dataframe

bootstrap_simulation_raw_summary_OLS<- bootstrapped_results_df

##Compile all data into 3 separate dataframes and clean them for plot----

###Arithmetic mean----

bootstrap_simulation_arith_mean_summary_SMA<- bootstrap_simulation_arith_mean_summary[8, ]

bootstrap_simulation_arith_mean_summary_SMA<- bootstrap_simulation_arith_mean_summary_SMA %>%
  mutate("Model" = "SMA")

bootstrap_simulation_arith_mean_summary_OLS<- bootstrap_simulation_arith_mean_summary_OLS %>%
  mutate("Model" = "OLS")

bootstrap_simulation_arith_mean_summary_OLS_SMA<- rbind(bootstrap_simulation_arith_mean_summary_SMA, bootstrap_simulation_arith_mean_summary_OLS)

###Hydraulic mean----

bootstrap_simulation_hydr_mean_summary_SMA<- bootstrap_simulation_hydr_mean_summary[8, ]

bootstrap_simulation_hydr_mean_summary_SMA<- bootstrap_simulation_hydr_mean_summary_SMA %>%
  mutate("Model" = "SMA")

bootstrap_simulation_hydr_mean_summary_OLS<- bootstrap_simulation_hydr_mean_summary_OLS %>%
  mutate("Model" = "OLS")

bootstrap_simulation_hydr_mean_summary_OLS_SMA<- rbind(bootstrap_simulation_hydr_mean_summary_SMA, bootstrap_simulation_hydr_mean_summary_OLS)

###Raw data----

bootstrap_simulation_raw_summary_SMA<- bootstrap_simulation_raw_summary[8, ]

bootstrap_simulation_raw_summary_SMA<- bootstrap_simulation_raw_summary_SMA %>%
  mutate("Model" = "SMA")

bootstrap_simulation_raw_summary_OLS<- bootstrap_simulation_raw_summary_OLS %>%
  mutate("Model" = "OLS")

bootstrap_simulation_raw_summary_OLS_SMA<- rbind(bootstrap_simulation_raw_summary_SMA, bootstrap_simulation_raw_summary_OLS)

##Plot the results----

##Build actual plot

bootstrap_simulation_model_plot<- ggplot()+
  geom_pointrange(data=bootstrap_simulation_arith_mean_summary_OLS_SMA, aes(y=mean_slope, x=Model, ymin = mean_slope_lowCI, ymax = mean_slope_highCI), stat = "identity", color="red", linewidth=1.5, size=1, position_nudge(x=-0.1)) +
  geom_line(data=bootstrap_simulation_arith_mean_summary_OLS_SMA, aes(y=mean_slope, x=Model, group=1), color="red", alpha=0.2,linewidth=1.5, position = position_nudge(x=-0.1))+
  geom_pointrange(data=bootstrap_simulation_hydr_mean_summary_OLS_SMA, aes(y=mean_slope, x=Model, ymin = mean_slope_lowCI, ymax = mean_slope_highCI), stat = "identity", color="yellow2",linewidth=1.5, size=1, position_nudge(x=0)) +
  geom_line(data=bootstrap_simulation_hydr_mean_summary_OLS_SMA, aes(y=mean_slope, x=Model, group=1), color="yellow2", alpha=0.2,linewidth=1.5, position = position_nudge(x=0))+
  geom_pointrange(data=bootstrap_simulation_raw_summary_OLS_SMA, aes(y=mean_slope, x=Model, ymin = mean_slope_lowCI, ymax = mean_slope_highCI), stat = "identity", color="blue3",linewidth=1.5,size=1, position_nudge(x=0.1)) +
  geom_line(data=bootstrap_simulation_raw_summary_OLS_SMA, aes(y=mean_slope, x=Model, group=1), color="blue3", alpha=0.2,linewidth=1.5, position = position_nudge(x=0.1))+
  xlab("Model type")+
  ylab(expression(paste("Fitted scaling exponent, ", italic("b"), " (dimensionless)")))+
  theme_classic()+
  theme(aspect.ratio = 1,
        axis.text.x = element_text(size=12, hjust = 0.5),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))+
  scale_x_discrete(expand = c(0.2, 0.2))

bootstrap_simulation_model_plot

bootstrap_simulation_model_plot_combined<- ggdraw() +
  draw_plot(bootstrap_simulation_model_plot, x = 0, y = 0, width = 1, height = 1)+
  draw_grob(data_fake_plot_legend, x = 0.37, y = 0.80, width = 0.1, height = 0.1)

bootstrap_simulation_model_plot_combined

#Save plot!

ggsave("bootstrap_simulation_plot_2.JPG", plot = bootstrap_simulation_model_plot_combined, width=6, height = 6, dpi = 600)
ggsave("bootstrap_simulation_plot_2.svg", plot = bootstrap_simulation_model_plot_combined, width=6, height = 6, dpi = 600)
