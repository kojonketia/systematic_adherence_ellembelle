suppressPackageStartupMessages({
  library(Rcpp)
  library(parallel)
  library(dplyr)
  library(ggplot2)
})

start.time <- Sys.time()

## MDA offset
convert_year_to_time <- function(year, start_year = 2001) {
  # convert year to time in months since the start of the simulation
  return((year - start_year) * 12)
}

sourceCpp(file.path("code", "Transfil.cpp")) # the TRANSFIL model

posterior_df <- read.csv(file.path("output", "posterior_samples.csv"))
posterior_k <- posterior_df$k
posterior_vth <- posterior_df$vth

n_runs <- 1000 # number of simulation runs

# intervention data
int_df <- read.csv(file.path("data", "interventions.csv")) %>%
  mutate(time = convert_year_to_time(year))

set.seed(123)  # for reproducibility
posterior_ks <- sample(posterior_k, n_runs, replace = FALSE)
posterior_vths <- sample(posterior_vth, n_runs, replace = FALSE)

output_path <- file.path("output", "scenarios")
if(!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
}

simulate_scenario <- function(k, vth, sys_compliance, int_df) {
  mda_settings <- list("time" = int_df$time + 36, ## lf prevalence run for 3 years
                       "coverage" = int_df$MDA_coverage )
  
  bednet_settings <- list("time" = int_df$time + 36,
                          "coverage" = int_df$bednet_coverage)
  
  x <- new(Data, 10000)
  x$set_parameter("k", k)
  x$set_parameter("vth", vth)  
  x$set_parameter("systematic_compliance", sys_compliance)
  x$set_parameter("importation_rate", 0.0003)
  x$burn_in()
  
  # set MDA and bednet coverage
  x$set_MDA(mda_settings)
  x$set_bednet(bednet_settings)
  
  # simulate through to 2040
  x$timestep(convert_year_to_time(2040) + 36) 
  
  out <- data.frame(sim_time = 0:(length(x$mf_prev)-1), # 2001 to 2040 in months
                    mf_prev = x$mf_prev,
                    ag_prev = x$worm_prev)
  
  return(out)
}

# assumed interventions from 2026 - 2040
forward_interventions <- data.frame(year = 2026:2040, 
                                    MDA_coverage = c(rep(0,15)),
                                    bednet_coverage = rep(0.64, 15)) %>% 
  mutate(time = convert_year_to_time(year))

all_int_df <- int_df %>%
  bind_rows(forward_interventions)


# modelling with observed (0.11) compliance rate -------------------------------
ag_prevalence_sim_data <- lapply(1:n_runs, function(i) {
  out <- simulate_scenario(k = posterior_ks[i], vth = posterior_vths[i],
                           sys_compliance = 0.11, all_int_df)
  out$sim_id <- i
  return(out)}) %>% bind_rows()

write.csv(ag_prevalence_sim_data, file.path(output_path, "ag_prevalence_sim_data.csv"), row.names = FALSE)

end.time <- Sys.time()
cat("Simulation completed in:", 
    round(difftime(end.time, start.time, units = "hours"), 2), "hours\n")

# plots ------------------------------------------------------------------------
prev_df <- read.csv(file.path("data", "prevalence.csv"))
sim_data <- read.csv(file.path("output", "scenarios", "ag_prevalence_sim_data.csv"))

summary_data <- sim_data %>% group_by(sim_time) %>%
  summarise(median = median(ag_prev),
            lower = quantile(ag_prev, 0.025), # 0.025
            upper = quantile(ag_prev, 0.975))

years <- c("1998", rep("", 2), "2001", rep("", 3), "2005", rep("", 4),
           "2010", rep("", 4), "2015", rep("", 4), "2020", rep("", 4), 
           "2025", rep("", 4), "2030", rep("", 4), "2035", rep("", 4), "2040")

ag_prev_plot <- ggplot(summary_data, aes(x = sim_time)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, color = "95% CI"), 
              fill = "lightblue", alpha = 0.8) +
  geom_line(aes(y = median, color = "Median"),
            linewidth = 0.5) +
  geom_hline(aes(yintercept = c(2), color = "2% threshold"), 
             lty="dashed", linewidth = 0.5) +
  geom_point(data = prev_df, aes(x = year2, y = worm_prevalence, 
                                 color = "Observed data"),
             size = 2, shape = 17) +
  labs(x = "Year", y = "Prevalence (%)") +
  theme_classic(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        #axis.title.x = element_text(face = "bold"), 
        #axis.title.y = element_text(face = "bold"),
        axis.line = element_line(linewidth = 0.2),
        axis.ticks.length = unit(0.2, "cm"),
        legend.position = c(0.8,0.9) ) +
  scale_x_continuous(breaks = seq(1, 505, 12), labels = years, expand = c(0,0)) +
  scale_y_continuous(breaks = seq(0,30,5), labels = seq(0,30,5), expand = c(0,0), 
                     limits = c(0,30)) +
  scale_color_manual(name = "", 
                     breaks = c("Median","95% CI","2% threshold", "Observed data"),
                     values = c("Median" = "blue",
                                "95% CI" = 'lightblue',
                                "2% threshold" = "red",
                                "Observed data" = "red"))
ag_prev_plot

output_file <- file.path("output", "figures", "ag_prevalence_plot.png")
ggsave(filename = output_file,  plot = ag_prev_plot, width = 6, height = 4, 
       bg = "white")
