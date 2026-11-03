library(foreach)
library(doParallel)
library(parallel)
library(Rcpp)
library(writexl)
library(readxl)
library(ggplot2)
library(dplyr)
library(ggthemes)

set.seed(123)
start.time <- Sys.time()


cpp_file <- normalizePath("code/Transfil.cpp", mustWork = TRUE)

output_dir <- "output"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# progess is saved as CSV during simulation for faster checkpointing.
# the final completed dataset is saved as Excel.
progress_file <- file.path(output_dir, "scenarios_data_progress.csv")
final_file <- file.path(output_dir, "scenarios_data.xlsx")

all_results <- list()

# scenario parameters
adherence <- c(0.11, 0.5, 0.99)
mda_rounds <- 1:10
mda_coverages <- c(0.7, 0.8, 0.9, 1.0)

n_simulations <- 1000
pop_size <- 1000


# set up parallel workers
available_cores <- parallel::detectCores(logical = TRUE)

if (is.na(available_cores)) {
  available_cores <- 2L
}

n_cores <- max(1L, available_cores - 1L)
cat("Available logical cores:", available_cores, "\n")
cat("Parallel workers:", n_cores, "\n")

cl <- parallel::makeCluster(n_cores)

# always stop the workers when the function exits.
run_simulations <- function(cl, cpp_file, output_dir, progress_file, final_file,
                            adherence, mda_rounds, mda_coverages, n_simulations,
                            pop_size) {
  
  on.exit({try(parallel::stopCluster(cl), silent = TRUE)}, add = TRUE)
  
  doParallel::registerDoParallel(cl)
  
  cat("Registered parallel workers:", foreach::getDoParWorkers(), "\n")
  

  parallel::clusterSetRNGStream(cl, iseed = 123)
  parallel::clusterExport(cl, varlist = "cpp_file", envir = environment())
  
  parallel::clusterEvalQ(cl, {library(Rcpp)
    Rcpp::sourceCpp(cpp_file, rebuild = FALSE)
    NULL})

  # run all simulations
  completed_scenarios <- 0L
  
  total_scenarios <- length(adherence) * length(mda_coverages) * length(mda_rounds)
  
  for (sys_compliance in adherence) {
    for (cov in mda_coverages) {
      for (j in seq_along(mda_rounds)) {
        
        n_rounds <- mda_rounds[j]
        
        # MDA and bednet coverage schedules
        mda_coverages_vec <- rep(cov, n_rounds)
        bednet_coverages <- rep(0.8, n_rounds)
        
        MDA_time <- seq(12, 12 * n_rounds, by = 12) + 100
        bednet_time <- MDA_time
        
        cat("\nScenario:", completed_scenarios + 1L, "of", total_scenarios,
            "| Compliance:", sys_compliance,
            "| Coverage:", cov,
            "| MDA rounds:", n_rounds, "\n")
        
        # run simulations in parallel
        sim_results <- foreach(i = seq_len(n_simulations),
                               .packages = "Rcpp",
                               .export = c("MDA_time", "bednet_time",
                                           "mda_coverages_vec", "bednet_coverages",
                                           "sys_compliance", "pop_size"),
                               .inorder = TRUE) %dopar% {
          
          # The Data class was compiled and loaded on each worker.
          bednet_settings <- list(time = bednet_time,
                                  coverage = bednet_coverages)
          
          MDA_settings <- list(time = MDA_time, coverage = mda_coverages_vec)
          
          x <- new(Data, pop_size)
          
          x$set_bednet(bednet_settings)
          x$set_MDA(MDA_settings)
          x$set_parameter("k", 0.008)
          x$set_parameter("vth", 35)
          x$set_parameter("systematic_compliance", sys_compliance)
          x$set_parameter("importation_rate", 0.0)
          
          x$burn_in()
          x$timestep(225)
          
          # return only worm prevalence, which is used below.
          x$worm_prev
        }
        
        # calculate elimination proportion
        antigenData <- as.data.frame(do.call(rbind, sim_results))
        
        eval_time <- tail(MDA_time, 1) + 2
        time_string <- sprintf("V%d", eval_time)
        
        if (!time_string %in% names(antigenData)) {
          stop(paste0("Time column ", time_string, 
                      " was not found in worm prevalence output. ",
                      "Check the structure and column names of x$worm_prev."))
        }
        
        elimination_prop <- mean(antigenData[[time_string]] < 2, na.rm = TRUE)
        
        # one result row for this completed scenario
        scenario_result <- data.frame(Rounds = n_rounds,
                                      Values = elimination_prop,
                                      Coverage = cov * 100,
                                      Adherence = sys_compliance)
        
        all_results[[length(all_results) + 1L]] <- scenario_result
        
      
        # save progress after every completed scenario
        current_results <- do.call(rbind, all_results)
        write.csv(current_results, progress_file, row.names = FALSE)
        
        completed_scenarios <- completed_scenarios + 1L
        
        cat("Elimination proportion:", round(elimination_prop, 4), "\n")
        cat("Progress:", round(100 * completed_scenarios / total_scenarios, 1),
            "% | Checkpoint saved.\n")
      }
    }
  }
  
  
  # save final data as excel
  final_data <- do.call(rbind, all_results)
  rownames(final_data) <- NULL
  
  writexl::write_xlsx(final_data, path = final_file)
  
  cat("\nAll simulations completed successfully.\n")
  cat("Total scenarios:", completed_scenarios, "\n")
  cat("Total simulations:", completed_scenarios * n_simulations, "\n")
  cat("Final Excel file:", normalizePath(final_file), "\n")
  cat("Progress file:", normalizePath(progress_file), "\n")
  
  return(final_data)
}


# exexute simulations 
final_data <- run_simulations(cl = cl,
                              cpp_file = cpp_file,
                              output_dir = output_dir,
                              progress_file = progress_file,
                              final_file = final_file,
                              adherence = adherence,
                              mda_rounds = mda_rounds,
                              mda_coverages = mda_coverages,
                              n_simulations = n_simulations,
                              pop_size = pop_size)

end.time <- Sys.time()

cat("\nTotal simulation time:", 
    round(as.numeric(difftime(end.time, start.time, units = "hours")), 2),
    "hours\n")


print(final_data)



# plotting simulation results --------------------------------------------------
final_data <- read_excel(file.path("output", "scenarios_data.xlsx"))


# Prepare data
plot_data <- final_data %>%
  mutate(Coverage = factor(Coverage, levels = c(70, 80, 90, 100), 
                           labels = c("70", "80", "90", "100")),
    Adherence = factor(Adherence, levels = c(0.11, 0.50, 0.99),
                        labels = c("ρ = 0.11", "ρ = 0.50", "ρ = 0.99")),
    Probability = Values * 100)

# scenarios plot
p <- ggplot(plot_data, aes(x = Rounds, y = Probability, colour = Coverage,
                           linetype = Adherence, 
                           group = interaction(Coverage, Adherence))) +
  geom_line(linewidth = 0.5) +
  geom_point(shape = 0, size = 1.5, stroke = 0.9) +
  scale_colour_manual(name = "Coverage:", 
                      values = c("70" = "#F12626", "80" = "#0879E8", 
                                 "90" = "#329640", "100" = "#8619B5")) +
  scale_linetype_manual(name = "Ahderence:",
                        values = c("ρ = 0.11" = "dashed",
                                   "ρ = 0.50" = "dotdash",
                                   "ρ = 0.99" = "solid")) +
  scale_x_continuous(breaks = 1:10, limits = c(1, 10)) +
  scale_y_continuous(breaks = seq(0, 100, by = 20),
                     minor_breaks = seq(0, 100, by = 10),
                     limits = c(0, 100)) +
  labs(x = "Rounds of mass drug administration (n)",
       y = "Probability of achieving <2% target (%)") +
  #theme_classic(base_size = 15) +
  theme_bw(base_size = 12) +
  theme(axis.title = element_text(size = 15, colour = "black"),
        axis.text = element_text(size = 12, colour = "black"),
        axis.line = element_line(colour = "black", linewidth = 0.6),
        axis.ticks = element_line(colour = "black", linewidth = 0.6),
        panel.grid = element_blank(),
        legend.position = c(0.05, 0.98),
        legend.justification = c(0, 1),
        legend.background = element_rect(fill = "white", colour = "black",
                                         linewidth = 0.5),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 11),
        legend.key.width = grid::unit(1.3, "cm"),
        legend.spacing.y = grid::unit(0.12, "cm")) 

print(p)

# Save high-resolution figure
ggsave("output/elimination_probability.png", plot = p, width = 7, height = 5,
       dpi = 600, bg = "white")

#ggsave("output/elimination_probability.pdf", plot = p, width = 8.5, height = 6)
