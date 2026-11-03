
rm(list = ls())

library(Rcpp)
library(parallel)

start.time <- Sys.time()

# settings ---------------------------------------------------------------------
out <- "output"
dir.create(file.path(out, "figures"), recursive = TRUE, showWarnings = FALSE)

N <- 1200
T <- 8
alpha <- 0.5
perturbation_scale <- 2
max_attempts <- 1000000

cores <- as.integer(Sys.getenv("N_WORKERS", "7"))

batch_size <- max(cores * 4L, 8L)

cat("Available cores:", detectCores(), "\n")
cat("Cores being used:", cores, "\n")
cat("Batch size:", batch_size, "\n\n")


# data and model settings ------------------------------------------------------
model_file <- normalizePath(file.path("code", "Transfil.cpp"), mustWork = TRUE)

prev <- read.csv(file.path("data", "prevalence.csv"))
s_obs <- as.numeric(prev$worm_prevalence[prev$year == 2001])

# parameter priors
lower <- c(k = 0.01, vth = 30)
upper <- c(k = 0.1, vth = 160)

p <- length(lower)

sample_prior <- function() {
  runif(p, min = lower, max = upper)
}

prior_density <- function(theta) {
  if (any(!is.finite(theta)) ||
      any(theta < lower) ||
      any(theta > upper)) {
    return(0)
  }
  
  prod(1 / (upper - lower))
}

# distance function ------------------------------------------------------------
distance_function <- function(simulated, observed) {
  simulated <- as.numeric(simulated)
  observed <- as.numeric(observed)
  
  if (length(simulated) != length(observed)) {
    stop("Simulator output and observed data have different lengths.")
  }
  
  if (any(!is.finite(simulated))) {
    return(Inf)
  }
  
  scale_obs <- if (length(observed) > 1) {
    sd(observed)
  } else {
    1
  }
  
  if (!is.finite(scale_obs) || scale_obs < 1e-8) {
    scale_obs <- 1
  }
  
  sqrt(mean(((simulated - observed) / scale_obs)^2))
}

# parallel TRANSFIL simulation -------------------------------------------------
# This function runs independently on each worker.
# TRANSFIL is compiled and initialized only once per worker.

parallel_simulate <- function(theta, seeds, cl) {
  tasks <- lapply(seq_len(nrow(theta)), function(i) {
    list(seed = seeds[i], k = theta[i, 1], vth = theta[i, 2])
  })
  
  parallel::parLapply(cl, tasks,
    function(task) {
      set.seed(task$seed)
      
      model <- new(Data, 10000)
      
      model$set_parameter("k", task$k)
      model$set_parameter("vth", task$vth)
      model$set_parameter("importation_rate", 0.0003)
      
      model$burn_in()
      
      as.numeric(model$worm_prev)
    }
  )
}


# importance weight functions --------------------------------------------------
log_sum_exp <- function(x) {
  m <- max(x)
  
  if (!is.finite(m)) {
    return(-Inf)
  }
  
  m + log(sum(exp(x - m)))
}

# Vectorized mixture proposal density.
# The covariance inverse and determinant are calculated only once
# per population rather than once per particle.

log_mixture_density <- function(theta, previous, w_prev, inv_cov, log_norm_const) {
  differences <- sweep(previous, 2, theta, "-")
  mahalanobis <- rowSums((differences %*% inv_cov) * differences)
  log_components <- log(w_prev) + log_norm_const - 0.5 * mahalanobis
  log_sum_exp(log_components)
}

# ABC-PMC algorithm ------------------------------------------------------------
abc_pmc <- function(N = 1200, T = 8, alpha = 0.5, perturbation_scale = 2,
                    max_attempts = 1000000, batch_size = 32, cl = NULL) {
  stopifnot(!is.null(cl))
  
  particles <- vector("list", T)
  weights <- vector("list", T)
  distances <- vector("list", T)
  epsilons <- numeric(T)
  
  nsim <- 0L
  
  # population 1: pilot simualtions
  cat("\nGenerating initial pilot simulations...\n")
  
  pilot_theta <- matrix(runif(N * p), nrow = N, ncol = p)
  pilot_theta <- sweep(pilot_theta, 2, upper - lower, "*")
  pilot_theta <- sweep(pilot_theta, 2, lower, "+")
  
  pilot_seeds <- sample.int(.Machine$integer.max, N)
  pilot_simulations <- parallel_simulate(pilot_theta, pilot_seeds, cl)
  pilot_distances <- vapply(pilot_simulations, distance_function, numeric(1),
                            observed = s_obs)
  
  nsim <- nsim + N
  epsilons[1] <- as.numeric(quantile(pilot_distances, probs = alpha,
                                     names = FALSE))
  
  cat("Initial tolerance:", round(epsilons[1], 6), "\n")
  

  # initial accepted population
  theta_pop <- matrix(NA_real_, nrow = N, ncol = p)
  d_pop <- numeric(N)
  
  i <- 1L
  attempts <- 0L
  
  while (i <= N && attempts < max_attempts) {
    n_batch <- min(batch_size, max_attempts - attempts)
    theta_batch <- matrix(runif(n_batch * p), nrow = n_batch, ncol = p)
    theta_batch <- sweep(theta_batch, 2, upper - lower, "*")
    theta_batch <- sweep(theta_batch, 2, lower, "+")
    
    seeds <- sample.int(.Machine$integer.max, n_batch)
    
    simulations <- parallel_simulate(theta_batch, seeds, cl)
    
    batch_distances <- vapply(simulations, distance_function, numeric(1),
                              observed = s_obs)
    
    attempts <- attempts + n_batch
    nsim <- nsim + n_batch
    
    # Preserve proposal order when accepting particles.
    
    for (j in seq_len(n_batch)) {
      if (batch_distances[j] <= epsilons[1]) {
        
        theta_pop[i, ] <- theta_batch[j, ]
        d_pop[i] <- batch_distances[j]
        
        i <- i + 1L
        
        if (i > N) {
          break
        }
      }
    }
  }
  
  if (i <= N) {
    stop("Could not fill initial population. Increase max_attempts.")
  }
  
  colnames(theta_pop) <- names(lower)
  
  particles[[1]] <- theta_pop
  weights[[1]] <- rep(1 / N, N)
  distances[[1]] <- d_pop
  
  cat("Population 1 completed |", "Tolerance:", round(epsilons[1], 6),
      "| Simulations:", nsim, "\n")
  
  # populations 2 to T
  if (T >= 2) {
    for (t in 2:T) {
      
      cat("\nGenerating population", t, "...\n")
      previous <- particles[[t - 1]]
      w_prev <- weights[[t - 1]]
      d_prev <- distances[[t - 1]]
      
      # tolerance reduction
      epsilons[t] <- min(epsilons[t - 1],
                         as.numeric(quantile(d_prev,probs = alpha, names = FALSE)))
      
      # weighted covariance
      center <- colSums(sweep(previous, 1, w_prev, "*"))
      centered <- sweep(previous, 2, center, "-")
      weighted_centered <- sweep(centered, 1, sqrt(w_prev), "*")
      
      cov_prev <- crossprod(weighted_centered)
      cov_prop <- perturbation_scale * cov_prev
      cov_prop <- cov_prop + diag(1e-8, p)
      
      # calculate covariance decomposition only once per population.
      chol_prop <- chol(cov_prop)
      inv_cov <- chol2inv(chol_prop)
      log_det <- 2 * sum(log(diag(chol_prop)))
      log_norm_const <- -0.5 * (p * log(2 * pi) + log_det)
      
      # generate particle in parallel batches
      theta_new <- matrix(NA_real_, nrow = N, ncol = p)
      
      d_new <- numeric(N)
      log_w_new <- numeric(N)
      
      i <- 1L
      attempts <- 0L
      
      while (i <= N && attempts < max_attempts) {
        n_batch <- min(batch_size, max_attempts - attempts)
        theta_batch <- matrix(NA_real_, nrow = n_batch, ncol = p)
        valid <- logical(n_batch)
        
        for (j in seq_len(n_batch)) {
          ancestor_id <- sample.int(N, size = 1, prob = w_prev)
          ancestor <- previous[ancestor_id, ]
          
          theta <- as.numeric(ancestor + t(chol_prop) %*% rnorm(p))
          
          attempts <- attempts + 1L
          
          if (prior_density(theta) <= 0) {
            next
          }
          
          theta_batch[j, ] <- theta
          valid[j] <- TRUE
        }
        
        valid_ids <- which(valid)
        
        if (length(valid_ids) == 0) {
          next
        }
        
        valid_theta <- theta_batch[valid_ids, ]
        
        seeds <- sample.int(.Machine$integer.max, length(valid_ids))
        
        # run TRANSFIL simulations simultaneously.
        simulations <- parallel_simulate(valid_theta, seeds, cl)
        batch_distances <- vapply(simulations, distance_function, numeric(1),
                                  observed = s_obs)
        
        nsim <- nsim + length(valid_ids)
        

        # accept particles and calculate importance weights
        for (j in seq_along(valid_ids)) {
          if (batch_distances[j] > epsilons[t]) {
            next
          }
          
          theta <- valid_theta[j, ]
          log_q <- log_mixture_density(theta = theta, previous = previous,
                                       w_prev = w_prev, inv_cov = inv_cov,
                                       log_norm_const = log_norm_const)
          
          if (!is.finite(log_q)) {
            next
          }
          
          prior_val <- prior_density(theta)
          
          log_w <- log(prior_val) - log_q
          
          theta_new[i, ] <- theta
          d_new[i] <- batch_distances[j]
          log_w_new[i] <- log_w
          
          i <- i + 1L
          
          if (i > N) {
            break
          }
        }
      }
      
      if (i <= N) {
        stop("Population ", t, " could not be completed. Increase max_attempts ",
             "or relax the tolerance.")
      }

      # normalize importance weights
      log_w_new <- log_w_new - max(log_w_new)
      w_new <- exp(log_w_new)
      w_new <- w_new / sum(w_new)
      colnames(theta_new) <- names(lower)
      
      particles[[t]] <- theta_new
      weights[[t]] <- w_new
      distances[[t]] <- d_new
      
      ess <- 1 / sum(w_new^2)
      
      cat("Population:", t, 
          "| Epsilon:", round(epsilons[t], 6),
          "| ESS:", round(ess, 1),
          "| Mean distance:", round(mean(d_new), 6),
          "| Simulations:", nsim, "\n")
    }
  }
  
  # return results
  list(particles = particles, weights = weights, distances = distances,
       epsilons = epsilons, nsim = nsim, N = N, T = T, alpha = alpha,
       observed_summary = s_obs)
}


# start parallel cluster and run AMC-PMC ---------------------------------------
cl <- parallel::makeCluster(cores, type = "PSOCK")

tryCatch({
  # export model file location to workers.
  parallel::clusterExport(cl, "model_file", envir = environment())
  parallel::clusterEvalQ(cl, {
    library(Rcpp)
    
    Rcpp::sourceCpp(model_file, rebuild = FALSE, showOutput = FALSE)
    
    NULL
  })
  
  cat("\nTRANSFIL successfully initialized on workers.\n")
  
  # run ABC-PMC.
  abc_fit <- abc_pmc(N = N, T = T, alpha = alpha,
                     perturbation_scale = perturbation_scale,
                     max_attempts = max_attempts,
                     batch_size = batch_size, cl = cl)
  
  # Save ABC-PMC object.
  
  saveRDS(abc_fit, file.path(out, "abc_fit.rds"))
  
  # extract final posterior
  final_particles <- abc_fit$particles[[T]]
  final_weights <- abc_fit$weights[[T]]
  final_distances <- abc_fit$distances[[T]]
  
  posterior <- data.frame(k = final_particles[, 1],
                          vth = final_particles[, 2],
                          weight = final_weights, 
                          distance = final_distances)
  
  write.csv(posterior, file.path(out, "posterior_samples.csv"), row.names = FALSE)
  

  # weighted posterior summaries
  weighted_mean <- function(x, w) {
    sum(x * w) / sum(w)
  }
  
  weighted_quantile <- function(x, w, probs = c(0.025, 0.5, 0.975)) {
    ord <- order(x)
    
    x <- x[ord]
    w <- w[ord] / sum(w)
    
    x[vapply(probs, function(prob) {
      which(cumsum(w) >= prob)[1]
    }, integer(1))]
  }
  
  posterior_summary <- data.frame(parameter = c("k", "vth"),
      mean = c(weighted_mean(posterior$k, final_weights),
               weighted_mean(posterior$vth, final_weights)),
      lower_95 = c(weighted_quantile(posterior$k, final_weights)[1],
                   weighted_quantile(posterior$vth, final_weights)[1]),
      median = c(weighted_quantile(posterior$k, final_weights)[2],
                 weighted_quantile(posterior$vth, final_weights)[2]),
      upper_95 = c(weighted_quantile(posterior$k, final_weights)[3],
                   weighted_quantile(posterior$vth, final_weights)[3]))
  
  print(posterior_summary)
  
  write.csv(posterior_summary, file.path(out, "posterior_summary.csv"), 
            row.names = FALSE)
  
  # parallel posterior predictive simulations
  cat("\nRunning posterior predictive simulations...\n")
  
  pp_seeds <- sample.int(.Machine$integer.max, nrow(posterior))
  pp_simulations <- parallel_simulate(as.matrix(posterior[, c("k", "vth")]),
                                      pp_seeds, cl)

  pp_worm_prev <- vapply(pp_simulations,
    function(x) {
      if (length(x) != 1 || !is.finite(x)) {
        stop("Invalid posterior predictive output.")
      }
      x
    },
    numeric(1)
  )
  
  pp <- data.frame(k = posterior$k,
                   vth = posterior$vth,
                   worm_prev = pp_worm_prev,
                   weight = posterior$weight)
  
  write.csv(pp, file.path(out, "posterior_predictive.csv"), row.names = FALSE)
  
  # ABC-PMC diagnostics --------------------------------------------------------
  ess_by_population <- vapply(abc_fit$weights,
                              function(w) 1 / sum(w^2),
                              numeric(1))
  
  diagnostics <- data.frame(population = seq_len(T),
                            epsilon = abc_fit$epsilons,
                            mean_distance = vapply(abc_fit$distances, mean,
                                                   numeric(1)),
                            ESS = ess_by_population)
  
  print(diagnostics)
  
  write.csv(diagnostics, file.path(out, "abc_diagnostics.csv"), 
            row.names = FALSE)

  # tolerance plot 
  png(file.path(out, "figures", "abc_pmc_tolerance.png"), 
      width = 900, height = 650)
  
  plot(seq_len(T), abc_fit$epsilons, type = "b", pch = 19, 
       xlab = "Population", ylab = expression(epsilon),
       main = "ABC-PMC tolerance schedule")
  
  dev.off()
  
  # posterior particle plot 
  png(file.path(out, "figures", "abc_pmc_posterior.png"),
      width = 900, height = 650)
  
  plot(posterior$k, posterior$vth, pch = 19, col = rgb(0, 0, 1, 0.3), 
       xlab = "k", ylab = "vth",
       main = "Final ABC-PMC posterior particles")
  
  dev.off()
  
  # final output 
  message("\nFinal tolerance: ", signif(tail(abc_fit$epsilons, 1), 6))
  message("Total simulations: ", abc_fit$nsim)
  message("Final effective sample size: ", round(tail(ess_by_population, 1), 1))
  message("Runtime: ", round(as.numeric(difftime(Sys.time(), start.time,
                                                 units = "hours")), 2), " hours")
  
}, finally = {
  parallel::stopCluster(cl)
  
})
