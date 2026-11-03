library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(scales)
library(ggpubr)

out <- "output"
fig_dir <- file.path(out, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

posterior <- read.csv(file.path(out, "posterior_samples.csv"))
pp <- read.csv(file.path(out, "posterior_predictive.csv"))

prev <- read.csv(file.path("data", "prevalence.csv"))
obs <- with(prev[prev$year == 2001, ], unique(worm_prevalence))
obs <- as.numeric(obs)

stopifnot(all(c("k", "vth") %in% names(posterior)))
stopifnot(all(c("k", "vth", "worm_prev") %in% names(pp)))

if (nrow(posterior) < 2 || nrow(pp) < 2) {
  stop("Insufficient valid posterior or predictive samples.")
}

# visual themes colors
blue <- "#3B75AF"
dark_blue <- "#214E78"
light_blue <- "#B9D5EA"
orange <- "#D97732"
gray <- "#777777"
light_gray <- "#E6E6E6"

theme_abc <- theme_classic(base_size = 12, base_family = "sans") +
  theme(plot.title = element_text(face = "bold", size = 13, color = "#222222"),
        plot.subtitle = element_text(size = 10, color = "#555555"),
        #axis.title = element_text(face = "bold", color = "#333333"),
        #axis.text = element_text(color = "#333333"),
        axis.line = element_line(color = "grey30"),
        legend.position = c(0.8,0.9),
        legend.title = element_blank(),
        strip.background = element_rect(fill = "#F0F0F0", color = NA),
        strip.text = element_text(face = "bold"),
        panel.grid = element_blank()) 

#scale_x_touch <- scale_x_continuous(expand = expansion(mult = 0, add = 0))
#scale_y_touch <- scale_y_continuous(expand = expansion(mult = 0, add = 0))

save_plot <- function(plot, name, width = 5, height = 4) {
  ggsave(file.path(fig_dir, paste0(name, ".png")), plot = plot, width = width,
         height = height, dpi = 320, bg = "white")
  
  #ggsave(file.path(fig_dir, paste0(name, ".pdf")), plot = plot, width = width,
  #       height = height, device = grDevices::cairo_pdf, bg = "white")
}

# posterior density plots ------------------------------------------------------
p_k <- ggplot(posterior, aes(x = k)) +
  geom_density(fill = light_blue, color = blue, linewidth = 0.5, alpha = 0.7) +
  geom_vline(aes(xintercept = median(k), linetype = "Median"), 
             color = orange, linewidth = 0.7) +
  scale_linetype_manual(values = c("Median" = "dashed")) +
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(limits = c(0, max(density(posterior$k)$y) * 1.115),
                     expand = c(0,0), breaks = breaks_extended(n = 5)) +
  labs(x = "k", y = "Posterior density") +
  theme_abc


p_vth <- ggplot(posterior, aes(x = vth)) +
  geom_density(fill = light_blue, color = blue, linewidth = 0.5, alpha = 0.7) +
  geom_vline(aes(xintercept = median(vth), linetype = "Median"),
             color = orange, linewidth = 0.7) +
  scale_linetype_manual(values = c("Median" = "dashed")) +
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(limits = c(0, max(density(posterior$vth)$y) * 1.1),
                     expand = c(0,0), breaks = breaks_extended(n = 5)) +
  labs(x = "VTH", y = "Posterior density") +
  theme_abc

save_plot(p_k, "01_posterior_k")
save_plot(p_vth, "02_posterior_vth")
save_plot(p_k + p_vth, "03_posterior_marginals", width = 10, height = 4.5)


# joint posterior distribution
# points are ABC posterior samples; contours show concentration
p_joint <- ggplot(posterior, aes(x = k, y = vth)) +
  geom_point(color = blue, alpha = 0.30, size = 1.5) +
  geom_density_2d(color = dark_blue, linewidth = 0.7, bins = 7) +
  scale_x_continuous(expand = c(0,0), limits = c(min(posterior$k),0.05)) + 
  scale_y_continuous(limits = c(30, 160),
                     expand = c(0,0), breaks = breaks_extended(n = 5)) +
  labs(x = "k", y = "VTH") +
  theme_abc

save_plot(p_joint, "04_joint_posterior")

pp_joint_plot <- (p_k | plot_spacer()) / (p_joint | p_vth) +
  plot_layout(widths = c(1, 1), heights = c(1, 1))

save_plot(pp_joint_plot, "all_posteriors", width = 8, height = 8)

# posterior pair plot ---------------------------------------------------------
p_pair <- (p_k | p_joint) +
  plot_annotation(title = "Posterior parameter diagnostics")

save_plot(p_pair, "05_parameter_diagnostics", width = 12, height = 5)


# prior vs posterior
# Uniform priors specified in the original model:
# k   ~ Uniform(0.001, 0.3)
# vth ~ Uniform(30, 150)

set.seed(123)
n_prior <- 1200
prior_draws <- data.frame(k = runif(n_prior, 0.001, 0.3),
                          vth = runif(n_prior, 30, 150))

prior_long <- prior_draws |>
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") |>
  mutate(distribution = "Prior")

post_long <- posterior |>
  pivot_longer(cols = c(k, vth), names_to = "parameter", values_to = "value") |>
  mutate(distribution = "Posterior")

comparison <- bind_rows(prior_long, post_long) |>
  mutate(parameter = factor(parameter, levels = c("k", "vth"), labels = c("k", "vth")),
         distribution = factor(distribution, levels = c("Prior", "Posterior")) )

# prior and posterior distribution
p_prior_post <- ggplot(comparison, aes(x = value, color = distribution, 
                                       fill = distribution)) +
  geom_density(alpha = 0.20, linewidth = 1, adjust = 1) +
  facet_wrap(~parameter, scales = "free", labeller = label_parsed) +
  scale_color_manual(values = c("Prior" = gray, "Posterior" = blue)) +
  scale_fill_manual(values = c("Prior" = light_gray, "Posterior" = light_blue)) +
  labs(x = "Parameter value", y = "Density") +
  theme_abc +
  theme(legend.position = "top")

save_plot(p_prior_post, "06_prior_vs_posterior", width = 9, height = 4.5)


# posterior predictive distribution --------------------------------------------
## posterior predictive check
## observed worm prevalence in 2001 versus simulated prevalence
p_ppc <- ggplot(pp, aes(x = worm_prev)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, fill = light_blue,
                 color = "white", linewidth = 0.3) +
  geom_density(color = blue, linewidth = 0.7, adjust = 1) +
  geom_vline(xintercept = obs, color = orange, linewidth = 0.7) +
  annotate("text", x = obs, y = Inf, label = "Observed", color = orange, 
           angle = 90, vjust = 1.4, hjust = 1.1, size = 3.5) +
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(limits = c(0, max(density(pp$worm_prev)$y) * 1.245),
                     expand = c(0,0), breaks = breaks_extended(n = 7)) +
  labs(x = "Simulated worm prevalence", y = "Density") +
  theme_abc
p_ppc

save_plot(p_ppc, "07_posterior_predictive")


# predictive interval plot -----------------------------------------------------
pp_summary <- data.frame(median = median(pp$worm_prev),
                         lower95 = unname(quantile(pp$worm_prev, 0.025)),
                         upper95 = unname(quantile(pp$worm_prev, 0.975)),
                         lower50 = unname(quantile(pp$worm_prev, 0.25)),
                         upper50 = unname(quantile(pp$worm_prev, 0.75)))

# observed prevalence and posterior predictive intervals
p_interval <- ggplot() +
  geom_linerange(data = pp_summary, 
                 aes(x = "Posterior predictive", ymin = lower95, ymax = upper95,
                     color = "95% predictive interval"), linewidth = 1.5) +
  geom_linerange(data = pp_summary,
                 aes( x = "Posterior predictive", ymin = lower50, ymax = upper50,
                      color = "50% predictive interval"), linewidth = 4) +
  geom_point(data = pp_summary, 
             aes( x = "Posterior predictive", y = median,
                  color = "Predictive median"), size = 5) +
  geom_point(aes(x = "Posterior predictive", y = obs, color = "Observed prevalence"),
             shape = 4, size = 4, stroke = 1.5) +
  scale_color_manual(name = "Legend",
                     values = c("Observed prevalence" = orange, 
                                "Predictive median" = dark_blue,
                                "50% predictive interval" = blue,
                                "95% predictive interval" = light_blue)) +
  labs(x = NULL, y = "Antigen prevalence") +
  coord_flip() +
  theme_abc

save_plot(p_interval, "08_predictive_intervals", 7, 4)


# parameter uncertainty --------------------------------------------------------
summarise_parameter <- function(x, parameter_name) {
  data.frame(parameter = parameter_name, median = median(x), 
             mean = mean(x), sd = sd(x),
             lower95 = unname(quantile(x, 0.025)),
             upper95 = unname(quantile(x, 0.975)) )
}

parameter_summary <- bind_rows(summarise_parameter(posterior$k, "k"),
                               summarise_parameter(posterior$vth, "vth"))

write.csv(parameter_summary, file.path(out, "parameter_summary.csv"),
          row.names = FALSE)

forest_data <- parameter_summary |>
  mutate(parameter = factor(parameter, levels = c("vth", "k")))

p_forest <- ggplot(forest_data, aes(y = parameter, x = median)) +
  geom_segment(aes(x = lower95, xend = upper95, yend = parameter),
               linewidth = 1.2, color = light_blue) +
  geom_point(color = blue, size = 3) +
  labs(title = "Posterior parameter estimates",
       subtitle = "Points: posterior medians; lines: 95% credible intervals",
       x = "Parameter value", y = NULL) +
  theme_abc

save_plot(p_forest, "09_parameter_forest")


#### save for paper
top <- ggarrange(p_prior_post, labels = "A", ncol = 1)
bottom <- ggarrange(p_joint, p_ppc, labels = c("B", "C"), ncol = 2)
all_com <- ggarrange(top, bottom, ncol = 1, nrow = 2, heights = c(1, 1))

save_plot(all_com, "parameter_stats", width = 8, height = 7)


# posterior predictive summary -------------------------------------------------
predictive_summary <- data.frame(observed = obs,
                                 pred_mean = mean(pp$worm_prev),
                                 pred_median = median(pp$worm_prev),
                                 pred_sd = sd(pp$worm_prev),
                                 pred_lower95 = unname(quantile(pp$worm_prev, 0.025)),
                                 pred_upper95 = unname(quantile(pp$worm_prev, 0.975)),
                                 prop_below_observed = mean(pp$worm_prev <= obs))

write.csv(predictive_summary, file.path(out, "predictive_summary.csv"), 
          row.names = FALSE)
