#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Estimating the systematic adherence parameter rho
#
# Based on the adherence mechanism described in:
# "Modelling strategies to break transmission of lymphatic
# filariasis - aggregation, adherence and vector competence
# greatly alter elimination"
#
# Interpretation:
#   rho = 0: adherence is completely random between rounds
#   rho = 1: every individual always adheres or never adheres
#
# In this example:
#   - rho_true is used to generate dummy data
#   - rho_hat is estimated from the simulated data
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rm(list=ls())

pacman::p_load(dplyr, tidyr, ggh4x, stats, ggpubr, readxl, writexl, ggplot2,
               officer, ggthemes, ggridges, reshape2, graphics, gtsummary, gt,
               flextable, hrbrthemes, purrr, binom, tribble, janitor, nnet,
               sandwich, brms, posterior, aod, knitr, kableExtra)


# settings ---------------------------------------------------------------------
out_dir <- "output"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
word_file <- file.path(out_dir, "Adherence results.docx")

# function to format and add a table to Word
add_word_table <- function(doc, title, data) {
  data <- as.data.frame(data)
  
  ft <- flextable(data)
  ft <- theme_booktabs(ft)
  ft <- autofit(ft)
  ft <- fontsize(ft, size = 10, part = "all")
  ft <- bold(ft, part = "header")
  ft <- align(ft, align = "left", part = "all")
  ft <- valign(ft, valign = "center", part = "all")
  
  doc <- body_add_par(doc, title, style = "heading 2")
  doc <- body_add_flextable(doc, value = ft)
  doc <- body_add_par(doc, "")
  
  doc
}

doc <- read_docx() # word document

# data -------------------------------------------------------------------------

study_data <- read_excel("Ellembelle MDA data.xlsx")
dat <- study_data |> 
  select(Community, Age, Sex, `2019 A`, `2021 A`, `2022 A`)

# selecting eligible persons in 2019 (baseline)
dat <- dat |> filter(Age > 4, `2019 A` %in% c(1, "AB", "RE")) # proportion at baseline

## dropping records of ineligible individuals in 2021 and 2022
dat <- dat |> filter(!(`2021 A` %in% c("BF", "DD", "PR", "SS")),
                     !(`2022 A` %in% c("BF", "DD", "PR", "SS")))

dat <- dat |> rename(status_2019 = `2019 A`, 
                     status_2021 = `2021 A`, 
                     status_2022 = `2022 A`)

## fill NA values in 2021 and 2022 as "Lost to follow-up"
status_cols <- c("status_2021", "status_2022")
dat[status_cols] <- lapply(dat[status_cols], function(x) {
  x[is.na(x)] <- "Lost to follow up" 
  x})

## creating different columns for each years reasons for missing MDA
dat <- dat |>
  mutate(reason_2019 = case_when(status_2019 != "1" ~ status_2019),
         reason_2021 = case_when(status_2021 != "1" ~ status_2021),
         reason_2022 = case_when(status_2022 != "1" ~ status_2022)) |>
  mutate(status_2019 = case_when(status_2019 != "1" ~ 0, status_2019 == "1" ~ 1),
         status_2021 = case_when(status_2021 != "1" ~ 0, status_2021 == "1" ~ 1),
         status_2022 = case_when(status_2022 != "1" ~ 0, status_2022 == "1" ~ 1))


## MDA data ---------------------------------------------------------------------
n_rounds <- 3
n_people <- nrow(dat)

adherence <- dat |> select(status_2019, status_2021, status_2022)
mda_data <- data.frame(id = 1:nrow(dat), adherence)

names(mda_data)[-1] <- paste0("round_", 1:n_rounds)
head(mda_data)

# now we combine each person's adherence outcomes into one string. 
# hence, for three rounds: "111" = adhered in every round
#                          "000" = did not adhere in any round
#                          "101"/"110"/"011" = intermittent adherence pattern

round_variables <- paste0("round_", 1:n_rounds)
mda_data$history <- apply(mda_data[, round_variables], MARGIN = 1,
                          FUN = paste0, collapse = "")

mda_data$systematic <- mda_data$history %in% c("111", "000")
head(mda_data)


### estimated rho --------------------------------------------------------------
# adding a column "systematic" where 111 or 000 MDA history is assigned TRUE and
# otherwise FALSE
# The mean then gives the proportion of individuals with systematic behavior.
#
# Mathematically:
# rho = number of systematic individuals / total number of individuals
rho <- mean(mda_data$systematic)
rho

history_frequency <- sort(table(mda_data$history), decreasing = TRUE)
history_frequency

n_systematic_observed <- sum(mda_data$systematic)
rho_conf_interval <- binom.test(x = n_systematic_observed, 
                                      n = n_people,
                                      conf.level = 0.95)
rho_conf_interval

rho_results <- data.frame(estimate_rho = unname(rho_conf_interval$estimate),
                          lower_95_CI = rho_conf_interval$conf.int[1],
                          upper_95_CI = rho_conf_interval$conf.int[2],
                          n_systematic = n_systematic_observed,
                          n_total = n_people)
rho_results


## rho function ----------------------------------------------------------------
### function to estimate rho given the MDA data

estimate_rho <- function(data, adherence_columns) {
  adherence_data <- data[, adherence_columns]
  
  # confirm that all adherence values are 0 or 1
  if (!all(unlist(adherence_data) %in% c(0, 1))) {
    stop("All adherence variables must be coded as 0 or 1.") }
  
  histories <- apply(adherence_data, MARGIN = 1, FUN = paste0, collapse = "")
  
  always_adherent <- paste(rep(1, ncol(adherence_data)), collapse = "")
  never_adherent <- paste(rep(0, ncol(adherence_data)), collapse = "")
  
  # individuals with systematic behaviour
  systematic <- histories %in% c(always_adherent, never_adherent)
  rho <- mean(systematic)
  
  list(rho = rho,
       n_systematic = sum(systematic),
       n_total = nrow(adherence_data),
       histories = histories,
       systematic_indicator = systematic)
}

## apply the function using the study data
rho_results <- estimate_rho(data = mda_data, adherence_columns = round_variables)
est_rho <- rho_results$rho
est_rho 



## mda patterns ----------------------------------------------------------------
n <- nrow(adherence)
all_patterns <- c("000", "001", "010", "011", "100", "101", "110", "111")
pattern_counts <- table(factor(mda_data$history, levels = all_patterns))
patterns_data <- data.frame(pattern = all_patterns,
                            count = as.integer(pattern_counts))

patterns_data$proportion <- patterns_data$count / n

cat("\nObserved pattern frequencies:\n")
print(patterns_data, row.names = FALSE)


##### testing assigning treatment prob from adherence parameter ----------------
## simulating the MDA patterns in the TRANSFIL model using the estimated rho value
## and coverage values from the MDA data

set.seed(123)
pop_size <- 100000
systematic_compliance <- rho    # rho value

# re-estimated coverage from analytical cohort\
coverage <- colMeans(mda_data[c("round_1", "round_2", "round_3")])
person_id <- seq_len(pop_size)
treatment_prob <- vector("list", length(coverage))

# treatment status matrix
treat <- matrix(0L, nrow = pop_size, ncol = length(coverage),
                dimnames = list(NULL, names(coverage)))

# simulate treatment for each year
for (yr in seq_along(coverage)) {
  treatment_coverage <- coverage[yr]
  beta_shape1 <- treatment_coverage *
    (1 - systematic_compliance) / systematic_compliance
  
  beta_shape2 <- (1 - treatment_coverage) *
    (1 - systematic_compliance) / systematic_compliance
  
  # individual treatment probabilities
  treatment_prob[[yr]] <- rbeta(pop_size, beta_shape1, beta_shape2)
  
  # simulate treatment status
  U <- runif(pop_size, 0, 1)
  
  treat[, yr] <- as.integer(U < treatment_prob[[yr]])
}

# individual data set
sim_mda_data <- data.frame(id = person_id, treat)
head(sim_mda_data)
names(sim_mda_data)[-1] <-  paste0("round_", 1:n_rounds)
round_variables <- paste0("round_", 1:n_rounds)
sim_mda_data$history <- apply(sim_mda_data[, round_variables], MARGIN = 1,
                          FUN = paste0, collapse = "")

sim_pattern_counts <- table(factor(sim_mda_data$history, levels = all_patterns))
sim_patterns_data <- data.frame(pattern = all_patterns,
                            count = as.integer(sim_pattern_counts))

# compare observed and simulated pattern frequencies
comparison <- data.frame(pattern = all_patterns,
                         observed_count = as.integer(pattern_counts),
                         simulated_count = as.integer(sim_pattern_counts))

comparison$observed_proportion <- comparison$observed_count / nrow(adherence)
comparison$simulated_proportion <- comparison$simulated_count / pop_size

cat("\nObserved versus simulated patterns:\n")
print(comparison, row.names = FALSE)

# plot patterns proportions
png("obs_vs_sim_MDA_patterns.png", width = 1200, height = 800, res = 150)

barplot(rbind(comparison$observed_proportion, comparison$simulated_proportion),
        beside = TRUE, names.arg = comparison$pattern,
        col = c("steelblue", "orange"),
        xlab = "Three-year attendance pattern",
        ylab = "Proportion",
        main = "Observed vs simulated MDA patterns",
        legend.text = c("Observed", "Simulated"),
        args.legend = list(x = "topleft", bty = "n"))
dev.off()


# flow summary -----------------------------------------------------------------
## all cohort ------------------------------------------------------------------
### frequency and counts of all data variables 

freq2019 <- data.frame(table(study_data$`2019 A`, useNA = "ifany")) %>%
  rename("2019" = "Freq") %>%
  mutate(`2019` = sprintf("%d (%.2f)", `2019`, `2019` / sum(`2019`) * 100))

freq2021 <- data.frame(table(study_data$`2021 A`, useNA = "ifany")) %>%
  rename("2021" = "Freq") %>%
  mutate(`2021` = sprintf("%d (%.2f)", `2021`, `2021` / sum(`2021`) * 100))

freq2022 <- data.frame(table(study_data$`2022 A`, useNA = "ifany")) %>%
  rename("2022" = "Freq") %>%
  mutate(`2022` = sprintf("%d (%.2f)", `2022`, `2022` / sum(`2022`) * 100))

count_table <- bind_rows(
  data.frame(Year = "2019",
             Reason = names(table(study_data$`2019 A`, useNA = "ifany")),
             Count = as.numeric(table(study_data$`2019 A`, useNA = "ifany"))),
  
  data.frame(Year = "2021",
             Reason = names(table(study_data$`2021 A`, useNA = "ifany")),
             Count = as.numeric(table(study_data$`2021 A`, useNA = "ifany"))),
  
  data.frame(Year = "2022",
             Reason = names(table(study_data$`2022 A`, useNA = "ifany")),
             Count = as.numeric(table(study_data$`2022 A`, useNA = "ifany")))
)

# create contigency table
contingency_table <- xtabs(Count ~ Reason + Year, data = count_table)
chi_test <- chisq.test(contingency_table)
p_value <- chi_test$p.value

summary_df <- freq2019 %>%
  right_join(freq2021, by = "Var1") %>%
  left_join(freq2022, by = "Var1") %>%
  rename(Reasons = Var1)


summary_df <- summary_df %>%
  mutate(`p-value` = ifelse(row_number() == 1,
                            ifelse(p_value < 0.001, "<0.001",
                                   sprintf("%.3f", p_value)), ""))

print(summary_df)


summary_df_transposed <- summary_df %>%
  t() %>%
  as.data.frame() %>%
  setNames(NULL) %>%
  row_to_names(row_number = 1)

print(summary_df_transposed)

doc <- add_word_table(doc, "Summary of all data variables", summary_df)
#doc <- add_word_table(doc, "Summary of all data variables", summary_df_transposed)

## final analytical cohort -----------------------------------------------------
### variable frequency ---------------------------------------------------------
reasons_year19 <- data.frame(table(dat$reason_2019, useNA = "ifany")) |> 
  rename("2019" = "Freq") |>
  mutate(`2019` = sprintf("%d (%.2f)", `2019`, `2019` / sum(`2019`) * 100))

reasons_year21 <- data.frame(table(dat$reason_2021, useNA = "ifany")) |> 
  rename("2021" = "Freq") |>
  mutate(`2021` = sprintf("%d (%.2f)", `2021`, `2021` / sum(`2021`) * 100))

reasons_year22 <- data.frame(table(dat$reason_2022, useNA = "ifany")) |> 
  rename("2022" = "Freq") |>
  mutate(`2022` = sprintf("%d (%.2f)", `2022`, `2022` / sum(`2022`) * 100))
#result_df <- left_join(reasons_year21, reasons_year22, by = "Var1")

reasons_for_missing_df <- reasons_year19 %>%
  right_join(reasons_year21, by = "Var1") %>%
  left_join(reasons_year22, by = "Var1") %>%
  rename(Reasons = Var1)
print(reasons_for_missing_df)

doc <- add_word_table(doc, "Reasons for missing MDA", reasons_for_missing_df)


# complete case selection bias -------------------------------------------------
### 0 = Excluded/Incomplete
### 1 = Complete Case

dat <- dat %>%
  mutate(is_complete_case = if_else(status_2019 == 1 & status_2021 == 1 & 
                                      status_2022 == 1, 1, 0))
bias_comparison <- dat %>%
  group_by(is_complete_case) %>%
  summarise(N = n(), 
            mean_age = mean(Age, na.rm = TRUE),
            sd_age = sd(Age, na.rm = TRUE),
            pct_female = mean(Sex == "Female", na.rm = TRUE) * 100)
print(bias_comparison)

## test for bias ---------------------------------------------------------------
t_test_result <- t.test(Age ~ is_complete_case, data = dat)
chi_sex <- chisq.test(table(dat$is_complete_case, dat$Sex))

print(t_test_result)
print(chi_sex)

  
# socio-demographic characteristics of the study population ------------------
# calculate age-group frequencies and percentages
## age groups
dat <- dat %>%
  mutate(age_group = cut(Age,
                         breaks = c(4, 10, 15, 20, 30, 40, 50, 60, Inf),
                         labels = c("5-10", "11-15", "16-20", "21-30",
                                    "31-40", "41-50", "51-60", ">60"),
                         right = TRUE) )


age_table <- dat %>%
  filter(!is.na(age_group)) %>%
  count(age_group, name = "Frequency") %>%
  mutate(Percentage = round(Frequency / sum(Frequency) * 100, 1),
         Variable = as.character(age_group)) %>%
  select(Variable, Frequency, Percentage)

sex_table <- dat %>%
  filter(!is.na(Sex)) %>%
  count(Sex, name = "Frequency") %>%
  mutate(Percentage = round(Frequency / sum(Frequency) * 100, 1),
         Variable = as.character(Sex)) %>%
  select(Variable, Frequency, Percentage)

demographic_table <- bind_rows(
  data.frame(Variable = "Age group", Frequency = NA_integer_, Percentage = NA_real_),
  age_table,
  data.frame(Variable = "Sex", Frequency = NA_integer_, Percentage = NA_real_),
  sex_table)

demographic_table
doc <- add_word_table(doc, "Demographic table of participants", demographic_table)

# frequency of participants taking part in MDA ---------------------------------
N <- nrow(mda_data)

# function to calculate frequency, percentage and 95% CI
calculate_stats <- function(x, label, N) {
  n <- sum(x, na.rm = TRUE)
  percentage <- 100 * n / N
  ci <- binom.confint(n, N, methods = "wilson")
  
  data.frame(`Participated in MDA` = label,
             Frequency = n,
             `Percentage % (95% CI)` = sprintf("%.1f (%.1f-%.1f)", percentage,
                                              ci$lower * 100, ci$upper * 100),
             check.names = FALSE)
}

# calculate participation patterns
mda_table <- bind_rows(
  calculate_stats(mda_data$history == "100", "Round 1 only", N),
  calculate_stats(mda_data$history == "010", "Round 2 only", N),
  calculate_stats(mda_data$history == "001", "Round 3 only", N),
  calculate_stats(mda_data$history == "110", "Round 1 and 2", N),
  calculate_stats(mda_data$history == "101", "Round 1 and 3", N),
  calculate_stats(mda_data$history == "011", "Round 2 and 3", N),
  calculate_stats(mda_data$history == "111", "Round 1, 2 and 3", N),
  calculate_stats(mda_data$history == "000", "Never", N) )

mda_table
doc <- add_word_table(doc, "MDA participation patterns", mda_table)

# levels of adherence ---------------------------------------------------------
dat <- dat |>
  mutate(adherence_levels = case_when(status_2019 == 1 & status_2021 == 1 & 
                                        status_2022 == 1 ~ "Systematically treated",
                                      status_2019 == 0 & status_2021 == 0 & 
                                        status_2022 == 0 ~ "Never treated",
                                      TRUE ~ "Non-consistently treated"))

group_order <- c("Systematically treated", "Non-consistently treated", "Never treated")

make_summary <- function(data, variable, section_label) {
  data %>%
    filter(!is.na(.data[[variable]]), !is.na(adherence_levels)) %>%
    count(category = .data[[variable]], adherence_levels, name = "n") %>%
    complete(category, adherence_levels = group_order,
      fill = list(n = 0)
    ) %>%
    group_by(category) %>%
    mutate(row_total = sum(n), percent = if_else(row_total > 0, 100 * n / row_total, NA_real_),
      result = sprintf("%d (%.1f%%)", n, percent)) %>%
    ungroup() %>%
    select(category, adherence_levels, result) %>%
    pivot_wider(names_from = adherence_levels, values_from = result) %>%
    mutate(category = as.character(category),
           section = section_label)
}

age_summary <- make_summary(dat, "age_group", "Age group")
sex_summary <- make_summary(dat, "Sex", "Sex")

final_table <- bind_rows(age_summary, sex_summary) %>%
  select(section, category, all_of(group_order))

final_table
doc <- add_word_table(doc, "Demographics against adherence levels", final_table)


### significance
age_adherence_chi_res <- chisq.test(table(dat$adherence_levels, dat$age_group))
age_adherence_chi_res

sex_adherence_chi_res <- chisq.test(table(dat$adherence_levels, dat$Sex))
sex_adherence_chi_res


####
####
# multinomial regression -------------------------------------------------------
dat <- dat %>%
  mutate(adherence_levels = factor(adherence_levels,
                                   levels = c("Systematically treated",
                                              "Non-consistently treated",
                                              "Never treated")),
         Sex = factor(Sex),
         Community = factor(Community),
         Age = as.numeric(Age))


print(table(dat$adherence_levels))
print(table(dat$Community, dat$adherence_levels))
print(table(dat$Sex, dat$adherence_levels))
print(summary(dat$Age))


# reference categories and center age
dat$adherence_levels <- relevel(dat$adherence_levels, ref = "Non-consistently treated")
dat$Sex <- relevel(dat$Sex, ref = "Female")

mean_age <- mean(dat$Age, na.rm = TRUE)

dat <- dat %>%
  mutate(age_c = Age - mean_age)

cat("Mean age:", mean_age, "\n")

# state age as at 2019
# Adherence levels by sex
sex_summary <- dat %>%
  count(Sex, adherence_levels, name = "n") %>%
  group_by(Sex) %>%
  mutate(denominator = sum(n),
         percentage = 100 * n / denominator) %>%
  ungroup()

print(sex_summary)

# Age distribution by adherence levels
age_summary <- dat %>%
  group_by(adherence_levels) %>%
  summarise(n = n(),
            mean_age = mean(Age),
            sd_age = sd(Age),
            median_age = median(Age),
            Q1 = quantile(Age, 0.25),
            Q3 = quantile(Age, 0.75),
            .groups = "drop")

print(age_summary)


## multinomial regression adjusted for age and sex ------------------------------
fit <- nnet::multinom(adherence_levels ~ age_c + Sex, data = dat, trace = FALSE,
                      Hess = TRUE)

print(summary(fit))
fit$convergence


## community-clustered covariance and adjusted RRRs -----------------------------
## cluster-robust covariance matrix with 53 communities
V_cluster <- sandwich::vcovCL(fit, cluster = dat$Community, type = "HC1")

# Align the coefficients with the covariance matrix
coef_mat <- coef(fit)
coef_vec <- as.vector(t(coef_mat))
names(coef_vec) <- rownames(V_cluster)

se_cluster <- sqrt(diag(V_cluster))

# RRRs and 95% Wald confidence intervals
z_975 <- qnorm(0.975)

regression_results <- data.frame(term = names(coef_vec),
                                 estimate = as.numeric(coef_vec),
                                 SE = as.numeric(se_cluster),
                                 row.names = NULL) %>%
  mutate(z = estimate / SE,
         p_value = 2 * pnorm(-abs(z)),
         CI_lower = estimate - z_975 * SE,
         CI_upper = estimate + z_975 * SE,
         RRR = exp(estimate),
         RRR_lower = exp(CI_lower),
         RRR_upper = exp(CI_upper))

print(regression_results)


### joint Wald tests for age and sex ---------------------------------------------
# age
wald.test(Sigma = V_cluster, b = coef_vec, 
          Terms = grep("age_c$", names(coef_vec)))

# sex
wald.test(Sigma = V_cluster, b = coef_vec,
          Terms = grep("Sex", names(coef_vec)))

# joint Wald test results
joint_tests <- data.frame(Predictor = c("Age (per additional year)", "Male (vs. female)"),
                          Wald_chisq = c(0.97, 23.6),
                          df = c(2, 2), 
                          joint_p = c(0.61, 7.7e-06))

### multi regression results ---------------------------------------------------
multi_res_table <- regression_results %>%
  mutate(Predictor = case_when(grepl("age_c", term) ~ "Age",
                               grepl("SexMale", term) ~ "Male (vs. female)",
                               TRUE ~ NA_character_),
         
         Outcome = case_when(grepl("^Systematically treated:", term) ~ "Systematically treated",
                             grepl("^Never treated:", term) ~ "Never treated",
                             TRUE ~ NA_character_),
         RRR_CI = sprintf("%.3f (%.3f–%.3f)", RRR, RRR_lower, RRR_upper),
         p_formatted = ifelse(p_value < 0.001, "<0.001", sprintf("%.3f", p_value))) %>%
  filter(!is.na(Predictor)) %>%
  select(Predictor, Outcome, RRR_CI, p_formatted) %>%
  pivot_wider(names_from = Outcome, 
              values_from = c(RRR_CI, p_formatted),
              names_glue = "{Outcome}_{.value}")

# adding joint Wald tests
multi_res_table <- multi_res_table %>%
  left_join(joint_tests, by = "Predictor") %>%
  mutate(`Joint Wald χ² (df)` = sprintf("%.2f (%d)", Wald_chisq, df),
         `Joint p-value` = ifelse(joint_p < 0.001, "<0.001", sprintf("%.3f", joint_p))) %>%
  select(Predictor, 
         `Systematically treated_RRR_CI`,
         `Systematically treated_p_formatted`,
         `Never treated_RRR_CI`,
         `Never treated_p_formatted`,
         `Joint Wald χ² (df)`,
         `Joint p-value`)

# Rename columns
names(multi_res_table) <- c("Predictor", "st RRR (95% CI)", "st p-value", 
                            "nt RRR (95% CI)", "nt p-value", "Joint Wald χ² (df)", 
                            "Joint p-value")

multi_res_table %>%
  kable(format = "html",
        align = c("l", "c", "c", "c", "c", "c", "c"),
        escape = TRUE) %>%
  kable_styling(bootstrap_options = c("striped", "condensed"), 
                full_width = TRUE, position = "center", font_size = 11) %>%
  add_header_above(c(" " = 1, "Systematically treated" = 2, 
                     "Never treated" = 2, "Joint association" = 2))

doc <- add_word_table(doc, "Multinomial regression results", multi_res_table)

print(doc, target = word_file)

