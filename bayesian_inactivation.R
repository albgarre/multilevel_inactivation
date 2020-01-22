
library(rethinking)
library(tidyverse)
library(readxl)
library(RColorBrewer)
library(broom)
library(ellipse)
library(bioinactivation)

## Set the seed

set.seed(1241)

## Load the data

data_55 <- c("ScottA", "LO28", "EGDe", "AOPM3", "F2365", "C5", "H7764", "H7962", "F12", "F13",
  "F14", "F15", "F16", "F17", "F18", "F19", "F20", "F21", "F33", "L6") %>%
    set_names(., .) %>%
    map(., ~ read_excel("./data/Data Diah/Inactivation_55C.xlsx",
                        sheet = ., na = c("", "NA"))
    ) %>%
    imap_dfr(., ~ mutate(.x, strain = .y)) %>%
    mutate(bio_c = paste(strain, day, sep = "-"),
           rep_c = paste(strain, day, rep, sep = "-")
    ) %>%
    mutate(temp = 55,
           strain = ifelse(strain == "ScottA", "Scott A", strain)
           )

data_60 <- c("Scott A", "LO28", "EGDe", "AOPM3", "F2365", "C5", "H7764", "H7962", "F12", "F13",
             "F14", "F15", "F16", "F17", "F18", "F19", "F20", "F21", "F33", "L6") %>%
    set_names(., .) %>%
    map(., ~ read_excel("./data/Data Diah/Inactivation_60C.xlsx",
                        sheet = ., na = c("", "NA"))
    ) %>%
    imap_dfr(., ~ mutate(.x, strain = .y)) %>%
    mutate(bio_c = paste(strain, day, sep = "-"),
           rep_c = paste(strain, day, rep, sep = "-")
    ) %>%
    mutate(temp = 60, time = time/60)

data_65 <- c("Scott A", "LO28", "EGDe", "AOPM3", "F2365", "C5", "H7764", 
             "H7962", "F12", "F13",
             "F14", "F15", "F16", "F17", "F18", "F19", "F20", "F21", "F33", "L6") %>%
    set_names(., .) %>%
    map(., ~ read_excel("./data/Data Diah/Inactivation_65C.xlsx",
                        sheet = ., na = c("", "NA"))
    ) %>%
    imap_dfr(., ~ mutate(.x, strain = .y)) %>%
    mutate(bio_c = paste(strain, day, sep = "-"),
           rep_c = paste(strain, day, rep, sep = "-")
    ) %>%
    mutate(temp = 65, time = time/60) 

## Figure 1

bind_rows(data_55, data_60, data_65) %>%
  mutate(temp = paste(temp, "ºC")) %>%
  mutate(strain = gsub("F1", "FBR1", strain)) %>%
  mutate(strain = ifelse(strain == "F20", "FBR20", strain)) %>%
  mutate(strain = ifelse(strain == "F21", "FBR21", strain)) %>%
  mutate(strain = ifelse(strain == "F33", "FBR33", strain)) %>%
  ggplot() +
    geom_point(aes(x = time, y = logN, colour = strain), shape = 1, size = 2) +
    facet_wrap("temp", scales = "free_x") +
    theme_bw() +
    theme(legend.position = "right",
          legend.title = element_blank()) +
    xlab("Time (min)") + ylab("Microbial count (log CFU/ml)") +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          strip.text = element_text(size = 12)
          )

bind_rows(data_55, data_60, data_65) %>%
    ggplot() +
    geom_point(aes(x = time, y = logN, colour = factor(temp))) +
    facet_wrap("strain", scales = "free_x") + scale_x_log10()

## Prepare the data for the bayesian model

dd_55 <- data_55 %>% 
    filter(!is.na(day)) %>%
    select(time, strain, bio_c, rep_c, rep, logN) %>%
    mutate(strain = gsub("F1", "FBR1", strain)) %>%
    mutate(strain = ifelse(strain == "F20", "FBR20", strain)) %>%
    mutate(strain = ifelse(strain == "F21", "FBR21", strain)) %>%
    mutate(strain = ifelse(strain == "F33", "FBR33", strain)) %>%
    na.omit() %>%
    mutate(strain_i = coerce_index(strain),
           bio_i = coerce_index(bio_c)) %>%
    separate(bio_c, into = c("foo", "bio"), sep = "-") %>% 
    mutate(bio = as.integer(bio)) %>%
    mutate(bio_1 = as.integer(bio == 1),
           bio_2 = as.integer(bio == 2),
           bio_3 = as.integer(bio == 3))

dd_60 <- data_60 %>% 
    filter(!is.na(day)) %>%
    mutate(strain = gsub("F1", "FBR1", strain)) %>%
    mutate(strain = ifelse(strain == "F20", "FBR20", strain)) %>%
    mutate(strain = ifelse(strain == "F21", "FBR21", strain)) %>%
    mutate(strain = ifelse(strain == "F33", "FBR33", strain)) %>%
    select(time, strain, bio_c, rep_c, rep, logN) %>%
    na.omit() %>%
    mutate(strain_i = coerce_index(strain),
           bio_i = coerce_index(bio_c)) %>%
    separate(bio_c, into = c("foo", "bio"), sep = "-") %>% 
    mutate(bio = as.integer(bio)) %>%
    mutate(bio_1 = as.integer(bio == 1),
           bio_2 = as.integer(bio == 2),
           bio_3 = as.integer(bio == 3))

dd_65 <- data_65 %>% 
    mutate(strain = gsub("F1", "FBR1", strain)) %>%
    mutate(strain = ifelse(strain == "F20", "FBR20", strain)) %>%
    mutate(strain = ifelse(strain == "F21", "FBR21", strain)) %>%
    mutate(strain = ifelse(strain == "F33", "FBR33", strain)) %>%
    filter(!is.na(day)) %>%
    select(time, strain, bio_c, rep_c, rep, logN) %>%
    na.omit() %>%
    mutate(strain_i = coerce_index(strain),
           bio_i = coerce_index(bio_c)) %>%
    separate(bio_c, into = c("foo", "bio"), sep = "-") %>% 
    mutate(bio = as.integer(bio)) %>%
    mutate(bio_1 = as.integer(bio == 1),
           bio_2 = as.integer(bio == 2),
           bio_3 = as.integer(bio == 3))

## Fit the models with random logN0

random_55 <- ulam(
    alist(
        
        # Fixed likelihood
        
        logN ~ dnorm(mu, sigma),
        mu <- logN0[strain_i] - 6*(t/6/exp(logD))^exp(logp),
        
        # Priors for the model parameters
        
        logD <- a_c1[strain_i]*bio_1 + a_c2[strain_i]*bio_2 + a_c3[strain_i]*bio_3,
        logp <- b_c1[strain_i]*bio_1 + b_c2[strain_i]*bio_2 + b_c3[strain_i]*bio_3,
        logN0[strain_i] ~ dnorm(8, sigma_logN0),
        
        # Priors for hyperparameters
        
        c(a_c1, a_c2, a_c3)[strain_i] ~ multi_normal( c(3,3,3) , Rho_ac , sigma_ac ),
        c(b_c1, b_c2, b_c3)[strain_i] ~ multi_normal( c(0,0,0) , Rho_bc , sigma_bc ),
        
        # Priors for variances
        
        sigma_ac ~ dexp(1),
        sigma_bc ~ dexp(1),
        sigma_logN0 ~ dexp(1),
        sigma ~ dexp(1),
        Rho_ac ~ lkj_corr(2),
        Rho_bc ~ lkj_corr(2)
    ),
    data = list(logN = dd_55$logN,
                t = dd_55$time,
                strain_i = dd_55$strain_i,
                bio_1 = dd_55$bio_1,
                bio_2 = dd_55$bio_2,
                bio_3 = dd_55$bio_3),
    iter = 4000, warmup = 1000, chains = 1, cores = 1, log_lik = TRUE
)

# traceplot(random_55)
par(mfrow = c(1,1))
precis(random_55, depth = 3) %>% plot()

random_60 <- ulam(
    alist(
        
        # Fixed likelihood
        
        logN ~ dnorm(mu, sigma),
        mu <- logN0[strain_i] - 6*(t/6/exp(logD))^exp(logp),
        
        # Priors for the model parameters
        
        logD <- a_c1[strain_i]*bio_1 + a_c2[strain_i]*bio_2 + a_c3[strain_i]*bio_3,
        logp <- b_c1[strain_i]*bio_1 + b_c2[strain_i]*bio_2 + b_c3[strain_i]*bio_3,
        logN0[strain_i] ~ dnorm(8, sigma_logN0),
        
        # Priors for hyperparameters
        
        c(a_c1, a_c2, a_c3)[strain_i] ~ multi_normal( c(-1,-1,-1) , Rho_ac , sigma_ac ),
        c(b_c1, b_c2, b_c3)[strain_i] ~ multi_normal( c(0,0,0) , Rho_bc , sigma_bc ),
        
        # Priors for variances
        
        sigma_ac ~ dexp(1),
        sigma_bc ~ dexp(1),
        sigma_logN0 ~ dexp(1),
        sigma ~ dexp(1),
        Rho_ac ~ lkj_corr(2),
        Rho_bc ~ lkj_corr(2)
    ),
    data = list(logN = dd_60$logN,
                t = dd_60$time,
                strain_i = dd_60$strain_i,
                bio_1 = dd_60$bio_1,
                bio_2 = dd_60$bio_2,
                bio_3 = dd_60$bio_3),
    iter = 4000, warmup = 1000, chains = 1, cores = 1, log_lik = TRUE
)

# traceplot(random_60)
par(mfrow = c(1,1))
precis(random_60, depth = 3) %>% plot()

random_65 <- ulam(
    alist(
        
        # Fixed likelihood
        
        logN ~ dnorm(mu, sigma),
        mu <- logN0[strain_i] - 6*(t/6/exp(logD))^exp(logp),
        
        # Priors for the model parameters
        
        logD <- a_c1[strain_i]*bio_1 + a_c2[strain_i]*bio_2 + a_c3[strain_i]*bio_3,
        logp <- b_c1[strain_i]*bio_1 + b_c2[strain_i]*bio_2 + b_c3[strain_i]*bio_3,
        logN0[strain_i] ~ dnorm(8, sigma_logN0),
        
        # Priors for hyperparameters
        
        c(a_c1, a_c2, a_c3)[strain_i] ~ multi_normal( c(-4,-4,-4) , Rho_ac , sigma_ac ),
        c(b_c1, b_c2, b_c3)[strain_i] ~ multi_normal( c(0,0,0) , Rho_bc , sigma_bc ),
        
        # Priors for variances
        
        sigma_ac ~ dexp(1),
        sigma_bc ~ dexp(1),
        sigma_logN0 ~ dexp(1),
        sigma ~ dexp(1),
        Rho_ac ~ lkj_corr(2),
        Rho_bc ~ lkj_corr(2)
    ),
    data = list(logN = dd_65$logN,
                t = dd_65$time,
                strain_i = dd_65$strain_i,
                bio_1 = dd_65$bio_1,
                bio_2 = dd_65$bio_2,
                bio_3 = dd_65$bio_3),
    iter = 2000, warmup = 1000, chains = 1, cores = 1, log_lik = TRUE
)

# traceplot(random_65)
par(mfrow = c(1,1))
precis(random_65, depth = 3) %>% plot()

## Fit the models with fix logN0

fixed_55 <- ulam(
    alist(
        
        # Fixed likelihood
        
        logN ~ dnorm(mu, sigma),
        mu <- logN0 - 6*(t/6/exp(logD))^exp(logp),
        
        # Priors for the model parameters
        
        logD <- a_c1[strain_i]*bio_1 + a_c2[strain_i]*bio_2 + a_c3[strain_i]*bio_3,
        logp <- b_c1[strain_i]*bio_1 + b_c2[strain_i]*bio_2 + b_c3[strain_i]*bio_3,
        logN0 ~ dnorm(8, 1),
        
        # Priors for hyperparameters
        
        c(a_c1, a_c2, a_c3)[strain_i] ~ multi_normal( c(3,3,3) , Rho_ac , sigma_ac ),
        c(b_c1, b_c2, b_c3)[strain_i] ~ multi_normal( c(0,0,0) , Rho_bc , sigma_bc ),
        
        # Priors for variances
        
        sigma_ac ~ dexp(1),
        sigma_bc ~ dexp(1),
        sigma ~ dexp(1),
        Rho_ac ~ lkj_corr(2),
        Rho_bc ~ lkj_corr(2)
    ),
    data = list(logN = dd_55$logN,
                t = dd_55$time,
                strain_i = dd_55$strain_i,
                bio_1 = dd_55$bio_1,
                bio_2 = dd_55$bio_2,
                bio_3 = dd_55$bio_3),
    iter = 2000, warmup = 1000, chains = 1, cores = 1, log_lik = TRUE
)

# traceplot(fixed_55)
par(mfrow = c(1,1))
precis(fixed_55, depth = 3) %>% plot()

fixed_60 <- ulam(
    alist(
        
        # Fixed likelihood
        
        logN ~ dnorm(mu, sigma),
        mu <- logN0 - 6*(t/6/exp(logD))^exp(logp),
        
        # Priors for the model parameters
        
        logD <- a_c1[strain_i]*bio_1 + a_c2[strain_i]*bio_2 + a_c3[strain_i]*bio_3,
        logp <- b_c1[strain_i]*bio_1 + b_c2[strain_i]*bio_2 + b_c3[strain_i]*bio_3,
        logN0 ~ dnorm(8, 1),
        
        # Priors for hyperparameters
        
        c(a_c1, a_c2, a_c3)[strain_i] ~ multi_normal( c(-1,-1,-1) , Rho_ac , sigma_ac ),
        c(b_c1, b_c2, b_c3)[strain_i] ~ multi_normal( c(0,0,0) , Rho_bc , sigma_bc ),
        
        # Priors for variances
        
        sigma_ac ~ dexp(1),
        sigma_bc ~ dexp(1),
        sigma ~ dexp(1),
        Rho_ac ~ lkj_corr(2),
        Rho_bc ~ lkj_corr(2)
    ),
    data = list(logN = dd_60$logN,
                t = dd_60$time,
                strain_i = dd_60$strain_i,
                bio_1 = dd_60$bio_1,
                bio_2 = dd_60$bio_2,
                bio_3 = dd_60$bio_3),
    iter = 2000, warmup = 1000, chains = 1, cores = 1, log_lik = TRUE
)

# traceplot(fixed_60)
par(mfrow = c(1,1))
precis(fixed_60, depth = 3) %>% plot()

fixed_65 <- ulam(
    alist(
        
        # Fixed likelihood
        
        logN ~ dnorm(mu, sigma),
        mu <- logN0 - 6*(t/6/exp(logD))^exp(logp),
        
        # Priors for the model parameters
        
        logD <- a_c1[strain_i]*bio_1 + a_c2[strain_i]*bio_2 + a_c3[strain_i]*bio_3,
        logp <- b_c1[strain_i]*bio_1 + b_c2[strain_i]*bio_2 + b_c3[strain_i]*bio_3,
        logN0 ~ dnorm(8, 1),
        
        # Priors for hyperparameters
        
        c(a_c1, a_c2, a_c3)[strain_i] ~ multi_normal( c(-4,-4,-4) , Rho_ac , sigma_ac ),
        c(b_c1, b_c2, b_c3)[strain_i] ~ multi_normal( c(0,0,0) , Rho_bc , sigma_bc ),
        
        # Priors for variances
        
        sigma_ac ~ dexp(1),
        sigma_bc ~ dexp(1),
        sigma ~ dexp(1),
        Rho_ac ~ lkj_corr(2),
        Rho_bc ~ lkj_corr(2)
    ),
    data = list(logN = dd_65$logN,
                t = dd_65$time,
                strain_i = dd_65$strain_i,
                bio_1 = dd_65$bio_1,
                bio_2 = dd_65$bio_2,
                bio_3 = dd_65$bio_3),
    iter = 2000, warmup = 1000, chains = 1, cores = 1, log_lik = TRUE
)

# traceplot(fixed_65)
par(mfrow = c(1,1))
precis(fixed_65, depth = 3) %>% plot()

## Comparison of information between models

compare(random_65, fixed_65) %>% plot()
compare(random_60, fixed_60) %>% plot()
compare(random_55, fixed_55) %>% plot()

## Compare model parameters

all_pars <- list(# random_55 = random_55,
     # random_60 = random_60,
     # random_65 = random_65,
     fixed_55 = fixed_55,
     fixed_60 = fixed_60,
     fixed_65 = fixed_65
     ) %>%
    map(precis, depth = 3) %>%
    map(as.data.frame) %>%
    map(., ~ rownames_to_column(., "par")) %>%
    imap_dfr(., ~ mutate(.x, data = .y)) %>%
    separate(data, into = c("logN0", "temp"), sep = "_")

precis(fixed_55, depth = 3) %>%
  as.data.frame() %>%
  rownames_to_column("par") %>%
  filter(grepl("sigma", par)) %>%
  ggplot(aes(x = par, y = mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = `5.5%`, ymax = `94.5%`), width = .5) +
  # facet_wrap("temp", nrow = 3) +
  coord_flip() +
  ylab("Estimate") + xlab("")

all_pars %>%
  filter(grepl("sigma", par)) %>%
  filter(logN0 == "fixed") %>%
  split(.$temp)
  

all_pars %>%
    filter(grepl("sigma", par)) %>%
    ggplot(aes(x = par, y = mean, colour = logN0)) +
        geom_point() +
        geom_errorbar(aes(ymin = `5.5%`, ymax = `94.5%`), width = .5) +
        facet_wrap("temp", nrow = 3) +
        coord_flip() +
        ylab("Estimate") + xlab("")

all_pars %>%
    filter(grepl("Rho_", par)) %>%
    ggplot(aes(x = par, y = mean, colour = logN0)) +
    geom_point() +
    geom_errorbar(aes(ymin = `5.5%`, ymax = `94.5%`), width = .5) +
    facet_wrap("temp", nrow = 3) +
    coord_flip() +
    ylab("Estimate") + xlab("") +
    ylim(-1, 1)

all_pars %>%
    filter(grepl("a_c", par)) %>%
    ggplot(aes(x = par, y = mean, colour = logN0)) +
    geom_point() +
    geom_errorbar(aes(ymin = `5.5%`, ymax = `94.5%`), width = .5) +
    facet_wrap("temp", nrow = 3) +
    coord_flip() +
    ylab("Estimate") + xlab("")


all_pars %>%
    filter(grepl("b_c", par)) %>%
    ggplot(aes(x = par, y = mean, colour = logN0)) +
    geom_point() +
    geom_errorbar(aes(ymin = `5.5%`, ymax = `94.5%`), width = .5) +
    facet_wrap("temp", nrow = 3) +
    coord_flip() +
    ylab("Estimate") + xlab("")

## Fit the fixed effects models

fit_meetselaar <- function(my_data, n_D) {
    nls(logN ~ logN0 - n_D*(time/D/n_D)^beta,
        start = list(logN0 = 8, D = 1/200, beta = 1),
        data = my_data)
}

n_D <- 6

by_bio_55 <- data_55 %>%
    filter(!is.na(day)) %>%
    # mutate(time = ifelse(time  0, 1e-6, time)) %>%
    split(.$bio_c) %>%
    map(., ~ fit_meetselaar(., n_D))

by_bio_60 <- data_60 %>%
    filter(!is.na(day)) %>%
    # mutate(time = ifelse(time == 0, 1e-6, time)) %>%
    split(.$bio_c) %>%
    map(., ~ fit_meetselaar(., n_D))

by_bio_65 <- data_65 %>%
    filter(!is.na(day)) %>%
    # mutate(time = ifelse(time == 0, 1e-6, time)) %>%
    split(.$bio_c) %>%
    map(., ~ fit_meetselaar(., n_D))


freq_pars <- list(T_55 = by_bio_55,
                  T_60 = by_bio_60,
                  T_65 = by_bio_65) %>%
    map(., ~ map(., summary)) %>%
    map(., 
        ~ map(., ~.$coefficients)
    ) %>%
    map(., 
        ~ map(., as.data.frame)
    ) %>%
    map(., 
        ~ map(., ~ rownames_to_column(., "par"))
    ) %>%
    map(., 
        ~ imap_dfr(., ~ mutate(.x, strain = .y))
    ) %>%
    imap_dfr(., ~ mutate(.x, temp = .y)) %>%
    mutate(temp = gsub("T_", "", temp)) %>%
    separate(strain, into = c("strain", "rep"), sep = "-") %>%
    mutate(strain = ifelse(strain == "ScottA", "Scott A", strain)) %>%
    mutate(strain = gsub("F1", "FBR1", strain)) %>%
    mutate(strain = ifelse(strain == "F20", "FBR20", strain)) %>%
    mutate(strain = ifelse(strain == "F21", "FBR21", strain)) %>%
    mutate(strain = ifelse(strain == "F33", "FBR33", strain))

## Extract the posteriors

all_samples <- list(fixed_55 = fixed_55,
                    fixed_60 = fixed_60,
                    fixed_65 = fixed_65
                    #random_55 = random_55,
                    # random_60 = random_60,
                    # random_65 = random_65,
) %>%
  map(extract.samples)

## D-value

# list(dd_55, dd_60, dd_65) %>%
#     map(., ~ group_by(., strain_i, strain)) %>%
#     map(summarize)  # They do match between strains

bayesian_Ds <- all_samples %>%
    map(., ~ list(rep1 = .$a_c1,
                  rep2 = .$a_c2,
                  rep3 = .$a_c3)
    ) %>%
    map(., ~ map(., as.tibble)) %>%
    map(., ~ map(.,
        ~ set_names(., paste0("strain_", 1:ncol(.)))
    )) %>%
    map(., 
        ~ imap_dfr(., ~ mutate(.x, rep = .y))
        ) %>%
    imap_dfr(., ~ mutate(.x, data = .y)) %>%
    separate(data, into = c("model", "temp"), sep = "_") %>%
    mutate(sim = row_number()) %>%
    gather(strain_i, lnD, -sim, -rep, -model, -temp) %>%
    mutate(strain_i = gsub("strain_", "", strain_i),
           strain_i = as.integer(strain_i)) %>%
    full_join(., select(dd_55, strain_i, strain), by = "strain_i") %>%
    mutate(D = exp(lnD))
    
# p1 <- bayesian_Ds %>%
#     filter(model == "fixed") %>%
#     sample_n(10000) %>%
#     ggplot() +
#         geom_boxplot(aes(x = strain, y = D, colour = rep)) +
#         facet_wrap("temp", scales = "free_y") +
#         xlab("") + ylab("D(min)") +
#         theme(axis.text.x = element_text(angle = 45, hjust = 1),
#               legend.position = "none")

# p2 <- freq_pars %>%
#     filter(par == "D") %>%
#     rename(std = `Std. Error`) %>%
#     ggplot(aes(x = strain, y = Estimate, colour = rep)) +
#         geom_point(position = position_dodge(1)) +
#         geom_errorbar(aes(ymin = Estimate - 1*std, ymax = Estimate + 1*std), 
#                       width = .2, position = position_dodge(1)) +
#         facet_wrap("temp", scales = "free_y") +
#         xlab("") + ylab("D(min)") +
#         theme(axis.text.x = element_text(angle = 45, hjust = 1),
#               legend.position = "none")
# 
# cowplot::plot_grid(p1, p2, nrow = 2)


p0 <- bayesian_Ds %>%
  sample_n(10000) %>%
  group_by(strain, rep, temp) %>%
  summarize(m_D = mean(D), q_up = quantile(D, probs = .84), 
            q_down = quantile(D, probs = .16)) %>%
  ungroup() %>%
  mutate(temp = paste0(temp, "ºC")) %>%
  ggplot(aes(x = strain, y = m_D, colour = rep)) +
  geom_point(position = position_dodge(1)) +
  geom_errorbar(aes(ymin = q_down, ymax = q_up), 
                width = .2, position = position_dodge(1))  +
  facet_wrap("temp", scales = "free_y")

p1_points <- freq_pars %>%
  filter(par == "D") %>%
  rename(std = `Std. Error`) %>%
  ungroup() %>%
  mutate(temp = paste0(temp, "ºC")) %>%
  geom_point(aes(x = strain, y = Estimate, colour = rep), shape = 1, 
             position = position_dodge(1), data = ., inherit.aes = FALSE)

p1_bars <- freq_pars %>%
  filter(par == "D") %>%
  rename(std = `Std. Error`) %>%
  ungroup() %>%
  mutate(temp = paste0(temp, "ºC")) %>%
  geom_errorbar(aes(x = strain, ymin = Estimate - 1*std, ymax = Estimate + 1*std,
                    colour = rep), linetype = 2,
                width = .2, position = position_dodge(1),
                data = ., inherit.aes = FALSE)

my_colors <- c(brewer.pal(9, "Blues")[7:9], brewer.pal(9, "Oranges")[7:9])

## Figure 3

p0 + p1_points + p1_bars +
  xlab("") + ylab("D(min)") +
  scale_color_manual(values = my_colors) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 12)
  )


## p-value

bayesian_ps <- all_samples %>%
    map(., ~ list(rep1 = .$b_c1,
                  rep2 = .$b_c2,
                  rep3 = .$b_c3)
    ) %>%
    map(., ~ map(., as.tibble)) %>%
    map(., ~ map(.,
                 ~ set_names(., paste0("strain_", 1:ncol(.)))
    )) %>%
    map(., 
        ~ imap_dfr(., ~ mutate(.x, rep = .y))
    ) %>%
    imap_dfr(., ~ mutate(.x, data = .y)) %>%
    separate(data, into = c("model", "temp"), sep = "_") %>%
    mutate(sim = row_number()) %>%
    gather(strain_i, lnp, -sim, -rep, -model, -temp) %>%
    mutate(strain_i = gsub("strain_", "", strain_i),
           strain_i = as.integer(strain_i)) %>%
    full_join(., select(dd_55, strain_i, strain), by = "strain_i") %>%
    mutate(p = exp(lnp))

# p1 <- bayesian_ps %>%
#     filter(model == "fixed") %>%
#     sample_n(10000) %>%
#     ggplot() +
#     geom_boxplot(aes(x = strain, y = p, colour = rep)) +
#     facet_wrap("temp") +
#     xlab("") + ylab("p(-)") +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1),
#           legend.position = "none") +
#     geom_hline(yintercept = 1, linetype = 2) +
#     ylim(0.25, 2.5)
# 
# 
# p2 <- freq_pars %>%
#     filter(par == "beta") %>%
#     rename(std = `Std. Error`) %>%
#     ggplot(aes(x = strain, y = Estimate, colour = rep)) +
#     geom_point() +
#     geom_errorbar(aes(ymin = Estimate - std, ymax = Estimate + std), width = .2) +
#     facet_wrap("temp") +
#     xlab("") + ylab("p(-)") +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1),
#           legend.position = "none") +
#     geom_hline(yintercept = 1, linetype = 2) +
#     ylim(0.25, 2.5)
# 
# cowplot::plot_grid(p1, p2, nrow = 2)

p0 <- bayesian_ps %>%
  sample_n(10000) %>%
  group_by(strain, rep, temp) %>%
  summarize(m_p = mean(p), q_up = quantile(p, probs = .84), 
            q_down = quantile(p, probs = .16)) %>%
  ungroup() %>%
  mutate(temp = paste0(temp, "ºC")) %>%
  ggplot(aes(x = strain, y = m_p, colour = rep)) +
  geom_point(position = position_dodge(1)) +
  geom_errorbar(aes(ymin = q_down, ymax = q_up), 
                width = .2, position = position_dodge(1))  +
  facet_wrap("temp")
  # facet_wrap("temp", scales = "free_y")

p1_points <- freq_pars %>%
  filter(par == "beta") %>%
  rename(std = `Std. Error`) %>%
  ungroup() %>%
  mutate(temp = paste0(temp, "ºC")) %>%
  geom_point(aes(x = strain, y = Estimate, colour = rep), shape = 1, 
             position = position_dodge(1), data = ., inherit.aes = FALSE)

p1_bars <- freq_pars %>%
  filter(par == "beta") %>%
  rename(std = `Std. Error`) %>%
  ungroup() %>%
  mutate(temp = paste0(temp, "ºC")) %>%
  geom_errorbar(aes(x = strain, ymin = Estimate - 1*std, ymax = Estimate + 1*std,
                    colour = rep), linetype = 2,
                width = .2, position = position_dodge(1),
                data = ., inherit.aes = FALSE)

my_colors <- c(brewer.pal(9, "Blues")[7:9], brewer.pal(9, "Oranges")[7:9])

p0 + p1_points + p1_bars +
  xlab("") + ylab(paste(expression(beta), "(·)")) +
  scale_color_manual(values = my_colors) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  geom_hline(yintercept = 1, linetype = 2)  +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 12)
  )

## A closer look at F12 at 55ºC

p0 <- data_55 %>%
  filter(strain == "F12") %>%
  ggplot(aes(x = time, y = logN, colour = factor(day))) +
    geom_point(aes(shape = factor(day)), size = 2) +
    cowplot::theme_cowplot() +
    theme(legend.position = "none") +
    xlab("Time (min)") + ylab("Microbial count (log CFU/ml)") 
    # geom_line()

F12_55_freq <- freq_pars %>%
  filter(temp == "55", strain == "FBR12") %>%
  select(par, Estimate, rep) %>%
  spread(par, Estimate)

F12_55_Ds <- bayesian_Ds %>%
  filter(temp == "55", strain == "FBR12") %>%
  group_by(rep) %>%
  summarize(m_lnD = median(lnD)) %>%
  mutate(D = exp(m_lnD))

F12_55_ps <- bayesian_ps %>%
  filter(temp == "55", strain == "FBR12") %>%
  group_by(rep) %>%
  summarize(m_lnp = median(lnp)) %>%
  mutate(p = exp(m_lnp))

p1 <- apply(full_join(F12_55_Ds, F12_55_ps), 1, function(x) {
    
    Ds <- as.numeric(x["D"])
    ps <- as.numeric(x["p"])
    tibble(time = seq(0, 90, length = 100),
           logS = -6*(time/6/Ds)^as.numeric(ps),
           logN = logS + 7.66
           )
  }) %>%
  set_names(., paste0("day_", 1:3)) %>%
  imap_dfr(., ~ mutate(.x, day = .y)) %>%
  separate(day, into = c("foo", "day"), sep = "_") %>%
  # ggplot() +
  geom_line(aes(x = time, y = logN, colour = day), 
            size = 1, data = ., inherit.aes = FALSE, linetype = 2)


p2 <- apply(F12_55_freq, 1, function(x) {
    
    Ds <- as.numeric(x["D"])
    ps <- as.numeric(x["beta"])
    tibble(time = seq(0, 90, length = 100),
           logS = -6*(time/6/Ds)^as.numeric(ps),
           logN = logS + as.numeric(x["logN0"])
    )
  }) %>%
  set_names(., paste0("day_", 1:3)) %>%
  imap_dfr(., ~ mutate(.x, day = .y)) %>%
  separate(day, into = c("foo", "day"), sep = "_") %>%
  # ggplot() +
  geom_line(aes(x = time, y = logN, colour = day), 
            size = 1, data = ., inherit.aes = FALSE, linetype = 2)

## Figure 6

# p0 + p1 + p2
cowplot::plot_grid(p0 + p2, p0 + p1, labels = NULL)

## A closer look at F12 at 55ºC

p0 <- data_55 %>%
  filter(strain == "F12") %>%
  ggplot(aes(x = time, y = logN, colour = factor(day))) +
  geom_point(aes(shape = factor(day))) +
  cowplot::theme_cowplot() +
  theme(legend.position = "none") +
  xlab("Time (min)") + ylab("Microbial count (log CFU/ml)")
# geom_line()

F12_55_freq <- freq_pars %>%
  filter(temp == "55", strain == "FBR12") %>%
  select(par, Estimate, rep) %>%
  spread(par, Estimate)

F12_55_Ds <- bayesian_Ds %>%
  filter(temp == "55", strain == "FBR12") %>%
  group_by(rep) %>%
  summarize(m_lnD = median(lnD)) %>%
  mutate(D = exp(m_lnD))

F12_55_ps <- bayesian_ps %>%
  filter(temp == "55", strain == "FBR12") %>%
  group_by(rep) %>%
  summarize(m_lnp = median(lnp)) %>%
  mutate(p = exp(m_lnp))

p1 <- apply(full_join(F12_55_Ds, F12_55_ps), 1, function(x) {
  
  Ds <- as.numeric(x["D"])
  ps <- as.numeric(x["p"])
  tibble(time = seq(0, 90, length = 100),
         logS = -6*(time/6/Ds)^as.numeric(ps),
         logN = logS + 7.66
  )
}) %>%
  set_names(., paste0("day_", 1:3)) %>%
  imap_dfr(., ~ mutate(.x, day = .y)) %>%
  separate(day, into = c("foo", "day"), sep = "_") %>%
  # ggplot() +
  geom_line(aes(x = time, y = logN, colour = day), data = ., inherit.aes = FALSE, linetype = 2)


p2 <- apply(F12_55_freq, 1, function(x) {
  
  Ds <- as.numeric(x["D"])
  ps <- as.numeric(x["beta"])
  tibble(time = seq(0, 90, length = 100),
         logS = -6*(time/6/Ds)^as.numeric(ps),
         logN = logS + as.numeric(x["logN0"])
  )
}) %>%
  set_names(., paste0("day_", 1:3)) %>%
  imap_dfr(., ~ mutate(.x, day = .y)) %>%
  separate(day, into = c("foo", "day"), sep = "_") %>%
  # ggplot() +
  geom_line(aes(x = time, y = logN, colour = day), data = ., inherit.aes = FALSE, linetype = 2)

# p0 + p1 + p2
cowplot::plot_grid(p0 + p2, p0 + p1, labels = "AUTO")


## A closer look at L6 at 60ºC

p0 <- data_60 %>%
  filter(!is.na(day)) %>%
  filter(strain == "L6") %>%
  ggplot(aes(x = time, y = logN, colour = factor(day))) +
  geom_point() # +
# geom_line()

L6_60_freq <- freq_pars %>%
  filter(temp == "60", strain == "L6") %>%
  select(par, Estimate, rep) %>%
  spread(par, Estimate)

L6_60_Ds <- bayesian_Ds %>%
  filter(temp == "60", strain == "L6") %>%
  group_by(rep) %>%
  summarize(m_lnD = median(lnD)) %>%
  mutate(D = exp(m_lnD))

L6_60_ps <- bayesian_ps %>%
  filter(temp == "60", strain == "L6") %>%
  group_by(rep) %>%
  summarize(m_lnp = median(lnp)) %>%
  mutate(p = exp(m_lnp))

p1 <- apply(full_join(L6_60_Ds, L6_60_ps), 1, function(x) {
  
  Ds <- as.numeric(x["D"])
  ps <- as.numeric(x["p"])
  tibble(time = seq(0, 4, length = 100),
         logS = -6*(time/6/Ds)^as.numeric(ps),
         logN = logS + 7.66
  )
}) %>%
  set_names(., paste0("day_", 1:3)) %>%
  imap_dfr(., ~ mutate(.x, day = .y)) %>%
  separate(day, into = c("foo", "day"), sep = "_") %>%
  # ggplot() +
  geom_line(aes(x = time, y = logN, colour = day), data = ., inherit.aes = FALSE, linetype = 2)


p2 <- apply(L6_60_freq, 1, function(x) {
  
  Ds <- as.numeric(x["D"])
  ps <- as.numeric(x["beta"])
  tibble(time = seq(0, 4, length = 100),
         logS = -6*(time/6/Ds)^as.numeric(ps),
         logN = logS + as.numeric(x["logN0"])
  )
}) %>%
  set_names(., paste0("day_", 1:3)) %>%
  imap_dfr(., ~ mutate(.x, day = .y)) %>%
  separate(day, into = c("foo", "day"), sep = "_") %>%
  # ggplot() +
  geom_line(aes(x = time, y = logN, colour = day), data = ., inherit.aes = FALSE, linetype = 2)

p0 + p1 + p2
p0 + p1
p0 + p2

## A closer look at L6 at 65ºC

p0 <- data_65 %>%
  filter(!is.na(day)) %>%
  filter(strain == "L6") %>%
  ggplot(aes(x = time, y = logN, colour = factor(day))) +
  geom_point()# +
# geom_line()

L6_65_freq <- freq_pars %>%
  filter(temp == "65", strain == "L6") %>%
  select(par, Estimate, rep) %>%
  spread(par, Estimate)

L6_65_Ds <- bayesian_Ds %>%
  sample_n(10000) %>%
  filter(temp == "65", strain == "L6") %>%
  group_by(rep) %>%
  summarize(m_lnD = median(lnD)) %>%
  mutate(D = exp(m_lnD))

L6_65_ps <- bayesian_ps %>%
  sample_n(10000) %>%
  filter(temp == "65", strain == "L6") %>%
  group_by(rep) %>%
  summarize(m_lnp = median(lnp)) %>%
  mutate(p = exp(m_lnp))

p1 <- apply(full_join(L6_65_Ds, L6_65_ps), 1, function(x) {
  
  Ds <- as.numeric(x["D"])
  ps <- as.numeric(x["p"])
  tibble(time = seq(0, .5, length = 100),
         logS = -6*(time/6/Ds)^as.numeric(ps),
         logN = logS + 7.66
  )
}) %>%
  set_names(., paste0("day_", 1:3)) %>%
  imap_dfr(., ~ mutate(.x, day = .y)) %>%
  separate(day, into = c("foo", "day"), sep = "_") %>%
  # ggplot() +
  geom_line(aes(x = time, y = logN, colour = day), data = ., inherit.aes = FALSE, linetype = 2)


p2 <- apply(L6_65_freq, 1, function(x) {
  
  Ds <- as.numeric(x["D"])
  ps <- as.numeric(x["beta"])
  tibble(time = seq(0, .5, length = 100),
         logS = -6*(time/6/Ds)^as.numeric(ps),
         logN = logS + as.numeric(x["logN0"])
  )
}) %>%
  set_names(., paste0("day_", 1:3)) %>%
  imap_dfr(., ~ mutate(.x, day = .y)) %>%
  separate(day, into = c("foo", "day"), sep = "_") %>%
  # ggplot() +
  geom_line(aes(x = time, y = logN, colour = day), data = ., inherit.aes = FALSE, linetype = 2)

p0 + p1 + p2
p0 + p1
p0 + p2

## Shrinkage of both the p and D-values

my_sims <- sample(1:max(bayesian_Ds$sim), 1000, replace = TRUE)

med_bayes_D <- bayesian_Ds %>%
  select(rep, model, temp, strain, lnD) %>%
  sample_n(50000) %>%
  filter(model == "fixed") %>%
  group_by(temp, rep, strain) %>%
  summarize(med_lnD = median(lnD)) %>% 
  mutate(D = exp(med_lnD))

med_bayes_p <- bayesian_ps %>%
  select(rep, model, temp, strain, lnp) %>%
  sample_n(50000) %>%
  filter(model == "fixed") %>%
  group_by(temp, rep, strain) %>%
  summarize(med_lnp = median(lnp)) %>% 
  mutate(p = exp(med_lnp))

full_join(med_bayes_D, med_bayes_p) %>%
  ggplot() +
    geom_point(aes(x = D, y = p, colour = strain)) +
    facet_wrap("temp", scales = "free") +
    theme(legend.position = "none")

# freq_pars %>%
#   select(par, Estimate, strain, rep, temp) %>%
#   filter(par %in% c("D", "beta")) %>%
#   spread(par, Estimate) %>%
#   rename(p = beta) %>%
#   ggplot() +
#     geom_point(aes(x = D, y = p, colour = strain)) +
#   facet_wrap("temp", scales = "free") +
#   theme(legend.position = "none")

aa <- freq_pars %>%
  select(par, Estimate, strain, rep, temp) %>%
  filter(par %in% c("D", "beta")) %>%
  spread(par, Estimate) %>%
  rename(p = beta)

pp <- freq_pars %>%
  select(par, std = `Std. Error`, strain, rep, temp) %>%
  filter(par %in% c("D", "beta")) %>%
  spread(par, std) %>%
  rename(p_std = beta, D_std = D) %>%
  full_join(., aa) %>% 
  # ggplot() +
  geom_point(aes(x = p, y = D, colour = strain, size = p_std), data = ., inherit.aes = FALSE) 
  # facet_wrap("temp", scales = "free") +
  # theme(legend.position = "none")

aa <- freq_pars %>%
  select(par, Estimate, strain, rep, temp) %>%
  filter(par %in% c("D", "beta")) %>%
  spread(par, Estimate) %>%
  mutate(rep = paste0("rep", rep)) %>%
  ungroup() %>%
  mutate(temp = paste0(temp, "ºC"))


p <- full_join(med_bayes_D, med_bayes_p)  %>%
  ungroup() %>%
  mutate(temp = paste0(temp, "ºC")) %>%
  select(temp, rep, strain, D_bayes = D, p_bayes = p) %>%
  full_join(., aa) %>% 
  ggplot(aes(x = beta, y = D, colour = strain)) +
  geom_point() +
  geom_segment(aes(xend = p_bayes, yend = D_bayes), 
               arrow = arrow(type = "closed", length = unit(0.1, "inches"))
               ) +
  facet_wrap("temp", scales = "free")

fix_summary <- aa %>%
  group_by(temp) %>%
  summarize(c = cor(beta, D), 
            m_beta = mean(beta), 
            m_D = mean(D),
            sd_beta = sd(beta),
            sd_D = sd(D))

p1 <- apply(fix_summary, 1, function(x) {
  
  x <- as.numeric(x)
  ellipse(x[2], c(x[5], x[6]), c(x[3], x[4]), level = 0.5) %>%
    as.data.frame() 
  
}) %>%
  map2_dfr(., c(55, 60, 65), 
           ~ mutate(.x, temp = .y)) %>%
  ungroup() %>%
  mutate(temp = paste0(temp, "ºC")) %>%
  # ggplot() +
  geom_path(aes(x, y), data = ., linetype = 2, inherit.aes = FALSE)
# facet_wrap("temp", scales = "free")


p2 <- apply(fix_summary, 1, function(x) {
  
  x <- as.numeric(x)
  ellipse(x[2], c(x[5], x[6]), c(x[3], x[4]), level = 0.9) %>%
    as.data.frame() 
  
}) %>%
  map2_dfr(., c(55, 60, 65), 
           ~ mutate(.x, temp = .y)) %>%
  ungroup() %>%
  mutate(temp = paste0(temp, "ºC")) %>%
  # ggplot() +
  geom_path(aes(x, y), data = ., linetype = 2, inherit.aes = FALSE)
# facet_wrap("temp", scales = "free")

p3 <- apply(fix_summary, 1, function(x) {
  
  x <- as.numeric(x)
  ellipse(x[2], c(x[5], x[6]), c(x[3], x[4]), level = 0.99) %>%
    as.data.frame() 
  
}) %>%
  map2_dfr(., c(55, 60, 65), 
           ~ mutate(.x, temp = .y)) %>%
  ungroup() %>%
  mutate(temp = paste0(temp, "ºC")) %>%
  # ggplot() +
  geom_path(aes(x, y), data = ., linetype = 2, inherit.aes = FALSE)
# facet_wrap("temp", scales = "free")

## Figure 5

p + p1 + p2 + p3 +
  theme_bw() +
  xlab(paste(expression(beta), "(·)")) + ylab("D (min)") +
  theme(legend.title = element_blank()) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 12)
  )

## D-values of unobserved strains

summary_lnD_vals <- freq_pars %>%
  filter(par == "D") %>%
  mutate(lnEstim = log(Estimate)) %>%
  group_by(temp) %>%
  summarize(m_lnD = mean(lnEstim), sd_lnD = sd(lnEstim))

p_55 <- all_samples %>%
  map(.,
      ~ list(a_c = cbind(rowMeans(.$a_c1), rowMeans(.$a_c2), rowMeans(.$a_c3)),
             sigma_ac = .$sigma_ac)
  ) %>%
  map(., 
      ~ map(., as.data.frame)
  ) %>%
  map(., 
      ~ imap_dfc(., ~ set_names(.x, paste(.y, names(.x), sep = "_")))
  ) %>%
  imap_dfr(., ~ mutate(.x, data = .y)) %>%
  sample_n(100) %>%
  filter(grepl("fixed_55", data)) %>%
  mutate(sim = row_number()) %>%
  split(.$sim) %>%
  map_dfr(.,
          ~ tibble(x = seq(1, 4.5, length = 500),
                   y1 = dnorm(x, .$a_c_V1, .$sigma_ac_V1),
                   y2 = dnorm(x, .$a_c_V2, .$sigma_ac_V2),
                   y3 = dnorm(x, .$a_c_V3, .$sigma_ac_V3)
          )
  ) %>%
  gather(foo, y, -x) %>%
  mutate(x = x/log(10)) %>%
  group_by(x) %>%
  summarize(med = median(y),
            q10 = quantile(y, probs = .1),
            q90 = quantile(y, probs = .9)
            ) %>%
  # filter(med > 1e-4) %>%
  ggplot() +
  geom_ribbon(aes(x = x, ymin = q10, ymax = q90), alpha = 0.5) 

# p1_55 <- summary_lnD_vals %>%
#   filter(temp == "55") %>%
#   tibble(x = seq(1, 4.5, length = 100),
#          y = dnorm(x, .$m_lnD, .$sd_lnD)) %>%
#   geom_line(aes(x, y), data = ., inherit.aes = FALSE, linetype = 2)
# 
# p2_55 <- freq_pars %>%
#   filter(par == "D") %>%
#   mutate(lnEstim = log(Estimate)) %>%
#   filter(temp == "55") %>%
#   # ggplot() + 
#   geom_density(aes(lnEstim), data = ., inherit.aes = FALSE, fill = "lightblue", alpha = 0.5)

p_60 <- all_samples %>%
  map(.,
      ~ list(a_c = cbind(rowMeans(.$a_c1), rowMeans(.$a_c2), rowMeans(.$a_c3)),
             sigma_ac = .$sigma_ac)
  ) %>%
  map(., 
      ~ map(., as.data.frame)
  ) %>%
  map(., 
      ~ imap_dfc(., ~ set_names(.x, paste(.y, names(.x), sep = "_")))
  ) %>%
  imap_dfr(., ~ mutate(.x, data = .y)) %>%
  sample_n(100) %>%
  filter(grepl("fixed_60", data)) %>%
  mutate(sim = row_number()) %>%
  split(.$sim) %>%
  map_dfr(.,
          ~ tibble(x = seq(-5, 1, length = 500),
                   y1 = dnorm(x, .$a_c_V1, .$sigma_ac_V1),
                   y2 = dnorm(x, .$a_c_V2, .$sigma_ac_V2),
                   y3 = dnorm(x, .$a_c_V3, .$sigma_ac_V3)
          )
  ) %>%
  gather(foo, y, -x) %>%
  mutate(x = x/log(10)) %>%
  group_by(x) %>%
  summarize(med = median(y),
            q10 = quantile(y, probs = .1),
            q90 = quantile(y, probs = .9)
  ) %>%
  # filter(med > 1e-4) %>%
  ggplot() +
  geom_ribbon(aes(x = x, ymin = q10, ymax = q90), alpha = 0.5)

# p1_60 <- summary_lnD_vals %>%
#   filter(temp == "60") %>%
#   tibble(x = seq(-5, 1, length = 100),
#          y = dnorm(x, .$m_lnD, .$sd_lnD)) %>%
#   geom_line(aes(x, y), data = ., inherit.aes = FALSE, linetype = 2)
# 
# 
# p2_60 <- freq_pars %>%
#   filter(par == "D") %>%
#   mutate(lnEstim = log(Estimate)) %>%
#   filter(temp == "60") %>%
#   # ggplot() + 
#   geom_density(aes(lnEstim), data = ., inherit.aes = FALSE, fill = "lightblue", alpha = 0.5)

p_65 <- all_samples %>%
  map(.,
      ~ list(a_c = cbind(rowMeans(.$a_c1), rowMeans(.$a_c2), rowMeans(.$a_c3)),
             sigma_ac = .$sigma_ac)
  ) %>%
  map(., 
      ~ map(., as.data.frame)
  ) %>%
  map(., 
      ~ imap_dfc(., ~ set_names(.x, paste(.y, names(.x), sep = "_")))
  ) %>%
  imap_dfr(., ~ mutate(.x, data = .y)) %>%
  sample_n(100) %>%
  filter(grepl("fixed_65", data)) %>%
  mutate(sim = row_number()) %>%
  split(.$sim) %>%
  map_dfr(.,
          ~ tibble(x = seq(-6, -1, length = 500),
                   y1 = dnorm(x, .$a_c_V1, .$sigma_ac_V1),
                   y2 = dnorm(x, .$a_c_V2, .$sigma_ac_V2),
                   y3 = dnorm(x, .$a_c_V3, .$sigma_ac_V3)
          )
  ) %>%
  gather(foo, y, -x) %>%
  mutate(x = x/log(10)) %>%
  group_by(x) %>%
  summarize(med = median(y),
            q10 = quantile(y, probs = .1),
            q90 = quantile(y, probs = .9)
  ) %>%
  # filter(med > 1e-4) %>%
  ggplot() +
  geom_ribbon(aes(x = x, ymin = q10, ymax = q90), alpha = 0.5)

# p1_65 <- summary_lnD_vals %>%
#   filter(temp == "65") %>%
#   tibble(x = seq(-6, -1, length = 100),
#          y = dnorm(x, .$m_lnD, .$sd_lnD)) %>%
#   geom_line(aes(x, y), data = ., inherit.aes = FALSE, linetype = 2)
# 
# 
# p2_65 <- freq_pars %>%
#   filter(par == "D") %>%
#   mutate(lnEstim = log(Estimate)) %>%
#   filter(temp == "65") %>%
#   # ggplot() + 
#   geom_density(aes(lnEstim), data = ., inherit.aes = FALSE, fill = "lightblue", alpha = 0.5)

# cowplot::plot_grid(
#   p_55 + p1_55 + p2_55 + xlab("ln of D-value") + ylab(""),
#   p_60 + p1_60 + p2_60 + xlab("ln of D-value") + ylab(""),
#   p_65 + p1_65 + p2_65 + xlab("ln of D-value") + ylab(""),
#   nrow = 1
# )

cowplot::plot_grid(
  p_55  + xlab("log10 of D-value") + ylab(""),
  p_60  + xlab("log10 of D-value") + ylab(""),
  p_65  + xlab("log10 of D-value") + ylab(""),
  nrow = 1
)
  
## Variance of fixed effects vs rarity

by_bio_65 %>%
  map(summary) %>%
  map(., 
      ~ tibble(sigma = .$sigma,
               df = .$df[2],
               D = .$coefficients[2,1],
               sd_D = .$coefficients[2,2],
               beta = .$coefficients[3,1],
               sd_beta = .$coefficients[3,2]
               )
      ) %>%
  imap_dfr(.,
           ~ mutate(.x, strain_day = .y)
           ) %>%
  separate(strain_day, into = c("strain", "day"), sep = "-") %>%
  ggplot(aes(x = strain, y = D, colour = day)) +
    geom_point(aes(size = sigma)) +
    geom_errorbar(aes(ymin = D - sd_D, ymax = D + sd_D), width = 0.2)

## Comparison D-values


#- MC of one strain

which_strain <- 14

make_MC <- function(which_strain) {
  tibble(logD_rep1 = all_samples$fixed_55$a_c1[,which_strain],
         logD_rep2 = all_samples$fixed_55$a_c2[,which_strain],
         logD_rep3 = all_samples$fixed_55$a_c3[,which_strain],
         logp_rep1 = all_samples$fixed_55$b_c1[,which_strain],
         logp_rep2 = all_samples$fixed_55$b_c2[,which_strain],
         logp_rep3 = all_samples$fixed_55$b_c3[,which_strain]) %>%
    mutate(sim = 1:nrow(.)) %>%
    gather(foo, value, -sim) %>%
    separate(foo, into = c("par", "rep")) %>%
    spread(par, value) %>%
    mutate(D = exp(logD), p = exp(logp)) %>%
    mutate(sim_rep = paste(sim, rep, sep = "_")) %>%
    sample_n(200) %>%
    split(.$sim_rep) %>%
    map(., ~ tibble(time = seq(0, 100, length = 100),
                    logN = 8 - 6*(time/6/.$D)^.$p
    )
    ) %>%
    imap_dfr(., ~ mutate(.x, sim = .y)) %>%
    group_by(time) %>%
    summarize(m_logN = median(logN), q_05 = quantile(logN, probs = 0.05),
              q_95 = quantile(logN, probs = 0.95)) %>%
    mutate(q_05 = ifelse(q_05 < 0, 0, q_05))
  
}

p_single <- c(4, 5, 7) %>%
  set_names(., .) %>%
  map(make_MC) %>%
  imap_dfr(., ~ mutate(.x, strain_i = .y)) %>%
  mutate(strain_i = as.integer(strain_i)) %>%
  left_join(., select(dd_55, strain_i, strain))  %>%
  ggplot(aes(x = time)) +
  geom_ribbon(aes(ymin = q_05, ymax = q_95, fill = strain,
                  linetype = strain, colour = strain), alpha = 0.5, size = 1) +
  cowplot::theme_cowplot() +
  xlab("Time (min)") + ylab("Microbial count (log CFU/ml)") +
  theme(legend.title = element_blank(),
        legend.position = "top")

#- MC with uncertainty

sigma_e <- mean(all_samples$fixed_55$sigma)

tibble(logD_rep1 = all_samples$fixed_55$a_c1[,14],
       logD_rep2 = all_samples$fixed_55$a_c2[,14],
       logD_rep3 = all_samples$fixed_55$a_c3[,14],
       logp_rep1 = all_samples$fixed_55$b_c1[,14],
       logp_rep2 = all_samples$fixed_55$b_c2[,14],
       logp_rep3 = all_samples$fixed_55$b_c3[,14]) %>%
  mutate(sim = 1:nrow(.)) %>%
  gather(foo, value, -sim) %>%
  separate(foo, into = c("par", "rep")) %>%
  spread(par, value) %>%
  mutate(D = exp(logD), p = exp(logp)) %>%
  mutate(sim_rep = paste(sim, rep, sep = "_")) %>%
  sample_n(200) %>%
  split(.$sim_rep) %>%
  map(., ~ tibble(time = seq(0, 100, length = 100),
                  logN = 8 - 6*(time/6/.$D)^.$p
  )
  ) %>%
  imap_dfr(., ~ mutate(.x, sim = .y)) %>%
  mutate(logN_unc = logN + rnorm(nrow(.), mean = 0, sd = sigma_e)) %>%
  group_by(time) %>%
  summarize(m_logN = median(logN), 
            q_05 = quantile(logN, probs = 0.05),
            q_95 = quantile(logN, probs = 0.95),
            # q_05_u = quantile(logN_unc, probs = 0.05),
            # q_95_u = quantile(logN_unc, probs = 0.95),
            ) %>%
  mutate(q_05 = ifelse(q_05 < 0, 0, q_05),
         # q_05_u = ifelse(q_05_u < 0, 0, q_05_u),
         ) %>%
  mutate(q_95_u = q_95 + 1.96*sigma_e,
         q_05_u = q_05 - 1.96*sigma_e) %>%
  ggplot() +
    geom_ribbon(aes(x = time, ymin = q_05, ymax = q_95), alpha = .5) +
    geom_ribbon(aes(x = time, ymin = q_05_u, ymax = q_95_u), alpha = .5) +
    xlab("Time (min)") + ylab("Microbial count (log CFU/ml)") +
    cowplot::theme_cowplot()

#- Monte Carlo for unobserved strains

bind_rows(
  as.data.frame(all_samples$fixed_55$a_c1),
  as.data.frame(all_samples$fixed_55$a_c2),
  as.data.frame(all_samples$fixed_55$a_c3)
) %>%
  gather(strain, a_c) %>%
  summarize(mean(a_c))

bind_rows(
  as.data.frame(all_samples$fixed_55$b_c1),
  as.data.frame(all_samples$fixed_55$b_c2),
  as.data.frame(all_samples$fixed_55$b_c3)
) %>%
  gather(strain, b_c) %>%
  summarize(mean(b_c))

p_unobs <- tibble(
  sigma_logD = all_samples$fixed_55$sigma_ac[,1],
  sigma_logp = all_samples$fixed_55$sigma_bc[,1]
) %>%
  sample_n(50) %>%
  mutate(iter = 1:nrow(.)) %>%
  split(.$iter) %>%
  map_dfr(., ~ tibble(logD = rnorm(20, 2.73, .$sigma_logD),
                      logp = rnorm(20, -0.0042, .$sigma_logp)
  )
  ) %>%
  mutate(D = exp(logD), p = exp(logp)) %>%
  mutate(sim = 1:nrow(.)) %>%
  split(.$sim) %>%
  map(., ~ tibble(time = seq(0, 100, length = 100),
                  logN = 8 - 6*(time/6/.$D)^.$p
  )
  ) %>%
  imap_dfr(., ~ mutate(.x, sim = .y)) %>%
  group_by(time) %>%
  summarize(m_logN = median(logN), q_05 = quantile(logN, probs = 0.05),
            q_95 = quantile(logN, probs = 0.95)) %>%
  mutate(q_05 = ifelse(q_05 < 0, 0, q_05)) %>%
  ggplot(aes(x = time)) +
  geom_ribbon(aes(ymin = q_05, ymax = q_95), alpha = 0.5, colour = "black") +
  cowplot::theme_cowplot() +
  xlab("Time (min)") + ylab("Microbial count (log CFU/ml)") +
  theme(legend.title = element_blank(),
        legend.position = "top")

p_beta_55 <- all_samples %>%
  map(.,
      ~ list(a_c = cbind(rowMeans(.$b_c1), rowMeans(.$b_c2), rowMeans(.$b_c3)),
             sigma_ac = .$sigma_bc)
  ) %>%
  map(., 
      ~ map(., as.data.frame)
  ) %>%
  map(., 
      ~ imap_dfc(., ~ set_names(.x, paste(.y, names(.x), sep = "_")))
  ) %>%
  imap_dfr(., ~ mutate(.x, data = .y)) %>%
  sample_n(100) %>%
  filter(grepl("fixed_55", data)) %>%
  mutate(sim = row_number()) %>%
  split(.$sim) %>%
  map_dfr(.,
          ~ tibble(x = seq(-1, 1, length = 500),
                   y1 = dnorm(x, .$a_c_V1, .$sigma_ac_V1),
                   y2 = dnorm(x, .$a_c_V2, .$sigma_ac_V2),
                   y3 = dnorm(x, .$a_c_V3, .$sigma_ac_V3)
          )
  ) %>%
  gather(foo, y, -x) %>%
  mutate(x = x/log(10)) %>%
  group_by(x) %>%
  summarize(med = median(y),
            q10 = quantile(y, probs = .1),
            q90 = quantile(y, probs = .9)
  ) %>%
  # filter(med > 1e-4) %>%
  ggplot() +
  geom_ribbon(aes(x = x, ymin = q10, ymax = q90), alpha = 0.5) 

# cowplot::plot_grid(p_single,
#                    p_55 + xlab("log10 of D-value (min)") + ylab("Probability density") + theme_minimal(),
#                    p_beta_55 + xlab(paste("log10 of", expression(beta))) + ylab("Probability density") + theme_minimal(),
#                    p_unobs,
#                    labels = "AUTO")

cowplot::plot_grid(
                   p_55 + xlab("log10 of D-value (min)") + ylab("Probability density") + theme_minimal(),
                   p_beta_55 + xlab(paste("log10 of", expression(beta))) + ylab("Probability density") + theme_minimal(),
                   p_single,
                   p_unobs,
                   labels = "AUTO")

## Output table of model parameters

bayesian_Ds %>%
  select(rep, model, temp, strain, lnD) %>%
  sample_n(50000) %>%
  filter(model == "fixed") %>%
  mutate(logD = lnD/log(10)) %>%
  group_by(temp, strain, rep) %>%
  summarize(med_logD = mean(logD),
            sd_logD = sd(logD)) %>%
  ungroup() %>%
  group_by(temp, strain) %>%
  mutate(med_logD_rep = mean(med_logD),
            sd_logD_rep = sd(med_logD)) %>%
  ungroup() %>%
  group_by(temp) %>%
  mutate(med_logD_strain = mean(med_logD_rep),
         sd_logD_strain = sd(med_logD_rep)) %>%
  split(.$temp) %>%
  map2(., c(55, 60, 65),
       ~ write_excel_csv(., path = paste0("out_Dvals_", .y, "C.csv"))
       )
  
bayesian_ps %>%
  select(rep, model, temp, strain, p) %>%
  sample_n(50000) %>%
  filter(model == "fixed") %>%
  group_by(temp, strain, rep) %>%
  summarize(med_p = mean(p),
            sd_p = sd(p)) %>%
  ungroup() %>%
  group_by(temp, strain) %>%
  mutate(med_p_rep = mean(med_p),
         sd_p_rep = sd(sd_p)) %>%
  ungroup() %>%
  group_by(temp) %>%
  mutate(med_p_strain = mean(med_p_rep),
         sd_p_strain = sd(med_p_rep)) %>%
  split(.$temp) %>%
  map2(., c(55, 60, 65),
       ~ write_excel_csv(., path = paste0("out_betas_", .y, "C.csv"))
  )

## Export the posteriors

# saveRDS(all_samples$fixed_55, "posteriors_55.rds")
# saveRDS(all_samples$fixed_60, "posteriors_60.rds")
# saveRDS(all_samples$fixed_65, "posteriors_65.rds")

# write_tsv(dd_55, "dd_55.csv")
# write_tsv(dd_60, "dd_60.csv")
# write_tsv(dd_65, "dd_65.csv")













































