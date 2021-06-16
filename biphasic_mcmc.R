#' ------------------------------------------------------------------------
#' ------------------------------------------------------------------------
#' 
#'         Reversible jump Markov Chain Monte Carlo (RJMCMC) 
#' for investigating responses of cetaceans to military sonar exposure
#' 
#'              --- Biphasic dose-response model --- 
#' 
#' ------------------------------------------------------------------------
#' ------------------------------------------------------------------------ 
#' 
#' Author: Phil Bouchet
#' Last update: 2021-06-16
#' Project: LMR Dose-Response
#' --------------------------------------------------------------------

# Initialisation -------------------------------------------------

# Load required libraries
# install.packages("pacman")
pacman::p_load(tidyverse,     # R packages for data science
               here,          # A Simpler Way to Find Your Files
               truncnorm,     # Truncated Normal Distribution
               magrittr,      # Pipes
               progress,      # ASCII progress bars
               rjags,         # Bayesian graphical models using MCMC
               Hmisc)         # Miscellaneous functions

# Set tibble options
options(tibble.width = Inf) # All tibble columns shown
options(pillar.neg = FALSE) # No colouring negative numbers
options(pillar.subtle = TRUE)
options(pillar.sigfig = 4)

# Load functions
source("biphasic_functions.R")

# Simulated data ----------------------------------------------------------

seed <- 122
set.seed(seed)

n.whales <- 10
max.trials <- 4

# Trials and individual IDs
n.trials.per.whale <- sample(x = 1:max.trials, size = n.whales, replace = TRUE)
n.trials <- sum(n.trials.per.whale)
whale.id <- purrr::map(.x = seq_len(n.whales), 
                       .f = ~rep(.x, each = n.trials.per.whale[.x])) %>% do.call(c, .)

# Parameters
nu <- c(105, 149)
tau <- c(20, 20)
psi <- 0.5
omega <- 1 
L <- 60
U <- 215
alpha <- 125

# Right-censoring ----
Rc <- c(150, 165) # Some right-censoring
# Rc <- c(210, 211) # No right-censoring

# Parameter bounds
alpha.bounds <- c(110, 160)
omega.bounds <- c(0, 5)
tau.bounds <- c(0, 20)
obs.sd <- 2.5

# Probability of whale i displaying the context-dependent response (on probit scale)
psi_i <- rnorm(n = n.whales, mean = psi, sd = omega)

# Priors
priors <- list(psi = c(0, 10))

# Covariate effects
covariates <- list(exposed = 5, sonar = c(0, -10), behaviour = c(0, 3, -9), range = 1)
covariate.types <- list(exposed = "i", sonar = "f", behaviour = "f",range = "d")

# covariate.names <- c("exposed", "sonar", "range")
covariate.names <- NULL # No covariates in this example
covariates <- covariates[covariate.names]
covariate.types <- covariate.types[covariate.names] %>% unlist()
n.covariates <- length(covariate.names)

if(n.covariates > 0){
  
  for (j in 1:n.covariates) priors[[covariate.names[j]]] <- c(0, 30) 
  
  covariates.list <- list(exposed = purrr::map(.x = n.trials.per.whale, 
                                               .f = ~c(0, rep(1, each = .x - 1))) %>% do.call(c, .),
                          sonar = sample(x = c("LFAS", "MFAS"), size = n.trials, replace = TRUE),
                          behaviour = sapply(seq_len(n.trials), sample, x = c("Feed", "Migrate", "Rest"), size = 1),
                          range = runif(n = n.trials, min = 0, max = 20))
  
  covariates.df <- as.data.frame(do.call(cbind, covariates.list[covariate.names]))
  colnames(covariates.df) <- covariate.names
  
  for(ii in 1:length(covariate.names)){
    if(covariate.types[ii] == "f") covariates.df[, covariate.names[ii]] <- as.factor(covariates.df[, covariate.names[ii]])
    if(covariate.types[ii] %in% c("d", "i")) covariates.df[, covariate.names[ii]] <- as.numeric(covariates.df[, covariate.names[ii]])
  }
  
  # Dummy coding
  dummy.cov <- purrr::map(.x = seq_len(n.covariates),
                          .f = ~{
                            if(covariate.types[.x] == "f"){
                              
                              fastDummies::dummy_cols(.data = covariates.df, 
                                                      select_columns = covariate.names[.x],
                                                      remove_first_dummy = FALSE, 
                                                      remove_selected_columns = TRUE) %>% 
                                dplyr::select(-tidyselect::any_of(covariate.names))
                              
                            } else { covariates.df[, .x, drop = FALSE] }}) %>% 
    purrr::set_names(., covariate.names)
  
  fL <- sapply(covariate.names, FUN = function(covname){
    
    # Number of factor levels
    nL <- nlevels(covariates.df[, covname])
    
    # Number of levels other than baseline
    if(nL == 0) nparam <- 1 else nparam <- nL - 1
    
    index <- lapply(X = nL, FUN = function(x){
      if(x == 0) res <- 1
      if(x == 2) res <- 2
      if(x > 2) res <- seq(from = x - 1, to = x)
      return(res)})
    
    return(list(nL = nL, nparam = nparam, index = index[[1]]))
    
  }, simplify = FALSE, USE.NAMES = TRUE)
  
  I.covariates <- purrr::map(.x = covariate.names, 
                             .f = ~ dummy.cov[[.x]][, fL[[.x]]$index, drop = FALSE])
  
  I.covariates <- purrr::map(.x = I.covariates,
                             .f = ~apply(.x, 2, list) %>% purrr::flatten(.)) %>% 
    purrr::flatten(.)
  
  for (j in 1:n.covariates) priors[[names(I.covariates)[j]]] <- c(0, 30) 
  
  # Translate to probit scale
  covariates.biphasic <- 
    purrr::map(.x = covariate.names,
               .f = ~{sapply(X = seq_along(covariates[[.x]]),
                             FUN = function(ind){
                               if(ind %in% fL[[.x]]$index){
                                 qnorm(pnorm(covariates[[.x]][ind], mean = priors[[.x]][1], sd = priors[[.x]][2]), mean = 0, sd = 1) } else { 0 }} )}) %>% 
    purrr::set_names(., covariate.names)
  
  psi_ij <- psi_i[whale.id]
  
  if(any(covariate.names == "exposed")) psi_ij <-  psi_ij + covariates.df[, "exposed"] * covariates.biphasic[["exposed"]] 
  
  if(any(covariate.names == "sonar")) psi_ij <-  psi_ij + Reduce(f = "+", lapply(X = seq_len(nlevels(covariates.df[, "sonar"])), FUN = function(fl){dummy.cov[["sonar"]][, fl] * covariates.biphasic[["sonar"]][fl]}))
  
  if(any(covariate.names == "behaviour")) psi_ij <-  psi_ij + Reduce(f = "+", lapply(X = seq_len(nlevels(covariates.df[, "behaviour"])), FUN = function(fl){dummy.cov[["behaviour"]][, fl] * covariates.biphasic[["behaviour"]][fl]}))
  
  if(any(covariate.names == "range")) psi_ij <-  psi_ij + covariates.df[, "range"] * covariates.biphasic[["range"]] 
  
  pi_ij <- pnorm(q = psi_ij)
  
} else {
  
  pi_ij <- pnorm(q = psi_i[whale.id])  
  
}


# Choose which part of the DR curve each whale responds to, with prob pi_ij
k_ij <- (1 - rbinom(n = n.trials, size = 1, prob = pi_ij)) + 1

# Simulate for each trial the value of the lower and upper mixture components
mu_ij <- matrix(data = 0, nrow = n.trials, ncol = 2)

mu_ij[, 1] <- rtnorm(n = n.trials, location = nu[1], scale = tau[1], L = L, U = alpha)
mu_ij[, 2] <- rtnorm(n = n.trials, location = nu[2], scale = tau[2], L = alpha, U = U)

t_ij <- numeric(n.trials)
for (u in 1:n.trials) t_ij[u] <- mu_ij[u, k_ij[u]]

y_ij <- rnorm(n = n.trials, mean = t_ij, sd = obs.sd)

# Right-censoring
max.dose <- runif(n = n.trials, min = Rc[1], Rc[2])
is.censored <- ifelse(y_ij > max.dose, 1, 0)
y_ij[is.censored == 1] <- NA
max.dose[is.censored == 0] <- L
rcensoring <- sum(is.censored) > 0

# MCMC ----------------------------------------------------------

n.burn.in <- 15000
n.iter <- 50000
tot.iter <- n.iter + n.burn.in
n.thin <- 1
n.chains <- 1

# Metropolis Hastings -----------------------------------------------------

# | Set up ----

# Standard deviations for the proposal distributions used in the Metropolis algorithm
prop.sd <- list(nu = 5, tau = 5, alpha = 2, mu.ij = 5, psi = 1, omega = 2, psi.i = 2)
if(n.covariates > 0){ 
  for (j in 1:n.covariates) prop.sd[[covariate.names[j]]] <- 10
  prop.sd[["range"]] <- 3}

# Sample matrices and vectors
mcmat <- list()

mcmat$t.ij <- mcmat$y.ij <- matrix(data = 0, nrow = (n.iter + n.burn.in), ncol = n.trials)
colnames(mcmat$t.ij) <- paste0("t.ij.", seq_len(n.trials))
colnames(mcmat$y.ij) <- paste0("y.ij.", seq_len(n.trials))

# Could use an array but cannot set different column names - which are useful for plotting
mcmat$mu.ij <- list('1' = matrix(data = 0, nrow = (n.iter + n.burn.in), ncol = n.trials),
                    '2' = matrix(data = 0, nrow = (n.iter + n.burn.in), ncol = n.trials))
colnames(mcmat$mu.ij$`1`) <- paste0("mu.ij.", seq_len(n.trials), ".1.")
colnames(mcmat$mu.ij$`2`) <- paste0("mu.ij.", seq_len(n.trials), ".2.")

mcmat$psi.i <- matrix(data = 0, nrow = (n.iter + n.burn.in), ncol = n.whales)
colnames(mcmat$psi.i) <- paste0("psi.i.", seq_len(n.whales))

mcmat$k.ij <- matrix(data = 0, nrow = (n.iter + n.burn.in), ncol = n.trials)
colnames(mcmat$k.ij) <- paste0("k.ij.", seq_len(n.trials))

mcmat$pi.ij <- matrix(data = 0, nrow = (n.iter + n.burn.in), ncol = n.trials)
colnames(mcmat$pi.ij) <- paste0("pi.ij.", seq_len(n.trials))

mcmat$nu <- matrix(data = 0, nrow = (n.iter + n.burn.in), ncol = 2)
mcmat$tau <- matrix(data = 0, nrow = (n.iter + n.burn.in), ncol = 2)
colnames(mcmat$nu) <- c("nu.1", "nu.2")
colnames(mcmat$tau) <- c("tau.1", "tau.2")

mcmat$alpha <- mcmat$psi <- mcmat$omega <- numeric((n.iter + n.burn.in))

if(n.covariates > 0){
  for (j in covariate.names) {
    tmp <- matrix(data = 0, nrow = n.iter + n.burn.in, ncol = ifelse(fL[[j]]$nL == 0, 1, fL[[j]]$nL))
    colnames(tmp) <- colnames(dummy.cov[[j]])
    mcmat[[j]] <- tmp}
} 

mcmat$L1 <- mcmat$L2 <-  matrix(data = 0, nrow = (n.iter + n.burn.in), ncol = n.trials)

# Starting values

mcmat$L1[1, ] <- L
mcmat$L2[1, ] <- alpha

mcmat$psi.i[1, ] <- psi_i
mcmat$k.ij[1, ] <- k_ij
mcmat$pi.ij[1, ] <- pi_ij

mcmat$nu[1, ] <- nu
mcmat$tau[1, ] <- tau
mcmat$alpha[1] <- alpha
mcmat$psi[1] <- psi
mcmat$omega[1] <- omega

mcmat$mu.ij$`1`[1, ] <- mu_ij[, 1]
mcmat$mu.ij$`2`[1, ] <- mu_ij[, 2]


if(sum(is.censored) > 0){
  for(k in 1:n.trials){
    mcmat$mu.ij[[mcmat$k.ij[1, k]]][1, k] <-
      rtnorm(n = 1,
             location = mcmat$nu[1, mcmat$k.ij[1, k]],
             scale = mcmat$tau[1, mcmat$k.ij[1, k]], # Set to a high value so that no Inf are returned
             L = ifelse(is.censored[k] == 1, max.dose[k], L),
             U = ifelse(mcmat$k.ij[1, k] == 1, mcmat$alpha[1], U))
  }
}

if(sum(is.na(mcmat$mu.ij$`1`[1, ])) > 0) stop("NAs found in starting values used for mu.ij[1]")
if(sum(is.na(mcmat$mu.ij$`2`[1, ])) > 0) stop("NAs found in starting values used for mu.ij[2]")

mcmat$t.ij[1, ] <- (2 - mcmat$k.ij[1, ]) * mcmat$mu.ij$`1`[1, ] + ((1 - (2 - mcmat$k.ij[1, ])) * mcmat$mu.ij$`2`[1, ])
mcmat$y.ij[1, ] <- y_ij

if(n.covariates > 0){
  for(j in covariate.names) mcmat[[j]][1, ] <- covariates[[j]]
}

if(sum(is.censored) > 0){ mcmat$y.ij[1, is.censored == 1] <-
  rnorm(n = sum(is.censored), 
        mean = mcmat$t.ij[1, is.censored == 1],
        sd = obs.sd) }


# | Algorithm ----

pb <- progress::progress_bar$new(
  format = "[:bar] :percent eta: :eta",
  total = (n.iter + n.burn.in), clear = FALSE, width = getOption("width") + 2)

for(i in 2:(n.iter + n.burn.in)){
  
  pb$tick()
  
  #' --------------------------------------------------------
  # | -- y.ij ----
  #' --------------------------------------------------------
  
  mcmat$y.ij[i, ] <- mcmat$y.ij[i - 1, ]
  
  if (sum(is.censored) > 0) {
    mcmat$y.ij[i, is.censored == 1] <- rnorm(
      n = sum(is.censored),
      mean = mcmat$t.ij[i - 1, is.censored == 1],
      sd = obs.sd)}
  
  #' --------------------------------------------------------
  # | -- alpha ----
  #' --------------------------------------------------------
  
  proposed.alpha <- rnorm(n = 1, mean = mcmat$alpha[i - 1], sd = prop.sd$alpha)
  
  if(proposed.alpha < alpha.bounds[1] | 
     proposed.alpha > alpha.bounds[2])
    
    # Should additional conditions be used here?
    # any(c(mcmat$mu.ij$`1`[i - 1, ], mcmat$nu[i - 1, 1]) > proposed.alpha) |
    # any(c(mcmat$mu.ij$`2`[i - 1, ], mcmat$nu[i - 1, 2]) < proposed.alpha))
{
    
    accept.prob <- 0
    
  } else {
    
    loglik.proposed <- sum(dtnorm(x = mcmat$mu.ij$`1`[i - 1, ], 
                                  location = mcmat$nu[i - 1, 1], 
                                  scale = mcmat$tau[i - 1, 1], 
                                  L = L, U = proposed.alpha, log = TRUE)) + 
      
      sum(dtnorm(x = mcmat$mu.ij$`2`[i - 1, ], 
                 location = mcmat$nu[i - 1, 2], 
                 scale = mcmat$tau[i - 1, 2], 
                 L = proposed.alpha, U = U, log = TRUE))
    
    loglik.current <- sum(dtnorm(x = mcmat$mu.ij$`1`[i - 1, ], 
                                 location = mcmat$nu[i - 1, 1], 
                                 scale = mcmat$tau[i - 1, 1], 
                                 L = L, U = mcmat$alpha[i - 1], 
                                 log = TRUE)) + 
      
      sum(dtnorm(x = mcmat$mu.ij$`2`[i - 1, ], 
                 location = mcmat$nu[i - 1, 2], 
                 scale = mcmat$tau[i - 1, 2], 
                 L = mcmat$alpha[i - 1], 
                 U = U, 
                 log = TRUE))
    
    accept.prob <- exp(loglik.proposed - loglik.current)
    
  }
  
  if (runif(1) < accept.prob) mcmat$alpha[i] <- proposed.alpha else mcmat$alpha[i] <- mcmat$alpha[i - 1]
  
  #' --------------------------------------------------------
  # | -- nu1 ----
  #' --------------------------------------------------------
  
  proposed.nu1 <- rnorm(n = 1, mean = mcmat$nu[i - 1, 1], sd = prop.sd$nu)
  
  if(proposed.nu1 > mcmat$alpha[i] | proposed.nu1 < L){
    
    accept.prob <- 0
    
  } else {
    
    loglik.proposed <- sum(dtnorm(x = mcmat$mu.ij$`1`[i - 1, ],
                                  location = proposed.nu1,
                                  scale = mcmat$tau[i - 1, 1],
                                  L = L, U = mcmat$alpha[i], log = TRUE))
    
    loglik.current <- sum(dtnorm(x = mcmat$mu.ij$`1`[i - 1, ],
                                 location = mcmat$nu[i - 1, 1],
                                 scale = mcmat$tau[i - 1, 1],
                                 L = L, U = mcmat$alpha[i], log = TRUE))
    
    accept.prob <- exp(loglik.proposed - loglik.current)
    
  }
  if (runif(1) < accept.prob) mcmat$nu[i, 1] <- proposed.nu1 else mcmat$nu[i, 1] <- mcmat$nu[i - 1, 1]
  
  #' --------------------------------------------------------
  # | -- nu2 ----
  #' --------------------------------------------------------
  
  proposed.nu2 <- rnorm(n = 1, mean = mcmat$nu[i - 1, 2], sd = prop.sd$nu)
  
  if(proposed.nu2 < mcmat$alpha[i] | proposed.nu2 > U){
    
    accept.prob <- 0
    
  } else {
    
    loglik.proposed <- sum(dtnorm(x = mcmat$mu.ij$`2`[i - 1, ],
                                  location = proposed.nu2,
                                  scale = mcmat$tau[i - 1, 2],
                                  L = mcmat$alpha[i], U = U, log = TRUE))
    
    loglik.current <- sum(dtnorm(x = mcmat$mu.ij$`2`[i - 1, ],
                                 location = mcmat$nu[i - 1, 2],
                                 scale = mcmat$tau[i - 1, 2],
                                 L = mcmat$alpha[i], U = U, log = TRUE))
    
    accept.prob <- exp(loglik.proposed - loglik.current)
    
  }
  if (runif(1) < accept.prob) mcmat$nu[i, 2] <- proposed.nu2 else mcmat$nu[i, 2] <- mcmat$nu[i - 1, 2]
  
  #' --------------------------------------------------------
  # | -- tau 1 ----
  #' --------------------------------------------------------
  
  proposed.tau1 <- rnorm(n = 1, mean = mcmat$tau[i - 1, 1], sd = prop.sd$tau)
  
  loglik.proposed <- sum(dtnorm(x = mcmat$mu.ij$`1`[i - 1, ],
                                location = mcmat$nu[i, 1],
                                scale = proposed.tau1,
                                L = L, U = mcmat$alpha[i], log = TRUE))
  
  loglik.current <- sum(dtnorm(x = mcmat$mu.ij$`1`[i - 1, ],
                               location = mcmat$nu[i, 1],
                               scale = mcmat$tau[i - 1, 1],
                               L = L, U = mcmat$alpha[i], log = TRUE))
  
  logprior.proposed <- uniform_prior_bi(proposed.tau1, lower = tau.bounds[1], upper = tau.bounds[2])
  logprior.current <- uniform_prior_bi(mcmat$tau[i - 1, 1], lower = tau.bounds[1], upper = tau.bounds[2])
  
  accept.prob <- exp((loglik.proposed + logprior.proposed) - (loglik.current + logprior.current))
  if (runif(1) < accept.prob) mcmat$tau[i, 1] <- proposed.tau1 else mcmat$tau[i, 1] <- mcmat$tau[i - 1, 1]
  
  #' --------------------------------------------------------
  # | -- tau 2 ----
  #' --------------------------------------------------------
  
  proposed.tau2 <- rnorm(n = 1, mean = mcmat$tau[i - 1, 2], sd = prop.sd$tau)
  
  loglik.proposed <- sum(dtnorm(x = mcmat$mu.ij$`2`[i - 1, ],
                                location = mcmat$nu[i, 2],
                                scale = proposed.tau2,
                                L = mcmat$alpha[i], U = U, log = TRUE))
  
  loglik.current <- sum(dtnorm(x = mcmat$mu.ij$`2`[i - 1, ],
                               location = mcmat$nu[i, 2],
                               scale = mcmat$tau[i - 1, 2],
                               L = mcmat$alpha[i], U = U, log = TRUE))
  
  logprior.proposed <- uniform_prior_bi(proposed.tau2, tau.bounds[1], tau.bounds[2])
  logprior.current <- uniform_prior_bi(mcmat$tau[i - 1, 2], tau.bounds[1], tau.bounds[2])
  
  accept.prob <- exp((loglik.proposed + logprior.proposed) - (loglik.current + logprior.current))
  if (runif(1) < accept.prob) mcmat$tau[i, 2] <- proposed.tau2 else mcmat$tau[i, 2] <- mcmat$tau[i - 1, 2]
  
  #' --------------------------------------------------------
  # | -- omega ----
  #' --------------------------------------------------------
  
  proposed.omega <- rnorm(n = 1, mean = mcmat$omega[i - 1], sd = prop.sd$omega)
  
  if(proposed.omega <  omega.bounds[1] | proposed.omega > omega.bounds[2]){
    
    accept.prob <- 0
    
  } else {
    
    loglik.proposed <- sum(dnorm(x = mcmat$psi.i[i - 1, ], mean = mcmat$psi[i - 1],
                                 sd = proposed.omega, log = TRUE))
    
    loglik.current <- sum(dnorm(x = mcmat$psi.i[i - 1, ], mean = mcmat$psi[i - 1],
                                sd = mcmat$omega[i - 1], log = TRUE))
    
    accept.prob <- exp(loglik.proposed - loglik.current)
  }
  
  if (runif(1) < accept.prob) mcmat$omega[i] <- proposed.omega else mcmat$omega[i] <- mcmat$omega[i - 1]
  
  #' --------------------------------------------------------
  # | -- psi ----
  #' --------------------------------------------------------
  
  proposed.psi <- rnorm(n = 1, mean = mcmat$psi[i - 1], sd = prop.sd$psi)
  
  loglik.proposed <- sum(dnorm(x = mcmat$psi.i[i - 1, ], mean = proposed.psi, 
                               sd = mcmat$omega[i], log = TRUE))
  
  loglik.current <- sum(dnorm(x = mcmat$psi.i[i - 1, ], mean = mcmat$psi[i - 1], 
                              sd = mcmat$omega[i], log = TRUE))
  
  logprior.proposed <- sum(dnorm(x = proposed.psi, mean = priors$psi[1], sd = priors$psi[2], log = TRUE))
  logprior.current <- sum(dnorm(x = mcmat$psi[i - 1], mean = priors$psi[1], sd = priors$psi[2], log = TRUE))
  
  accept.prob <- exp((loglik.proposed + logprior.proposed) - (loglik.current + logprior.current))
  if (runif(1) < accept.prob) mcmat$psi[i] <- proposed.psi else mcmat$psi[i] <- mcmat$psi[i - 1]
  
  #' --------------------------------------------------------
  # | -- mu.ij ----
  #' --------------------------------------------------------
  
  L.1 <- rep(L, n.trials)
  L.2 <- rep(mcmat$alpha[i], n.trials)
  
  # When right-censoring occcurs, the lower bounds of the proposal distribution must be adjusted
  
  ind.1 <- which(is.censored == 1)[which(mcmat$k.ij[i - 1, which(is.censored == 1)] == 1)]
  ind.2 <- which(is.censored == 1)[which(mcmat$k.ij[i - 1, which(is.censored == 1)] == 2)]

  L.1[ind.1] <- purrr::map_dbl(.x = max.dose[ind.1], .f = ~max(.x, L))
  L.2[ind.2] <- purrr::map_dbl(.x = max.dose[ind.2], .f = ~max(.x, mcmat$alpha[i]))
  
  # L.1[ind.1] <- max.dose[ind.1]
  # L.2[ind.2] <- max.dose[ind.2]
  
  mcmat$L1[i, ] <- L.1
  mcmat$L2[i, ] <- L.2
  
  mu_ij_1 <- rtnorm(n = n.trials,
                    location = mcmat$mu.ij$`1`[i - 1, ],
                    scale = prop.sd$mu.ij,
                    L = L.1, 
                    U = mcmat$alpha[i])
  
  mu_ij_2 <- rtnorm(n = n.trials,
                    location = mcmat$mu.ij$`2`[i - 1, ],
                    scale = prop.sd$mu.ij,
                    L = L.2, 
                    U = U)
  
  if(sum(is.na(mu_ij_1)) > 0) stop("NA in mu_ij1")
  if(sum(is.na(mu_ij_2)) > 0) stop("NA in mu_ij2")
  
  proposed.t.ij <- (2 - mcmat$k.ij[i - 1, ]) * mu_ij_1 + ((1 - (2 - mcmat$k.ij[i - 1, ])) * mu_ij_2)
  
  loglik.proposed <- loglik.current <- vector(mode = "list", length = 2)
  
  # The likelihood is evaluated using the standard bounds
  loglik.proposed[[1]] <- unname(dtnorm(x = mu_ij_1, 
                                        location = mcmat$nu[i, 1], 
                                        scale = mcmat$tau[i, 1],
                                        L = L, 
                                        U = mcmat$alpha[i], 
                                        log = TRUE))
  
  loglik.proposed[[2]] <- unname(dtnorm(x = mu_ij_2,
                                        location = mcmat$nu[i, 2],
                                        scale = mcmat$tau[i, 2],
                                        L = mcmat$alpha[i],
                                        U = U, 
                                        log = TRUE))
  
  loglik.proposed <- lapply(loglik.proposed, function(x) x + dnorm(x = mcmat$y.ij[i, ],
                                                                   mean = proposed.t.ij,
                                                                   sd = obs.sd,
                                                                   log = TRUE))
  
  loglik.current[[1]] <- unname(dtnorm(x = mcmat$mu.ij$`1`[i - 1, ], 
                                       location = mcmat$nu[i, 1], 
                                       scale = mcmat$tau[i, 1],
                                       L = L, 
                                       U = mcmat$alpha[i], 
                                       log = TRUE))
  
  loglik.current[[2]] <- unname(dtnorm(x = mcmat$mu.ij$`2`[i - 1, ],
                                       location = mcmat$nu[i, 2],
                                       scale = mcmat$tau[i, 2],
                                       L = mcmat$alpha[i],
                                       U = U, 
                                       log = TRUE))
  
  
  loglik.current <- lapply(loglik.current, function(x) x + dnorm(x = mcmat$y.ij[i, ],
                                                                 mean = mcmat$t.ij[i - 1, ],
                                                                 sd = obs.sd,
                                                                 log = TRUE))
  
  
  logprop.forward <- list()
  logprop.forward[[1]] <- dtnorm(x = mu_ij_1, 
                                 location = mcmat$mu.ij$`1`[i - 1, ],
                                 scale = prop.sd$mu.ij, 
                                 L = L.1, 
                                 U = mcmat$alpha[i], log = TRUE)
  
  logprop.forward[[2]] <- dtnorm(x = mu_ij_2, 
                                 location = mcmat$mu.ij$`2`[i - 1, ],
                                 scale = prop.sd$mu.ij, 
                                 L = L.2, 
                                 U = U, log = TRUE)
  
  logprop.backward <- list()
  logprop.backward[[1]] <- dtnorm(x = mcmat$mu.ij$`1`[i - 1, ], 
                                  location = mu_ij_1,
                                  scale = prop.sd$mu.ij, 
                                  L = L.1, 
                                  U = mcmat$alpha[i], log = TRUE)
  
  logprop.backward[[2]] <- dtnorm(x = mcmat$mu.ij$`2`[i - 1, ], 
                                  location = mu_ij_2,
                                  scale = prop.sd$mu.ij, 
                                  L = L.2, 
                                  U = U, log = TRUE)
  
  accept.prob <- runif(1) < exp((loglik.proposed[[1]] + logprop.backward[[1]]) - 
                                  (loglik.current[[1]] + logprop.forward[[1]]))
  
  mcmat$mu.ij$`1`[i, ] <- mcmat$mu.ij$`1`[i - 1, ]
  mcmat$mu.ij$`1`[i, accept.prob] <- mu_ij_1[accept.prob]
  
  accept.prob <- runif(1) < exp((loglik.proposed[[2]] + logprop.backward[[2]]) - 
                                  (loglik.current[[2]] + logprop.forward[[2]]))
  
  mcmat$mu.ij$`2`[i, ] <- mcmat$mu.ij$`2`[i - 1, ]
  mcmat$mu.ij$`2`[i, accept.prob] <- mu_ij_2[accept.prob]
  
  mcmat$t.ij[i, ] <- ((2 - mcmat$k.ij[i - 1, ]) * mcmat$mu.ij$`1`[i, ]) +
    ((1 - (2 - mcmat$k.ij[i - 1, ])) * mcmat$mu.ij$`2`[i, ])
  
  #' --------------------------------------------------------
  # | -- psi.i ----
  #' --------------------------------------------------------
  
  proposed.psi.i <- rnorm(n = n.whales, mean = mcmat$psi.i[i - 1, ], sd = prop.sd$psi.i)
  
  proposed.psi.ij <- Reduce("+", append(list(proposed.psi.i[whale.id]),
                                        lapply(X = covariate.names, 
                                               FUN = function(k){
                                                 apply(t(t(dummy.cov[[k]]) * 
                                                           qnorm(pnorm(mcmat[[k]][i - 1, ], 
                                                                       mean = priors[[k]][1],
                                                                       sd = priors[[k]][2]))), 1, sum)})))
  
  current.psi.ij <- Reduce("+", append(list(mcmat$psi.i[i - 1, whale.id]),
                                       lapply(X = covariate.names, 
                                              FUN = function(k){
                                                apply(t(t(dummy.cov[[k]]) * 
                                                          qnorm(pnorm(mcmat[[k]][i - 1, ], 
                                                                      mean = priors[[k]][1],
                                                                      sd = priors[[k]][2]))), 1, sum)})))
  
  loglik.proposed <- dnorm(x = proposed.psi.i,
                           mean = mcmat$psi[i],
                           sd = mcmat$omega[i],
                           log = TRUE) +
    
    purrr::map_dbl(.x = seq_len(n.whales),
                   .f = ~ sum(d_binom(x = 2 - mcmat$k.ij[i - 1, ], size = 1,
                                      prob = pnorm(q = proposed.psi.ij), log = TRUE)[whale.id == .x]))
  
  loglik.current <- dnorm(x = mcmat$psi.i[i - 1, ],
                          mean = mcmat$psi[i],
                          sd = mcmat$omega[i],
                          log = TRUE) +
    
    purrr::map_dbl(.x = seq_len(n.whales),
                   .f = ~ sum(d_binom(x = 2 - mcmat$k.ij[i - 1, ], size = 1,
                                      prob = pnorm(q = current.psi.ij), log = TRUE)[whale.id == .x]))
  
  accept.prob <- runif(1) < exp(loglik.proposed - loglik.current)
  mcmat$psi.i[i, ] <- mcmat$psi.i[i - 1, ]
  mcmat$psi.i[i, accept.prob] <- proposed.psi.i[accept.prob]
  
  mcmat$pi.ij[i, ] <- pnorm(Reduce("+", append(list(mcmat$psi.i[i, whale.id]),
                                               lapply(X = covariate.names, 
                                                      FUN = function(k){
                                                        apply(t(t(dummy.cov[[k]]) * 
                                                                  qnorm(pnorm(mcmat[[k]][i - 1, ], 
                                                                              mean = priors[[k]][1],
                                                                              sd = priors[[k]][2]))), 1, sum)}))))
  
  #' --------------------------------------------------------
  # | -- covariates ----
  #' --------------------------------------------------------
  
  if(n.covariates > 0){
    
    for (a in covariate.names) {
      
      covariate.proposed <- covariate.current <- 
        sapply(X = covariate.names, FUN = function(x) mcmat[[x]][i - 1, ], USE.NAMES = TRUE)
      
      proposed.value <- rnorm(n = fL[[a]]$nparam, 
                              mean = mcmat[[a]][i - 1, fL[[a]]$index], 
                              sd = prop.sd[[a]])
      
      if(fL[[a]]$nL >= 2) proposed.value <- c(0, proposed.value)
      covariate.proposed[[a]] <- proposed.value
      
      proposed.psi.ij <- Reduce("+", append(list(mcmat$psi.i[i, whale.id]),
                                            lapply(X = covariate.names, FUN = function(k){
                                              apply(t(t(dummy.cov[[k]]) * 
                                                        qnorm(pnorm(covariate.proposed[[k]], 
                                                                    mean = priors[[k]][1], 
                                                                    sd = priors[[k]][2]))), 1, sum)})))
      
      current.psi.ij <- Reduce("+", append(list(mcmat$psi.i[i, whale.id]),
                                           lapply(X = covariate.names, FUN = function(k){
                                             apply(t(t(dummy.cov[[k]]) * 
                                                       qnorm(pnorm(covariate.current[[k]], 
                                                                   mean = priors[[k]][1], 
                                                                   sd = priors[[k]][2]))), 1, sum)})))
      
      
      loglik.proposed <- sum(d_binom(x = 2 - mcmat$k.ij[i - 1, ], size = 1,
                                     prob = pnorm(q = proposed.psi.ij), log = TRUE))
      
      loglik.current <- sum(d_binom(x = 2 - mcmat$k.ij[i - 1, ], size = 1,
                                    prob = pnorm(q = current.psi.ij), log = TRUE))
      
      logprior.proposed <- sum(dnorm(x = proposed.value, 
                                     mean = priors[[a]][1], 
                                     sd = priors[[a]][2],
                                     log = TRUE))
      
      logprior.current <- sum(dnorm(x = mcmat[[a]][i - 1, ], 
                                    mean = priors[[a]][1], 
                                    sd = priors[[a]][2],
                                    log = TRUE))
      
      accept.prob <- exp((loglik.proposed + logprior.proposed) - (loglik.current + logprior.current))
      
      if (runif(1) < accept.prob){
        mcmat[[a]][i, ] <- proposed.value 
      } else {mcmat[[a]][i, ] <- mcmat[[a]][i - 1, ] }
    }
    
    mcmat$pi.ij[i, ] <- pnorm(Reduce("+", append(list(mcmat$psi.i[i, whale.id]),
                                                 lapply(X = covariate.names, 
                                                        FUN = function(k){
                                                          apply(t(t(dummy.cov[[k]]) * 
                                                                    qnorm(pnorm(mcmat[[k]][i, ], 
                                                                                mean = priors[[k]][1],
                                                                                sd = priors[[k]][2]))), 1, sum)}
                                                 ))))
    
    
  }
  
  #' --------------------------------------------------------
  # | -- k.ij ----
  #' --------------------------------------------------------
  
  #(1) Implementation with an arbitrary prob of 0.5 as an equivalent "proposal SD"
  
  proposed.k.ij <- (1 - rbinom(n = n.trials, size = 1, prob = 0.5)) + 1
  
  proposed.t.ij <- purrr::map_dbl(.x = seq_len(n.trials),
                                  .f = ~ mcmat$mu.ij[[proposed.k.ij[.x]]][i, .x])
  
  loglik.proposed <- dnorm(x = mcmat$y.ij[i, ], mean = proposed.t.ij, sd = obs.sd, log = TRUE) +
    d_binom(x = 2 - proposed.k.ij, size = 1, prob = mcmat$pi.ij[i, ], log = TRUE)
  
  loglik.current <- dnorm(x = mcmat$y.ij[i, ], mean = mcmat$t.ij[i, ], sd = obs.sd, log = TRUE) +
    d_binom(x = 2 - mcmat$k.ij[i - 1, ], size = 1, prob = mcmat$pi.ij[i, ], log = TRUE)
  
  logprop.forward <- d_binom(x = 2 - proposed.k.ij, size = 1, prob = 0.5, log = TRUE)
  logprop.backward <- d_binom(x = 2 - mcmat$k.ij[i - 1, ], size = 1, prob = 0.5, log = TRUE)
  
  accept.prob <- runif(1) < exp((loglik.proposed + logprop.backward) - (loglik.current + logprop.forward))
  

  # Original implementation
  # proposed.k.ij <- (1 - rbinom(n = n.trials, size = 1, prob = mcmat$pi.ij[i, ])) + 1
  # 
  # proposed.t.ij <- purrr::map_dbl(.x = seq_len(n.trials),
  #                                 .f = ~ mcmat$mu.ij[[proposed.k.ij[.x]]][i, .x])
  # 
  # loglik.proposed <- dnorm(x = mcmat$y.ij[i, ], mean = proposed.t.ij, sd = obs.sd, log = TRUE) +
  #   d_binom(x = 2 - proposed.k.ij, size = 1, prob = mcmat$pi.ij[i, ], log = TRUE)
  # 
  # loglik.current <- dnorm(x = mcmat$y.ij[i, ], mean = mcmat$t.ij[i, ], sd = obs.sd, log = TRUE) +
  #   d_binom(x = 2 - mcmat$k.ij[i - 1, ], size = 1, prob = mcmat$pi.ij[i, ], log = TRUE)
  # 
  # logprop.forward <- d_binom(x = 2 - proposed.k.ij, size = 1, prob = mcmat$pi.ij[i, ], log = TRUE)
  # logprop.backward <- d_binom(x = 2 - mcmat$k.ij[i - 1, ], size = 1, prob = mcmat$pi.ij[i, ], log = TRUE)
  # 
  # accept.prob <- runif(1) < exp((loglik.proposed + logprop.backward) - (loglik.current + logprop.forward))
  
  # This check is done to prevent the sampler from accepting values of k_ij that would lead to
  # numerical problems with proposals for mu_ij at next iteration when right-censoring occurs
  # e.g., a right-censored observation assigned to the context-dependent phase, leading to 
  # L.1 > mcmat$alpha[i] and generating NAs
  rc.check <- proposed.t.ij[is.censored == 1] < max.dose[is.censored == 1]
  accept.prob[is.censored == 1][which(rc.check == 1)] <- FALSE
  
  mcmat$k.ij[i, ] <- mcmat$k.ij[i - 1, ]
  mcmat$k.ij[i, accept.prob] <- proposed.k.ij[accept.prob]
  
  mcmat$t.ij[i, ] <- purrr::map_dbl(.x = seq_len(n.trials),
                                    .f = ~ mcmat$mu.ij[[mcmat$k.ij[i, .x]]][i, .x])
  
  if(!sum(mcmat$t.ij[i, ]) == sum((2 - mcmat$k.ij[i, ]) * mcmat$mu.ij$`1`[i, ] + ((1 - (2 - mcmat$k.ij[i, ])) * mcmat$mu.ij$`2`[i, ]))) stop("Sums not equal")
  
}

# rjags ----------------------------------------------------------

# Define model

cat("\n")
modelstring <- "model{

#-------------------------------------
# PRIORS
#-------------------------------------

# Mean/variance of lower mixture component
nu[1] ~ dunif(L, alpha)
tau[1] ~ dunif(tau.bounds[1], tau.bounds[2])
inv_tau2[1] <- pow(tau[1], -2)

# Mean/variance of upper mixture component
nu[2] ~ dunif(alpha, U)
tau[2] ~ dunif(tau.bounds[1], tau.bounds[2])
inv_tau2[2] <- pow(tau[2], -2)

# Mixture probability (aka probability of exhibiting context-dependent response)
# Note: on probit scale
psi ~ dnorm(psi.prior[1], psi.prior[2])

# Between-whale variance in context-dependent prob
omega ~ dunif(omega.bounds[1], omega.bounds[2])
inv_omega2 <- pow(omega, -2)

# Breakpoint between lower and upper mixture components
alpha ~ dunif(alpha.bounds[1], alpha.bounds[2])

# Priors on covariates

#-------------------------------------
# PROCESS MODEL
#-------------------------------------

# Between whales ------------------

for(i in 1:n.whales){
  # Context-dependent probabilities
  psi_i[i] ~ dnorm(psi, inv_omega2)
}

# Within whales, between trials ------------------

for(j in 1:n.trials){

  mu_ij[j,1] ~ dnorm(nu[1], inv_tau2[1]) T(L, alpha)
  mu_ij[j,2] ~ dnorm(nu[2], inv_tau2[2]) T(alpha, U)
  
  psi_ij[j] <- psi_i[whale.id[j]] #/#
  
  pi_ij[j, 1] <- pnorm(psi_ij[j], 0, 1)
  pi_ij[j, 2] <- 1 - pi_ij[j, 1]
  pi[j] <- pi_ij[j, 1] # Simply for density plot comparisons
  
  k_ij[j] ~ dcat(pi_ij[j, ])
  I.censored[j] ~ dinterval(mu_ij[j, k_ij[j]], Rc[j])
  
}

#-------------------------------------
# OBSERVATION MODEL
#-------------------------------------
for(k in 1:n.trials){
  t_ij[k] <- mu_ij[k, k_ij[k]] 
  y_ij[k] ~ dnorm(t_ij[k], obs.precision)
}

}"

if(n.covariates > 0){
  
  modelstring <- gsub(
    pattern = "# Priors on covariates",
    replacement = paste0(paste0(names(I.covariates),
                                " ~ dnorm(", paste0(names(I.covariates), ".prior.mean"),
                                ", 1/(", paste0(names(I.covariates), ".prior.sd"), "^2))",
                                sep = "\n "),
                         paste0(paste0("beta_", names(I.covariates)),
                                " <- qnorm(pnorm(", 
                                paste0(names(I.covariates)), ", ", 
                                paste0(names(I.covariates), ".prior.mean, "), 
                                "1/(", paste0(names(I.covariates), ".prior.sd"), "^2))", 
                                ", 0, 1)", sep = "\n "), collapse = "\n"),
    x = modelstring)
  
  modelstring <- gsub(
    pattern = "#/#",
    replacement = paste0("+ ", paste0("beta_", names(I.covariates), " * I.", 
                                      names(I.covariates), "[j]", collapse = " + ")), x = modelstring)
}

# dose.jags <- max.dose
# dose.jags[is.censored == 0] <- U

dose.jags <- rep(U, n.trials)
dose.jags[is.censored == 1] <- max.dose[is.censored == 1]

# Data list
data <- list(
  n.trials = n.trials, 
  n.whales = n.whales,
  L = L, 
  U = U,
  tau.bounds = tau.bounds,
  alpha.bounds = alpha.bounds,
  omega.bounds = omega.bounds,
  whale.id = whale.id,
  y_ij = y_ij, 
  Rc = dose.jags,
  I.censored = is.censored,
  psi.prior = c(priors$psi[1], 1/(priors$psi[2]^2)),
  obs.precision = 1 / (obs.sd^2))

if(n.covariates > 0){
  for(nc in names(I.covariates)){
    data[[paste0(nc, ".prior.mean")]] <- priors[[nc]][1]
    data[[paste0(nc, ".prior.sd")]] <- priors[[nc]][2]
    data[[paste0("I.", nc)]] <- I.covariates[[nc]]}}

# Initial values - k_ij needed when censoring occurs
inits <- list(alpha = alpha, nu = nu, tau = tau, psi = psi, mu_ij = mu_ij, psi_i = psi_i, k_ij = k_ij)

if(sum(is.censored) > 0){
  for(cc in which(is.censored == 1)){
    inits$mu_ij[cc, k_ij[cc]] <- rtnorm(n = 1, location = nu[k_ij[cc]], scale = tau[k_ij[cc]], 
                                        L = dose.jags[cc],
                                        U = ifelse(k_ij[cc] == 1, alpha, U))
  } 
}

# Get the model loaded
m <- jags.model(file = textConnection(modelstring), 
                data = data, inits = inits, n.chains = n.chains, n.adapt = 1000)

# Burn-in
update(m, n.burn.in)

cov.to.monitor <- c("alpha", "nu", "tau", "omega", "psi", "mu_ij", "psi_i", "k_ij", "y_ij", "t_ij", "pi")
if(n.covariates > 0) cov.to.monitor <- c(cov.to.monitor, names(I.covariates))

# Samples for inference
samples <- coda.samples(model = m, variable.names = cov.to.monitor, n.iter = n.iter, thin = n.thin)

# Density plots -----------------------------------------------------------

densplot(params = c("alpha", 
                    "nu", 
                    "tau", 
                    "omega",
                    "psi",
                    "mu_ij", 
                    "psi_i",
                    "pi_ij",
                    "k_ij",
                    "y_ij",
                    "t_ij", 
                    covariate.names),
         show.legend = TRUE,
         out.dir = getwd(),
         save.pdf = FALSE,
         file.name = "biphasic_mcmc_001")
