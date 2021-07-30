glance <- function(dat, f = "head", start = 1, end = 6, reset = TRUE) {
  options(max.print = 999999)
  purrr::map(.x = dat, .f = ~ {
    if(is.list(.x)) {
      lapply(X = .x, FUN = function(x) {tail(get(f)(x, end), (end - start) + 1)})
    } else {
      tail(get(f)(.x, end), (end - start) + 1)
    }
})}
    

densplot <- function(params = c("beta", "nu", "omega", "psi", "tau", "psi_i"), 
                     show.legend = FALSE, 
                     save.pdf = FALSE,
                     out.dir,
                     file.name = NULL){
  
  params <- gsub(pattern = "_", replacement = ".", x = params)
  
  # params <- gsub(pattern = "mu_ij", replacement = "mu.ij", x = params)
  # params <- gsub(pattern = "psi_i", replacement = "psi.i", x = params)
  # params <- gsub(pattern = "k_ij", replacement = "k.ij", x = params)
  # params <- gsub(pattern = "y_ij", replacement = "y.ij", x = params)
  
  col.1 <- "gray"
  col.2 <- "blue"
  
  jagsdat <- data.frame(samples[[1]]) 
  
  if(any(grepl(pattern = "beta_", x = names(jagsdat)))){
    names(jagsdat)[which(grepl(pattern = "beta_", x = names(jagsdat)))] <- 
      gsub(pattern = "beta_", replacement = "", x = names(jagsdat)[which(grepl(pattern = "beta_", x = names(jagsdat)))])
  }
  
  # print(nrow(jagsdat))
  
  ll <- purrr::map_lgl(.x = mcmat[params], .f = ~is.list(.x))
  params.list <- params[ll]
  params.nonlist <- params[!ll]
  
  mcmcdat <- do.call(cbind, mcmat[params.nonlist])
  
  for(p in params.list){
    mcmcdat <- cbind(mcmcdat, do.call(cbind, mcmat[[p]]))
  }
  
  mcmcdat <- mcmcdat[(n.burn.in + 1):(n.iter + n.burn.in), , drop = FALSE]
  # print(nrow(mcmcdat))
  
  # n <- do.call(c, purrr::map(.x = params, 
  #                            .f = ~{
  #                              if(.x %in% c("nu", "tau")) 
  #                                paste0(.x, ".", 1:2, ".") else .x}))
  # colnames(mcmcdat) <- n
  colnames(mcmcdat) <- gsub(pattern = "_", replacement = ".", x = colnames(mcmcdat))
  colnames(jagsdat) <- gsub(pattern = "_", replacement = ".", x = colnames(jagsdat))
  # colnames(jagsdat) <- purrr::map(.x = 1:length(n),
  #            .f = ~{
  #              last_char <- substr(n[.x], nchar(n[.x])-1+1, nchar(n[.x]))
  #              if(last_char == ".") substr(n[.x], 1, nchar(n[.x]) - 1) else n[.x]
  #            }) %>% unlist()
  
  
  colnames(jagsdat) <- stringr::str_remove(colnames(jagsdat), "\\.$")
  colnames(mcmcdat) <- stringr::str_remove(colnames(mcmcdat), "\\.$")
  
  if("pi.ij" %in% params) colnames(jagsdat)[which(grepl(pattern = "pi", x = colnames(jagsdat)))] <- 
    gsub(pattern = "pi.", replacement = "pi.ij.",
                               x = colnames(jagsdat)[which(grepl(pattern = "pi", x = colnames(jagsdat)))])
 
  params.col <- intersect(sort(colnames(jagsdat)), sort(colnames(mcmcdat)))
  
  mcmcdat <- mcmcdat[, params.col, drop = FALSE]
  jagsdat <- jagsdat[, params.col, drop = FALSE]
  
  if(is.null(file.name)) file.name <- "biphasic_densityplots.pdf" 
  if(save.pdf){
    
    pdf(file.path(out.dir, paste0(file.name, ".pdf")))
    
    par(mfrow = c(1, 1))
    plot.new()
    y.seq <- seq(1, 0.04, length.out = 19)
    text(x = 0.1, y = y.seq[1], paste0("Random seed: ", seed), adj = 0)
    text(x = 0.1, y = y.seq[2], paste0("Number of whales: ", n.whales), adj = 0)
    text(x = 0.1, y = y.seq[3], paste0("Max trials per whale: ", max.trials), adj = 0)
    text(x = 0.1, y = y.seq[4], paste0("Number of trials: ", n.trials), adj = 0)
    text(x = 0.1, y = y.seq[5], paste0("Means of mixtures (nu): ", paste0(nu, collapse = "; ")), adj = 0)
    text(x = 0.1, y = y.seq[6], paste0("SD of mixtures (tau): ", paste0(tau, collapse = "; ")), adj = 0)
    text(x = 0.1, y = y.seq[7], paste0("tau bounds: ", paste0(tau.bounds, collapse = "; ")), adj = 0)
    text(x = 0.1, y = y.seq[8], paste0("psi: ", psi), adj = 0)
    text(x = 0.1, y = y.seq[9], paste0("omega: ", omega), adj = 0)
    text(x = 0.1, y = y.seq[10], paste0("omega bounds: ", paste0(omega.bounds, collapse = "; ")), adj = 0)
    text(x = 0.1, y = y.seq[11], paste0("[L, U]: ", paste0("[", L, "-", U, "]")), adj = 0)
    text(x = 0.1, y = y.seq[12], paste0("alpha: ", alpha), adj = 0)
    text(x = 0.1, y = y.seq[13], paste0("alpha bounds: ", paste0(alpha.bounds, collapse = "; ")), adj = 0)
    text(x = 0.1, y = y.seq[14], paste0("Rc: ", paste0(Rc, collapse = "; ")), adj = 0)
    text(x = 0.1, y = y.seq[15], paste0("Covariates: ", paste0(covariate.names, collapse = "; ")), adj = 0)
    text(x = 0.1, y = y.seq[16], paste0("Censored: ", sum(is.censored)), adj = 0)
    text(x = 0.1, y = y.seq[17], paste0("--------------------------- "), adj = 0)
    text(x = 0.1, y = y.seq[18], paste0("Chain length: ", format(n.iter, big.mark   = ",")), adj = 0)
    text(x = 0.1, y = y.seq[19], paste0("Burn in: ", format(n.burn.in, big.mark   = ",")), adj = 0)
  }
  
  par(mfrow = c(3, 4))
  
  for(pp in params.col){
    if(grepl(pattern = "k.ij", x = pp)){
      jags.tb <- c(sum(jagsdat[, pp]==1), sum(jagsdat[, pp]==2))
      mcmc.tb <- c(sum(mcmcdat[, pp]==1), sum(mcmcdat[, pp]==2))
      barplot(t(cbind(jags.tb, mcmc.tb)), beside = TRUE, main = pp,
              col = c(col.1, col.2), names.arg = 1:2)
    } else {
      d1 <- density(jagsdat[, pp], adjust = 2)
      d2 <- density(mcmcdat[, pp], adjust = 2)
      xlimits <- c(min(c(d1$x, d2$x)), max(c(d1$x, d2$x)))
      ylimits <- c(0, max(c(d1$y, d2$y)))
      plot(d1, main = pp, 
           xlim = xlimits, ylim = ylimits, col = col.1)
      lines(d2, main = pp, xlim = xlimits, col = col.2)
    }

  }
  
  if(show.legend) {plot.new(); legend("topleft", c("jags", "R"), lty = c(1, 1), col = c(col.1, col.2))}
  par(mfrow = c(1, 1))
  if(save.pdf) dev.off()
}

p.biphasic <- function(dose, pi, nu, tau, L, U, B) {
  p1 <- truncnorm::ptruncnorm(q = dose, a = L, b = B, mean = nu[1], sd = tau[1])
  p2 <- truncnorm::ptruncnorm(q = dose, a = B, b = U, mean = nu[2], sd = tau[2])
  return(pi * p1 + (1 - pi) * p2)
}

dtnorm <- function(x, location = 0, scale = 1, log = FALSE, L = -Inf, U = Inf) {
  
  suppressWarnings(d <- dnorm(x, location, scale, log = TRUE))
  if(all(is.na(d))) {
    d <- rep(-100000, length(d))
  } else {
  denom <- log(pnorm(U, location, scale) - pnorm(L, location, scale))
  d <- d - denom
  d[(x < L) | (x > U)] <- -100000 # When input quantile is outside bounds
  d[is.infinite(d)] <- -100000 # When input location is outside bounds
  if(!log) d <- exp(d)
  }
  return(d)}

rtnorm <- function(n, location, scale, L, U){

  res <- location + scale * qnorm(pnorm(L, location, scale) + runif(n)*(pnorm(U, location, scale) - pnorm(L, location, scale)))
  res[L > U] <- NA
  return(res)
  }

d_binom <- function(x, size, prob, log){
  d <- dbinom(x = x, size = size, prob = prob, log = log)
  d[is.infinite(d)] <- -100000 # To avoid Inf that cause numerical issues
  return(d)
}

uniform_prior_bi <- function(param = NULL, lower, upper) {
  loglik.unif <- dunif(x = param, min = lower, max = upper, log = TRUE)
  if (any(abs(loglik.unif) == Inf)) loglik.unif <- -100000
  return(sum(loglik.unif))
}


extract_values <- function(dat1, dat2, k){
  res <- cbind(dat1, dat2)
  purrr::map_dbl(.x = seq_len(n.trials), .f = ~res[.x, k[.x]])
}




#'--------------------------------------------------------------------
# Function to find the HEX colour code corresponding to an input colour 
# with a set opacity level (i.e. emulate transparency)
#'--------------------------------------------------------------------

hexa2hex <- function(input.colour, 
                     opacity, 
                     bg.colour = "white"){
  
  #'---------------------------------------------
  # PARAMETERS
  #'---------------------------------------------
  #' @param input.colour Initial colour.
  #' @param opacity Desired level of transparency (number between 0 and 1).
  #' @param bg.colour Colour of the background. Defaults to 'white'.
  #'---------------------------------------------
  
  # White background
  
  bg <- grDevices::col2rgb(bg.colour, alpha = FALSE)
  
  # Convert input colour to RGB
  
  rgbcol <- grDevices::col2rgb(input.colour, alpha = FALSE)
  
  # Calculate red, green, blue values corresponding to input colour at chosen transparency level
  
  rc <- (1 - opacity) * bg[1,] + opacity * rgbcol[1,]
  gc <- (1 - opacity) * bg[2,] + opacity * rgbcol[2,]
  bc <- (1 - opacity) * bg[3,] + opacity * rgbcol[3,]
  
  # Convert back to hex
  
  rgb2hex <- function(r,g,b) rgb(r, g, b, maxColorValue = 255)
  return(rgb2hex(r = rc, g = gc, b = bc))
}

# Plot dose-response functions from the biphasic model
plot.DR <- function(jagsdata = NULL, rjdata = NULL, simulation = TRUE, quantiles = seq(95, 1, by = -5)) {
  
  if(is.null(rjdata) & is.null(jagsdata)) stop("Cannot find input data")
  if(!is.null(rjdata) & !is.null(jagsdata)) stop("Multiple input datasets")
  
  # Set up plot
  par(mfrow = c(1, 1))
  
  summary.jags <- summary(jagsdata)
  n.thin <- summary.jags$thin
  n.chains <- summary.jags$nchain
  n.iter <- summary.jags$end - summary.jags$start + 1

  
  # Dose range
  dose <- seq(from = L, to = U, by = 1)
  
  # Set up matrix for probabilities of response
  p.response <- matrix(data = 0, nrow = (n.iter / n.thin) * n.chains, 
                       ncol = length(dose))
  
  ioffset <- (0:(n.chains - 1)) * (n.iter / n.thin)
  
  npts <- 20 # only 20 quadrature points for now in integrating out the rf
  p.response.individ <- matrix(data = 0, nrow = npts, ncol = length(dose))
  
  pb <- dplyr::progress_estimated(n = (n.iter / n.thin) * n.chains)
  
  for (i in 1:(n.iter / n.thin)) {
    
    for (k in 1:n.chains) {
      
      pb$tick()$print()
      
      # Integrate out the random effect
      pi.individ <- pnorm(q = qnorm(p = seq(0, 1, length = (npts + 2))[-c(1, (npts + 2))],
                                    mean = jagsdata[[k]][i, "psi"], 
                                    sd = jagsdata[[k]][i, "omega"]))
      
      # Dose-response curve corresponding to lower mixture component
      p.response.individ.lower <- truncnorm::ptruncnorm(q = dose, 
                                                        a = L, 
                                                        b = jagsdata[[k]][i, "alpha"],
                                                        mean = jagsdata[[k]][i, "nu[1]"], 
                                                        sd = jagsdata[[k]][i, "tau[1]"])
      
      # Dose-response curve corresponding to higher mixture component
      p.response.individ.upper <- truncnorm::ptruncnorm(q = dose, 
                                                        a = jagsdata[[k]][i, "alpha"], 
                                                        b = U, 
                                                        mean = jagsdata[[k]][i, "nu[2]"], 
                                                        sd = jagsdata[[k]][i, "tau[2]"])
      
      for (j in 1:npts) {
        p.response.individ[j, ] <- pi.individ[j] * p.response.individ.lower +
          (1 - pi.individ[j]) * p.response.individ.upper}
      
      p.response[i + ioffset[k], ] <- apply(X = p.response.individ, MARGIN = 2, FUN = mean)
      
    } # End for k loop
  } # End for i loop
  
  # Calculate median
  p.response.med <- apply(X = p.response, MARGIN = 2, FUN = median)
  
  # Calculate mean 
  p.response.mean <- apply(X = p.response, MARGIN = 2, FUN = mean)
  
  # Produce plot
  plot(dose, p.response.med, type = "n", ylab = "p(response)", xlab = "dose (dB re 1uPa)",
       main = "", ylim = c(0, 1), cex.lab = 1, cex.main = 1, cex.axis = 1)
  
  # Add quantiles
  purrr::walk(.x = quantiles, .f = ~{
    yy <- apply(p.response, 2, quantile, c((50 - .x / 2) / 100, (50 + .x / 2) / 100))
    y1 <- c(yy[1, ], yy[2, seq(from = length(dose), to = 1, by = -1)])
    x1 <- c(dose, dose[seq(from = length(dose), to = 1, by = -1)])
    polygon(x1, y1, col =  hexa2hex(input.colour = "#066e9e", opacity = 1-(.x/100)), border = NA)
  })
  
  line.colors <- c("orange", "deepskyblue4")
  
  # Add posterior median
  lines(dose, p.response.med, lwd = 2, col = line.colors[1])
  
  # Add the true function if it's a simulation
  if (simulation) {
    
    # p.response.norf <- p.biphasic(dose, pnorm(simdata$psi), simdata$nu, simdata$tau, data$L, data$U, simdata$beta)
    
    # integrate out random effect
    pi.individ <- pnorm(qnorm(seq(0, 1, length = (npts + 2))[-c(1, (npts + 2))], psi, omega))
    
    for (i in 1:npts){
      p.response.individ[i, ] <- p.biphasic(dose = dose, pi.individ[i], nu, tau, L, U, alpha)}
    
    p.response.withrf <- apply(X = p.response.individ, MARGIN = 2, FUN = mean)
    
    lines(dose, p.response.withrf, lwd = 2, col = line.colors[2])
    # lines(dose, p.response.norf, col = "lightblue1", lwd = 2)
    legend(x = "topleft", legend = c("posterior", "simulated"), col = line.colors, lty = 1, lwd = 1)
  }
  
}


