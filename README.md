# Biphasic Bayesian dose-response model

Comparison of JAGS and R implementations [currently not working]

Current observations (`biphasic_mcmc.R`):

-   Parameter estimates appear to match when there is no right-censoring.
-   The previous issue encountered with individual \#22 may have been the result of an error in the calculation of the acceptance ratio calculation for parameter `k_ij` .
-   Mismatch still occurs when right-censored observations are included.
-   The current example `seed = 122` seems to run ok but other simulated datasets do not.
-   Issues are most likely related to `alpha`, `mu_ij`, and `k_ij` , as far as I can tell.
