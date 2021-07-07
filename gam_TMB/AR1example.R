#AR1 example
library(TMB)
library(tidyverse)
library(latex2exp)
library(ggpubr)
setwd("~/Documents/yield-analysis-2021/gam_TMB")

## Simulate data
set.seed(123)
n <- 200
sigma <- .1
phi <- .95
simdata <- function(){
  u <- numeric(n)
  u[1] = rnorm(1)
  if(n>=2)
    for(i in 2:n){
      u[i] = phi * u[i-1] + rnorm(1) #AR1 process
    }
  u <- u * sigma
  x <- as.numeric( rbinom(n, 1, plogis(u)) )
  data <- list(obs=x,latent=u)
  return(data)
}
##
data <- simdata()

par(mfrow=c(2,1))
plot(data$obs,pch=19,ylab='Observed')
plot(data$latent,type='l',pch=19,ylab='Latent process'); abline(h=0,lty='dashed',col='red')

compile("AR1example.cpp") # does not work
dyn.load(dynlib("AR1example"))

# Fixed effects
parameters <- list(phi=phi, logSigma=log(sigma))

# Random effects
parameters$u <- rep(0,n)

# Create TMB object
obj <- MakeADFun(data, parameters, random="u", DLL="AR1example")
obj$fn()

opt <- nlminb(obj$par, obj$fn, obj$gr, control=list(trace=0))
opt$convergence == 0 #Did we converge?

rep <- sdreport(obj)
summary(rep, select=c("fixed", "report"))[c("phi","logSigma","sigma"),] # sigma is found by delta method on logSigma

#Estimates of latent variables
options(repr.plot.width = 15, repr.plot.height = 7)

# Information from report, fixed effects
srep_fixed <- summary(rep, select = c("fixed"), p.value = TRUE) %>% 
  tibble::as_tibble(rownames = NA) %>%
  tibble::rownames_to_column() %>% 
  dplyr::rename(parameter = rowname, 
                estimate = Estimate, 
                std_error = `Std. Error`,
                z_value = `z value`,
                p_value = `Pr(>|z^2|)`) %>% 
  mutate(type = "fixed")

# Information from report, random effects
srep_random <- summary(rep, select = c("random"), p.value = TRUE) %>% 
  tibble::as_tibble(rownames = NA) %>%
  tibble::rownames_to_column() %>% 
  dplyr::rename(parameter = rowname, 
                estimate = Estimate, 
                std_error = `Std. Error`,
                z_value = `z value`,
                p_value = `Pr(>|z^2|)`) %>% 
  mutate(type = "random") %>% 
  dplyr::mutate(parameter = ifelse(parameter == "h", paste0("h", 1:n()), parameter))

# Combine information
srep <- dplyr::bind_rows(srep_fixed, srep_random)

# Add data information and confidence estimates
report <- srep %>% 
  filter(parameter == "u") %>% 
  mutate(
    obs = data$obs,
    time = 1:n(),
    u_upper = estimate + 2 * std_error,
    u_lower = estimate - 2 * std_error,
    actual = data$latent)

(p1 <- report %>%
  ggplot() + 
  geom_ribbon(aes(time, ymax = u_upper, ymin = u_lower), alpha= 0.10) + 
  geom_line(aes(time,estimate), size=1) +
    geom_line(aes(time,actual),col='red')+
  ggtitle(TeX("$\\log U_t$ process: Estimates and error bounds")) + 
  xlab("i") + 
  ylab("value") + 
  theme_bw())

## PROBABILITIES
# Information from report, random effects
srep_report <- summary(rep, select = c("report"), p.value = TRUE) %>% 
  tibble::as_tibble(rownames = NA) %>%
  tibble::rownames_to_column() %>% 
  dplyr::rename(parameter = rowname, 
                estimate = Estimate, 
                std_error = `Std. Error`,
                z_value = `z value`,
                p_value = `Pr(>|z^2|)`) %>% 
  mutate(type = "random") %>% 
  dplyr::mutate(parameter = ifelse(parameter == "h", paste0("h", 1:n()), parameter))

# Add data information and confidence estimates
report <- srep_report %>% 
  filter(parameter == "p") %>% 
  mutate(
    obs = data$obs,
    time = 1:n(),
    p_upper = estimate + 2 * std_error,
    p_lower = estimate - 2 * std_error)

p2 <- report %>%
  ggplot() + 
  geom_point(aes(time, obs)) + 
  geom_line(aes(time,estimate), size=1) + 
  geom_ribbon(aes(time, ymax = p_upper, ymin = p_lower), alpha= 0.10) + 
  ggtitle(TeX("Underlying probability process, $\\P_t$, and observations 0 and 1")) + 
  xlab("i") + 
  ylab("value") + 
  theme_bw()

ggarrange(p1, p2, ncol=2)



