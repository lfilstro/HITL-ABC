###### This script generates the misspecification results

# Load functions from Turin_generateData.R
library(abc)
library(rmatio)
library(pracma)
library(ggplot2)
require(gridExtra)

# Loading simulated statistics and parameters from pre-saved data
# To generate new data, run Turin_generateData.R first
data_param <- readRDS("Turin_param_miss")[,1:3]
sumStats <- readRDS("Turin_ssim_miss")[,1:6]

# Plotting functions
ggplot_expert <- function(x,data_param, index, prior_bounds, parameter_name, param_true){
  ggplot() + 
    geom_density(aes(x[,index], fill="Posterior"), alpha = 0.5) +
    geom_density(aes(data_param[,index], fill="Prior"), alpha = 0.5) +
    geom_vline(aes(xintercept = param_true[index], colour = "a"), linetype = "dashed")+
    xlim(prior_bounds[index, 1], prior_bounds[index, 2])+
    xlab(parameter_name[index]) +
    scale_fill_manual(values = c("blue", "orange")) +
    scale_colour_manual(values = c( "a" = "green4"), labels = c("True value")) +
    theme(legend.title=element_blank()) +
    theme(legend.position = "none") +
    theme(axis.title.y = element_blank())
}

c1 <- ggplot_expert(samples_hitl$adj.values, data_param, 1, prior_bounds, param_name, param_true)+ ggtitle(expression(zeta~"= 5")) 
c2 <- ggplot_expert(samples_hitl$adj.values, data_param, 2, prior_bounds, param_name, param_true)
c3 <- ggplot_expert(samples_hitl$adj.values, data_param, 3, prior_bounds, param_name, param_true)

d1 <- ggplot_expert(samples_hitl$adj.values, data_param, 1, prior_bounds, param_name, param_true)+ ggtitle(expression(zeta~"= 10"))+
  theme(legend.position = "right") + theme(axis.text.x = element_blank())
d2 <- ggplot_expert(samples_hitl$adj.values, data_param, 2, prior_bounds, param_name, param_true)
d3 <- ggplot_expert(samples_hitl$adj.values, data_param, 3, prior_bounds, param_name, param_true)
grid.arrange(c1, d1, c2, d2, c3,d3, ncol = 2)

ggplot_measured <- function(x,index, prior_bounds, parameter_name, true_parameter, fill_col){
  
  library(tidyverse)
  theme_set(theme_classic())
  
  # distribution <- prior_bounds[[1]][1]
  x_start <- prior_bounds[index, 1]
  x_end <- prior_bounds[index, 2]
  
  
  p1 <- ggplot() +
    geom_density(aes(x[,index], fill="Posterior\nestimate"), alpha = 0.75) +
    # stat_density(aes(y[,index], fill="Posterior\nestimate"), alpha = 0.2) +
    geom_vline(aes(xintercept = true_parameter[index], colour = "a"), linetype = "dashed")+
    # xlim(x_start, x_end)+
    xlab(parameter_name[index]) +
    scale_linetype_manual(values = c("dashed" = "dashed"), labels = "Prior") +
    scale_colour_manual(values = c( "a" = "green4"), labels = c("True value")) +
    scale_x_continuous(limits=c(x_start, x_end), breaks = scales::pretty_breaks(n = 2)) +
    theme(legend.title=element_blank()) +
    scale_fill_manual(values = fill_col) +
    theme(axis.title.y = element_blank(),
          axis.text = element_text(size=8),
          # axis.text.x = element_text(angle = 30),
          # axis.ticks.y = element_blank(),
          axis.title=element_text(size=14, family = "mono")) +
    theme(legend.position = "none") #+theme(text = element_text(size=15)) 
}

######################----------------MAIN-------------######################

# Settings
B <- 4e9
Ns = 801
nReal <- 300 # Number of data realizations
tol = 0.05 # tolerance threshold in ABC

# Prior ranges
G0_start <- 1e-10
G0_end <- 3e-9
T_start <- 1e-9
T_end <- 2e-8
lambda_start <- 1e8
lambda_end <- 4e9

prior_bounds <- rbind(c(G0_start, G0_end), c(T_start, T_end),c(lambda_start, lambda_end))
param_name = c(expression(G[0]), "T", expression(lambda))

#################--------- Generating observed data------------#############
# True parameter values
G0_true = 1e-9
T_true = 1e-8
lambda_true <- 1e9
sigma2_W_true <- 0
param_true <- c(G0_true, T_true, lambda_true)

# Generate observed data
Y_obs <- TurinModel(G0_true, T_true, lambda_true, sigma2_W_true, B = B, Ns = Ns, N = nReal, tau0 = 0)

# Compute observed temporal moments
m_obs <- log(temporalMomentsGeneral(Y_obs, 3))
# m_obs <- readRDS("Turin_mobs") # Uncomment this to load pre-saved data used in the paper

# Computing summary statistics
sumStats_obs <- summaryMoments(m_obs)[1:6]


# Corrupting data
zeta <- 1 # Value of misspecification level zeta. Set to either 1, 5, or 10 in the paper
m_obs_miss <- array(m_obs, dim = dim(m_obs))
set.seed(1)
m_obs_miss[, 1] <- m_obs_miss[, 1] + rnorm(nReal, mean = 0, sd = sqrt(zeta)) # Adding noise to the first temporal moment
sumStats_obs <- summaryMoments(m_obs_miss)[1:6] # Misspecified statistics

# sumStats_obs <- readRDS("Turin_sobs_zeta1")

########----------- Running ABC ------------###############
# Run the HITL-ABC method from TurinExample.ipynb file using the generated datasets to get the misspecified statistic.
# In this code, the misspecified statistic is set to the variance of m_0 (4th statistic).
samples_hitl <- abc::abc(sumStats_obs[c(-4)], data_param, sumStats[,c(-4)], tol=tol, method = "loclinear",  transf = 'logit', logit.bounds = prior_bounds)

# linear-ABC
samples_reg <- abc::abc(sumStats_obs, data_param, sumStats, tol=tol, method = "loclinear",  transf = 'logit', logit.bounds = prior_bounds)


# Plotting
p1 <- ggplot_measured(samples_reg$adj.values, 1, prior_bounds, param_name, param_true, c("grey")) +  ggtitle("linear-ABC")
p2 <- ggplot_measured(samples_reg$adj.values, 2, prior_bounds, param_name, param_true, c("grey"))
p3 <- ggplot_measured(samples_reg$adj.values, 3, prior_bounds, param_name, param_true, c("grey"))

# p4 <- ggplot_measured(samples_reg$adj.values, 3, prior_bounds, param_name, param_true)

h1 <- ggplot_measured(samples_hitl$adj.values, 1, prior_bounds, param_name, param_true, c("red"))+ ggtitle("HITL-ABC")
h2 <- ggplot_measured(samples_hitl$adj.values, 2, prior_bounds, param_name, param_true, c("red"))
h3 <- ggplot_measured(samples_hitl$adj.values, 3, prior_bounds, param_name, param_true, c("red"))

grid.arrange(p1, h1, p2, h2, p3,h3, ncol = 2)




##########---Additional results: Plotting neural-ABC and ridge-ABC results -----#################

samples_neural <- abc::abc(sumStats_obs, data_param, sumStats, tol=tol, method = "neuralnet",
                        transf = 'logit', logit.bounds = prior_bounds)

# Plotting
n1 <- ggplot_measured(samples_neural$adj.values, 1, prior_bounds, param_name, param_true, c("grey")) +  ggtitle("neural-ABC")
n2 <- ggplot_measured(samples_neural$adj.values, 2, prior_bounds, param_name, param_true, c("grey"))
n3 <- ggplot_measured(samples_neural$adj.values, 3, prior_bounds, param_name, param_true, c("grey"))

samples_ridge <- abc::abc(sumStats_obs, data_param, sumStats, tol=tol, method = "ridge",
                           transf = 'logit', logit.bounds = prior_bounds)

# Plotting
r1 <- ggplot_measured(samples_ridge$adj.values, 1, prior_bounds, param_name, param_true, c("grey")) +  ggtitle("ridge-ABC")
r2 <- ggplot_measured(samples_ridge$adj.values, 2, prior_bounds, param_name, param_true, c("grey"))
r3 <- ggplot_measured(samples_ridge$adj.values, 3, prior_bounds, param_name, param_true, c("grey"))

grid.arrange(n1, r1, n2, r2, n3,r3, ncol = 2)
