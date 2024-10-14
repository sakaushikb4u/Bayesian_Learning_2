# ===================================================================
# Init Library
# ===================================================================
knitr::opts_chunk$set(echo = TRUE)
library(ggplot2)
library(gridExtra)
library(readxl)
library(LaplacesDemon)
library(MASS)
library(reshape2)
library(matrixStats)
library(mvtnorm)


# ===================================================================
# 1. Linear and polynomial regression
# ===================================================================

#--------------------------------------------------------------------
# Init code for question 1
#--------------------------------------------------------------------
rm(list = ls())

# ===================================================================
# 1.a. Use the conjugate prior for the linear regression model
# ===================================================================

# given values
mu_0 <- c(0, 100, -100)
omega_0 <- 0.01 * diag(3)
v_0 <- 1
sigma2_0 <- 0.1
days_in_year <- 365
first_day_of_year <- as.Date("2022-01-01")
#--------------------------------------------------------------------
# load data and plot
#--------------------------------------------------------------------
# load the data
temperature_data <- read_xlsx(path="Linkoping2022.xlsx")
temperature_data$time <- as.numeric(as.Date(temperature_data$datetime) - first_day_of_year) /
                         days_in_year
temperature_data$time_square <- temperature_data$time ** 2 

# plot points to visualize 
plot_1a <- ggplot(temperature_data, aes(x = time, y=temp)) +
    geom_point(aes(color = "red")) +
    ylab("Temperature") + xlab("Time") + 
    ggtitle("Time vs Temperature") +
    scale_color_identity(name = "Colour", 
                       breaks = c("red"),
                       labels = c("Data"),
                       guide = "legend") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position="bottom")

plot_1a

#--------------------------------------------------------------------
# make a linear regression model
#--------------------------------------------------------------------
linRegModel <- lm(temp ~ time + time_square, data = temperature_data)
summary(linRegModel)
#--------------------------------------------------------------------
# function to simulate draws from the joint prior distribution
#--------------------------------------------------------------------
joint_conjugate_prior <- function(nr_of_iterations, v_0, sigma2_0, mu_0, omega_0, days_in_year=365){  
  # data frame to save prior coef data
  prior <- data.frame(matrix(ncol = 3, nrow = nr_of_iterations)) # for every beta one column

  # data frame to save linear regression line
  regline <- data.frame(matrix(nrow = days_in_year, ncol = nr_of_iterations)) # every col for one day
  
  # joint conjugate prior
  for (i in 1:nr_of_iterations) {
    
    #  calculate var
    var <- LaplacesDemon::rinvchisq(1,v_0,sigma2_0)
    
    # calculate beta
    beta <- MASS::mvrnorm(1, mu_0, var*solve(omega_0))
    
    # set values
    prior[i,1:3] <- beta
    regline[,i] <- beta[1] + beta[2] * temperature_data$time + beta[3] * temperature_data$time_square
  }
  
  # add meaningful column and row names
  for (i in 1:nr_of_iterations) {
    colnames(regline)[i] <- paste("pred","", i)
    rownames(regline)[i] <- paste("day","",i)
  }
  
  # return the prior and regline
  return(list(p = prior, r = regline))
}
#--------------------------------------------------------------------
# draw from joint_conjugate_prior and plot
#--------------------------------------------------------------------
# set random seed
set.seed(123456)

# joint conjugate prior
N <- 50

result <- joint_conjugate_prior(N, v_0, sigma2_0, mu_0, omega_0)
prior <- result$p
regline <- result$r

colors <- sample(colors(), N)

plot_1a_lines <- plot_1a

plot_1a_lines +
  mapply(function(i, col) {
    geom_line(aes(x = time, y = regline[, i]), col = col)
  }, 1:N, colors)
#--------------------------------------------------------------------
# draw from joint_conjugate_prior and plot(new parameters)
#--------------------------------------------------------------------
# new given parameters
mu_0 <- c(-10, 100, -100)
omega_0 <- 0.03 * diag(3)
v_0 <- 2
sigma2_0 <- 0.1
# set random seed
set.seed(123456)

# joint conjugate prior
N <- 50

result <- joint_conjugate_prior(N, v_0, sigma2_0, mu_0, omega_0)
prior <- result$p
regline <- result$r

colors <- sample(colors(), N)

plot_1a_lines <- plot_1a

plot_1a_lines +
  mapply(function(i, col) {
    geom_line(aes(x = time, y = regline[, i]), col = col)
  }, 1:N, colors)

# ===================================================================
# 1.b. Write a function that simulates draws from the joint posterior distribution of $\beta_{0},\beta_{1},\beta_{2} \text{ and } \sigma^2$
# ===================================================================
#--------------------------------------------------------------------
# function to simulate draws from the joint posterior distribution
#--------------------------------------------------------------------
joint_conjugate_posterior <- function(nr_of_iterations, v_0, sigma2_0, mu_0, omega_0, days_in_year=365){  
  # data frame to save prior coef data
  prior <- data.frame(matrix(ncol = 3, nrow = nr_of_iterations)) # for every beta one column

  # data frame to save linear regression line
  regline <- data.frame(matrix(nrow = days_in_year, ncol = nr_of_iterations)) # every col for one day
  
  # joint conjugate prior
  for (i in 1:nr_of_iterations) {
    
    #  calculate var
    var <- LaplacesDemon::rinvchisq(1,v_0,sigma2_0)
    
    # calculate beta
    beta <- MASS::mvrnorm(1, mu_0, var*solve(omega_0))
    
    # set values
    prior[i,1:3] <- beta
    regline[,i] <- beta[1] + beta[2] * temperature_data$time + beta[3] * temperature_data$time_square
  }
  
  # add meaningful column and row names
  for (i in 1:nr_of_iterations) {
    colnames(regline)[i] <- paste("pred","", i)
    rownames(regline)[i] <- paste("day","",i)
  }
  
  # return the prior and regline
  return(list(p = prior, r = regline))
}

# ===================================================================
# 1.b.i Plot a histogram for each marginal posterior of the parameters
# ===================================================================
#--------------------------------------------------------------------
# calculate the posterior distribution parameters
#--------------------------------------------------------------------
k <- 3
X <- cbind(1, temperature_data$time, temperature_data$time_square) 
y <- temperature_data$temp
n <- nrow(temperature_data)
beta_hat <- solve(t(X) %*% X) %*% t(X) %*% y
mu_n <- solve(t(X) %*% X + omega_0) %*% (t(X) %*% X %*% beta_hat + omega_0 %*% mu_0)
omega_n <- t(X) %*% X + omega_0
v_n <- v_0 + n
sigma2_n <- (v_0 * sigma2_0 + t(y) %*% y + t(mu_0) %*% omega_0 %*% 
            mu_0 - t(mu_n) %*% omega_n %*% mu_n) / v_n

#--------------------------------------------------------------------
# simulate the beta parameter and draw
#--------------------------------------------------------------------
# simulate the beta parameter
data_hist <- as.data.frame(
  mvtnorm::rmvt(n = 1000, delta = mu_n, df = n-k,
                sigma = as.numeric(sigma2_n) * solve(t(X) %*% X)))

# simulate the sigma^2 parameter 
var_hist <- LaplacesDemon::rinvchisq(n = 1000, v_n, sigma2_n)
data_hist <- cbind(data_hist, var_hist)
cnames <- c("beta0", "beta1", "beta2", "sigma2")
colnames(data_hist) <- cnames

# plot the histogram of 4 parameters
plot_func <- function(cname){
  ggplot(data_hist, aes(x = !!sym(cname))) +
    geom_histogram(aes(y = after_stat(density)),
                   colour = "black",
                   fill   = "white",
                   bins   = 30) +
    geom_density(alpha = .2, fill = "blue") 
}
plot(arrangeGrob(grobs = lapply(cnames, plot_func)))

# ===================================================================
# 1.b.ii Make a plot, and overlay curves for the 90% equal tail posterior probability intervals of f(time), then comment
# ===================================================================

#--------------------------------------------------------------------
# calculate the median of beta and f(time)
#--------------------------------------------------------------------
data_hist = data_hist[,1:3] 

beta_median = matrixStats::colMedians(as.matrix(data_hist))
pred.1b <- beta_median %*% t(X)

#--------------------------------------------------------------------
# calculate the 90% credible interval
#--------------------------------------------------------------------
preds <- as.matrix(data_hist) %*% t(X) # regression function
pred_interval <- data.frame(nrow = n, nrow = 2)
colnames(pred_interval) <- c("lower","upper")

for(i in 1:days_in_year){
  data_t <- preds[,i]
  pred_interval[i,] <- quantile(data_t, probs = c(0.05,0.95))
}
# plot data 
data.1b <- cbind(temperature_data, t(pred.1b),pred_interval)

ggplot(data.1b) +
  geom_point(aes(x = time, y = temp,color = "red")) +
  geom_line(aes(x = time, y = t(pred.1b),color = "blue"),linewidth=1) +
  geom_line(aes(x = time, y = lower,color = "green")) +
  geom_line(aes(x = time, y = upper,color = "green")) + 
  ggtitle("Time vs Temperature with pred curve and pred interval") +    
  scale_color_identity(name = "Colour", 
                       breaks = c("red", "blue", "green"),
                       labels = c("Data", "Pred Curve", "Posterior Probability Interval"),
                       guide = "legend") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position="bottom")

# ===================================================================
# 1.c. Use the simulated draws in (b) to simulate the posterior distribution of  $\widetilde{x}$
# ===================================================================

#--------------------------------------------------------------------
# calculate the highest prediction for every day
#--------------------------------------------------------------------

pred_highest <- c()
for (i in 1:days_in_year) {
  pred_highest[i] <- max(preds[,i])
}

data.1c <- cbind(temperature_data, t(pred.1b), pred_interval, pred_highest)
ggplot(data.1b) +
  geom_point(aes(x = time, y = temp,color = "red")) +
  geom_line(aes(x = time, y = t(pred.1b),color = "blue"),linewidth=1) +
  geom_line(aes(x = time, y = lower,color = "green")) +
  geom_line(aes(x = time, y = upper,color = "green")) + 
  geom_line(aes(x = time, y = pred_highest,color = "skyblue")) + 
  ggtitle("Time vs Temperature with pred curve, pred_interval and pred_highest") +
  scale_color_identity(name = "Colour", 
                       breaks = c("red", "blue", "green", "skyblue"),
                       labels = c("Data", "Pred Curve", "Posterior Probability Interval", "Pred Highest"),
                       guide = "legend") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position="bottom")


# ===================================================================
# 2. Posterior approximation for classification with logistic regression
# ===================================================================

# clean up environment
rm(list=ls())

# ===================================================================
# 2.a. Approximate the posterior distribution of the parameter vector $\beta$ with a multivariate normal distribution and comment. compute 95% equal tail posterior probability interval and comment.
# ===================================================================


#--------------------------------------------------------------------
# given values 
#--------------------------------------------------------------------
mu <- 0
tau <- 2
#--------------------------------------------------------------------
# load the data
#--------------------------------------------------------------------
WomenAtWork = read.table("WomenAtWork.dat", header = TRUE)
y = WomenAtWork$Work
X = WomenAtWork[,2:8]
X = as.matrix(X)
# Xnames will be used as the names of the output of beta and sigma
Xnames <- colnames(X)
#--------------------------------------------------------------------
# Functions that returns the log posterior for the logistic
#--------------------------------------------------------------------
LogPostLogistic <- function(betas,y,X,mu,sigma){
  linPred <- X%*%betas;
  logLik <- sum( linPred*y - log(1 + exp(linPred)) );  
  # Likelihood is not finite, stear the optimizer away from here!
  # Idea from cource code by teacher
  if (abs(logLik) == Inf)
    logLik = -20000;  
  # prior follows multi-normal distribution
  logPrior <- dmvnorm(betas, mu, sigma, log=TRUE);
  return(logLik + logPrior)
}
#--------------------------------------------------------------------
# get the optimized beta and inverse jacobian
#--------------------------------------------------------------------
# number of features
n <- dim(X)[2]
# setting up the prior
mu    <- as.vector(rep(mu,n))  # Prior mean vector
sigma <- tau^2 * diag(n)  # Prior variance matrix

# use random initial values
init_val <- as.vector(rnorm(dim(X)[2]))
# optimize the log posterior
OptimRes <- optim(init_val,
                  LogPostLogistic,
                  y = y,
                  X = X,
                  mu = mu,
                  sigma = sigma,
                  method=c("BFGS"),
                  control=list(fnscale=-1),
                  hessian=TRUE)
# set values to print out
posterior_mode  <- OptimRes$par

# hessian is the negative of the observed Hessian evaluated at the posterior mode
# Jacobian = (-hessian)
beta_jacobian <- -OptimRes$hessian
beta_inverse_jacobian <- solve(beta_jacobian) 
# Naming the coefficient by Xnames
names(posterior_mode) <- Xnames 
# Naming the coefficient by Xnames
colnames(beta_inverse_jacobian) <- colnames(X) 

print('The posterior beta is:')
print(posterior_mode)
print('The Inverse Jacobian of beta is:')
print(beta_inverse_jacobian)
#--------------------------------------------------------------------
# Compute an approximate 95% equal tail posterior probability interval
#--------------------------------------------------------------------
beta_sim <- mvtnorm::rmvnorm(n = 1000, mean = posterior_mode, sigma = beta_inverse_jacobian)
data.NSmallChild <- beta_sim[,6]

# calculate the 95% credible interval
pred_interval <- quantile(data.NSmallChild, probs = c(0.025,0.975))
pred_interval
#--------------------------------------------------------------------
# logistic regression
#--------------------------------------------------------------------
glmModel <- glm(Work ~ 0 + ., data = WomenAtWork, family = binomial)
summary(glmModel)

# ===================================================================
# 2.b. Write a function that simulates draws from the posterior predictive distribution and plot
# ===================================================================

#--------------------------------------------------------------------
# Simulates draws from the posterior predictive distribution
#--------------------------------------------------------------------
post_pred <- function(X, beta) {  
  pred <- beta %*% X #linear predictions
  pred_logit <- 1- exp(pred)/(1 + exp(pred)) #sigmoid function to model probabilities
  return(data.frame(Probability = pred_logit))
}
#Create matrix of sample features
X_data <- as.matrix(c(0, 18, 11, 7, 40, 1, 1))
#function call on sample data
post_pred_sim <- post_pred(X_data, beta_sim)
#plotting the posterior prediction density & histogram
ggplot(data = post_pred_sim, aes(x = Probability)) + xlim(c(0,1)) +
  geom_histogram(aes(y = after_stat(density)), 
                 bins = 100, 
                 color = "black", 
                 fill = "grey") +
  geom_density(alpha = 0.2, fill = "blue") + 
  geom_vline(xintercept = 0.5, alpha = 0.8,  color = "red")

# ===================================================================
# 2.c. Rewrite your function and plot
# ===================================================================

#--------------------------------------------------------------------
# Get the trial number and plot
#--------------------------------------------------------------------
n <- 13
post_pred_binom <- function(X, beta) {
  
  pred <- beta %*% X #linear predictions
  pred_logit <- 1- exp(pred)/(1 + exp(pred)) #sigmoid function to model probabilities  
  trial <- c()
  for (i in 1:length(pred_logit)) {
    trial[i] <- rbinom(n = 1, size = n, prob = pred_logit[i])
  }  
  return(trial)
}

test_trials <- post_pred_binom(X_data, beta_sim)
hist(test_trials)

