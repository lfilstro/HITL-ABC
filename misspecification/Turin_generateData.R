library(pracma)

##########-----------FUNCTIONS----------#######
# Function to generate samples from the Turin Model
TurinModel <- function(G0, T, lambda_0, sigma2_N,B = 4e9, Ns = 801, N = 100, tau0 = 0){
  
  nRx <- N
  
  delta_f <- B/(Ns-1)   # Frequency step size
  t_max <- 1 / delta_f
  
  tau <- seq(from = 0, to = t_max, length.out = Ns)
  
  H <- array(0, dim = c( nRx, Ns))
  
  mu_poisson <- lambda_0*t_max # Mean of Poisson process
  
  for (jR in 1:nRx) {
    
    repeat{
      n_points <- rpois(1, mu_poisson) # Number of delay points sampled from Poisson process
      if (n_points > 0) {
        break
      }
    }
      
    delays <- runif(n_points,0,t_max) # Delays sampled from a 1-dimensional Poisson point process
      
    delays <- delays[order(delays)]
      
    alpha <- array(0, dim = c(n_points,1))     # Initialising vector of gains of length equal to the number of delay points
      
    sigma2 <- G0*exp(-delays/T)/lambda_0 * B
      
    for (l in 1:n_points) {

      alpha[l] <- (rnorm(1,mean = 0,sd = sqrt(sigma2[l]/2)) + rnorm(1,mean = 0,sd = sqrt(sigma2[l]/2))*1i)

    }
      
    H[jR, ] <- exp(-1i*2*pi*delta_f * (0:(Ns-1))%*%t(delays)) %*% alpha
  }

  # Noise power by setting SNR
  Noise <- array(data = 0, dim = c(nRx, Ns))
  
  for (j in 1:nRx) {
    Noise[j,] <- rnorm(Ns, mean = 0, sd= sqrt(sigma2_N / 2)) + rnorm(Ns, mean = 0,  sd= sqrt(sigma2_N / 2)) * 1i
  }
  
  # Received signal in frequency domain
  
  Y <- H + Noise
  
  # Y <- pracma::Reshape(Y, N, Ns)
  
  return(Y)
  
}

# Function to compute the means and covariances of the temporal moments
summaryMoments <- function(x){
  library(matrixcalc)
  z <-  c(apply(x, 2, mean), upper.triangle(cov(x)))
  z <- z[z!=0]
  y <- array(z, dim = c(1, length(z)))
  # colnames(y) <- c("Mean(m0)", "Mean(m1)", "Mean(m2)", "Var(m0)", "Cov(m0,m1)", "Var(m1)", "Cov(m0,m2)", "Cov(m1,m2)", "Var(m2)")
  return(y)
}

# Function to compute the temporal moments 
temporalMomentsGeneral <- function(Y, K=3, B = 4e9){
  library(pracma)
  
  N <- dim(Y)[1]
  Ns <- dim(Y)[2]
  
  delta_f <- B / (Ns-1)
  t_max <- 1 / delta_f
  
  tau <- seq(from = 0, to = t_max, length.out = Ns)
  
  m <- array(0, dim = c(N, K))
  
  for (k in 1:K) {
    for (i in 1:N) {
      y <- ifft(Y[i,])
      m[i,k] <- trapz(tau, tau^(k-1) * abs(y)^2)
    }
  }
  
  return(m)
}

# Function to compute the means and covariances of the temporal moments given the data from the Turin model
temporalMomentsStats <- function(x, order = 2, B = 4e9){
  
  library(pracma)
  
  if (is.matrix(x) == TRUE) {
    nRealizations <- dim(x)[1]
    Ns <- dim(x)[2]
    
    delta_f <- B / Ns
    t_max <- 1 / delta_f
    
    tau <- seq(from = 0, to = t_max, length.out = Ns)
    
    m0 <- array(0, dim = c(nRealizations, 1))
    m1 <- array(0, dim = c(nRealizations, 1))
    m2 <- array(0, dim = c(nRealizations, 1))
    
    for (i in 1:nRealizations) {
      y <- ifft(x[i,])
      m0[i] <- trapz(tau,  abs(y)^order)
      m1[i] <- trapz(tau, tau * abs(y)^order)
      m2[i] <- trapz(tau, tau^2 * abs(y)^order)
    }
    
    if (length(x[,1]) > 1) {
      a <- summaryMoments(cbind(m0, m1, m2))
    } else if (length(x[,1]) == 1) {
      a <- cbind(m0, m1, m2)
      colnames(a) <- c("m0", "m1", "m2")
    }
    
  } else if (is.matrix(x) == FALSE) {
    N <- 1
    Ns <- length(x)
    
    delta_f <- B / Ns
    t_max <- 1 / delta_f
    
    tau <- seq(from = 0, to = t_max, length.out = Ns)
    
    y <- ifft(x)
    
    m0 <- trapz(tau,  abs(y)^2)
    m1 <- trapz(tau, tau * abs(y)^2)
    m2 <- trapz(tau, tau^2 * abs(y)^2)
    
    a <- cbind(m0, m1, m2)
    colnames(a) <- c("m0", "m1", "m2")
  }
  
  return(a)
  
}

# Function to find the peaks in the data
peaks <- function(X){
  
  a <- array(dim(findpeaks(pdp(X)))[1], dim = c(1,1))
  colnames(a) <- c("Number of Peaks")
  return(a)
  
}

# Function to compute all the summary statistics of the Turin model. 
# Note that we only use the first 6 statistics for the misspecification experiment.
computeTurinStatistics <- function(Y){
  
  sumStats <- array(0, dim = c(20,1))
  
  sumStats[1:9] <- temporalMomentsStats(Y, order = 2)
  sumStats[10:18] <- temporalMomentsStats(Y, order = 4)
  sumStats[19] <- max(findpeaks(pdp(Y))[,1]) - min(findpeaks(pdp(Y))[,1])
  sumStats[20] <- mean(pdp(Y))
  
  return(sumStats)
}

# Function to transform x from linear scale to decibel
db <- function(x){
  return(10*log10(x))
}


####-------Parameters--------------#########
B <- 4e9  # Bandwidth
Ns <- 801 # Number of samples
delta_f <- B/(Ns-1)   # Frequency step size
t_max <- 1 / delta_f

tau <- seq(from = 0, to = t_max, length.out = Ns) # Time axis

nReal <- 300 # Number of realizations


#############---------Generating simulated statistics----------############

start.time <- Sys.time()

M <- 2000 # Number of samples from the prior
nStats <- 9 # Number of statistics (only the first 6 are used in the misspecification experiment)

# Sampling from the uniform prior
G0_start <- 1e-10
G0_end <- 3e-9
G0 <- runif(M, G0_start, G0_end)

T_start <- 1e-9
T_end <- 2e-8
T <- runif(M, T_start, T_end)

lambda_start <- 1e8
lambda_end <- 4e9
lambda <- runif(M, lambda_start, lambda_end)

sigma_start <- 1e-13
sigma_end <- 1e-12
sigma2_W <- runif(M, sigma_start, sigma_end)
sigma2_W <- 0

s_sim <- array(0, dim = c(M, nStats))

# Running model evaluations
for (i in 1:M) {
  Y_sim <- TurinModel(G0[i], T[i], lambda[i], sigma2_W, B = B, Ns = Ns, N = nReal, tau0 = 0)
  m_sim <- log(temporalMomentsGeneral(Y_sim, 3))
  s_sim[i, ] <- summaryMoments(m_sim) 
  print(i)
}

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

saveRDS(cbind(G0, T, lambda, sigma2_W), file = "Turin_param_miss")

saveRDS(s_sim, file = "Turin_ssim_miss")
