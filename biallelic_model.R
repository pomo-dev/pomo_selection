# normalization constant
# normalizes the stationary distribution
# N is the virtual population size
# pi is the proportion of GC alleles
# rho is the echangeability rate between GC and AT alleles
# sigma is the selection coefficion of GC alleles

normalization_constant <- function(N,pi,rho,sigma) {

  n  <- (N-1):1
  nc <- sum(pi*(1+sigma)^(N-1),1-pi,pi*(1-pi)*rho*(1+sigma)^(n-1)*N/(n*(N-n)))
  return(nc)
  
}

normalization_constant(10,0.5,0.001,0.01)

# calculates the stationary distribution for the biallelic case
# N is the virtual population size
# pi is the proportion of reference allele
# rho is the echangeability rate between alleles
# sigma is the selection coefficion of the reference allele

stationary_distribution <- function(N,pi,rho,sigma){
  
  # polymorphic sites
  n <- (N-1):1
  poly <- (1+sigma)^(n-1)*N/(n*(N-n))
  poly <- pi*(1-pi)*rho*poly

  # include monomorphic sites
  sd <- c(pi*(1+sigma)^(N-1),poly,1-pi)
  print(sd)
  
  # normalization
  sd <- sd/(normalization_constant(N,pi,rho,sigma))
  return(sd)
}

sd <- stationary_distribution(100,0.5,0.001,0.01)
plot(sd)


# calculates usefull quantities from counts
# these quantities are needed to calculate the 
# likelihood function
counts_to_data <- function(counts){
  
  N <- length(counts)-1
  S <- sum(counts)
  number_Ni <- counts[1] 
  number_Nj <- counts[length(counts)]
  number_poly <- S-number_Ni-number_Nj
  sum_poly <- sum(counts[-c(1,length(counts))]*((N-1):1))
  
  return(list(N=N,S=S,number_Ni=number_Ni,number_Nj=number_Nj,number_poly=number_poly ,sum_poly=sum_poly))
  
} 

pi     <- 0.5
rho    <- 0.001
sigma  <- 0
counts <- c(10,1,1,1,1,1,1,1,1,1,1,1,1,10)
data   <- counts_to_data(counts)


# log posterior of pi
log_posterior_pi <- function(data,pi,rho,sigma){
  lp <- -data$S*log(normalization_constant(data$N,pi,rho,sigma)) + (data$S-data$number_Nj)*log(pi)+(data$number_Nj+data$number_poly)*log(1-pi)
  return(lp)
}
log_posterior_pi(data,pi,rho,sigma)

# log posterior of rho
#-data$S*log(normalization_constant(data$N,pi,rho,sigma))
log_posterior_rho <- function(data,pi,rho,sigma){
  lp <- -data$S*log(normalization_constant(data$N,pi,rho,sigma)) + (data$number_poly)*log(rho)

  return(lp)
}
log_posterior_rho(data,pi,rho,sigma)

# log posterior of sigma
log_posterior_sigma <- function(data,pi,rho,sigma){
  lp <- -data$S*log(normalization_constant(data$N,pi,rho,sigma)) + ((data$N-1)*data$number_Ni+data$sum_poly)*log(1+sigma)
  return(lp)
}
log_posterior_rho(data,pi,rho,sigma)

#likelihood
likelihood <- function(data,pi,rho,sigma){
  lp <- -data$S*log(normalization_constant(data$N,pi,rho,sigma))+((data$N-1)*data$number_Ni+data$sum_poly)*log(1+sigma)+
                                                                 (data$number_poly)*log(rho)+
                                                                 (data$S-data$number_Nj)*log(pi)+(data$number_Nj+data$number_poly)*log(1-pi)
  return(lp)
}

# metropolis hastings step pi with beta proposal
# proposal_parameter is the variance of a beta distribution
mh_pi <- function(data,pi,rho,sigma,proposal_parameter_p){
  
  lp1   <- log_posterior_pi(data,pi,rho,sigma)
  pi2   <- rbeta(1,shape1=pi*(1-pi-proposal_parameter_p)/proposal_parameter_p,shape2=(1-pi)*(1-pi-proposal_parameter_p)/proposal_parameter_p)
  lp2   <- log_posterior_pi(data,pi2,rho,sigma)
  alpha <- exp(lp2-lp1)
  
  if (alpha > runif(1,0,1)){
    return(pi2)
  } else {
    return(pi)
  }
}
proposal_parameter_p <- 0.001
mh_pi(data,pi,rho,sigma,proposal_parameter_p)

# metropolis hastings step rho with multiplier proposal
# proposal_parameter is a fine tunning parameter
mh_rho <- function(data,pi,rho,sigma,proposal_parameter_r){
  
  lp1   <- log_posterior_rho(data,pi,rho,sigma)
  rho2  <- rho*exp(proposal_parameter_r*(runif(1,0,1)-0.5))
  lp2   <- log_posterior_rho(data,pi,rho2,sigma)
  alpha <- exp(lp2-lp1)
  
  if (alpha > runif(1,0,1)){
    return(rho2)
  } else {
    return(rho)
  }
}
proposal_parameter_r <- 0.1
mh_rho(data,pi,rho,sigma,proposal_parameter_r)

# metropolis hastings step sigma with multiplier proposal
# proposal_parameter is a fine tunning parameter
mh_sigma <- function(data,pi,rho,sigma,proposal_parameter_s){
  
  lp1    <- log_posterior_sigma(data,pi,rho,sigma)
  sigma2 <- (1+sigma)*exp(proposal_parameter_s*(runif(1,0,1)-0.5))-1
  lp2    <- log_posterior_sigma(data,pi,rho,sigma2)
  alpha  <- exp(lp2-lp1)
  
  if (alpha > runif(1,0,1)){
    return(sigma2)
  } else {
    return(sigma)
  }
}
proposal_parameter_s <- 0.1
mh_sigma(data,pi,rho,sigma,proposal_parameter_s)

pi_sampler <- function(data,rho,sigma,r_pi){
  c_pi <- rep(NA,length(r_pi))
  for (i in 1:length(r_pi)){
    c_pi[i] <- log_posterior_pi(data,r_pi[i],rho,sigma)
  }
  n_pi <- sample(r_pi,1,prob=exp(c_pi-max(c_pi)))
  return(n_pi)
}

rho_sampler <- function(data,pi,sigma,r_rho){
  c_rho <- rep(NA,length(r_rho))
  for (i in 1:length(r_rho)){
    c_rho[i] <- log_posterior_rho(data,pi,r_rho[i],sigma)
  }
  n_rho <- sample(r_rho,1,prob=exp(c_rho-max(c_rho)))
  return(n_rho)
}
rho_sampler(data,pi,sigma,r_rho)
plot(r_rho,exp(c_rho))

sigma_sampler <- function(data,pi,rho,r_sigma){
  c_sigma <- rep(NA,length(r_sigma))
  for (i in 1:length(r_sigma)){
    c_sigma[i] <- log_posterior_sigma(data,pi,rho,r_sigma[i])
  }
  n_sigma <- sample(r_sigma,1,prob=exp(c_sigma-max(c_sigma)))
  return(n_sigma)
}
sigma_sampler(data,pi,rho,r_sigma)

# mcmc scheme 
# generations
mcmc <- function(generations,counts,proposal_parameter_p,proposal_parameter_r,proposal_parameter_s,initial_parameters){
  
  p_mcmc <- matrix(NA,ncol=4,nrow=generations)
  
  data <- counts_to_data(counts)
  
  pi    <- initial_parameters[1]
  rho   <- initial_parameters[2]
  sigma <- initial_parameters[3]
  
  for (i in 1:generations){
    
    pi    <- mh_pi(data,pi,rho,sigma,proposal_parameter_p)
    rho   <- mh_rho(data,pi,rho,sigma,proposal_parameter_r)
    sigma <- mh_sigma(data,pi,rho,sigma,proposal_parameter_s)
    print(paste0("pi: ",round(pi,2)," rho: ",round(rho,5)," sigma: ",round(sigma,2)))
    
    p_mcmc[i,] <- c(pi,rho,sigma,lp_total(data,pi,rho,sigma))
  
  }
  return(p_mcmc)
}

# gibbs sampler 
# generations
mcmc_gibbs <- function(generations,counts,r_pi,r_rho,r_sigma){
  
  p_mcmc <- matrix(NA,ncol=4,nrow=generations)
  
  data <- counts_to_data(counts)
  
  pi    <- 0.5
  rho   <- 0.1
  sigma <- 0.5

  for (i in 1:generations){
    
    pi    <- pi_sampler(data,rho,sigma,r_pi)
    rho   <- rho_sampler(data,pi,sigma,r_rho)
    sigma <- sigma_sampler(data,pi,rho,r_sigma)

    p_mcmc[i,] <- c(pi,rho,sigma,lp_total(data,pi,rho,sigma))
    
  }
  return(p_mcmc)
}

# maximum likelihood estimates
# r_pi, r_rho and r_sigma are the coordintates (pi,rho,sigma) 
# the likelihood will be assessed
# return the likelihood surface evalutated on r_parameter points

likelihood_surface <- function(r_pi,r_rho,r_sigma,counts) {
  
  data <- counts_to_data(counts)
  lik  <- matrix(NA,ncol=4,nrow=length(r_pi)*length(r_rho)*length(r_sigma))
  p    <- 1
  for ( i in r_pi){
    for (j in r_rho){
      for(k in r_sigma) {
        lik[p,] <- c(i,j,k,likelihood(data,i,j,k))
        p <- p+1
      }
    }
  }
  return(lik)
}

r_pi <- seq(1/100,99/100,length.out = 75)
r_rho <- exp(-seq(1/100,15,length.out = 50))
r_sigma <- seq(-0.5,0.5,0.01)

counts <- 1000000*stationary_distribution(100,pi0,rho0,sigma0)
likelihood_surface(r_pi,r_rho,r_sigma,counts)
