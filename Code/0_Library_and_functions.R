library(tidyverse)
library(Hmisc)
library(odin)
library(doSNOW)
library(doParallel)
library(data.table)
library(cowplot)
library(sf)
library(ggmap)
library(ggbeeswarm)
library(ggthemes)
library(ggridges)
library(sp)
library(maptools)
library(spatstat)
library(packcircles)
library(FRK)
library(ggraph)
library(ggrepel)
library(ranger)



#===============================================================
#                        odin model 
#===============================================================

SFTSV_multihost <- odin::odin({
  
  n_host <- user()
  
  ## Params
  # estimate excessive moratality rate based on case fatality rate
  rho_H[] <- v_H[i]/(1 - v_H[i]) * (mu_H[i] + sigma_H[i])
  
  # estimate beta_TH[i] based on beta_HT[i] and the number of each species
  alpha[] <- beta[i] * chi_H[i]
  lambda_H[] <- alpha[i] * I_H[i]
  lambda <- sum(lambda_H[])
  
  ## Derivatives
  # hosts
  deriv(S_H[]) <- mu_H[i] - beta[i]*S_H[i]*I_T - mu_H[i]*S_H[i]
  deriv(E_H[]) <- beta[i]*S_H[i]*I_T - (gamma_H[i] + mu_H[i])*E_H[i]
  deriv(I_H[]) <- gamma_H[i]*E_H[i] - (sigma_H[i] + rho_H[i] + mu_H[i])*I_H[i]
  deriv(R_H[]) <- sigma_H[i]*I_H[i] - mu_H[i]*R_H[i]
  
  # ticks
  deriv(S_T) <- mu_T*S_T + mu_T*(1-phi)*I_T - lambda*S_T - mu_T*S_T
  deriv(I_T) <- mu_T*I_T*phi + lambda*S_T - mu_T*I_T
  
  ## Initial conditions
  # hosts
  initial(S_H[]) <- 1 - 1E-6
  initial(E_H[]) <- 0
  initial(I_H[]) <- 1E-6
  initial(R_H[]) <- 0
  
  # ticks
  initial(S_T) <- 1
  initial(I_T) <- 0
  
  ## parameters
  chi_H[] <- user()         # number of individuals of host species i/number of ticks
  
  beta[] <- user()   # on average, each host species i contact beta_HT[i] ticks
  mu_H[] <- user()      # natural mortality rate of host species i
  gamma_H[] <- user()   # incubation period of host species i
  sigma_H[] <- user()   # recovery rate of host species i
  v_H[] <- user()     # case fatality rate due to SFTSV of host species i
  
  mu_T <- user()        # mortality rate of ticks
  phi <- user()           # a proportion k of the eggs laid by infectious ticks are infectious through vertical transmission
  
  
  
  ## dimensions
  dim(chi_H) <- n_host
  dim(beta) <- n_host
  dim(alpha) <- n_host
  dim(mu_H) <- n_host
  dim(gamma_H) <- n_host
  dim(sigma_H) <- n_host
  dim(v_H) <- n_host
  dim(rho_H) <- n_host
  dim(lambda_H) <- n_host
  
  dim(S_H) <- n_host
  dim(E_H) <- n_host
  dim(I_H) <- n_host
  dim(R_H) <- n_host
})




#===============================================================
#                    simulate function
#===============================================================
simulate.fun <- function(n.species,     # number of species for the odin model
                         n.times,       # number of times to simulate
                         i,             # ID for the core
                         current.sampler,  # sampling function
                         keep.k)        # only output results when at least keep.k criteria were met; to save storage
{
  result.core <- NULL
  
  for(n.rep in 1:n.times)
  {
    if(n.rep %% 5000 == 1){
      system(paste("echo 'now processing survey ", survey.i, " rep ", n.rep, " on core", i, " ", Sys.time(),  "'"))
    }
    
    
    
    current.sample <- current.sampler(n.species)
    current.sample$alpha <- current.sample$beta * current.sample$chi_H
    current.sample$rho_H <- current.sample$v_H/(1 - current.sample$v_H) * (current.sample$mu_H + current.sample$sigma_H)
    
    
    K_mat <- matrix(0, n.species + 1, n.species + 1)
    colnames(K_mat) <- rownames(K_mat) <- c(paste("E_H", 1:n.species, sep = ""), "I_T")
    K_mat[n.species + 1, n.species + 1] <- current.sample$phi
    K_mat[1:n.species, n.species + 1] <- current.sample$beta * current.sample$chi_H / current.sample$mu_T
    K_mat[n.species + 1, 1:n.species] <- current.sample$alpha / current.sample$chi_H * current.sample$gamma_H / (current.sample$mu_H + current.sample$gamma_H) / (current.sample$mu_H + current.sample$rho_H + current.sample$sigma_H)
    
    R0 <- Re(eigen(K_mat)$values[1])
    
    if(R0 > 1)     # only run the model when R0 > 1
    {
      
      # plos comp bio formula
      R_i <- R_no_i <- array(0, n.species)
      for(k in 1:n.species)
      {
        P_V <- matrix(0, ncol = n.species + 1, nrow = n.species + 1)
        P_V[n.species + 1, n.species + 1] <- 1
        
        I_V <- diag(n.species + 1) - P_V
        
        P_A <- matrix(0, ncol = n.species + 1, nrow = n.species + 1)
        P_A[k, k] <- 1
        
        R_i[k] <- Re(eigen((P_V + P_A) %*% K_mat)$values[1])
        R_no_i[k] <- Re(eigen((P_V + I_V - P_A) %*% K_mat)$values[1])
      }
      
      mod <- SFTSV_multihost$new(n_host = n.species,       # number of host species
                                 chi_H = current.sample$chi_H,     # number of individuals of each species/number of ticks
                                 beta = current.sample$beta,     # biting rate ticks/host
                                 mu_H = current.sample$mu_H,     # mortality rate of each host species
                                 gamma_H = current.sample$gamma_H,    # incubation period of each host species
                                 sigma_H = current.sample$sigma_H,    # duration of viremia
                                 v_H = current.sample$v_H,  # case-fatality rate of each species
                                 mu_T = current.sample$mu_T,     # mortality rate of ticks
                                 phi = current.sample$phi)      # proportion of eggs infected transovarially
      
      t <- seq(0, 10000, length.out = 10000)
      out <- mod$run(t)
      out <- mod$transform_variables(out) 
      
      
      current.result <- c(unlist(current.sample), R = tail(out$R_H, 1), I = tail(out$I_H, 1), R0 = R0, R_i = R_i, R_no_i = R_no_i)
      result.core <- rbind(result.core, current.result)
    }
  }
  
  R.matrix <- select(as.data.frame(result.core), matches("R[1-9]"))
  # pass.ID <- apply(R.matrix, 1, function(x){sum(x < Animal_prev$prev_high & x > Animal_prev$prev_low) >= 2})
  
  pass.ID <- apply(R.matrix, 1, function(x){sum(x < Animal_prev$prev_high & x > Animal_prev$prev_low) >= keep.k})
  
  result.core[pass.ID,]
}




#===============================================================
#             estimate ranges of beta for R2
#===============================================================
beta_for_R2 <- function(R1,                 # results from R1
                        Animal_prev,        # Animal prevalence
                        beta.min.value,     # minimum values for betas
                        beta.max.value)     # maximum values for betas
{
  R.final <- select(R1, matches("R[1-9]"))
  n.species <- ncol(R.final)
  
  pass.mat <- t(apply(R.final, 1, function(x) x <= Animal_prev$prev_high & x >= Animal_prev$prev_low))
  colnames(pass.mat) <- paste("R", 1:n.species, ".pass", sep = "")
  
  
  R1 <- cbind(R1, pass.mat) %>% as.data.table()
  R1[, pass.no := rowSums(.SD), .SDcols = grep(".pass", names(R1), value = TRUE)]
  
  # print(table(R1$pass.no))
  
  
  beta.range <- data.frame(betai = 1:n.species, beta.min = rep(beta.min.value, n.species), beta.max = rep(beta.max.value, n.species))
  for(pass.i in 2:n.species)
  {
    cat(pass.i, "\n")
    all.comb <- combn(n.species, pass.i)
    
    for(comb.n in 1:ncol(all.comb))
    {
      current.beta <- R1[, current.pass := rowSums(.SD), .SDcols = paste("R", all.comb[,comb.n], ".pass", sep = "")][current.pass == pass.i, .SD, .SDcols = grep("beta", names(R1), value = TRUE)]
      
      if(nrow(current.beta) > 100)
      {
        current.beta.min <- apply(current.beta, 2, min)
        current.beta.max <- apply(current.beta, 2, max)
        
        beta.range$beta.min <- ifelse(beta.range$beta.min <= current.beta.min, current.beta.min, beta.range$beta.min)
        beta.range$beta.max <- ifelse(beta.range$beta.max >= current.beta.max, current.beta.max, beta.range$beta.max)
      }
    }
  }
  
  beta.range %>%
    mutate(beta.min = beta.min * 0.8,
           beta.max = pmin(beta.max * 1.2, beta.max.value))
}




#===============================================================
#             estimate chi_max * beta for R3
#===============================================================
CB_high_R3 <- function(R2)    # results from R2
{
  R.final <- select(R2, matches("R[1-9]"))
  n.species <- ncol(R.final)
  
  pass.mat <- t(apply(R.final, 1, function(x) x <= Animal_prev$prev_high & x >= Animal_prev$prev_low))
  colnames(pass.mat) <- paste("R", 1:n.species, ".pass", sep = "")
  
  
  R2 <- cbind(R2, pass.mat) %>% as.data.table()
  R2[, pass.no := rowSums(.SD), .SDcols = grep(".pass", names(R2), value = TRUE)]
  
  # print(table(R2$pass.no))
  
  host.lowest.prev <- which.min(Animal_prev$prev_high)
  
  R2 %>%
    filter(pass.no == nrow(Animal_prev)) %>%
    select(paste("beta", host.lowest.prev, sep = ""), paste("chi_H", which.max(Animal_prev$N), sep = "")) %>%
    mutate(product = get(paste("chi_H", which.max(Animal_prev$N), sep = ""))*get(paste("beta", host.lowest.prev, sep = ""))) %>%
    pull(product) %>%
    max()
}



#===============================================================
#           get passes  
#===============================================================
get_pass <- function(R.result, Animal_prev)    # Results from any round
{
  R.final <- select(R.result, matches("R[1-9]"))
  n.species <- ncol(R.final)
  
  pass.mat <- t(apply(R.final, 1, function(x) x <= Animal_prev$prev_high & x >= Animal_prev$prev_low))
  colnames(pass.mat) <- paste("R", 1:n.species, ".pass", sep = "")
  
  c(1:nrow(R.result))[rowSums(pass.mat) == n.species]
}



#===============================================================
#           get the number of passes
#===============================================================
get_pass_no <- function(R.result)    # Results from any round
{
  R.final <- select(R.result, matches("R[1-9]"))
  n.species <- ncol(R.final)
  
  pass.mat <- t(apply(R.final, 1, function(x) x <= Animal_prev$prev_high & x >= Animal_prev$prev_low))
  colnames(pass.mat) <- paste("R", 1:n.species, ".pass", sep = "")
  
  sum(rowSums(pass.mat) == n.species)
}






#===============================================================
#              Sample function for round 1
#===============================================================
sample.function1 <- function(n)   # n - number of species
{
  beta <- runif(n, 0, 2)
  
  current_animal_par <- Animal_par %>%
    filter(Animal %in% Animal_prev$Host.species)
  
  chi_max <- runif(1, 0, 1)      # sample the reference X_max
  N_max <- max(Animal_prev$N)
  
  chi_H <- chi_max/N_max*Animal_prev$N
  
  
  mu_H <- 1/(current_animal_par$life_span_mean*365)
  gamma_H <- runif(n, 1/current_animal_par$latent_period_U, 1/current_animal_par$latent_period_L)
  sigma_H <- runif(n, 1/current_animal_par$infectious_period_U, 1/current_animal_par$infectious_period_L)
  v_H <- runif(n, current_animal_par$CFR_L, current_animal_par$CFR_U)
  
  mu_T <- runif(1, 0.000913242, 0.1)   # 10 days to 3 years
  phi <- runif(1, 0, 1)
  
  list(chi_H = chi_H,
       beta = beta,
       mu_H = mu_H,
       gamma_H = gamma_H,
       sigma_H = sigma_H,
       v_H = v_H,
       mu_T = mu_T,
       phi = phi)
}







#===============================================================
#              Sample function for round 2
#===============================================================
sample.function2 <- function(n)    # n - number of species
{
  beta <- runif(n, beta.range$beta.min, beta.range$beta.max)
  
  current_animal_par <- Animal_par %>%
    filter(Animal %in% Animal_prev$Host.species)
  
  chi_max <- runif(1, 0, 1)      # sample the reference X_max
  N_max <- max(Animal_prev$N)
  
  chi_H <- chi_max/N_max*Animal_prev$N
  
  
  mu_H <- 1/(current_animal_par$life_span_mean*365)
  gamma_H <- runif(n, 1/current_animal_par$latent_period_U, 1/current_animal_par$latent_period_L)
  sigma_H <- runif(n, 1/current_animal_par$infectious_period_U, 1/current_animal_par$infectious_period_L)
  v_H <- runif(n, current_animal_par$CFR_L, current_animal_par$CFR_U)
  
  mu_T <- runif(1, 0.000913242, 0.1)   # 10 days to 3 years
  phi <- runif(1, 0, 1)
  
  list(chi_H = chi_H,
       beta = beta,
       mu_H = mu_H,
       gamma_H = gamma_H,
       sigma_H = sigma_H,
       v_H = v_H,
       mu_T = mu_T,
       phi = phi)
}








sample.function3 <- function(n)
{
  beta <- runif(n, beta.range$beta.min, beta.range$beta.max)
  
  current_animal_par <- Animal_par %>%
    filter(Animal %in% Animal_prev$Host.species)
  
  chi_max <- runif(1, 0, min(1, Chi_max_high/beta[which.min(Animal_prev$prev_high)]))      # sample the reference X_max
  N_max <- max(Animal_prev$N)
  
  chi_H <- chi_max/N_max*Animal_prev$N
  
  
  mu_H <- 1/(current_animal_par$life_span_mean*365)
  gamma_H <- runif(n, 1/current_animal_par$latent_period_U, 1/current_animal_par$latent_period_L)
  sigma_H <- runif(n, 1/current_animal_par$infectious_period_U, 1/current_animal_par$infectious_period_L)
  v_H <- runif(n, current_animal_par$CFR_L, current_animal_par$CFR_U)
  
  mu_T <- runif(1, 0.000913242, 0.1)   # 10 days to 3 years
  phi <- runif(1, 0, 1)
  
  list(chi_H = chi_H,
       beta = beta,
       mu_H = mu_H,
       gamma_H = gamma_H,
       sigma_H = sigma_H,
       v_H = v_H,
       mu_T = mu_T,
       phi = phi)
}
