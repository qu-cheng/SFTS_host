source("0_Library_and_functions.R")


cores = detectCores()
cl <- makeCluster(cores) #register cores
registerDoParallel(cl, cores = cores)

for(survey.i in 1:9)
{
  sens.par <- read.csv(paste("3_Parameter_permutation/3_2_Perturbed_parameters/Survey_", survey.i, "_Sensitivity.csv", sep = "")) 
  
  result <- foreach(i = 1:nrow(sens.par), .combine = rbind, .packages = c("tidyverse"), .errorhandling = "pass") %dopar% {
    
    if(i %% 5000 == 1){
      system(paste("echo 'now processing survey", survey.i, " rep ", i, " ", Sys.time(),  "'"))
    }
    
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
    
    current.sample <- sens.par[i,] %>%
      pivot_longer(cols = !c(SurveyID, ID, mu_T, phi, change.par, change.value),
                   names_to = c(".value", "Host"),
                   names_pattern = "(.*)_(.*)") 
    
    current.sample$alpha <- current.sample$beta * current.sample$chi_H
    current.sample$rho_H <- current.sample$v_H/(1 - current.sample$v_H) * (current.sample$mu_H + current.sample$sigma_H)
    
    tryCatch({  mod <- SFTSV_multihost$new(n_host = nrow(current.sample),       # number of host species
                                           chi_H = current.sample$chi_H,     # number of individuals of each species/number of ticks
                                           beta = current.sample$beta,     # biting rate ticks/host
                                           mu_H = current.sample$mu_H,     # mortality rate of each host species
                                           gamma_H = current.sample$gamma_H,    # incubation period of each host species
                                           sigma_H = current.sample$sigma_H,    # duration of viremia
                                           v_H = current.sample$v_H,  # case-fatality rate of each species
                                           mu_T = unique(current.sample$mu_T),     # mortality rate of ticks
                                           phi = unique(current.sample$phi))      # proportion of eggs infected transovarially
    
    t <- seq(0, 10000, length.out = 10000)
    out <- mod$run(t)
    out <- mod$transform_variables(out) 
    
    c(R = tail(out$R_H, 1))},
    error = function(err){
      rep(NA, nrow(current.sample))
    })
  }
  
  write.csv(result, paste("3_Parameter_permutation/3_4_Perturbed_results/Survey_", survey.i, "_Sensitivity_results.csv", sep = ""), row.names = FALSE)
}

stopCluster(cl)












