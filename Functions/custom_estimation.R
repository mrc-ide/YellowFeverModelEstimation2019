
GLM_MCMC_step = function(param, x, y, chain_cov, adapt,  accCurrent) {
  
  ### propose new param ###
  param_prop = YFestimation::GLMproposal(param, chain_cov, adapt)
  
  ### priors ###
  prior_prop = sum( GLMprior(param_prop) )
  
  ### if prior finite, evaluate likelihood ###
  if (is.finite(prior_prop)) {
    
    like_prop = YFestimation::GLMlike(param_prop, x, y)
    
    ### accept prob ###
    accProp = like_prop + prior_prop
    
    p_accept= accProp - accCurrent
    if(is.na(p_accept) ){
      p_accept = -Inf
    }
  } else {
    p_accept = log(0)
  }
  
  ## accept/reject step:
  tmp = (runif(1))
  if( tmp< min( exp(p_accept), 1)  ) { # accept:
    param = param_prop
    accCurrent = accProp
    accept = 1
  } else { # reject:
    accept = 0
  }
  
  return(list(param = param,  accCurrent = accCurrent, accept = accept)  )
}



GLM_MCMC = function(Niter, name_dir, pars_ini, x, y, plot_chain, run_id = 1){
  
  burnin = min(2*length(pars_ini), Niter)
  
  ### find posterior probability at start ###
  out = GLM_MCMC_step(pars_ini, x, y, chain_cov=1, adapt=0, accCurrent=-Inf)
  
  chain = posteriorProb = acceptRate = NULL
  
  fileIndex = 0
  iter = 1
  
  for (iter in iter:Niter){
    #save current step
    param = out$param;  names(param) = names(pars_ini);  accCurrent = out$accCurrent;  accept = out$accept
    
    chain = rbind(chain, param);
    posteriorProb = rbind(posteriorProb, accCurrent)
    acceptRate = rbind(acceptRate, accept)
    
    if(iter>100 & plot_chain) plot(chain[,1] , type= "l")
    
    if (iter %% 1000 == 0){
      
      colnames(acceptRate) = "acceptRate"
      colnames(posteriorProb) = "posteriorProb"
      if (iter %% 10000 == 0){
        fileIndex  = iter/10000
      }
      
      out_state = cbind(chain,  posteriorProb, acceptRate)[min((fileIndex * 10000+1),iter):iter,]
      
      write.csv(out_state, paste0(name_dir,"/","GLM_chain",run_id,"_output_",fileIndex, " .csv") )
    }
    
    #adapt?
    if (iter>burnin & runif(1)<0.9){ #adapt
      adapt = 1
      chain_cov  = cov(chain[max(nrow(chain)-10000, 1):nrow(chain),])
    } else {
      adapt = 0
      chain_cov = 1
    }
    
    #new step
    out = GLM_MCMC_step(param, x, y, chain_cov, adapt, accCurrent)
    
  }
  
}


GLMprior = function(param) {
  
  Prior = rep(0,2)

  Prior[1] =  sum(dnorm(param[grep("predicted", names(param), invert = TRUE)],
                        mean = 0,
                        sd = 30,
                        log = TRUE))
  Prior[2] = dunif(param[grep("predicted", names(param))],
                   min = 0, max = 10, log = TRUE)
  
  # this term is for normally distributed non-country parameters : normal distrib with high sd (flat distrib)
  out = as.numeric( Prior )
  return( out )
}
