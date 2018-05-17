# NEE memory model for each Flux site
# data files for Pinus Edulis (PIED) growing near Montrose, CO:
# ring widths (mm): BUGS_data_widths_Montrose_PIED
# total monthly precipitation (mm): BUGS_data_PPT_montrose
# mean daily temperature (C): BUGS_data_Tave_montose
# sample sizes and time-scale info: BUGS_data_sample_sizes
# In the PPT and Tave data, the variables are tave[r,c], where r (row)
# represents a year (as indicated by Year.tave or Year.ppt), and c is the 
# month (site = 1 for montrose datasets).

# Model code. Note that this model set "antecedent importance" weights = 0 in year 1 
#for Oct-Dec (post-growing season), and uses varying time periods (for weigths) 
#for each year such that yrs 1-2 = 12 blocks of 1 month, year 3 = 6 blocks of 2 months, 
#yrs 4-5 = 4 blocks of 3 months. Total # weights = 2*12 + 1*6 + 2*4 = 38.

nee_SAM_model <- function(){
  # Monitor the model deviance:
  deviance <- -2 * sum(loglike[1:Nmem])
  
  # Likelihood and mean model, looping over the NEE records
  for(r in 1:Nmem){ # r is the t in the supplemental material model description
    # Monitor the log-likelihood:
    loglike[r] <- 0.5*log(tau_y/(2*3.141593))-((tau_y/2)*pow(NEE[Mem_records[r]]-muNEE[r],2))
    # Likelihood for NEE data:
    muNEE[r] <- sum(effect[, r])
    # Replicated data for evaluating model fit:
    NEE_pred[r] ~ ddexp(muNEE[r], tau_y)
    NEE_resi[r] <- NEE_pred[r] - NEE[Mem_records[r]]
    ag_part[r] <- (1 + phi0 * (NDVI[NDVI_index[Mem_records[r]]] - 1)) * NDVI[NDVI_index[Mem_records[r]]] 
    
    # the climate effects
    for(i in 1:Nsen){
      effect[i, r] <- (1 - ag_part[r]) * an[i] * antClim[i, r] + ag_part[r] * ag[i] * antClim[i, r]
    }
    antClim[1, r] <- 1
    antClim[22, r] <- antAP[r]
    for(e in 1:5){
      antClim[(e+1), r] <- antA[e, r]
      antClim[(e+6), r] <- antA[e, r] * antA[e, r] 
      Sen_main[e, r] <- (1 - ag_part[r]) * an[(e+1)]  + ag_part[r] * ag[(e+1)] 
      Sen_quad[e, r] <- (1 - ag_part[r]) * an[(e+6)] * antClim[(e+1), r] + ag_part[r] * ag[(e+6)] * antClim[(e+1), r]
      Sen_inter[e, r] <- (e==1)*((1 - ag_part[r]) * an[12] * antClim[3, r] + ag_part[r] * ag[12] * antClim[3, r] + 
                                   (1 - ag_part[r]) * an[13] * antClim[4, r] + ag_part[r] * ag[13] * antClim[4, r] +
                                   (1 - ag_part[r]) * an[14] * antClim[5, r] + ag_part[r] * ag[14] * antClim[5, r] +
                                   (1 - ag_part[r]) * an[15] * antClim[6, r] + ag_part[r] * ag[15] * antClim[6, r]) +
        (e==2)*((1 - ag_part[r]) * an[16] * antClim[4, r] + ag_part[r] * ag[16] * antClim[4, r] + 
                  (1 - ag_part[r]) * an[17] * antClim[5, r] + ag_part[r] * ag[17] * antClim[5, r] +
                  (1 - ag_part[r]) * an[18] * antClim[6, r] + ag_part[r] * ag[18] * antClim[6, r] +
                  (1 - ag_part[r]) * an[12] * antClim[2, r] + ag_part[r] * ag[12] * antClim[2, r]) +
        (e==3)*((1 - ag_part[r]) * an[19] * antClim[5, r] + ag_part[r] * ag[19] * antClim[5, r] + 
                  (1 - ag_part[r]) * an[20] * antClim[6, r] + ag_part[r] * ag[20] * antClim[6, r] +
                  (1 - ag_part[r]) * an[13] * antClim[2, r] + ag_part[r] * ag[13] * antClim[2, r] +
                  (1 - ag_part[r]) * an[16] * antClim[3, r] + ag_part[r] * ag[16] * antClim[3, r]) +
        (e==4)*((1 - ag_part[r]) * an[21] * antClim[6, r] + ag_part[r] * ag[21] * antClim[6, r] + 
                  (1 - ag_part[r]) * an[14] * antClim[2, r] + ag_part[r] * ag[14] * antClim[2, r] +
                  (1 - ag_part[r]) * an[17] * antClim[3, r] + ag_part[r] * ag[17] * antClim[3, r] +
                  (1 - ag_part[r]) * an[19] * antClim[4, r] + ag_part[r] * ag[19] * antClim[4, r]) +
        (e==5)*((1 - ag_part[r]) * an[15] * antClim[2, r] + ag_part[r] * ag[15] * antClim[2, r] + 
                  (1 - ag_part[r]) * an[18] * antClim[3, r] + ag_part[r] * ag[18] * antClim[3, r] +
                  (1 - ag_part[r]) * an[20] * antClim[4, r] + ag_part[r] * ag[20] * antClim[4, r] +
                  (1 - ag_part[r]) * an[21] * antClim[5, r] + ag_part[r] * ag[21] * antClim[5, r])
      
      Sen[e, r] <- sum(Sen_main[e, r] + Sen_quad[e, r] + Sen_inter[e, r])
      RSen[e, r] <- Sen[e, r] * abs(NEE[Mem_records[r]]) * Nmem/SUMNEE
    }
    Sen[6, r] <- (1 - ag_part[r]) * an[22] + ag_part[r] * ag[22]
    RSen[6, r] <- Sen[6, r] * abs(NEE[Mem_records[r]]) * Nmem/SUMNEE
    Sen[7, r] <- Sen[4, r] + Sen[5, r]
    RSen[7, r] <- RSen[4, r] + RSen[5, r]
    
    for(e in 1:4){
      for(f in (e+1):5){
        antClim[((e==1)*(f+10)+(e==2)*(f+13)+(e==3)*(f+15)+(e==4)*(f+16)), r] <- antA[e, r] * antA[f, r]}
    }
  }
  
  
  # Calculate the overall (with and without rescaling) climate effects (or the intercepts when e=1) 
  for(e in 1:7){
    ESen[e] <- mean(Sen[e,])
    ERSen[e] <- mean(RSen[e,])
  }
  
  # Calculate the antecedent climate covariates
  for(r in 1:Nmem){
    for(v in 1:Nv){
      for(d in 1:Nlag){antAD[d,v,r] <- weightA[v, d] * clim[(Mem_records[r]+1-d), v]}
      antA[v, r] <- sum(antAD[,v,r])
    }
    for(d in 1:NlagP){
      antAPD[d,r] <- weightAP[d] * ppt_multiscale[Mem_records[r], d]
    }
    antAP[r] <- sum(antAPD[,r])
  }
  
  for(d in 1:NlagP){
    # !!!!!! the use of delta and weight is flipped here, fix before archiving
    deltaAP[d] <- (d>1)*deltaXAP[d] + (d==1)*0
    weightAP[d] <- deltaAP[d]/sum(deltaAP[])
  }
  
  # define the weights and calculate cummulative weights
  for(v in 1:Nv){
    sumDA[v] <- sum(deltaA[v,])
    for(l in 1:Nlag){
      deltaA[v, l] <- (v<4)*deltaXA[v, block[l]]/BlockSize[l] + (v==4)*(l==1) + 
        (v==5)*(l>1)*deltaXA[v, block[l]]/BlockSize[l] #(l<8)*
      weightA[v, l] <- deltaA[v, l]/sumDA[v]
      cum_weightA[v, l] <- sum(weightA[v, 1:l])
    }
  }
  
  for(l in 1:NlagP){cum_weightAP[l] <- sum(weightAP[1:l])}
  
  # standard, relatively non-informative priors for climate weights, climate effects (an and ag)
  for(v in 1:Nv){deltaXA[v, 1:Nblocks] ~ ddirch(rep(1, Nblocks))}
  deltaXAP[1:NblocksP] ~ ddirch(rep(1, NblocksP))
    for(i in 1:Nsen){
      an[i] ~ dnorm(0, 0.001)
      ag[i] ~ dnorm(0, 0.001)
      # Monitor the Bayesian p-value of the climate effects:
      p_an[i] <- (an[i] > 0) * 1
      p_ag[i] <- (ag[i] > 0) * 1
    }
  sig_y ~ dunif(0, 100)
  tau_y <- pow(sig_y, -2)
  phi0 ~ dunif(-1, 1)
}