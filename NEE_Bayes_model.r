# Bayesian statistical model of daily NEE response to current and antecedent climate covariates, 
# including determining the antecedent climatic drivers following the stochastic antecedent
# modeling (SAM) framework (Ogle et al., 2015). Model is fit to one site at a time.

# ------------------ Inputs ------------------
# Nv = 5, for the five short-term climatic predictors in the order of: Tave, SW, VPD, SWCcurrent, and antSWC
# Nlag= 14 days into the past
# NlagP = 8, Tlag for precipitation (see Materials and Methods)
# Ns = 22, total number of climate covariates (see Materials and Methods)
# Nmem = A scalar that equals the record-length at a site subtracted by 365
#       (i.e. we can calculate memory terms starting at 366 days into the record)
# Mem_records = A vector (366, 367, ... and up to the record-length at the site), 
#       listing the indices of the days for which we calculate memory terms
# clim = A matrix of the daily climates (Tave, SW, VPD, SWC, and SWC). 
#       The dimension of the matrix is record-length by 5 (Nv). 
# ppt_multiscale = A matrix of precipitation at varying time periods into the past, see Methods for details. 
#       The dimension of the matrix is record-length by 8 (NlagP). 
# NEE = A vector, daily NEE values
# NDVI = A vector, bi-weekly NDVI values
# NDVI_index = A vector, looking up the corresponding bi-weekly NDVI for each day

# Varying time periods are used on the weights for the five short-term predictors: 
#       starting the current day up to 6 days previously = 7 blocks of 1 day each,
#       starting 7 days previously and up to 10 days previously = 2 blocks of 2 days each, 
#       starting 11 days previously and up to 13 days previously = 1 block of 3 days.
# Therefore below we have:
# block = c(1:7, rep(8:9, each = 2), rep(10, 3))
# BlockSize = c(rep(1, 7), rep(2, 4),rep(3, 3))

nee_SAM_model <- function(){
  # Likelihood and mean model, looping over daily NEE records at this site (starting at 366 days into the record)
  for(r in 1:Nmem){ # r is the t in the supplemental material model description
    # Likelihood for daily NEE data:
    NEE[Mem_records[r]] ~ ddexp(mu[r], tau_y)
    # Replicated data for evaluating model fit:
    NEE_pred[r] ~ ddexp(muNEE[r], tau_y)
    # Mean model: summing the effectClim (climate effects multiplied by the associated current 
    # and antecedent climate covariates)
    muNEE[r] <- sum(effectClim[, r])
    
    # phi: the mixture function for growing-season vs non-growing season effects, 
    # which changes over time depending on NDVI value; see below for more details
    phi[r] <- (1 + phi0 * (NDVI[NDVI_index[Mem_records[r]]] - 1)) * NDVI[NDVI_index[Mem_records[r]]] 
    
    # Compute the effectClim 
    for(i in 1:Ns){
      # For each climate covariate, the associated effect is modeled as a mixture of 
      #   growing-season effects (an) and non-growing season effects. Then the effectClim 
      #   is calcualted as climate effects multiplied by the associated current 
      #   and antecedent climate covariates
      effectClim[i, r] <- (1 - phi[r]) * an[i] * antClim[i, r] + phi[r] * ag[i] * antClim[i, r]
    }
    antClim[1, r] <- 1 # an[1] and ag[1] are intercepts
    antClim[22, r] <- antAP[r] # antecedent precipitation
    
    # Compute the 22 climate covariates, their quadratic terms, and their pairwise interactions
    for(e in 1:Nv){
      # main-effect terms of antTave, antSW, antVPD, SWCcurrent, and antSWC:
      antClim[(e+1), r] <- antA[e, r] 
      # quadratic terms of the aforementioned climatic predictors:
      antClim[(e+6), r] <- antA[e, r] * antA[e, r] 
    }  
    for(e in 1:(Nv-1)){
      for(f in (e+1):Nv){
        # interaction terms:
        antClim[((e==1)*(f+10)+(e==2)*(f+13)+(e==3)*(f+15)+(e==4)*(f+16)), r] <- antA[e, r] * antA[f, r]} 
    }
  }

  # Calculate the antecedent climate covariates
  for(r in 1:Nmem){
    # for the five short-term predictors:
    for(v in 1:Nv){
      for(d in 1:Nlag){
          antAD[d,v,r] <- weightA[v, d] * clim[(Mem_records[r]+1-d), v]
      }
      antA[v, r] <- sum(antAD[,v,r])
    }
    # for precipitation:
    for(d in 1:NlagP){
      antAPD[d,r] <- weightAP[d] * ppt_multiscale[Mem_records[r], d]
    }
    antAP[r] <- sum(antAPD[,r])
  }
  
  # Define the antecedent importance weights for the five short-term predictors: 
  #   antTave, antSW, antVPD, SWCcurrent, and antSWC.
  # Note: SWCcurrent is not an antecedent variable, however, it is included here for coding convenience, 
  #   no weights are estimated for it.
  # We used varying time periods/blocks (deltaA) for weights (weightA):
  #   i) seven daily blocks, starting current day to 6 days previously (current week, l = 1,..., and 7)
  #   ii) two 2-day blocks, starting 7 days previously to 10 days previously (l = 8,.., and 11)
  #   iii) one 3-day block, starting 11 days previously to 13 days previously (l = 12, 13, and 14).
  #   see comment section "Inputs" above for more details.
  for(v in 1:Nv){
    sumDA[v] <- sum(deltaA[v,])
    for(l in 1:Nlag){
      # weight[4,1] = 1 : current day weight of current SWC is set to 1.
      #     All other day weights of current SWC are set to 0. See below.
      # weight[5,1] = 0 : current day weight of past SWC is set to 0.
      #     Other day day weights of past SWC are estimated. See below.
      deltaA[v, l] <- (v<4)*deltaXA[v, block[l]]/BlockSize[l] + (v==4)*(l==1) + 
        (v==5)*(l>1)*deltaXA[v, block[l]]/BlockSize[l]
      weightA[v, l] <- deltaA[v, l]/sumDA[v]
      # Compute the cumulative weights:
      cum_weightA[v, l] <- sum(weightA[v, 1:l])
    }
  }
  
  # Define the antecedent importance weights for precipitation. 
  #   weightAP[1] = 0 (1 week previously) because SWC effects are accounted for up to 13 days previously
  # WeightsAP are on varying time scale: 
  #   i) weekly scale, starting 2 weeks previous and up to 3 weeks previous (i = 2 and 3), 
  #   ii) monthly scale, for 1 month previously (i = 4), 
  #   iii) bi-monthly scale, starting 2 months previously and up to 5 months previously (i = 5 and 6), 
  #   iv) seasonal scale, starting 6 months previous and over the past year (i = 7 and 8). 
  for(i in 1:NlagP){
    deltaAP[i] <- (i>1)*deltaXAP[d] + (i==1)*0  
    weightAP[i] <- deltaAP[i]/sum(deltaAP[])
    # Compute the cumulative weights:
    cum_weightAP[i] <- sum(weightAP[1:i])
  }
  
  # Standard, relatively non-informative priors for the antecedent climate weights 
  for(v in 1:Nv){
    deltaXA[v, 1:Nblocks] ~ ddirch(rep(1, Nblocks)) # for weights of the five short-term predictors
  }
  deltaXAP[1:NblocksP] ~ ddirch(rep(1, NblocksP)) # for precipitation weights

  # Standard, relatively non-informative priors for the climate effects (an and ag)
  for(i in 1:Ns){
    an[i] ~ dnorm(0, 0.001)
    ag[i] ~ dnorm(0, 0.001)
  }
  
  # Standard, relatively non-informative priors for the NEE standard deviation 
  sig_y ~ dunif(0, 100)
  tau_y <- pow(sig_y, -2)
  
  # Prior for the growing season vs. non-growing season partitioning term 
  phi0 ~ dunif(-1, 1)
  
  # ------ Note: Code below is for creating Fig. S1. And is not necessarily part of the Bayesian model. -----
  # Calculate the overall climate effects (combining growing season and non-growing season effects). 
  for(r in 1:Nmem){
    for(e in 1:Nv){
      # main effects:
      Sen_main[e, r] <- (1 - phi[r]) * an[(e+1)]  + phi[r] * ag[(e+1)]
      # quadratic effects:
      Sen_quad[e, r] <- (1 - phi[r]) * an[(e+6)] * antClim[(e+1), r] + phi[r] * ag[(e+6)] * antClim[(e+1), r]
      # interactions when other covariates take values of the site averages:
      Sen_inter[e, r] <- (e==1)*((1 - phi[r]) * an[12] * antClim[3, r] + phi[r] * ag[12] * antClim[3, r] + 
                                   (1 - phi[r]) * an[13] * antClim[4, r] + phi[r] * ag[13] * antClim[4, r] +
                                   (1 - phi[r]) * an[14] * antClim[5, r] + phi[r] * ag[14] * antClim[5, r] +
                                   (1 - phi[r]) * an[15] * antClim[6, r] + phi[r] * ag[15] * antClim[6, r]) +
        (e==2)*((1 - phi[r]) * an[16] * antClim[4, r] + phi[r] * ag[16] * antClim[4, r] + 
                  (1 - phi[r]) * an[17] * antClim[5, r] + phi[r] * ag[17] * antClim[5, r] +
                  (1 - phi[r]) * an[18] * antClim[6, r] + phi[r] * ag[18] * antClim[6, r] +
                  (1 - phi[r]) * an[12] * antClim[2, r] + phi[r] * ag[12] * antClim[2, r]) +
        (e==3)*((1 - phi[r]) * an[19] * antClim[5, r] + phi[r] * ag[19] * antClim[5, r] + 
                  (1 - phi[r]) * an[20] * antClim[6, r] + phi[r] * ag[20] * antClim[6, r] +
                  (1 - phi[r]) * an[13] * antClim[2, r] + phi[r] * ag[13] * antClim[2, r] +
                  (1 - phi[r]) * an[16] * antClim[3, r] + phi[r] * ag[16] * antClim[3, r]) +
        (e==4)*((1 - phi[r]) * an[21] * antClim[6, r] + phi[r] * ag[21] * antClim[6, r] + 
                  (1 - phi[r]) * an[14] * antClim[2, r] + phi[r] * ag[14] * antClim[2, r] +
                  (1 - phi[r]) * an[17] * antClim[3, r] + phi[r] * ag[17] * antClim[3, r] +
                  (1 - phi[r]) * an[19] * antClim[4, r] + phi[r] * ag[19] * antClim[4, r]) +
        (e==5)*((1 - phi[r]) * an[15] * antClim[2, r] + phi[r] * ag[15] * antClim[2, r] + 
                  (1 - phi[r]) * an[18] * antClim[3, r] + phi[r] * ag[18] * antClim[3, r] +
                  (1 - phi[r]) * an[20] * antClim[4, r] + phi[r] * ag[20] * antClim[4, r] +
                  (1 - phi[r]) * an[21] * antClim[5, r] + phi[r] * ag[21] * antClim[5, r])
      # Summing the main, quadratic, and interaction effects
      Sen[e, r] <- sum(Sen_main[e, r] + Sen_quad[e, r] + Sen_inter[e, r])
    }
    # precipitation effect:
    Sen[6, r] <- (1 - phi[r]) * an[22] + phi[r] * ag[22]
    # the effect of SWC (both current and antecedent):
    Sen[7, r] <- Sen[4, r] + Sen[5, r]
  }
  
  # Obtain the mean over the entire record
  for(e in 1:7){
    # ESen's, therefore, represent the overall effects (combining growing season and non-growing season)
    # of climate covariate e when all other climate covariates take values of the site averages.
    ESen[e] <- mean(Sen[e,])
  }
}
