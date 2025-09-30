#############################################
##
## dat_gen: function to generate data
##
#############################################

### scenario = 1, high-dimensional setting
### scenario = 2, observational study with 10 covariates for HAL implenetation (Ertefaie et al. 2023)
dat_gen <- function(n, ps, seed1, seed2){
  if(scenario == 1){### high-dimensional data
    nstudy<- n
    nvars<- 1000
    
    ## defining alpha and beta coefficients
    nc<- 100
    ni<- 0
    nr<- 0
    ns<- nvars-(nc+ni+nr)
    
    ## global parameters for sim (only want to set once so they are same for all sims)
    ## seed2 should remain same across simulation runs
    set.seed(seed2)
    coef_strength<- 0.693
    alpha_conf<- runif(nc, 0.0, coef_strength)
    beta_conf<- runif(nc, 0.0, coef_strength)
    alpha_temp<- c(alpha_conf)
    beta_temp<- c( beta_conf)
    random_neg<- sample(1:length(alpha_temp), 0.5*length(alpha_temp), replace=FALSE)
    alpha_temp[random_neg]<- -1*alpha_temp[random_neg]
    beta_temp[random_neg]<- -1*beta_temp[random_neg]
    alpha<- matrix(c(alpha_temp, rep(0, ns)), ncol=1)
    beta<- matrix(c(beta_temp, rep(0, ns)), ncol=1)
    betaE<- 0.0
    cprev<- runif(nvars, 0.01, 0.1)
    oprev<- 0.05
    tprev<- 0.30
    
    ## resetting seed so that it is unique for each simulation
    ## seed1 should change each run to get new data
    set.seed(seed1)
    
    # generate synthetic matrix of baseline covariates
    Xcovs<- matrix(NA, nrow=nstudy, ncol=nvars)
    for(pp in 1:nvars){
      Xcovs[,pp]<- rbinom(nstudy, 1, cprev[pp])
    }
    Xcovs<- as.data.frame(Xcovs)
    names(Xcovs)<- c(paste0("x", 1:nvars))
    W<- as.matrix(Xcovs)
    linear_pred_e<- W %*% alpha
    linear_pred_y<- W %*% beta
    
    ##function to find intercept to get specified treatment prevalence
    treatment_inc<- tprev
    fn <- function(c) mean(plogis(c + linear_pred_e)) - treatment_inc
    alpha0 <- uniroot(fn, lower = -20, upper = 20)$root
    Ee <- (1 + exp( -(alpha0 + linear_pred_e) ))^-1
    e<- rbinom(nstudy, 1, Ee)
    A<- e
    p<- Ee
    
    ##function to find intercept to get specified outcome incidence
    outcome_inc<- oprev
    fn <- function(c) mean(plogis(c + betaE*e + linear_pred_y )) - outcome_inc
    beta0 <- uniroot(fn, lower = -20, upper = 20)$root
    Ey <- (1 + exp( -( beta0 + betaE*e + linear_pred_y )))^-1
    Y <- rbinom(nstudy, 1, Ey)
    Y1<- (1 + exp( -( beta0 + betaE*1 + linear_pred_y )))^-1
    Y0<- (1 + exp( -( beta0 + betaE*0 + linear_pred_y )))^-1
  }
  if(scenario == 2){### Simulation in Supplemental Material of Ertefaie et al. (2023)
    set.seed(seed1)
    x1 = runif(n,-2,2)
    x2 = runif(n,-2,2)
    x3 = runif(n,-2,2)
    x4 = runif(n,-2,2)
    x5 = runif(n,-2,2)
    x6 = rbinom(n,1,0.6)
    x7 = rbinom(n,1,0.6)
    x8 = rbinom(n,1,0.6)
    x9 = rbinom(n,1,0.6)
    x10 = rbinom(n,1,0.6)
    
    ##function to find intercept to get specified treatment prevalence
    linear_pred_e = (1*x2^2-exp(x1/2)-x3+x4-exp(x5/2)+x6+x7)/2
    tprev<- 0.50
    treatment_inc<- tprev
    fn <- function(c) mean(plogis(c + linear_pred_e)) - treatment_inc
    alpha0 <- 0 #uniroot(fn, lower = -20, upper = 20)$root
    p <- (1 + exp( -(alpha0 + linear_pred_e) ))^-1
    A = 1-rbinom(n, size = 1, prob = p)
    Y_error<- rnorm(n, mean = 0, sd = 0.1)
    Y1 = +0*1 +(-2*x2^2+2*x1+2*mean(x2^2)+x2+x1*x2+x3+x4+2*x5^2-2*mean(x5^2))/1 + Y_error
    Y0 = +0*0 +(-2*x2^2+2*x1+2*mean(x2^2)+x2+x1*x2+x3+x4+2*x5^2-2*mean(x5^2))/1 + Y_error
    Y = A*Y1 + (1-A)*Y0
  }
  if(ps == 1) {return(data.frame(Y, Y1, Y0, A, p, Xcovs))}
  if(ps == 2) {return(data.frame(Y, Y1, Y0, A, p, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10))}
}
