### LOAD THE LIBRARIES ###
library(R2OpenBUGS)
library(coda)
library(mcmcplots)

### LOAD THE DATASET ###
dtf_input <- read.table("day.csv",
                        head=T,
                        sep=",")

dtf_input$instant <- NULL
dtf_input$dteday <- NULL
dtf_input$casual <- NULL
dtf_input$registered <- NULL
# Number of observations
int_n <- nrow(dtf_input)
# Number of variables

int_p = ncol(dtf_input)
dtf_X_t <- matrix(unlist(dtf_input[, -int_p]),
                nrow=dim(dtf_input[, -int_p])[1],
                byrow=FALSE)

dtf_y_t = unlist(dtf_input$cnt)

pos_scaler <- function(lst_el) {
  return((lst_el - mean(lst_el))/sd(lst_el))
}

dtf_X <- matrix(nrow = dim(dtf_X_t)[1],
                ncol = dim(dtf_X_t)[2])

for (int_num in (1:dim(dtf_X_t)[2])) {
  dtf_X[, int_num] <- pos_scaler(lst_el = dtf_X_t[, int_num])
}

dtf_y <- pos_scaler(lst_el = dtf_y_t)

### GIBBS SAMPLING ###
# Function that computes log of marginal distribution
cls_log_marg <- function(dtf_y,
                         dtf_X,
                         int_g=length(dtf_y),
                         flt_nu0=1,
                         flt_s20=try(summary(lm(dtf_y~dtf_X-1))$sigma^2,
                                     silent=TRUE)) {
  
  int_n <- dim(dtf_X)[1]
  int_p <- dim(dtf_X)[2]
  # Compute SSR_z^g = y^T(I-H_{g,z})y
  if (int_p == 0) {
    lst_Hg <- 0;
    flt_s20 <- mean(dtf_y^2)
  }
  
  if (int_p > 0) {
    lst_Hg <- (int_g/(int_g+1))*dtf_X%*%solve(t(dtf_X)%*%dtf_X)%*%t(dtf_X)
  }
  
  flt_SSRg <- t(dtf_y)%*%(diag(x=1 ,
                               nrow=int_n)-lst_Hg)%*%dtf_y
  
  # Return the log of the distribution of dtf_y conditional X and z (pp 165)
  return(-0.5*(int_n*log(pi)+int_p*log(1+int_g)+(flt_nu0+int_n)*log(flt_nu0*flt_s20+flt_SSRg)-flt_nu0*log(flt_nu0*flt_s20))+lgamma((flt_nu0+int_n)/2)-lgamma (flt_nu0/2))
}


# GIBBS SAMPLER
# Function that develops the Gibbs sampler
# I sample z (the vector of 0/1's for feature importance)
# At the same time, I sample Beta and sigma2 for the regression
# PRIORS: parameters are g, nu0, s20
# sigma2: 1/gamma(nu0, sigma02)
# beta: g-prior
cls_gibbs <- function(int_G,
                      int_burnin,
                      int_thin,
                      dtf_X,
                      dtf_y,
                      lst_z0,
                      flt_s20=try(summary(lm(dtf_y~dtf_X-1))$sigma^2,
                                  silent=TRUE),
                      
                      int_g=length(dtf_y),
                      flt_nu0=1) {
  
  # We start computing the number of iteration of the algorithm
  int_iterations <- int_burnin + int_thin * int_G
  cat("Number of iterations:", int_iterations)
  int_acc <- 1
  
  # Define the output vector for z
  lst_Z <- matrix(data=NA,
                  nrow=int_G,
                  ncol=dim(dtf_X)[2])
  
  # Define the output vector for s2
  lst_s2 <- matrix(data=NA,
                   nrow=int_G,
                   ncol=1)
  
  # Define the output vector for beta
  lst_beta <- matrix(data=NA,
                     nrow=int_G,
                     ncol=dim(dtf_X)[2])
  
  # The current state of the chains
  lst_z_curr <- lst_z0
  lst_s2_curr <- flt_s20
  
  # Sample z
  log_distr_curr <- cls_log_marg(dtf_y = dtf_y,
                                 dtf_X = dtf_X[, lst_z_curr==1, drop=FALSE])
  
  for (int_i in 1:int_iterations) {
    if (int_i %% 50 == 0) {print(int_i)}
    for (int_j in sample(1:dim(dtf_X)[2])) {
      lst_zp <- lst_z_curr;
      lst_zp[int_j] <- 1 - lst_zp[int_j]
      log_distr_next <- cls_log_marg(dtf_y = dtf_y,
                                     dtf_X = dtf_X[, lst_zp==1, drop=FALSE])
      
      flt_log_odds <- (log_distr_next - log_distr_curr)*(-1)^(lst_zp[int_j]==0)
      lst_z_curr[int_j] <- rbinom(n = 1,
                                  size = 1,
                                  prob = 1/(1+exp(-flt_log_odds)))
      
      if (lst_z_curr[int_j] == lst_zp[int_j]) {
        log_distr_curr <- log_distr_next
        }
    }
    # NOW I WILL SIMULATE SIGMA, THEN BETA
    dtf_Xz <- dtf_X[, lst_z_curr==1, drop=FALSE]
    lst_hgz <- (int_g/(int_g+1))*dtf_Xz%*%solve(t(dtf_Xz)%*%dtf_Xz)%*%t(dtf_Xz)
    flt_SSRgz <- t(dtf_y)%*%(diag(x=1,
                                  nrow=dim(dtf_Xz)[1])-lst_hgz)%*%dtf_y
    
    flt_s2_curr <- 1/rgamma(n=1,
                            shape=(flt_nu0+int_n)/2,
                            rate=(flt_nu0*flt_s20 + flt_SSRgz)/2)
    
    lst_var_beta <- int_g*solve(t(dtf_Xz)%*%dtf_Xz)/(int_g+1)
    lst_mean_beta <- lst_var_beta%*%t(dtf_Xz)%*%dtf_y
    lst_mean <- matrix(data=rnorm(n=dim(dtf_Xz)[2],
                                  mean=0,
                                  sd=sqrt(flt_s2_curr)),
                       
                       nrow=1,
                       ncol=dim(dtf_Xz)[2])
    
    # I've compute s2 and beta using X_z, only the relevant features according to z
    # The simulation will be this one and zero when z=0.
    lst_beta_z <- t(t(lst_mean%*%chol(x=lst_var_beta))+c(lst_mean_beta))
    lst_beta_curr <- rep(0, dim(dtf_X)[2])
    lst_beta_curr[which(lst_z_curr==1)] <- lst_beta_z
    ## In the output matrix we save only the states
    ## visited after the burn-in period and each "thin" iterations
    if((int_i > int_burnin) & (int_i %% int_thin == 0)){
      lst_Z[int_acc, ] <- lst_z_curr
      lst_s2[int_acc, ] <- flt_s2_curr
      lst_beta[int_acc, ] <- lst_beta_curr
      int_acc = int_acc + 1
    }
  }
  return(list("z" = lst_Z,
              "beta" = lst_beta,
              "sigma2" = lst_s2))
}

cls_diagnostic <- function(lst_beta_chain, int_dim) {
  layout(mat=matrix(data=c(1,2),
                    nrow=2,
                    ncol=1,
                    byrow=T))
  
  plot(x=lst_beta_chain,
       type="l",
       xlab = "SIM",
       ylab = "VALUE",
       main = paste("TRACEPLOT of beta[", int_dim, "]",
                    sep=""))
  
  plot(acf(lst_beta_chain,
           plot=FALSE),
       
       xlab = "LAGS",
       ylab = "ACF",
       main = paste("ACF of beta[", int_dim, "]",
                    sep=""))
}

int_G <- 5000
int_burnin <- 1000
int_thin <- 5

# Let's start with a z vector all equal to 1
lst_z0 <- rep(1, dim(dtf_X)[2])
set.seed(567)

lst_post <- cls_gibbs(int_G = int_G,
                      int_burnin = int_burnin,
                      int_thin = int_thin,
                      dtf_X = dtf_X,
                      dtf_y = dtf_y,
                      lst_z0 = lst_z0)

lst_z_post <- lst_post$z
lst_beta_post <- lst_post$beta
lst_sigma2_post <- lst_post$sigma2

dim(lst_z_post)
head(lst_z_post)

### CHECK FOR CONVERGENCE AND MIXING ###
# A. FOR sigma2
# Traceplot and ACF
# But the main results are:
# B. FOR BETA
# I choose some dimensions of beta and look at
# 1. Trace plot
# 2. ACF plot
# 3. The histogram of the (simulated) posterior density, along with
#   - The estimated density (via built-in kernel density estimation function of R)
#   - The 5%, 50% and 95% quantiles of the (simulated) posterior distribution
# Finally: geweke test for convergence of beta's

set.seed(345)
lst_dim_checks <- sample(x = (1:dim(lst_z_post)[2]),
                         size = 3,
                         replace = FALSE)

print(lst_dim_checks)
################## DIAGNOSTIC FOR SIGMA2 ##################
layout(mat=matrix(data=c(1,2),
                  nrow=2,
                  ncol=1,
                  byrow=T))

plot(x=lst_sigma2_post[, 1],
     type="l",
     xlab = "SIM",
     ylab = "VALUE",
     main = "TRACEPLOT of sigma2")

plot(acf(lst_sigma2_post[, 1],
         plot=FALSE),
     
     xlab = "LAGS",
     ylab = "ACF",
     main = "ACF of sigma2")

geweke.diag(mcmc(lst_sigma2_post[, 1]))

################## DIAGNOSTIC FOR BETA ##################
# 1. IN ONE PLOT: trace plot and ACF
cls_diagnostic(lst_beta_chain = lst_beta_post[, lst_dim_checks[1]],
               int_dim = lst_dim_checks[1])

cls_diagnostic(lst_beta_chain = lst_beta_post[, lst_dim_checks[2]],
               int_dim = lst_dim_checks[2])

cls_diagnostic(lst_beta_chain = lst_beta_post[, lst_dim_checks[3]],
               int_dim = lst_dim_checks[3])

# 2. GEWEKE DIAGNOSTIC
geweke.diag(mcmc(lst_beta_post[, lst_dim_checks[1]]))
geweke.diag(mcmc(lst_beta_post[, lst_dim_checks[2]]))
geweke.diag(mcmc(lst_beta_post[, lst_dim_checks[3]]))

################## GELMAN AND RUBIN DIAGNOSTIC ##################
lst_post_gr_1 <- cls_gibbs(int_G = int_G,
                           int_burnin = int_burnin,
                           int_thin = int_thin,
                           dtf_X = dtf_X,
                           dtf_y = dtf_y,
                           lst_z0 = lst_z0,
                           flt_s20=try(10*summary(lm(dtf_y~dtf_X-1))$sigma^2,
                                       silent=TRUE),
                           
                           int_g=10*length(dtf_y),
                           flt_nu0=1*2)

lst_post_gr_2 <- cls_gibbs(int_G = int_G,
                           int_burnin = int_burnin,
                           int_thin = int_thin,
                           dtf_X = dtf_X,
                           dtf_y = dtf_y,
                           lst_z0 = lst_z0,
                           flt_s20=try(0.1*summary(lm(dtf_y~dtf_X-1))$sigma^2,
                                       silent=TRUE),
                           
                           int_g=0.1*length(dtf_y),
                           flt_nu0=1/2)

lst_beta_post_1 <- lst_post_gr_1$beta
lst_beta_post_2 <- lst_post_gr_2$beta
lst_sigma2_post_1 <- lst_post_gr_1$sigma2
lst_sigma2_post_2 <- lst_post_gr_2$sigma2

lst_gr_3 <- mcmc.list(mcmc(lst_beta_post[, lst_dim_checks]),
                      mcmc(lst_beta_post_1[, lst_dim_checks]),
                      mcmc(lst_beta_post_2[, lst_dim_checks]))

lst_gr_sigma2 <- mcmc.list(mcmc(lst_sigma2_post[, 1]),
                           mcmc(lst_sigma2_post_1[, 1]),
                           mcmc(lst_sigma2_post_2[, 1]))

gelman.diag(lst_gr_3)
gelman.diag(mcmc(lst_gr_sigma2))

#4. EFFECTIVE SAMPLE SIZE (mixing)
effectiveSize(mcmc(lst_beta_post[, lst_dim_checks])) # single chain
effectiveSize(lst_gr_3) # multiple chain

effectiveSize(mcmc(lst_sigma2_post[, 1])) # single chain
effectiveSize(lst_gr_sigma2) # multiple chain
 #########################################################################

# Evaluation
summary(object=lst_beta_post)
summary(object=lst_sigma2_post)

pdf("output_2.pdf",family="Times",height=1.75,width=5)
par(mar=c(2.75,2.75,.5,.5),
    mgp=c(1.7,.7,0))

plot(apply(X = lst_z_post,
           MARGIN = 2,
           FUN = mean,
           na.rm = TRUE),
     
     xlab="Regressor index",
     ylab=expression(paste( "Pr(",italic(z[j] == 1),"|", italic(y),", X)",
                            sep="")),
     
     type="h",
     lwd=2) ; abline(h=0.5,
                     col="red",
                     lty=2)

dev.off()


