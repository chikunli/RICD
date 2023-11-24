#Installing required packages and libraries
library(MASS)
library(robustbase)
library(rlist)
#---------------------------------------------------------------------------------------------------------
# Function RICD (min the ridgelized covariance determinant & min the modified mah.distance)
#---------------------------------------------------------------------------------------------------------
RICD <- function(x, alpha=NULL, h=NULL, lambda = NULL, m=NULL, l=NULL, h_ini=NULL){
  x <- data.frame(as.matrix(x, nrow = nrow(x), ncol = ncol(x)))
  p <- ncol(x) #number of varaibles in data
  n <- nrow(x) #number of observations in data
  #------------------------------------------------------------------------------------------
  # Before we start we might check whether some conditions for inputes are satisified
  #------------------------------------------------------------------------------------------
  if(n < 3){
    stop("Too small sample!")
  }
  if(is.null(m)){
    m <- 100
  }
  if(is.null(l)){
    l <- 10
  }
  if(is.null(alpha)){
    alpha <- 0.05
  }
  if(is.null(h)){
    h <- floor(n/2)+1
  }
  if(is.null(h_ini)){
    h_ini <- h
  }
  c_h <- p/h
  delta <- alpha/2
  I_p <- diag(p)
  if(is.null(lambda)){
    lambda <- fnlambda(x)
  }
  determ <- as.vector(rep(0, m)) #initialization of determinant vector
  #------------------------------------------------------------------------------------------
  # Determine top l initial subsets
  #------------------------------------------------------------------------------------------
  #Initialization of all required lists in which we will store the generated samples, matrices and vectors 
  T_0 <- list()
  R_0 <- list()
  H_0 <- list()
  T_1 <- list()
  R_1 <- list()
  H_1 <- list()
  T_2 <- list()
  R_2 <- list()
  H_2 <- list()
  T_3 <- list()
  R_3 <- list()
  H_3 <- list()
  H_final <- list()
  list_sets <- list()
  R_final <- list()
  center <- list()
  
  for(i in 1:m)
    try({
      #-----------------------------------------------------------
      # C0-step looking at H0
      #-----------------------------------------------------------
      H_0[[i]] <- x[sample(1:n, h_ini), ] #randomly select h_ini observations from data
      T_0[[i]] <- apply(H_0[[i]], 2, mean)
      R_0[[i]] <- var(H_0[[i]])+lambda*I_p
      mah_dist0 <- mahalanobis(x, T_0[[i]], R_0[[i]], inverted = FALSE)
      d_sqr0 <- sort(mah_dist0, decreasing = FALSE, index.return = TRUE)
      H1_index <- d_sqr0$ix[1:h]
      H_1[[i]] <- x[H1_index, ] #the first subset with h components that will be used in C1-step
      #-----------------------------------------------------------
      # C1-step looking at H1
      #-----------------------------------------------------------
      T_1[[i]] <- apply(H_1[[i]], 2, mean) #center based on H1
      R_1[[i]] <- var(H_1[[i]])+lambda*I_p
      mah_dist1 <- mahalanobis(x, T_1[[i]], R_1[[i]], inverted = FALSE)
      d_sqr1 <- sort(mah_dist1, decreasing = FALSE, index.return = TRUE)
      H2_index <- d_sqr1$ix[1:h] #indeces of h small d^2's observations
      H_2[[i]] <- x[H2_index, ] #H2 subset
      #-----------------------------------------------------------
      # C2-step looking at H2
      #-----------------------------------------------------------
      T_2[[i]] <- apply(H_2[[i]], 2, mean) #center based on H2
      R_2[[i]] <- var(H_2[[i]])+lambda*I_p
      mah_dist2 <- mahalanobis(x, T_2[[i]], R_2[[i]], inverted = FALSE)
      d_sqr2 <- sort(mah_dist2, decreasing = FALSE, index.return = TRUE)
      H3_index <- d_sqr2$ix[1:h] #indeces of h small d^2's observations
      H_3[[i]] <- x[H3_index, ] #H3 subset
      T_3[[i]] <- apply(H_3[[i]], 2, mean) #center based on H3
      R_3[[i]] <- var(H_3[[i]])+lambda*I_p
      determ[i] <- det(R_3[[i]]) #take the det of the best subset
    }, silent = FALSE)
  
  D <- sort(determ, decreasing = FALSE, index.return = TRUE) #sort the deteminants of all m subsets and pick smallest l determinants
  detl <- D$ix[1:l] #get indexes of these l low det sets
  #-----------------------------------------------------------
  # Get the best subsets H and the corresponding R_RICD
  #-----------------------------------------------------------
  for(r in 1:l){
    list_sets[[r]] <- H_3[[detl[r]]]
    R_final[[r]] <- R_3[[detl[r]]]
  }
  objectRICD <- list(call = match.call(), H0sets = H_0, H1sets = H_1, H2sets = H_2, H3sets = H_3,
                     finall_set = list_sets, finalR = R_final, index_l = detl, determ = determ)
  #---------------------------------------------------------------------------------------
  # While loop for C-steps of convergence for top l subsets from initial step
  #---------------------------------------------------------------------------------------
  for(k in 1:l){
    H_new <- data.frame(as.matrix(objectRICD$finall_set[[k]], nrow = h, ncol = p)) #the first subset of l finial subsets
    term1 <- 1 #determinant of new one
    term2 <- 2 #determinant of old one
    iter <- 1
    while(term1 != term2 && iter < 10){ #continue C-steps until det(R_new) = det(R_old)
      T_old <- apply(H_new, 2, mean)
      R_old <- var(H_new)+lambda*I_p
      term2 <- det(R_old) #determinant of the initial subset
      mah_distnew <- mahalanobis(x, T_old, R_old, inverted = FALSE)
      d_loop <- sort(mah_distnew, decreasing = FALSE, index.return = TRUE)
      index_new <- d_loop$ix[1:h]
      #Update infromation based on new indices
      H_old <- H_new #new sumbset becomes old subset to compare it with new one
      H_new <- x[index_new, ] #update H_new using the indices determined by H_old
      T_new <- apply(H_new, 2, mean)
      R_new <- var(H_new)+lambda*I_p
      term1 <- det(R_new) #update new determinant
      iter <- iter + 1
    }
    center[[k]] <- T_new # store the final center estimate in list of this final l centers
    R_final[[k]] <- R_new
    H_final[[k]] <- H_new  #store the final l subsets
  }  
  #Decide which one has the min determinant of these final l sets
  d <- rep(0, l) 
  for(i in 1:l){
    d[i] <- det(R_final[[i]])
  }
  detf <- sort(d, decreasing = FALSE, index.return = TRUE) #detf is sorted vector of l final determinants
  mindet <- detf$ix[1] #the index of the best set corresponding to the smallest det
  #------------------------------------------------------------------------------------
  # Final estimators with index mindet 
  #------------------------------------------------------------------------------------
  T_RICD <- center[[mindet]] #center estimater of the minimal det set
  R_RICD <- R_final[[mindet]]
  H_RICD <- H_final[[mindet]] #final subset with smallest det
  indeces_H_RICD <- as.integer(rownames(H_RICD)) #indeces of the final subset
  inverse_R_RICD <- solve(R_RICD)
  raw.mfnp <- sum(diag(inverse_R_RICD))/p
  raw.mfnp_prime <- sum(diag(inverse_R_RICD%*%inverse_R_RICD))/p
  raw.Theta1 <- (1-lambda*raw.mfnp)/(1-c_h*(1-lambda*raw.mfnp))
  raw.Theta2 <- (1-lambda*raw.mfnp)/((1-c_h+c_h*lambda*raw.mfnp)^3)-lambda*(raw.mfnp-lambda*raw.mfnp_prime)/((1-c_h+c_h*lambda*raw.mfnp)^4)
  raw.cutoff <- (qnorm((1-delta))*sqrt(2*raw.Theta2/p)+raw.Theta1)*p
  result_raw <- list(raw.center=T_RICD, raw.R=R_RICD, raw.inv.R=inverse_R_RICD, result.best=indeces_H_RICD, lambda=lambda, H_RICD=H_RICD, determ=determ, raw.Theta1=raw.Theta1, raw.Theta2=raw.Theta2, raw.cutoff=raw.cutoff)
  #------------------------------------------------------------------------------------
  # Reweighting step
  #------------------------------------------------------------------------------------
  # We will use delta to expect that (1-delta)*100% of data will be included
  w <- rep(0, n) #initialization of weight vector
  d_raw <- mahalanobis(x, result_raw$raw.center, result_raw$raw.inv.R, inverted = TRUE)#mah.distances for all observations
  x <- as.matrix(x, nrow=n, ncol=p)
  for(i in 1:n){
    if(d_raw[i] <= result_raw$raw.cutoff){
      w[i] <- 1
    }
  }
  x_rw <- x[w==1, ]
  #--------------------------------------
  # Refined RICD estimator
  #--------------------------------------
  c_rw <- p/sum(w)
  T_rw <- apply(x_rw, 2, mean)
  S_tilde <- var(x_rw)
  delta_w <- 1-sum(w)/n
  k_ricd <- 1+2*dnorm(qnorm(1-delta_w))/sqrt(2*p*result_raw$raw.Theta2)*result_raw$raw.Theta1/(1-delta_w)
  R_refine <- k_ricd*S_tilde+lambda*I_p
  inv_R_refine <- solve(R_refine)
  mfnp <- sum(diag(inv_R_refine))/p
  mfnp_prime <- sum(diag(inv_R_refine%*%inv_R_refine))/p
  Theta1 <- (1-lambda*mfnp)/(1-c_rw*(1-lambda*mfnp))
  Theta2 <- (1-lambda*mfnp)/((1-c_rw+c_rw*lambda*mfnp)^3)-lambda*(mfnp-lambda*mfnp_prime)/((1-c_rw+c_rw*lambda*mfnp)^4)
  cutoff <- (qnorm((1-alpha))*sqrt(2*Theta2/p)+Theta1)*p
  #--------------------------------------
  # Identify outliers
  #--------------------------------------
  dist_RICD <- mahalanobis(x, T_rw, inv_R_refine, inverted = TRUE)
  ou_RICD <- rep(0, n)
  for(dd in 1:n){
    if(dist_RICD[dd] > cutoff ){
      ou_RICD[dd] <- 1
    }
  }
  #The output adding the new results to result_raw with raw estimators
  result_RICD <- list.append(result_raw, k_ricd=k_ricd, w=w, x_rw=x_rw, center=T_rw, R_refine = R_refine, Theta1 = Theta1, Theta2 = Theta2, cutoff = cutoff, dist_RICD = dist_RICD, ou_RICD=ou_RICD)
  return(result_RICD)
} ####end function RICD 
#-----------------------------------------------------------
# Calculate lambda
#-----------------------------------------------------------
fnlambda <- function(X)
{
  X <- as.matrix(X)
  n <- dim(X)[1]
  p <- dim(X)[2]
  alpha <- 0.05
  cov_X <- var(X)
  c_test <- p/n
  grid <- seq(0.05, 200, by=2)
  len_grid <- length(grid)
  mu_hat <- colMeans(X)
  d2 <- 2
  Theta1_test <- Theta2_test <- 0
  i <- 0
  while (abs(median(d2)-p*Theta1_test-qnorm(1-alpha)*sqrt(2*p*Theta2_test))>1 && i < len_grid){
    i <- i + 1
    lambda <- grid[i]
    R_test <- cov_X+lambda*diag(p)
    inverse_R_test <- solve(R_test)
    d2 <- mahalanobis(X, mu_hat, inverse_R_test, inverted = TRUE)
    mfnp_test <- (sum(diag(inverse_R_test)))/p
    mfnp_prime_test <- sum(diag(inverse_R_test%*%inverse_R_test))/p
    Theta1_test <- (1-lambda*mfnp_test)/(1-c_test*(1-lambda*mfnp_test))
    Theta2_test <- (1-lambda*mfnp_test)/((1-c_test+c_test*lambda*mfnp_test)^3)-lambda*(mfnp_test-lambda*mfnp_prime_test)/((1-c_test+c_test*lambda*mfnp_test)^4)
  }
  return(lambda)
} ####end function fnlambda 
