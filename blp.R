BerryLevinsohnPakes <- function(dat, mkt.id.fld = "mkt.id", prod.id.fld = "prod.id", prc.fld = "px", share.fld = "share", x.var.flds = c("x1", "x2", "x3"), prc.iv.flds = "z", n.sim = 1500, sigma.guess){
  # all data should appear in the data.frame "dat"
  # input variables ending in ".fld" are the names of the columns
  # n.sim = number of simulated "indviduals" per market 
  # sigma guess may be missing
  #
  # Required packages (installed as the function runs):
  # SQUAREM, BB, AER
  #
  # Based on code written by Aviv Nevo, May 1998.
  # Adapted by Michael Carniol, January 2015
  
  blp_inner <- function(delta.in, mu.in) {
    # Computes a single update of the BLP (1995) contraction mapping.
    # of market level predicted shares.
    # This single-update function is required by SQUAREM, see Varadhan and
    # Roland (SJS, 2008), and Roland and Varadhan (ANM, 2005)
    # INPUT
    # 	delta.in : current value of delta vector
    # 	mu.in: current mu matrix
    # Requires global variables: s.jt
    # OUTPUT
    # 	delta.out : delta vector that equates observed with predicted market shares
    pred.s <- rowMeans(ind_sh(delta.in, mu.in));
    delta.out <- delta.in + log(s.jt) - log(pred.s)
    return(delta.out)
  }
  ind_sh <- function(delta.in, mu.in){
    # This function computes the "individual" probabilities of choosing each brand
    # Requires global variables: mkt.id, X, v
    numer <- exp(mu.in) * matrix(rep(exp(delta.in), n.sim), ncol = n.sim);
    denom <- as.matrix(do.call("rbind", lapply(mkt.id, function(tt){
      1 + colSums(numer[mkt.id %in% tt, ])
    })))
    return(numer / denom);	
  }
  gmm_obj <- function(theta2){
    # This function computes the GMM objective function
    # Requires global variable inputs: X, v, delta, a, W
    # Outputs: theta1, xi.hat
    print(paste0("GMM Loop number: ", Sys.time()))
    print(a <<- a + 1)
    print("Updated theta2 estimate:")
    print(theta2)
    print("Change in theta2 estimate:")
    print(theta.chg <- as.numeric(theta2 - theta2.prev));
    if(sum(theta.chg != 0) <= 2){
      delta <- dat[, "delta"];
    } else {
      delta <- Y;
    }
    theta2.prev <<- theta2;
    
    mu <- X %*% diag(theta2) %*% v;
    
    print("Running SQUAREM contraction mapping")
    print(system.time(
      squarem.output <- squarem(par = delta, fixptfn = blp_inner, mu.in = mu, control = list(trace = TRUE))
    ));
    delta <- squarem.output$par
    print(summary(dat[, "delta"] - delta));
    dat[, "delta"] <<- delta;
    
    mo.ivreg <- ivreg(fm.ivreg, data = dat, x = TRUE)
    theta1 <<- coef(mo.ivreg);
    xi.hat <<- as.vector(mo.ivreg$resid);
    Z.hat <- Z * matrix(rep(xi.hat, ncol(Z)), ncol = ncol(Z))
    W.inv <- try(solve(t(Z.hat) %*% Z.hat), silent = FALSE)
    if("matrix" == class(W.inv)){
      PZ <<- Z %*% W.inv %*% t(Z);
      PX.inv <- solve(t(X) %*% PZ %*% X)
      theta1 <<- PX.inv %*% t(X) %*% PZ %*% delta
      xi.hat <<- delta - X %*% theta1
      X.hat <- (PZ %*% X) * matrix(rep(xi.hat, K), ncol = K)
      tsls.se <- sqrt(diag(PX.inv %*% t(X.hat) %*% X.hat %*% PX.inv))
      print("GMM step 2 updated theta1 estimate:")
      print(beta.est <<- data.frame(beta.est = theta1, beta.se = tsls.se, sigma.est = theta2))
    }
    dat[, "xi.hat"] <<- xi.hat
    f <- t(xi.hat) %*% PZ %*% xi.hat;
    print("Updated GMM objective:")
    print(f <- as.numeric(f));
    return(f)
  }
  jacobian <- function(delta.in, theta.in){
    print(paste0("Calculating Jacobian matrix, ", Sys.time()))
    #Requires global variables X, v, mkt.id
    mu1 <- X %*% diag(theta.in) %*% v;
    ind.shares <- ind_sh(delta.in, mu1);
    K <- ncol(X);
    print(paste0("Calculating dsigma matrix, ", Sys.time()))
    dsigma <- lapply(l.Xv, function(x){
      temp2 <- x * ind.shares;
      temp3 <- as.matrix(do.call("rbind", lapply(mkt.id, function(m){
        colSums(temp2[mkt.id %in% m, ])
      })));
      dsigma.res <- rowMeans(temp2 - ind.shares * temp3);
      return(dsigma.res)
    })
    dsigma <- as.matrix(do.call("cbind", dsigma))
    print(paste0("Calculating ddelta matrices, ", Sys.time()))
    ddelta <- list()
    for(m in mkt.id){
      if(m %in% names(ddelta)){next}
      temp1 <- as.matrix(ind.shares[mkt.id %in% m, ]);
      H1 <- temp1 %*% t(temp1);
      H2 <- diag(rowSums(temp1));
      H <- (H2 - H1) / n.sim;
      H.inv <- solve(H);
      ddelta[[as.character(m)]] <- H.inv %*% dsigma[mkt.id %in% m, ];
      rm(temp1, H1, H2, H, H.inv)
    }
    ddelta <- as.matrix(do.call("rbind", ddelta));
    return(ddelta)
  }
  gradient_obj <- function(theta2){
    #Requires global variables PZ, delta, xi.hat
    print(system.time(jacobian_res <<- jacobian(as.vector(dat[, "delta"]), theta2)))
    print(paste0("Updated gradient:", Sys.time()))
    print(f <- -2 * as.numeric(t(jacobian_res) %*% PZ %*% xi.hat));
    #######
    L <- ncol(Z)
    covg <- matrix(0, nrow = L, ncol = L)
    for(i in 1:JT){
      covg <- covg + (Z[i, ] %*% t(Z[i, ])) * xi.hat[i]^2
    }
    d.delta <- jacobian_res;
    Dg <- t(d.delta) %*% Z
    p.Dg <- try(solve(Dg %*% W.inv %*% t(Dg)))
    cov.mat <- p.Dg %*% (Dg %*% W.inv %*% covg %*% W.inv %*% t(Dg)) %*% p.Dg
    beta.est$sigma.se <<- sqrt(diag(cov.mat));
    print(paste0("Updated coefficients table:", Sys.time()))
    print(beta.est)
    write.csv(beta.est, file = paste0("BLP_beta_est_", Sys.Date(), ".csv"))
    #######
    return(as.numeric(f))
  }
  
  #Set up data
  dat <- dat[dat[, share.fld] > 0, ]
  dat <- dat[order(dat[, mkt.id.fld], dat[, prod.id.fld]), ]
  JT <- nrow(dat)
  #market identifier variable
  mkt.id <- dat[, mkt.id.fld];
  #Number of characteristics (including constant and price)
  X <- as.matrix(cbind(ones = rep(1, JT), dat[, c(x.var.flds, prc.fld)]));
  K <- ncol(X)
  #Compute the outside good market share by market
  s.jt <- as.vector(dat[, share.fld]);
  temp <- aggregate(s.jt, by = list(mkt.id = mkt.id), sum);
  sum1 <- temp$x[match(mkt.id, temp$mkt.id)];
  s.j0 <- as.vector(1 - sum1);
  rm(temp, sum1);
  dat[, "delta"] <- Y <- log(s.jt) - log(s.j0);
  iv <- dat[, prc.iv.flds]
  
  while(!require(AER)){install.packages("AER")}
  #Construct 2SLS regression specification
  str.ivreg.y <- "delta ~ "
  str.ivreg.x <- paste(x.var.flds, collapse = " + ")
  str.ivreg.prc <- paste(prc.fld, collapse = " + ")
  str.ivreg.iv <- paste(prc.iv.flds, collapse = " + ")
  print("2SLS specification:")
  print(fm.ivreg <- paste0(str.ivreg.y, str.ivreg.x, " + ", str.ivreg.prc, " | ", str.ivreg.x, " + ", str.ivreg.iv))
  rm(str.ivreg.y, str.ivreg.x, str.ivreg.prc, str.ivreg.iv)
  
  print("2SLS beta estimate:")
  print(summary(mo.ivreg <- ivreg(fm.ivreg, data = dat, x = TRUE)))
  beta.est <- summary(mo.ivreg)$coef[, 1:2]
  #Z = instrumental variable matrix include exogenous X's
  Z <- as.matrix(mo.ivreg$x$instruments)
  PZ <- Z %*% solve(t(Z) %*% Z) %*% t(Z);
  theta1 <- coef(mo.ivreg);
  xi.hat <- as.vector(mo.ivreg$resid);
  Z.hat <- Z * matrix(rep(xi.hat, ncol(Z)), ncol = ncol(Z))
  W.inv <- try(solve(t(Z.hat) %*% Z.hat), silent = FALSE)
  if("matrix" == class(W.inv)){
    PZ <- Z %*% W.inv %*% t(Z);
    PX.inv <- solve(t(X) %*% PZ %*% X)
    theta1 <- PX.inv %*% t(X) %*% PZ %*% Y
    xi.hat <- Y - X %*% theta1
    X.hat <- (PZ %*% X) * matrix(rep(xi.hat, K), ncol = K)
    tsls.se <- sqrt(diag(PX.inv %*% t(X.hat) %*% X.hat %*% PX.inv))
    print("GMM step 2 updated theta1 estimate:")
    print(beta.est <- data.frame(beta.est = theta1, se.est = tsls.se))
  }
  dat[, "xi.hat"] <- xi.hat
  
  #Starting point
  print("Sigma guess:")
  if(missing(sigma.guess)){
    tsls.se <- beta.est[, 2]
    print(theta2 <- 0.5 * tsls.se);
  } else {
    print(theta2 <- sigma.guess);
  }
  theta2.prev <- theta2;
  
  # Matrix of individuals' characteristics
  #	Standard normal distribution draws, one for each characteristic
  v <- matrix(rnorm(K * n.sim), nrow = K, ncol = n.sim)
  # Break X and v matrices into list variables 
  # in attempt to expedite calculation of the Jacobian matrix
  l.X <- lapply(1:K, function(k){
    return(X[, k])
  })
  l.v <- lapply(1:K, function(k){
    return(v[k, ])
  })
  l.Xv <- lapply(1:K, function(k){
    l.X[[k]] %*% t(l.v[[k]]);
  })
  
  print("Estimating random coefficients multinomial logit")
  a <- 0;
  beta.est <- NULL;
  while(!require(SQUAREM)){install.packages("SQUAREM")}
  while(!require(BB)){install.packages("BB")}
  print(system.time(
    theta.est <- multiStart(par = theta2, fn = gmm_obj, gr = gradient_obj, lower = 0, control = list(trace = TRUE), action = "optimize")
  ));
  save(theta.est, file = paste0("theta_est_", Sys.time(), ".RData"))
  
  print("Final coefficients estimate:")
  theta2 <- theta.est$par
  gmm.res <- gmm_obj(theta2)
  grad.res <- gradient_obj(theta2)	
  return(list(coef.mat = beta.est, gmm.obj.func = gmm.res, gmm_est = theta.est, final.data = dat))
}