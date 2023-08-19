
# estimate mortality adjusted by hiv prevalence and art proportion.
# hiv svd.comp model. Coded by Sarah Hertog based on model from Sam Clark (https://github.com/sinafala/svd-comp)

predictNQX <- function(sex, cm, am=NULL, hiv, art, adult, im=NULL) {
  # sex: "female" or "male"
  # cm: vector of 5q0 values
  # am: vector of 45q15 or 35q15 values
  # hiv: vector of HIV values (%)
  # art: vector of ART values (%)
  # adult: "q45" or "q35"
  # im: vector of 1q0 values
  
  library(rms)
  
  # Logit transform
  cml <- logit(cm)
  
  # predict 1q0
  if(missing(im)) {
    cmls <- cml^2
    preds.q0 <- data.frame(
      cml = as.numeric(cml),
      cmls = as.numeric(cmls)
    )
    iml <- predict(mods.r[[sex]]$q0, newdata=preds.q0)
  } else{
    iml <- logit(im)
  }
  
  # predict adult mx
  if(missing(am)) {
    preds.aml <- data.frame(
      cm = as.numeric(cm),
      cml = as.numeric(cml),
      hiv = as.numeric(hiv),
      art = as.numeric(art)
    )
    aml <- predict(mods.r[[sex]][[adult]]$aml, newdata=preds.aml)	
  } else {
    aml <- logit(am)
  }
  
  preds.vs <- data.frame(
    cml = as.numeric(cml),
    aml = as.numeric(aml),
    hiv = as.numeric(hiv),
    art = as.numeric(art)
  )
  
  v1 <- predict(mods.r[[sex]][[adult]]$v1, newdata=preds.vs)
  v2 <- predict(mods.r[[sex]][[adult]]$v2, newdata=preds.vs)
  v3 <- predict(mods.r[[sex]][[adult]]$v3, newdata=preds.vs)
  v4 <- predict(mods.r[[sex]][[adult]]$v4, newdata=preds.vs)
  
  # adjust weights
  input.list <- list()
  for (i in 1:length(v1)) {   
    input.list[[i]] <- list(ws = c(v1[i], v2[i], v3[i], v4[i]),
                            q5.ref = cm[i], sex = sex, qadult.ref = am[i], 
                            adult = adult, q0.ref = iml[i])
    if(missing(am)) input.list[[i]]$qadult.ref <- -9999
    #            if(missing(am)) input.list[[i]]$qadult.ref <- expit(aml[i])
  } 
  
  opt.out <- lapply(input.list, function (x) { optim(x$ws, fn = error.svd, 
                                                     q5.ref = x$q5.ref, sex = x$sex, qadult.ref = x$qadult.ref, 
                                                     adult = x$adult, q0.ref = x$q0.ref,
                                                     lower=c((x$ws[1] - abs(0.0001*x$ws[1])), (x$ws[2] - abs(1*x$ws[2])), (x$ws[3] - abs(.25*x$ws[3])), (x$ws[4] - abs(.25*x$ws[4]))),
                                                     upper=c((x$ws[1] + abs(0.0001*x$ws[1])), (x$ws[2] + abs(1*x$ws[2])), (x$ws[3] + abs(.25*x$ws[3])), (x$ws[4] + abs(.25*x$ws[4]))),
                                                     method="L-BFGS-B")$par})
  names(opt.out) <- 1:length(opt.out)
  
  # Predict qx values
  v <- matrix(unlist(opt.out), ncol = length(cml))
  r.p <- t(mods.r[[sex]]$components) %*% v + mods.r[[sex]][[adult]]$offset
  r.p <- data.frame(r.p)
  # splice in q0, predicted or original
  r.p[1,] <- iml
  
  return(r.p)
  
}

# compute error estimates fro optimization in svd.comp model

error.svd <- function(weights, sex, q5.ref, qadult.ref, adult, q0.ref) {
  
  adult_q <- "none"
  if (qadult.ref!=-9999 & adult == "q45") adult_q <- "q45"
  if (qadult.ref!=-9999 & adult == "q35") adult_q <- "q35"
  
  # Predict qx values
  r.p <- matrix(data = 0, ncol = 1, nrow = dim(mods.r[[sex]]$components)[2])
  for (z in 1:4) {
    r.p <- r.p + mods.r[[sex]]$components[z,] * weights[z]
    # r.p <- r.p + mods.r[[sex]]$components[z,] %*% t(weights[,z])
  }
  r.p <- r.p + mods.r[[sex]][[adult]]$offset
  
  r.p <- data.frame(r.p)
  
  # Splice in q0
  r.p[1,] <- q0.ref
  
  r.p.exp <- expit(r.p)
  
  # Predict q5
  q5.pred <- 1-sapply(1-r.p.exp[1:5,,drop=FALSE], prod)
  
  # Predict qAdult
  if (adult_q == "q45") {
    qadult.pred <- 1-sapply(1-r.p.exp[16:60,,drop=FALSE], prod)
  } else if (adult_q == "q35" ) {
    qadult.pred <- 1-sapply(1-r.p.exp[16:50,,drop=FALSE], prod)
  } else {
    qadult.pred <- qadult.ref <- 0
  }
  # Calculate and return the errors
  return((abs(q5.ref-q5.pred)*.1)+(abs(qadult.ref-qadult.pred)*.01))
}