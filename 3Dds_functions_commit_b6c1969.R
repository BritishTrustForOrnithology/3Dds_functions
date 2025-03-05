library(tidyverse)
library(truncdist)
library(fastDummies)


## detection function
if(p_meth=="hn"){
  if(vertdist=="norm"){
    full_det_func <- function(rho,phi,theta,sigma_det,mu_fh,sigma_fh)((exp(-(rho^2)/(2*(sigma_det^2))))*
                                                                        (dtrunc((rho*sin(phi)+camera_height), a = 0, b = Inf,
                                                                                mean = mu_fh, sd = sigma_fh,spec="norm"))*(rho^2*cos(phi)))


  } else if(vertdist=="cauchy"){
    full_det_func <- function(rho,phi,theta,sigma_det,mu_fh,sigma_fh)((exp(-(rho^2)/(2*(sigma_det^2))))*
                                                                        (dtrunc((rho*sin(phi)+camera_height), a = 0, b = Inf,
                                                                                location = mu_fh, scale = sigma_fh,spec="cauchy"))*(rho^2*cos(phi)))

  } else if (vertdist=="norm_mix"){
    full_det_func <- function(rho,phi,theta,sigma_det,mu_fh,sigma_fh,w.est,mu.diff,sigma_fh_2)((exp(-(rho^2)/(2*(sigma_det^2))))*
                                                                        ((w.est*dtrunc((rho*sin(phi)+camera_height), a = 0, b = Inf,
                                                                                mean = mu_fh, sd = sigma_fh,spec="norm"))+((1-w.est)*dtrunc((rho*sin(phi)+camera_height), a = 0, b = Inf,
                                                                                                                                          mean = mu_fh+mu.diff, sd = sigma_fh_2,spec="norm")))*(rho^2*cos(phi)))

  }
} else if (p_meth=="hr"){
  if(vertdist=="norm"){
    full_det_func <- function(rho,phi,theta,sigma_det,beta_det,mu_fh,sigma_fh)((1-exp(-(rho/sigma_det)^(-beta_det)))*
                                                                        (dtrunc((rho*sin(phi)+camera_height), a = 0, b = Inf,
                                                                                mean = mu_fh, sd = sigma_fh,spec="norm"))*(rho^2*cos(phi)))

  } else if(vertdist=="cauchy"){
    full_det_func <- function(rho,phi,theta,sigma_det,beta_det,mu_fh,sigma_fh)((1-exp(-(rho/sigma_det)^(-beta_det)))*
                                                                        (dtrunc((rho*sin(phi)+camera_height), a = 0, b = Inf,
                                                                                location = mu_fh, scale = sigma_fh,spec="cauchy"))*(rho^2*cos(phi)))

  }
}



## For the observation level
if (p_meth=="hn"){
  if(detcov==FALSE& fhcov==FALSE){
    glik_obs = function(parm,data){

      par.index <- 0

      sigma_det = parm[par.index+1];par.index <- par.index+1
      mu_bird = parm[par.index+1];par.index <- par.index+1
      sigma_bird = parm[par.index+1];par.index <- par.index+1

      if(vertdist=="norm_mix"){
        w.est =parm[par.index+1];par.index <- par.index+1
        mu.diff =parm[par.index+1];par.index <- par.index+1
        sigma_fh_2=parm[par.index+1]
      }


      if (vertdist=="norm_mix"){
        full_int = integral3(full_det_func, r_min,r_max,phi_min,phi_max,theta_min,theta_max, sigma_det = exp(sigma_det),mu_fh=exp(mu_bird),sigma_fh=exp(sigma_bird), w.est = expit(w.est),
                             mu.diff=exp(mu.diff),sigma_fh_2=exp(sigma_fh_2),
                             reltol=.Machine$double.eps^.05)

        obs_out = data %>%
          rowwise %>%
          mutate(p_g_pi= (full_det_func(range,phi,theta,sigma_det=exp(sigma_det),mu_fh=exp(mu_bird),sigma_fh=exp(sigma_bird),w.est=expit(w.est),mu.diff=exp(mu.diff),sigma_fh_2=exp(sigma_fh_2))/full_int))

      }else {
        full_int = integral3(full_det_func, r_min,r_max,phi_min,phi_max,theta_min,theta_max, sigma_det = exp(sigma_det),mu_fh=exp(mu_bird),sigma_fh=exp(sigma_bird),
                             reltol=.Machine$double.eps^.05)

        obs_out = data %>%
          rowwise %>%
          mutate(p_g_pi= (full_det_func(range,phi,theta,sigma_det=exp(sigma_det),mu_fh=exp(mu_bird),sigma_fh=exp(sigma_bird))/full_int))

      }



      nll = -1*sum(log(obs_out$p_g_pi))
      return(nll)

    }
  } else if(detcov!=FALSE&fhcov==FALSE){
    glik_obs = function(parm,data){
      if (class(data$var_detcov)!="numeric"){
        nfac <- length(unique(fac_choose_det))

        par.index <- 0
        beta_par <- vector()
        for (i in 1:nfac){
          beta_par[i] <- parm[par.index+1]; par.index <- par.index + 1
        }

        mu_bird = parm[par.index+1]; par.index <- par.index +1
        sigma_bird = parm[par.index+1]

        term <- vector()
        weath2 <- vector()
        term[1] <- "beta_par[1]"

        for (i in 2:nfac){
          weath2[i] = paste0(detcov,"_",fac_choose_det[i])
          term[i] = paste0("+","beta_par[",i,"]*as.vector(select(data,paste(weath2[",i,"])))[[1]]")
        }

        expr <- paste(term,collapse=" ")



        sigma_det = log(eval(parse(text=expr)))

        data$sigma_det = sigma_det
      } else {
        par.index <- 0

        beta0 = parm[par.index+1]; par.index <- par.index+1
        beta1 = parm[par.index+1]; par.index <- par.index+1
        mu_bird = parm[par.index+1]; par.index <- par.index+1
        sigma_bird = parm[par.index+1]

        sigma_det = log(beta0 + beta1*data$var_detcov)

        data$sigma_det = sigma_det

      }

      obs_out = data %>%
        rowwise %>%
        mutate(full_int=integral3(full_det_func, r_min,r_max,phi_min,phi_max,theta_min,theta_max, sigma_det = exp(sigma_det),mu_fh=exp(mu_bird),sigma_fh=exp(sigma_bird),
                                  reltol=.Machine$double.eps^.05),
               p_g_pi= (full_det_func(range,phi,theta,sigma_det=exp(sigma_det),mu_fh=exp(mu_bird),sigma_fh=exp(sigma_bird))/full_int))

      nll = -1*sum(log(obs_out$p_g_pi))
      return(nll)

    }
  } else if(detcov!=FALSE&fhcov!=FALSE){
    glik_obs = function(parm,data){
      if (class(data$var_detcov)!="numeric" & class(data$var_fhcov)!="numeric"){
        nfac <- length(unique(fac_choose_det))
        nfac_fh <- length(unique(fac_choose_fh))

        par.index <- 0
        beta_par <- vector()
        for (i in 1:nfac){
          beta_par[i] <- parm[par.index+1]; par.index <- par.index + 1
        }

        gamma_par <- vector()
        for (i in 1:nfac_fh){
          gamma_par[i] <- parm[par.index+1]; par.index <- par.index + 1
        }

        if(fhsigcov==FALSE){
          sigma_bird = parm[par.index+1]
        } else {
          eta_par <- vector()
          for (i in 1:nfac_fh){
            eta_par[i] <- parm[par.index+1]; par.index <- par.index + 1
          }
          term_fhsig <- vector()
          weath2 <- vector()
          term_fhsig[1] <- "eta_par[1]"

          for (i in 2:nfac_fh){
            weath2[i] = paste0(fhcov,"_",fac_choose_fh[i])
            term_fhsig[i] = paste0("+","eta_par[",i,"]*as.vector(select(data,paste(weath2[",i,"])))[[1]]")
          }

          expr_fhsig <- paste(term_fhsig,collapse=" ")


          sigma_bird = log(eval(parse(text=expr_fhsig)))

          data$sigma_bird = sigma_bird
          }


        term <- vector()
        weath2 <- vector()
        term[1] <- "beta_par[1]"

        for (i in 2:nfac){
          weath2[i] = paste0(detcov,"_",fac_choose_det[i])
          term[i] = paste0("+","beta_par[",i,"]*as.vector(select(data,paste(weath2[",i,"])))[[1]]")
        }

        expr <- paste(term,collapse=" ")


        sigma_det = log(eval(parse(text=expr)))

        term_fh <- vector()
        weath2 <- vector()
        term_fh[1] <- "gamma_par[1]"

        for (i in 2:nfac_fh){
          weath2[i] = paste0(fhcov,"_",fac_choose_fh[i])
          term_fh[i] = paste0("+","gamma_par[",i,"]*as.vector(select(data,paste(weath2[",i,"])))[[1]]")
        }

        expr_fh <- paste(term_fh,collapse=" ")


        mu_bird = log(eval(parse(text=expr_fh)))


        data$sigma_det = sigma_det
        data$mu_bird= mu_bird
      } else if (class(data$var_detcov)!="numeric" & class(data$var_fhcov)=="numeric"){
        nfac <- length(unique(fac_choose_det))

        par.index <- 0
        beta_par <- vector()
        for (i in 1:nfac){
          beta_par[i] <- parm[par.index+1]; par.index <- par.index + 1
        }


        gamma0 <- parm[par.index+1]; par.index <- par.index + 1
        gamma1 <- parm[par.index+1]; par.index <- par.index + 1

        sigma_bird = parm[par.index+1]

        term <- vector()
        weath2 <- vector()
        term[1] <- "beta_par[1]"

        for (i in 2:nfac){
          weath2[i] = paste0(detcov,"_",fac_choose_det[i])
          term[i] = paste0("+","beta_par[",i,"]*as.vector(select(data,paste(weath2[",i,"])))[[1]]")
        }

        expr <- paste(term,collapse=" ")


        sigma_det = log(eval(parse(text=expr)))

        mu_bird = log(gamma0 + gamma1*data$var_fhcov)

        data$sigma_det = sigma_det
        data$mu_bird= mu_bird

      } else if (class(data$var_detcov)=="numeric" & class(data$var_fhcov)=="numeric"){
        par.index <- 0

        beta0 <- parm[par.index+1]; par.index <- par.index + 1
        beta1 <- parm[par.index+1]; par.index <- par.index + 1

        gamma0 <- parm[par.index+1]; par.index <- par.index + 1
        gamma1 <- parm[par.index+1]; par.index <- par.index + 1

        sigma_bird = parm[par.index+1]


        sigma_det = log(beta0 + beta1*data$var_detcov)

        mu_bird = log(gamma0 + gamma1*data$var_fhcov)

        data$sigma_det = sigma_det
        data$mu_bird= mu_bird

      } else if (class(data$var_detcov)=="numeric" & class(data$var_fhcov)!="numeric"){

        nfac_fh <- length(unique(fac_choose_fh))

        par.index <- 0

        beta0 <- parm[par.index+1]; par.index <- par.index + 1
        beta1 <- parm[par.index+1]; par.index <- par.index + 1


        gamma_par <- vector()
        for (i in 1:nfac_fh){
          gamma_par[i] <- parm[par.index+1]; par.index <- par.index + 1
        }

        if (fhsigcov==FALSE){
          sigma_bird = parm[par.index+1]
        } else {
          eta_par <- vector()
          for (i in 1:nfac_fh){
            eta_par[i] <- parm[par.index+1]; par.index <- par.index + 1
          }
          term_fhsig <- vector()
          weath2 <- vector()
          term_fhsig[1] <- "eta_par[1]"

          for (i in 2:nfac_fh){
            weath2[i] = paste0(fhcov,"_",fac_choose_fh[i])
            term_fhsig[i] = paste0("+","eta_par[",i,"]*as.vector(select(data,paste(weath2[",i,"])))[[1]]")
          }

          expr_fhsig <- paste(term_fhsig,collapse=" ")


          sigma_bird = log(eval(parse(text=expr_fhsig)))

          data$sigma_bird <- sigma_bird
        }

        sigma_det = log(beta0 + beta1*data$var_detcov)

        term_fh <- vector()
        weath2 <- vector()
        term_fh[1] <- "gamma_par[1]"

        for (i in 2:nfac_fh){
          weath2[i] = paste0(fhcov,"_",fac_choose_fh[i])
          term_fh[i] = paste0("+","gamma_par[",i,"]*as.vector(select(data,paste(weath2[",i,"])))[[1]]")
        }

        expr_fh <- paste(term_fh,collapse=" ")


        mu_bird = log(eval(parse(text=expr_fh)))


        data$sigma_det = sigma_det
        data$mu_bird= mu_bird
      }

      obs_out = data %>%
        rowwise %>%
        mutate(full_int=integral3(full_det_func, r_min,r_max,phi_min,phi_max,theta_min,theta_max, sigma_det = exp(sigma_det),mu_fh=exp(mu_bird),sigma_fh=exp(sigma_bird),
                                  reltol=.Machine$double.eps^.05),
               p_g_pi= (full_det_func(range,phi,theta,sigma_det=exp(sigma_det),mu_fh=exp(mu_bird),sigma_fh=exp(sigma_bird))/full_int))

      nll = -1*sum(log(obs_out$p_g_pi))
      return(nll)

    }

  }

} else if (p_meth=="hr"){
  if (detcov==FALSE&fhcov==FALSE){
    glik_obs = function(parm,data){
    par.index <- 0

    sigma_det = parm[par.index+1];par.index <- par.index+1
    beta_det = parm[par.index+1];par.index <- par.index+1
    mu_bird = parm[par.index+1];par.index <- par.index+1
    sigma_bird = parm[par.index+1]


    full_int = integral3(full_det_func, r_min,r_max,phi_min,phi_max,theta_min,theta_max, sigma_det = exp(sigma_det),beta_det = exp(beta_det),mu_fh=exp(mu_bird),sigma_fh=exp(sigma_bird),
                         reltol=.Machine$double.eps^.05)

    obs_out = data %>%
      rowwise %>%
      mutate(p_g_pi= (full_det_func(range,phi,theta,sigma_det=exp(sigma_det),beta_det=exp(beta_det),mu_fh=exp(mu_bird),sigma_fh=exp(sigma_bird))/full_int))

    nll = -1*sum(log(obs_out$p_g_pi))
    return(nll)

  }
  } else if (detcov!=FALSE&fhcov==FALSE){
    glik_obs = function(parm,data){
      nfac <- length(unique(fac_choose_det))

      par.index <- 0

      beta_par <- vector()
      for (i in 1:nfac){
        beta_par[i] <- parm[par.index+1]; par.index <- par.index + 1
      }

      beta_det = parm[par.index+1]; par.index <- par.index + 1

      mu_bird = parm[par.index+1]; par.index <- par.index + 1
      sigma_bird = parm[par.index+1]


      term <- vector()
      weath2 <- vector()
      term[1] <- "beta_par[1]"

      for (i in 2:nfac){
        weath2[i] = paste0(detcov,"_",fac_choose_det[i])
        term[i] = paste0("+","beta_par[",i,"]*as.vector(select(data,paste(weath2[",i,"])))[[1]]")
      }

      expr <- paste(term,collapse=" ")


      sigma_det = log(eval(parse(text=expr)))



      data$sigma_det = sigma_det


      obs_out = data %>%
        rowwise %>%
        mutate(full_int=integral3(full_det_func, r_min,r_max,phi_min,phi_max,theta_min,theta_max, sigma_det = exp(sigma_det),beta_det=exp(beta_det),mu_fh=exp(mu_bird),sigma_fh=exp(sigma_bird),
                                  reltol=.Machine$double.eps^.05),
               p_g_pi= (full_det_func(range,phi,theta,sigma_det=exp(sigma_det),beta_det=exp(beta_det),mu_fh=exp(mu_bird),sigma_fh=exp(sigma_bird))/full_int))

      nll = -1*sum(log(obs_out$p_g_pi))
      return(nll)
    }
  } else if (detcov!=FALSE&fhcov!=FALSE){
    glik_obs = function(parm,data){
      nfac <- length(unique(fac_choose_det))
      nfac_fh <- length(unique(fac_choose_fh))

      par.index <- 0

      beta_par <- vector()
      for (i in 1:nfac){
        beta_par[i] <- parm[par.index+1]; par.index <- par.index + 1
      }

      beta_det = parm[par.index+1]; par.index <- par.index + 1

      gamma_par <- vector()
      for (i in 1:nfac_fh){
        gamma_par[i] <- parm[par.index+1]; par.index <- par.index + 1
      }

      sigma_bird = parm[par.index+1]


      term <- vector()
      weath2 <- vector()
      term[1] <- "beta_par[1]"

      for (i in 2:nfac){
        weath2[i] = paste0(detcov,"_",fac_choose_det[i])
        term[i] = paste0("+","beta_par[",i,"]*as.vector(select(data,paste(weath2[",i,"])))[[1]]")
      }

      expr <- paste(term,collapse=" ")


      sigma_det = log(eval(parse(text=expr)))


      term_fh <- vector()
      weath2 <- vector()
      term_fh[1] <- "gamma_par[1]"

      for (i in 2:nfac_fh){
        weath2[i] = paste0(fhcov,"_",fac_choose_fh[i])
        term_fh[i] = paste0("+","gamma_par[",i,"]*as.vector(select(data,paste(weath2[",i,"])))[[1]]")
      }

      expr_fh <- paste(term_fh,collapse=" ")


      mu_bird = log(eval(parse(text=expr_fh)))



      data$sigma_det = sigma_det
      data$mu_bird = mu_bird


      obs_out = data %>%
        rowwise %>%
        mutate(full_int=integral3(full_det_func, r_min,r_max,phi_min,phi_max,theta_min,theta_max, sigma_det = exp(sigma_det),beta_det=exp(beta_det),mu_fh=exp(mu_bird),sigma_fh=exp(sigma_bird),
                                  reltol=.Machine$double.eps^.05),
               p_g_pi= (full_det_func(range,phi,theta,sigma_det=exp(sigma_det),beta_det=exp(beta_det),mu_fh=exp(mu_bird),sigma_fh=exp(sigma_bird))/full_int))

      nll = -1*sum(log(obs_out$p_g_pi))
      return(nll)
    }
  }
}









hazard_rate <- function(data, sigma,beta){
 hr = 1-exp(-(data/sigma)^(-beta))
 return(hr)
}

expit <- function(x){ exp(x)/(1+exp(x)) }
logit <- function(x){ log(x/(1-x)) }


if (vertdist=="norm"){
  det_func_fh <- function(rho,phi,theta,mu_fh,sigma_fh){(dtrunc((rho*sin(phi)+camera_height), a = 0, b = Inf,
                                                                mean = mu_fh, sd = sigma_fh,spec="norm"))*(rho^2*cos(phi))}
} else if (vertdist=="cauchy"){
  det_func_fh <- function(rho,phi,theta,mu_fh,sigma_fh){(dtrunc((rho*sin(phi)+camera_height), a = 0, b = Inf,
                                                                mean = mu_fh, sd = sigma_fh,spec="cauchy"))*(rho^2*cos(phi))}
} else if (vertdist=="norm_mix"){
  det_func_fh <- function(rho,phi,theta,mu_fh,sigma_fh,w.est,mu.diff,sigma_fh_2){(w.est*(dtrunc((rho*sin(phi)+camera_height), a = 0, b = Inf,
                                                                mean = mu_fh, sd = sigma_fh,spec="norm"))+(1-w.est)*(dtrunc((rho*sin(phi)+camera_height), a = 0, b = Inf,
                                                                                                                            mean = mu_fh+mu.diff, sd = sigma_fh_2,spec="norm")))*(rho^2*cos(phi))}
}




















################################################################################
################################################################################
##model_fit function

fit_3Dds <- function(data,p_meth="hn", vertdist="norm", detcov=FALSE,fhcov=FALSE,fhsigcov=FALSE,dens_boot=FALSE, parm, control){
  if(detcov!=FALSE){data <-  data %>% mutate(var_detcov=select(data,paste(detcov))[[1]])}
  if(fhcov!=FALSE){data$var_fhcov = select(data,paste(fhcov))[[1]]}

  if(detcov!=FALSE& class(data$var_detcov)!="numeric"){fac_choose_det <<- unique(data$var_detcov)}
  if(fhcov!=FALSE & class(data$var_fhcov)!="numeric"){fac_choose_fh <<- unique(data$var_fhcov)}

  det_fit_obs = optim(par=parm, fn=glik_obs, data=data,method="Nelder-Mead",
                      hessian=TRUE,
                      control=control)


  ################################################################################
  ################################################################################
  ## Check identifiability
  ################################################################################

  Hess <- det_fit_obs$hessian
  eigout <- eigen(Hess)
  min_eig_scaled <- eigout$values[length(eigout$values)]/eigout$values[1]


  inv.Hess <- solve(Hess)
  res <- sqrt(diag(inv.Hess))



  ################################################################################
  ################################################################################
  ##Output values


  if (detcov!=FALSE & fhcov ==FALSE){
    if(class(data$var_detcov)!="numeric"){
      nfac <- length(unique(fac_choose_det))

      out.index <- 0

      sig_det <- vector()
      for (i in 1:nfac){
        sig_det[i] <- det_fit_obs$par[out.index+1]; out.index <- out.index+1
      }

      for (i in 2:nfac){
        sig_det[i] <- sig_det[1] + sig_det[i]
      }

      if(p_meth=="hr"){beta_det_out= exp(det_fit_obs$par[out.index+1]); out.index <- out.index + 1}

      mu_fh_out = exp(det_fit_obs$par[out.index+1]); out.index <- out.index+1
      sig_fh_out = exp(det_fit_obs$par[out.index+1])

      n_days <- vector()
      n_hours <- vector()
      for (i in 1:nfac){
        n_days[i] <- length(unique(data$date[data$var_detcov==unique(data$var_detcov)[i]]))
        n_hours[i] <- length(unique(data$seg_hour[data$var_detcov==unique(data$var_detcov)[i]]))
      }


      output=list(
        sig_det_out = sig_det,
        beta_det_out = if (p_meth=="hr"){beta_det_out},
        mu_fh_out = mu_fh_out,
        sig_fh_out = sig_fh_out,
        min_eig_scaled = min_eig_scaled,
        n_days = n_days,
        n_hours = n_hours,
        fac_cat = unique(data$var_detcov),
        Hess=Hess,
        integral = list(
          r_min = 0,
          r_max  =as.numeric(quantile(data$range,rigtru)),
          phi_min = min(data$phi),
          phi_max = max(data$phi),
          theta_min = min(data$theta),
          theta_max = max(data$theta)
        ),
        model=list(
          parm_raw_out = det_fit_obs$par
        ),
        ll.val = -det_fit_obs$value,
        aic =  (2*length(det_fit_obs$par))-(2*(-det_fit_obs$value))
      )
    } else {

      out.index <- 0


      sig_det <- det_fit_obs$par[out.index+1]; out.index <- out.index+1

      sig_det_slope <- det_fit_obs$par[out.index+1]; out.index <- out.index+1


      if(p_meth=="hr"){beta_det_out= exp(det_fit_obs$par[out.index+1]); out.index <- out.index + 1}

      mu_fh_out = exp(det_fit_obs$par[out.index+1]); out.index <- out.index+1
      sig_fh_out = exp(det_fit_obs$par[out.index+1])

      output=list(
        sig_det_out = sig_det,
        sig_det_slope = sig_det_slope,
        beta_det_out = if (p_meth=="hr"){beta_det_out},
        mu_fh_out = mu_fh_out,
        sig_fh_out = sig_fh_out,
        min_eig_scaled = min_eig_scaled,
        n_days = length(unique(data$date)),
        n_hours = length(unique(data$seg_hour)),
        fac_cat = unique(data$var_detcov),
        Hess=Hess,
        integral = list(
          r_min = 0,
          r_max  =as.numeric(quantile(data$range,rigtru)),
          phi_min = min(data$phi),
          phi_max = max(data$phi),
          theta_min = min(data$theta),
          theta_max = max(data$theta)
        ),
        model=list(
          parm_raw_out = det_fit_obs$par
        ),
        ll.val = -det_fit_obs$value,
        aic =  (2*length(det_fit_obs$par))-(2*(-det_fit_obs$value))
      )
    }

  } else if (detcov!=FALSE & fhcov !=FALSE){
    if (class(data$var_detcov)!="numeric" & class(data$var_fhcov)!="numeric"){
      nfac_det <- length(unique(fac_choose_det))
      nfac_fh <- length(unique(fac_choose_fh))

      out.index <- 0

      sig_det <- vector()
      for (i in 1:nfac_det){
        sig_det[i] <- det_fit_obs$par[out.index+1]; out.index <- out.index+1
      }

      for (i in 2:nfac_det){
        sig_det[i] <- sig_det[1] + sig_det[i]
      }

      if(p_meth=="hr"){beta_det_out= exp(det_fit_obs$par[out.index+1]); out.index <- out.index + 1}


      mu_fh <- vector()
      for (i in 1:nfac_fh){
        mu_fh[i] <- det_fit_obs$par[out.index+1]; out.index <- out.index+1
      }

      for (i in 2:nfac_fh){
        mu_fh[i] <- mu_fh[1] + mu_fh[i]
      }

      if (fhsigcov==FALSE){
        sig_fh_out = exp(det_fit_obs$par[out.index+1])
      } else {
        sig_fh <- vector()
        for (i in 1:nfac_fh){
          sig_fh[i] <- det_fit_obs$par[out.index+1]; out.index <- out.index+1
        }

        for (i in 2:nfac_fh){
          sig_fh[i] <- sig_fh[1] + sig_fh[i]
        }

        sig_fh_out <- sig_fh
      }


      n_days <- vector()
      n_hours <- vector()
      for (i in 1:nfac){
        n_days[i] <- length(unique(data$date[data$var_detcov==unique(data$var_detcov)[i]]))
        n_hours[i] <- length(unique(data$seg_hour[data$var_detcov==unique(data$var_detcov)[i]]))
      }

      output=list(
        sig_det_out = sig_det,
        beta_det_out = if (p_meth=="hr"){beta_det_out},
        mu_fh_out = mu_fh,
        sig_fh_out = sig_fh_out,
        min_eig_scaled = min_eig_scaled,
        n_days=n_days,
        n_hours = n_hours,
        fac_cat = unique(data$var_detcov),
        Hess=Hess,
        integral = list(
          r_min = 0,
          r_max  =as.numeric(quantile(data$range,rigtru)),
          phi_min = min(data$phi),
          phi_max = max(data$phi),
          theta_min = min(data$theta),
          theta_max = max(data$theta)
        ),
        model=list(
          parm_raw_out = det_fit_obs$par
        ),
        ll.val = -det_fit_obs$value,
        aic =  (2*length(det_fit_obs$par))-(2*(-det_fit_obs$value))
      )
    } else if (class(data$var_detcov)=="numeric" & class(data$var_fhcov)=="numeric") {
      out.index <- 0

        sig_det <- det_fit_obs$par[out.index+1]; out.index <- out.index+1
        sig_det_slope <- det_fit_obs$par[out.index+1]; out.index <- out.index+1

      if(p_meth=="hr"){beta_det_out= exp(det_fit_obs$par[out.index+1]); out.index <- out.index + 1}

      mu_fh <- det_fit_obs$par[out.index+1]; out.index <- out.index+1

      mu_fh_slope<- det_fit_obs$par[out.index+1]; out.index <- out.index+1

      sig_fh_out = exp(det_fit_obs$par[out.index+1])


      output=list(
        sig_det_out = sig_det,
        sig_det_slope=sig_det_slope,
        beta_det_out = if (p_meth=="hr"){beta_det_out},
        mu_fh_out = mu_fh,
        mu_fh_slope=mu_fh_slope,
        sig_fh_out = sig_fh_out,
        min_eig_scaled = min_eig_scaled,
        n_days = length(unique(data$date)),
        n_hours = length(unique(data$seg_hour)),
        fac_cat = unique(data$var_detcov),
        Hess=Hess,
        integral = list(
          r_min = 0,
          r_max  =as.numeric(quantile(data$range,rigtru)),
          phi_min = min(data$phi),
          phi_max = max(data$phi),
          theta_min = min(data$theta),
          theta_max = max(data$theta)
        ),
        model=list(
          parm_raw_out = det_fit_obs$par
        ),
        ll.val = -det_fit_obs$value,
        aic =  (2*length(det_fit_obs$par))-(2*(-det_fit_obs$value))
      )
    }else if (class(data$var_detcov)!="numeric" & class(data$var_fhcov)=="numeric"){
      nfac_det <- length(unique(fac_choose_det))


      out.index <- 0

      sig_det <- vector()
      for (i in 1:nfac_det){
        sig_det[i] <- det_fit_obs$par[out.index+1]; out.index <- out.index+1
      }

      for (i in 2:nfac_det){
        sig_det[i] <- sig_det[1] + sig_det[i]
      }

      if(p_meth=="hr"){beta_det_out= exp(det_fit_obs$par[out.index+1]); out.index <- out.index + 1}



        mu_fh <- det_fit_obs$par[out.index+1]; out.index <- out.index+1

        mu_fh_slope <- det_fit_obs$par[out.index+1]; out.index <- out.index+1



      sig_fh_out = exp(det_fit_obs$par[out.index+1])

      n_days <- vector()
      n_hours <- vector()
      for (i in 1:nfac){
        n_days[i] <- length(unique(data$date[data$var_detcov==unique(data$var_detcov)[i]]))
        n_hours[i] <- length(unique(data$seg_hour[data$var_detcov==unique(data$var_detcov)[i]]))
      }

      output=list(
        sig_det_out = sig_det,
        beta_det_out = if (p_meth=="hr"){beta_det_out},
        mu_fh_out = mu_fh,
        mu_fh_slope= mu_fh_slope,
        sig_fh_out = sig_fh_out,
        min_eig_scaled = min_eig_scaled,
        n_days=n_days,
        n_hours = n_hours,
        fac_cat = unique(data$var_detcov),
        Hess=Hess,
        integral = list(
          r_min = 0,
          r_max  =as.numeric(quantile(data$range,rigtru)),
          phi_min = min(data$phi),
          phi_max = max(data$phi),
          theta_min = min(data$theta),
          theta_max = max(data$theta)
        ),
        model=list(
          parm_raw_out = det_fit_obs$par
        ),
        ll.val = -det_fit_obs$value,
        aic =  (2*length(det_fit_obs$par))-(2*(-det_fit_obs$value))
      )
    }else if (class(data$var_detcov)=="numeric" & class(data$var_fhcov)!="numeric"){
      nfac_fh <- length(unique(fac_choose_fh))

      out.index <- 0


      sig_det <- det_fit_obs$par[out.index+1]; out.index <- out.index+1

      sig_det_slope <-  det_fit_obs$par[out.index+1]; out.index <- out.index+1

      if(p_meth=="hr"){beta_det_out= exp(det_fit_obs$par[out.index+1]); out.index <- out.index + 1}


      mu_fh <- vector()
      for (i in 1:nfac_fh){
        mu_fh[i] <- det_fit_obs$par[out.index+1]; out.index <- out.index+1
      }

      for (i in 2:nfac_fh){
        mu_fh[i] <- mu_fh[1] + mu_fh[i]
      }

      if (fhsigcov==FALSE){
        sig_fh_out = exp(det_fit_obs$par[out.index+1])
      } else {
        sig_fh <- vector()
        for (i in 1:nfac_fh){
          sig_fh[i] <- det_fit_obs$par[out.index+1]; out.index <- out.index+1
        }

        for (i in 2:nfac_fh){
          sig_fh[i] <- sig_fh[1] + sig_fh[i]
        }

        sig_fh_out <- sig_fh
      }


      n_days <- vector()
      n_hours <- vector()
      for (i in 1:nfac){
        n_days[i] <- length(unique(data$date[data$var_fhcov==unique(data$var_fhcov)[i]]))
        n_hours[i] <- length(unique(data$seg_hour[data$var_fhcov==unique(data$var_fhcov)[i]]))
      }

      output=list(
        sig_det_out = sig_det,
        sig_det_slope=sig_det_slope,
        beta_det_out = if (p_meth=="hr"){beta_det_out},
        mu_fh_out = mu_fh,
        sig_fh_out = sig_fh_out,
        min_eig_scaled = min_eig_scaled,
        n_days=n_days,
        n_hours = n_hours,
        fac_cat = unique(data$var_detcov),
        Hess=Hess,
        integral = list(
          r_min = 0,
          r_max  =as.numeric(quantile(data$range,rigtru)),
          phi_min = min(data$phi),
          phi_max = max(data$phi),
          theta_min = min(data$theta),
          theta_max = max(data$theta)
        ),
        model=list(
          parm_raw_out = det_fit_obs$par
        ),
        ll.val = -det_fit_obs$value,
        aic =  (2*length(det_fit_obs$par))-(2*(-det_fit_obs$value))
      )
    }

  } else if (detcov==FALSE & fhcov==FALSE) {

    out.index <- 0

    sig_det_out = exp(det_fit_obs$par[out.index+1]); out.index <- out.index + 1
    if(p_meth=="hr"){beta_det_out= exp(det_fit_obs$par[out.index+1]); out.index <- out.index + 1}
    mu_fh_out = exp(det_fit_obs$par[out.index+1]); out.index <- out.index + 1
    sig_fh_out = exp(det_fit_obs$par[out.index+1]); out.index <- out.index+1
    if (vertdist=="norm_mix"){
      w.out = expit(det_fit_obs$par[out.index+1]); out.index <- out.index+1
      mu_diff_out = exp(det_fit_obs$par[out.index+1]); out.index <- out.index+1
      sig_2_out = expit(det_fit_obs$par[out.index+1])

      mu_2_out = mu_diff_out + mu_fh_out
    }


    output=list(
      sig_det_out = sig_det_out,
      beta_det_out = if (p_meth=="hr"){beta_det_out},
      mu_fh_out = mu_fh_out,
      sig_fh_out = sig_fh_out,
      w.out=if(vertdist=="norm_mix"){w.out},
      mu_2_out=if(vertdist=="norm_mix"){mu_2_out},
      sig_2_out=if(vertdist=="norm_mix"){sig_2_out},
      min_eig_scaled = min_eig_scaled,
      n_days = length(unique(data$date)),
      n_hours = length(unique(data$seg_hour)),
      Hess=Hess,
      integral = list(
        r_min = 0,
        r_max  =as.numeric(quantile(data$range,rigtru)),
        phi_min = min(data$phi),
        phi_max = max(data$phi),
        theta_min = min(data$theta),
        theta_max = max(data$theta)
      ),
      model=list(
        parm_raw_out = det_fit_obs$par
      ),
      ll.val = -det_fit_obs$value,
      aic =  (2*length(det_fit_obs$par))-(2*(-det_fit_obs$value))
    )
  }


  ##############################################################################
  ##############################################################################
  ## CI calculations
  ##############################################################################


  if (detcov==FALSE&fhcov==FALSE){
    out.index <- 0

    se_sigdet <- output[["sig_det_out"]]*res[out.index+1];out.index <- out.index + 1
    if(p_meth=="hr"){se_betadet <- output[["beta_det_out"]]*res[out.index+1];out.index <- out.index + 1}
    se_mufh <- output[["mu_fh_out"]]*res[out.index+1];out.index <- out.index + 1
    se_sigfh <- output[["sig_fh_out"]]*res[out.index+1];out.index <- out.index + 1
    if (vertdist=="norm_mix"){
      se_w <- output[["w.out"]]*res[out.index+1];out.index <- out.index + 1
      se_mu_2 <- output[["mu_2_out"]]*res[out.index+1];out.index <- out.index + 1
      se_sig_2 <- output[["sig_2_out"]]*res[out.index+1]
    }


    sig_det_up <- output[["sig_det_out"]]+1.96*se_sigdet
    sig_det_down <- output[["sig_det_out"]]-1.96*se_sigdet

    ci_sig_det <- c(sig_det_down,sig_det_up)

    if(p_meth=="hr"){
      beta_det_up <- output[["beta_det_out"]]+1.96*se_betadet
      beta_det_down <- output[["beta_det_out"]]-1.96*se_betadet

      ci_beta_det <- c(beta_det_down,beta_det_up)
      output$ci_beta_det<- ci_beta_det
    }

    mu_fh_up <- output[["mu_fh_out"]]+1.96*se_mufh
    mu_fh_down <- output[["mu_fh_out"]]-1.96*se_mufh

    ci_mu_fh <- c(mu_fh_down,mu_fh_up)

    sig_fh_up <- output[["sig_fh_out"]]+1.96*se_sigfh
    sig_fh_down <- output[["sig_fh_out"]]-1.96*se_sigfh

    ci_sig_fh <- c(sig_fh_down,sig_fh_up)

    if (vertdist=="norm_mix"){
      mu_2_up <- output[["mu_2_out"]]+1.96*(se_mu_2+se_mufh)
      mu_2_down <- output[["mu_2_out"]]-1.96*(se_mu_2+se_mufh)

      ci_mu_2 <- c(mu_2_down,mu_2_up)


      sig_2_up <- output[["sig_2_out"]]+1.96*se_sig_2
      sig_2_down <- output[["sig_2_out"]]-1.96*se_sig_2

      ci_sig_2 <- c(sig_2_down,sig_2_up)
    }


    output$ci_sig_det <- ci_sig_det
    output$ci_mu_fh <- ci_mu_fh
    output$ci_sig_fh <- ci_sig_fh
    if(vertdist=="norm_mix"){
      output$ci_mu_2 <- ci_mu_2
      output$ci_sig_2 <- ci_sig_2
    }
  } else if (detcov!=FALSE & fhcov!=FALSE){
    if (class(data$var_detcov)!="numeric" & class(data$var_fhcov)!="numeric"){
      out.index <- 0
      nfac_det <- length(unique(fac_choose_det))
      nfac_fh <- length(unique(fac_choose_fh))

      ci_sig_det <- vector(mode='list', length=nfac_det)
      ci_mu_fh <- vector(mode='list', length=nfac_fh)


      for (i in 1:nfac_det){
        se_sig_det <- res[out.index+1]; out.index <- out.index+1
        se_sig_det_down <- output[["sig_det_out"]][i]-1.96*se_sig_det
        se_sig_det_up <- output[["sig_det_out"]][i]+1.96*se_sig_det
        ci_sig_det[[i]] <- c(se_sig_det_down,se_sig_det_up)
      }

      if(p_meth=="hr"){se_betadet <- output[["beta_det_out"]]*res[out.index+1];out.index <- out.index + 1}


      if(p_meth=="hr"){
        beta_det_up <- output[["beta_det_out"]]+1.96*se_betadet
        beta_det_down <- output[["beta_det_out"]]-1.96*se_betadet

        ci_beta_det <- c(beta_det_down,beta_det_up)
        output$ci_beta_det<- ci_beta_det
      }


      for (i in 1:nfac_fh){
        se_mu_fh <- res[out.index+1]; out.index <- out.index+1
        se_mu_fh_down <- output[["mu_fh_out"]][i]-1.96*se_mu_fh
        se_mu_fh_up <- output[["mu_fh_out"]][i]+1.96*se_mu_fh
        ci_mu_fh[[i]] <- c(se_mu_fh_down,se_mu_fh_up)
      }

      if (fhsigcov==FALSE){
        se_sigfh <- output[["sig_fh_out"]]*res[out.index+1]

        sig_fh_up <- output[["sig_fh_out"]]+1.96*se_sigfh
        sig_fh_down <- output[["sig_fh_out"]]-1.96*se_sigfh

        ci_sig_fh <- c(sig_fh_down,sig_fh_up)
      } else {
        ci_sig_fh <- vector(mode='list', length=nfac_fh)
        for (i in 1:nfac_fh){
          se_sig_fh <- res[out.index+1]; out.index <- out.index+1
          se_sig_fh_down <- output[["sig_fh_out"]][i]-1.96*se_sig_fh
          se_sig_fh_up <- output[["sig_fh_out"]][i]+1.96*se_sig_fh
          ci_sig_fh[[i]] <- c(se_sig_fh_down,se_sig_fh_up)
        }
      }


      output$ci_sig_det <- ci_sig_det
      output$ci_mu_fh <- ci_mu_fh
      output$ci_sig_fh <- ci_sig_fh
    } else if (class(data$var_detcov)!="numeric" & class(data$var_fhcov)=="numeric"){
      out.index <- 0
      nfac_det <- length(unique(fac_choose_det))

      ci_sig_det <- vector(mode='list', length=nfac_det)


      for (i in 1:nfac_det){
        se_sig_det <- res[out.index+1]; out.index <- out.index+1
        se_sig_det_down <- output[["sig_det_out"]][i]-1.96*se_sig_det
        se_sig_det_up <- output[["sig_det_out"]][i]+1.96*se_sig_det
        ci_sig_det[[i]] <- c(se_sig_det_down,se_sig_det_up)
      }

      if(p_meth=="hr"){se_betadet <- output[["beta_det_out"]]*res[out.index+1];out.index <- out.index + 1}


      if(p_meth=="hr"){
        beta_det_up <- output[["beta_det_out"]]+1.96*se_betadet
        beta_det_down <- output[["beta_det_out"]]-1.96*se_betadet

        ci_beta_det <- c(beta_det_down,beta_det_up)
        output$ci_beta_det<- ci_beta_det
      }



        se_mu_fh <- res[out.index+1]; out.index <- out.index+1
        se_mu_fh_down <- output[["mu_fh_out"]]-1.96*se_mu_fh
        se_mu_fh_up <- output[["mu_fh_out"]]+1.96*se_mu_fh
        ci_mu_fh <- c(se_mu_fh_down,se_mu_fh_up)

        se_mu_fh_slope <- res[out.index+1]; out.index <- out.index+1
        se_mu_fh_slope_down <- output[["mu_fh_slope"]]-1.96*se_mu_fh_slope
        se_mu_fh_slope_up <- output[["mu_fh_slope"]]+1.96*se_mu_fh_slope
        ci_mu_fh_slope <- c(se_mu_fh_slope_down,se_mu_fh_slope_up)


      se_sigfh <- output[["sig_fh_out"]]*res[out.index+1]

      sig_fh_up <- output[["sig_fh_out"]]+1.96*se_sigfh
      sig_fh_down <- output[["sig_fh_out"]]-1.96*se_sigfh

      ci_sig_fh <- c(sig_fh_down,sig_fh_up)


      output$ci_sig_det <- ci_sig_det
      output$ci_mu_fh <- ci_mu_fh
      output$ci_mu_fh_slope <- ci_mu_fh_slope
      output$ci_sig_fh <- ci_sig_fh
    } else if (class(data$var_detcov)=="numeric" & class(data$var_fhcov)!="numeric"){
      out.index <- 0

      nfac_fh <- length(unique(fac_choose_fh))

      ci_mu_fh <- vector(mode='list', length=nfac_fh)



        se_sig_det <- res[out.index+1]; out.index <- out.index+1
        se_sig_det_down <- output[["sig_det_out"]]-1.96*se_sig_det
        se_sig_det_up <- output[["sig_det_out"]]+1.96*se_sig_det
        ci_sig_det <- c(se_sig_det_down,se_sig_det_up)

        se_sig_det_slope <- res[out.index+1]; out.index <- out.index+1
        se_sig_det_slope_down <- output[["sig_det_slope"]]-1.96*se_sig_det_slope
        se_sig_det_slope_up <- output[["sig_det_slope"]]+1.96*se_sig_det_slope
        ci_sig_det_slope <- c(se_sig_det_slope_down,se_sig_det_slope_up)


      if(p_meth=="hr"){se_betadet <- output[["beta_det_out"]]*res[out.index+1];out.index <- out.index + 1}


      if(p_meth=="hr"){
        beta_det_up <- output[["beta_det_out"]]+1.96*se_betadet
        beta_det_down <- output[["beta_det_out"]]-1.96*se_betadet

        ci_beta_det <- c(beta_det_down,beta_det_up)
        output$ci_beta_det<- ci_beta_det
      }


      for (i in 1:nfac_fh){
        se_mu_fh <- res[out.index+1]; out.index <- out.index+1
        se_mu_fh_down <- output[["mu_fh_out"]][i]-1.96*se_mu_fh
        se_mu_fh_up <- output[["mu_fh_out"]][i]+1.96*se_mu_fh
        ci_mu_fh[[i]] <- c(se_mu_fh_down,se_mu_fh_up)
      }

      if (fhsigcov==FALSE){
        se_sigfh <- output[["sig_fh_out"]]*res[out.index+1]

        sig_fh_up <- output[["sig_fh_out"]]+1.96*se_sigfh
        sig_fh_down <- output[["sig_fh_out"]]-1.96*se_sigfh

        ci_sig_fh <- c(sig_fh_down,sig_fh_up)
      } else {
        ci_sig_fh <- vector(mode='list', length=nfac_fh)
        for (i in 1:nfac_fh){
          se_sig_fh <- res[out.index+1]; out.index <- out.index+1
          se_sig_fh_down <- output[["sig_fh_out"]][i]-1.96*se_sig_fh
          se_sig_fh_up <- output[["sig_fh_out"]][i]+1.96*se_sig_fh
          ci_sig_fh[[i]] <- c(se_sig_fh_down,se_sig_fh_up)
        }
      }



      output$ci_sig_det <- ci_sig_det
      output$ci_sig_det_slope <- ci_sig_det_slope
      output$ci_mu_fh <- ci_mu_fh
      output$ci_sig_fh <- ci_sig_fh
    } else if (class(data$var_detcov)=="numeric" & class(data$var_fhcov)=="numeric"){
      out.index <- 0




        se_sig_det <- res[out.index+1]; out.index <- out.index+1
        se_sig_det_down <- output[["sig_det_out"]]-1.96*se_sig_det
        se_sig_det_up <- output[["sig_det_out"]]+1.96*se_sig_det
        ci_sig_det <- c(se_sig_det_down,se_sig_det_up)

        se_sig_det_slope <- res[out.index+1]; out.index <- out.index+1
        se_sig_det_slope_down <- output[["sig_det_slope"]]-1.96*se_sig_det_slope
        se_sig_det_slope_up <- output[["sig_det_slope"]]+1.96*se_sig_det_slope
        ci_sig_det_slope <- c(se_sig_det_slope_down,se_sig_det_slope_up)


      if(p_meth=="hr"){se_betadet <- output[["beta_det_out"]]*res[out.index+1];out.index <- out.index + 1}


      if(p_meth=="hr"){
        beta_det_up <- output[["beta_det_out"]]+1.96*se_betadet
        beta_det_down <- output[["beta_det_out"]]-1.96*se_betadet

        ci_beta_det <- c(beta_det_down,beta_det_up)
        output$ci_beta_det<- ci_beta_det
      }



        se_mu_fh <- res[out.index+1]; out.index <- out.index+1
        se_mu_fh_down <- output[["mu_fh_out"]]-1.96*se_mu_fh
        se_mu_fh_up <- output[["mu_fh_out"]]+1.96*se_mu_fh
        ci_mu_fh <- c(se_mu_fh_down,se_mu_fh_up)

        se_mu_fh_slope <- res[out.index+1]; out.index <- out.index+1
        se_mu_fh_slope_down <- output[["mu_fh_slope"]]-1.96*se_mu_fh_slope
        se_mu_fh_slope_up <- output[["mu_fh_slope"]]+1.96*se_mu_fh_slope
        ci_mu_fh_slope <- c(se_mu_fh_slope_down,se_mu_fh_slope_up)



      se_sigfh <- output[["sig_fh_out"]]*res[out.index+1]

      sig_fh_up <- output[["sig_fh_out"]]+1.96*se_sigfh
      sig_fh_down <- output[["sig_fh_out"]]-1.96*se_sigfh

      ci_sig_fh <- c(sig_fh_down,sig_fh_up)


      output$ci_sig_det <- ci_sig_det
      output$ci_sig_det_slope <- ci_sig_det_slope
      output$ci_mu_fh <- ci_mu_fh
      output$mu_fh_slope <- ci_mu_fh_slope
      output$ci_sig_fh <- ci_sig_fh
    }


  }else if (detcov!=FALSE & fhcov==FALSE){
   if (class(data$var_detcov)!="numeric"){
     out.index <- 0
     nfac <- length(unique(fac_choose_det))

     ci_sig_det <- vector(mode='list', length=nfac)
     ci_mu_fh <- vector(mode='list', length=nfac)


     for (i in 1:nfac){
       se_sig_det <- res[out.index+1]; out.index <- out.index+1
       se_sig_det_down <- output[["sig_det_out"]][i]-1.96*se_sig_det
       se_sig_det_up <- output[["sig_det_out"]][i]+1.96*se_sig_det
       ci_sig_det[[i]] <- c(se_sig_det_down,se_sig_det_up)
     }

     if(p_meth=="hr"){se_betadet <- output[["beta_det_out"]]*res[out.index+1];out.index <- out.index + 1}


     if(p_meth=="hr"){
       beta_det_up <- output[["beta_det_out"]]+1.96*se_betadet
       beta_det_down <- output[["beta_det_out"]]-1.96*se_betadet

       ci_beta_det <- c(beta_det_down,beta_det_up)
       output$ci_beta_det <- ci_beta_det
     }


     se_mu_fh <- output[["mu_fh_out"]]*res[out.index+1]; out.index <- out.index+1

     mu_fh_up <- output[["mu_fh_out"]]+1.96*se_mu_fh
     mu_fh_down <- output[["mu_fh_out"]]-1.96*se_mu_fh

     ci_mu_fh <- c(mu_fh_down,mu_fh_up)



     se_sigfh <- output[["sig_fh_out"]]*res[out.index+1]

     sig_fh_up <- output[["sig_fh_out"]]+1.96*se_sigfh
     sig_fh_down <- output[["sig_fh_out"]]-1.96*se_sigfh

     ci_sig_fh <- c(sig_fh_down,sig_fh_up)


     output$ci_sig_det <- ci_sig_det
     output$ci_mu_fh <- ci_mu_fh
     output$ci_sig_fh <- ci_sig_fh
   } else {
     out.index <- 0

       se_sig_det <- res[out.index+1]; out.index <- out.index+1
       se_sig_det_down <- output[["sig_det_out"]]-1.96*se_sig_det
       se_sig_det_up <- output[["sig_det_out"]]+1.96*se_sig_det
       ci_sig_det <- c(se_sig_det_down,se_sig_det_up)


       se_sig_det_slope <- res[out.index+1]; out.index <- out.index+1
       se_sig_det_slope_down <- output[["sig_det_slope"]]-1.96*se_sig_det_slope
       se_sig_det_slope_up <- output[["sig_det_slope"]]+1.96*se_sig_det_slope
       ci_sig_det_slope <- c(se_sig_det_slope_down,se_sig_det_slope_up)


     if(p_meth=="hr"){se_betadet <- output[["beta_det_out"]]*res[out.index+1];out.index <- out.index + 1}


     if(p_meth=="hr"){
       beta_det_up <- output[["beta_det_out"]]+1.96*se_betadet
       beta_det_down <- output[["beta_det_out"]]-1.96*se_betadet

       ci_beta_det <- c(beta_det_down,beta_det_up)
       output$ci_beta_det <- ci_beta_det
     }


     se_mu_fh <- output[["mu_fh_out"]]*res[out.index+1]; out.index <- out.index+1

     mu_fh_up <- output[["mu_fh_out"]]+1.96*se_mu_fh
     mu_fh_down <- output[["mu_fh_out"]]-1.96*se_mu_fh

     ci_mu_fh <- c(mu_fh_down,mu_fh_up)



     se_sigfh <- output[["sig_fh_out"]]*res[out.index+1]

     sig_fh_up <- output[["sig_fh_out"]]+1.96*se_sigfh
     sig_fh_down <- output[["sig_fh_out"]]-1.96*se_sigfh

     ci_sig_fh <- c(sig_fh_down,sig_fh_up)


     output$ci_sig_det <- ci_sig_det
     output$ci_sig_det_slope <- ci_sig_det_slope
     output$ci_mu_fh <- ci_mu_fh
     output$ci_sig_fh <- ci_sig_fh
   }


  }

  ################################################################################
  ################################################################################
  ##Density calculation
  ################################################################################
  r_min <- output[["integral"]][["r_min"]]
  r_max <- output[["integral"]][["r_max"]]
  phi_min <- output[["integral"]][["phi_min"]]
  phi_max <- output[["integral"]][["phi_max"]]
  theta_min <- output[["integral"]][["theta_min"]]
  theta_max <- output[["integral"]][["theta_max"]]

  if(detcov!=FALSE & fhcov!=FALSE){
    if (class(data$var_detcov)!="numeric" & class(data$var_fhcov)!="numeric"){
      nfac_det <- length(unique(fac_choose_det))
      nfac_fh <- length(unique(fac_choose_fh))

      p_hat_cov <- vector()

      if(p_meth=="hn"){
        for (i in 1:nfac_det){
          if (fhsigcov==FALSE){
            p_hat_cov[i]<- integral3(full_det_func, r_min,r_max,phi_min,phi_max,theta_min,theta_max,
                                     sigma_det = output[["sig_det_out"]][i], mu_fh=output[["mu_fh_out"]][i], sigma_fh=output[["sig_fh_out"]],
                                     reltol=.Machine$double.eps^.05)
          }else{
            p_hat_cov[i]<- integral3(full_det_func, r_min,r_max,phi_min,phi_max,theta_min,theta_max,
                                     sigma_det = output[["sig_det_out"]][i], mu_fh=output[["mu_fh_out"]][i], sigma_fh=output[["sig_fh_out"]][i],
                                     reltol=.Machine$double.eps^.05)
          }

        }
      } else if (p_meth=="hr"){
        for (i in 1:nfac_det){
          if(fhsigcov==FALSE){
            p_hat_cov[i]<- integral3(full_det_func, r_min,r_max,phi_min,phi_max,theta_min,theta_max,
                                     sigma_det = output[["sig_det_out"]][i], beta_det = output[["beta_det_out"]], mu_fh=output[["mu_fh_out"]][i], sigma_fh=output[["sig_fh_out"]],
                                     reltol=.Machine$double.eps^.05)
          } else{
            p_hat_cov[i]<- integral3(full_det_func, r_min,r_max,phi_min,phi_max,theta_min,theta_max,
                                     sigma_det = output[["sig_det_out"]][i], beta_det = output[["beta_det_out"]], mu_fh=output[["mu_fh_out"]][i], sigma_fh=output[["sig_fh_out"]][i],
                                     reltol=.Machine$double.eps^.05)
          }

        }
      }


      p_hat_fh <- vector()
      for (i in 1:nfac_det){
        if (fhsigcov==FALSE){
          p_hat_fh[i] <- integral3(det_func_fh, r_min,r_max,phi_min,phi_max,theta_min,theta_max,
                                   mu_fh=output[["mu_fh_out"]][i], sigma_fh=output[["sig_fh_out"]],
                                   reltol=.Machine$double.eps^.05)
        }else{
          p_hat_fh[i] <- integral3(det_func_fh, r_min,r_max,phi_min,phi_max,theta_min,theta_max,
                                   mu_fh=output[["mu_fh_out"]][i], sigma_fh=output[["sig_fh_out"]][i],
                                   reltol=.Machine$double.eps^.05)
        }

      }


      vol_sph <- function(r,phi,theta) (r^2*cos(phi))
      volume_obs=integral3(vol_sph,r_min,r_max,phi_min,phi_max,theta_min,theta_max)

      count_fov <- vector()
      for (i in 1:nfac_det){
        count_fov[i] = sum(data$var_detcov==fac_choose_det[i])
      }

      dens_fov <- vector()
      for (i in 1:nfac_det){
        dens_fov[i] = ((count_fov[i]*p_hat_fh[i]*12)/(volume_obs*p_hat_cov[i]*output[["n_hours"]][i]))*10000000
      }

      output$dens_hat = dens_fov
    } else if (class(data$var_detcov)!="numeric" & class(data$var_fhcov)=="numeric"){
      nfac_det <- length(unique(fac_choose_det))

      p_hat_cov <- vector()

      if(p_meth=="hn"){
        for (i in 1:nfac_det){
          p_hat_cov[i]<- integral3(full_det_func, r_min,r_max,phi_min,phi_max,theta_min,theta_max,
                                   sigma_det = output[["sig_det_out"]][i], mu_fh=output[["mu_fh_out"]], sigma_fh=output[["sig_fh_out"]],
                                   reltol=.Machine$double.eps^.05)
        }
      } else if (p_meth=="hr"){
        for (i in 1:nfac_det){
          p_hat_cov[i]<- integral3(full_det_func, r_min,r_max,phi_min,phi_max,theta_min,theta_max,
                                   sigma_det = output[["sig_det_out"]][i], beta_det = output[["beta_det_out"]], mu_fh=output[["mu_fh_out"]], sigma_fh=output[["sig_fh_out"]],
                                   reltol=.Machine$double.eps^.05)
        }
      }



        p_hat_fh <- integral3(det_func_fh, r_min,r_max,phi_min,phi_max,theta_min,theta_max,
                                 mu_fh=output[["mu_fh_out"]], sigma_fh=output[["sig_fh_out"]],
                                 reltol=.Machine$double.eps^.05)



      vol_sph <- function(r,phi,theta) (r^2*cos(phi))
      volume_obs=integral3(vol_sph,r_min,r_max,phi_min,phi_max,theta_min,theta_max)

      count_fov <- vector()
      for (i in 1:nfac_det){
        count_fov[i] = sum(data$var_detcov==fac_choose_det[i])
      }

      dens_fov <- vector()
      for (i in 1:nfac_det){
        dens_fov[i] = ((count_fov[i]*p_hat_fh*12)/(volume_obs*p_hat_cov[i]*output[["n_hours"]][i]))*10000000
      }

      output$dens_hat = dens_fov
    }else if (class(data$var_detcov)=="numeric" & class(data$var_fhcov)!="numeric"){

      nfac_fh <- length(unique(fac_choose_fh))

      p_hat_cov <- vector()

      if(p_meth=="hn"){
        for (i in 1:nfac_fh){
          if(fhsigcov==FALSE){
            p_hat_cov[i]<- integral3(full_det_func, r_min,r_max,phi_min,phi_max,theta_min,theta_max,
                                     sigma_det = output[["sig_det_out"]], mu_fh=output[["mu_fh_out"]][i], sigma_fh=output[["sig_fh_out"]],
                                     reltol=.Machine$double.eps^.05)
          } else{
            p_hat_cov[i]<- integral3(full_det_func, r_min,r_max,phi_min,phi_max,theta_min,theta_max,
                                     sigma_det = output[["sig_det_out"]], mu_fh=output[["mu_fh_out"]][i], sigma_fh=output[["sig_fh_out"]][i],
                                     reltol=.Machine$double.eps^.05)
          }

        }
      } else if (p_meth=="hr"){
        for (i in 1:nfac_fh){
          if (fhsigcov==FALSE){
            p_hat_cov[i]<- integral3(full_det_func, r_min,r_max,phi_min,phi_max,theta_min,theta_max,
                                     sigma_det = output[["sig_det_out"]], beta_det = output[["beta_det_out"]], mu_fh=output[["mu_fh_out"]][i], sigma_fh=output[["sig_fh_out"]],
                                     reltol=.Machine$double.eps^.05)
          } else{
            p_hat_cov[i]<- integral3(full_det_func, r_min,r_max,phi_min,phi_max,theta_min,theta_max,
                                     sigma_det = output[["sig_det_out"]], beta_det = output[["beta_det_out"]], mu_fh=output[["mu_fh_out"]][i], sigma_fh=output[["sig_fh_out"]][i],
                                     reltol=.Machine$double.eps^.05)
          }

        }
      }


      p_hat_fh <- vector()
      for (i in 1:nfac_fh){
        if (fhsigcov==FALSE){
          p_hat_fh[i] <- integral3(det_func_fh, r_min,r_max,phi_min,phi_max,theta_min,theta_max,
                                   mu_fh=output[["mu_fh_out"]][i], sigma_fh=output[["sig_fh_out"]],
                                   reltol=.Machine$double.eps^.05)
        } else {
          p_hat_fh[i] <- integral3(det_func_fh, r_min,r_max,phi_min,phi_max,theta_min,theta_max,
                                   mu_fh=output[["mu_fh_out"]][i], sigma_fh=output[["sig_fh_out"]][i],
                                   reltol=.Machine$double.eps^.05)
        }

      }


      vol_sph <- function(r,phi,theta) (r^2*cos(phi))
      volume_obs=integral3(vol_sph,r_min,r_max,phi_min,phi_max,theta_min,theta_max)

      count_fov <- vector()
      for (i in 1:nfac_fh){
        count_fov[i] = sum(data$var_fhcov==fac_choose_fh[i])
      }

      dens_fov <- vector()
      for (i in 1:nfac_det){
        dens_fov[i] = ((count_fov[i]*p_hat_fh[i]*12)/(volume_obs*p_hat_cov[i]*output[["n_hours"]][i]))*10000000
      }

      output$dens_hat = dens_fov
    } else if (class(data$var_detcov)=="numeric" & class(data$var_fhcov)=="numeric"){
      p_hat_cov <- vector()

      if(p_meth=="hn"){

          p_hat_cov<- integral3(full_det_func, r_min,r_max,phi_min,phi_max,theta_min,theta_max,
                                   sigma_det = output[["sig_det_out"]], mu_fh=output[["mu_fh_out"]], sigma_fh=output[["sig_fh_out"]],
                                   reltol=.Machine$double.eps^.05)
      } else if (p_meth=="hr"){

          p_hat_cov<- integral3(full_det_func, r_min,r_max,phi_min,phi_max,theta_min,theta_max,
                                   sigma_det = output[["sig_det_out"]], beta_det = output[["beta_det_out"]], mu_fh=output[["mu_fh_out"]], sigma_fh=output[["sig_fh_out"]],
                                   reltol=.Machine$double.eps^.05)

      }



        p_hat_fh <- integral3(det_func_fh, r_min,r_max,phi_min,phi_max,theta_min,theta_max,
                                 mu_fh=output[["mu_fh_out"]], sigma_fh=output[["sig_fh_out"]],
                                 reltol=.Machine$double.eps^.05)



      vol_sph <- function(r,phi,theta) (r^2*cos(phi))
      volume_obs=integral3(vol_sph,r_min,r_max,phi_min,phi_max,theta_min,theta_max)


        count_fov = nrow(data)



        dens_fov = ((count_fov*p_hat_fh*12)/(volume_obs*p_hat_cov*output[["n_hours"]]))*10000000


      output$dens_hat = dens_fov
    }

  } else if(detcov!=FALSE & fhcov==FALSE){
    if (class(data$var_detcov)!="numeric"){
    nfac <- length(unique(fac_choose_det))

    p_hat_cov <- vector()
    if (p_meth=="hn"){
      for (i in 1:nfac){
        p_hat_cov[i]<- integral3(full_det_func, r_min,r_max,phi_min,phi_max,theta_min,theta_max,
                                 sigma_det = output[["sig_det_out"]][i], mu_fh=output[["mu_fh_out"]], sigma_fh=output[["sig_fh_out"]],
                                 reltol=.Machine$double.eps^.05)
      }
    } else if (p_meth=="hr"){
      for (i in 1:nfac){
        p_hat_cov[i]<- integral3(full_det_func, r_min,r_max,phi_min,phi_max,theta_min,theta_max,
                                 sigma_det = output[["sig_det_out"]][i], beta_det = output[["beta_det_out"]], mu_fh=output[["mu_fh_out"]], sigma_fh=output[["sig_fh_out"]],
                                 reltol=.Machine$double.eps^.05)
      }
    }



    p_hat_fh<- integral3(det_func_fh, r_min,r_max,phi_min,phi_max,theta_min,theta_max,
                         mu_fh=output[["mu_fh_out"]], sigma_fh=output[["sig_fh_out"]],
                         reltol=.Machine$double.eps^.05)



    vol_sph <- function(r,phi,theta) (r^2*cos(phi))
    volume_obs=integral3(vol_sph,r_min,r_max,phi_min,phi_max,theta_min,theta_max)

    count_fov <- vector()
    for (i in 1:nfac){
      count_fov[i] = sum(data$var_detcov==fac_choose_det[i])
    }

    dens_fov <- vector()
    for (i in 1:nfac){
      dens_fov[i] = ((count_fov[i]*p_hat_fh*12)/(volume_obs*p_hat_cov[i]*output[["n_hours"]][i]))*10000000
    }

    output$dens_hat = dens_fov
    }else{

      if (p_meth=="hn"){
          p_hat_cov <- integral3(full_det_func, r_min,r_max,phi_min,phi_max,theta_min,theta_max,
                                   sigma_det = output[["sig_det_out"]], mu_fh=output[["mu_fh_out"]], sigma_fh=output[["sig_fh_out"]],
                                   reltol=.Machine$double.eps^.05)

      } else if (p_meth=="hr"){
          p_hat_cov<- integral3(full_det_func, r_min,r_max,phi_min,phi_max,theta_min,theta_max,
                                   sigma_det = output[["sig_det_out"]], beta_det = output[["beta_det_out"]], mu_fh=output[["mu_fh_out"]], sigma_fh=output[["sig_fh_out"]],
                                   reltol=.Machine$double.eps^.05)
      }



      p_hat_fh<- integral3(det_func_fh, r_min,r_max,phi_min,phi_max,theta_min,theta_max,
                           mu_fh=output[["mu_fh_out"]], sigma_fh=output[["sig_fh_out"]],
                           reltol=.Machine$double.eps^.05)



      vol_sph <- function(r,phi,theta) (r^2*cos(phi))
      volume_obs=integral3(vol_sph,r_min,r_max,phi_min,phi_max,theta_min,theta_max)


        count_fov = nrow(data)



        dens_fov = ((count_fov*p_hat_fh*12)/(volume_obs*p_hat_cov*output[["n_hours"]]))*10000000


      output$dens_hat = dens_fov
    }
  } else if (detcov==FALSE & fhcov==FALSE) {
   if(vertdist!="norm_mix"){
     if(p_meth=="hr"){
       p_hat = integral3(full_det_func, r_min,r_max,phi_min,phi_max,theta_min,theta_max,
                         sigma_det = output[["sig_det_out"]], beta_det = output[["beta_det_out"]], mu_fh=output[["mu_fh_out"]], sigma_fh=output[["sig_fh_out"]],
                         reltol=.Machine$double.eps^.05)
     } else if (p_meth=="hn") {
       p_hat = integral3(full_det_func, r_min,r_max,phi_min,phi_max,theta_min,theta_max,
                         sigma_det = output[["sig_det_out"]], mu_fh=output[["mu_fh_out"]], sigma_fh=output[["sig_fh_out"]],
                         reltol=.Machine$double.eps^.05)
     }



     p_hat_fh <- integral3(det_func_fh, r_min,r_max,phi_min,phi_max,theta_min,theta_max,
                           mu_fh=output[["mu_fh_out"]], sigma_fh=output[["sig_fh_out"]],
                           reltol=.Machine$double.eps^.05)
   }else if (vertdist=="norm_mix"){
     if(p_meth=="hr"){
       p_hat = integral3(full_det_func, r_min,r_max,phi_min,phi_max,theta_min,theta_max,
                         sigma_det = output[["sig_det_out"]], beta_det = output[["beta_det_out"]], mu_fh=output[["mu_fh_out"]], sigma_fh=output[["sig_fh_out"]],
                         reltol=.Machine$double.eps^.05)
     } else if (p_meth=="hn") {
       p_hat = integral3(full_det_func, r_min,r_max,phi_min,phi_max,theta_min,theta_max, sigma_det = output[["sig_det_out"]],mu_fh=output[["mu_fh_out"]],sigma_fh=output[["sig_fh_out"]], w.est = output[["w.out"]],
                         mu.diff=exp(output[["model"]][["parm_raw_out"]][5]),sigma_fh_2=output[["sig_2_out"]],
                         reltol=.Machine$double.eps^.05)
     }



     p_hat_fh <- integral3(det_func_fh, r_min,r_max,phi_min,phi_max,theta_min,theta_max,
                           mu_fh=output[["mu_fh_out"]], sigma_fh=output[["sig_fh_out"]],w.est = output[["w.out"]],
                           mu.diff=exp(output[["model"]][["parm_raw_out"]][5]),sigma_fh_2=output[["sig_2_out"]],
                           reltol=.Machine$double.eps^.05)
   }

    vol_sph <- function(r,phi,theta) (r^2*cos(phi))
    volume_obs=integral3(vol_sph,r_min,r_max,phi_min,phi_max,theta_min,theta_max)

    count = nrow(data)
    n_days = output[["n_days"]]
    n_hours = output[["n_hours"]]

    dens_hat = (count*p_hat_fh*12)/(volume_obs*p_hat*n_hours)
    dens_hat=dens_hat*10000000

    output$dens_hat = dens_hat
  }

  ##############################################################################
  ##############################################################################
  ###density ci calculation
  ##############################################################################

    if (dens_boot!=FALSE){

      resample <- mvtnorm::rmvnorm(dens_boot, mean=output[["model"]][["parm_raw_out"]],sigma=solve(output[["Hess"]]))

      if (detcov!=FALSE & fhcov!=FALSE){
        nfac <- length(output[["fac_cat"]])

        dens_ci <- vector(mode='list', length=nfac)
        ci_dens <- vector(mode='list', length=nfac)

        for (i in 1:dens_boot){
          r_min <- output[["integral"]][["r_min"]]
          r_max <- output[["integral"]][["r_max"]]
          phi_min <- output[["integral"]][["phi_min"]]
          phi_max <- output[["integral"]][["phi_max"]]
          theta_min <- output[["integral"]][["theta_min"]]
          theta_max <- output[["integral"]][["theta_max"]]


          p_hat_cov_ci <- vector()
          p_hat_fh_ci <- vector()

          if (p_meth=="hn"){
            if(fhsigcov==FALSE){
              p_hat_cov_ci[1] <- integral3(full_det_func, r_min,r_max,phi_min,phi_max,theta_min,theta_max,
                                           sigma_det = resample[i,1], mu_fh=resample[i,1+nfac], sigma_fh=exp(resample[i,1+(2*nfac)]),
                                           reltol=.Machine$double.eps^.05)

              p_hat_fh_ci[1] <- integral3(det_func_fh, r_min,r_max,phi_min,phi_max,theta_min,theta_max,
                                          mu_fh=resample[i,1+nfac], sigma_fh=exp(resample[i,1+(2*nfac)]),
                                          reltol=.Machine$double.eps^.05)

              for (j in 2:nfac){
                p_hat_cov_ci[j]<- integral3(full_det_func, r_min,r_max,phi_min,phi_max,theta_min,theta_max,
                                            sigma_det = resample[i,1]+resample[i,j], mu_fh=resample[i,1+nfac]+resample[i,j+nfac], sigma_fh=exp(resample[i,1+2*nfac]),
                                            reltol=.Machine$double.eps^.05)

                p_hat_fh_ci[j]<- integral3(det_func_fh, r_min,r_max,phi_min,phi_max,theta_min,theta_max,
                                           mu_fh=resample[i,1+nfac]+resample[i,j+nfac], sigma_fh=exp(resample[i,1+2*nfac]),
                                           reltol=.Machine$double.eps^.05)
              }
            }else{
              p_hat_cov_ci[1] <- integral3(full_det_func, r_min,r_max,phi_min,phi_max,theta_min,theta_max,
                                           sigma_det = resample[i,1], mu_fh=resample[i,1+nfac], sigma_fh=resample[i,1+(2*nfac)],
                                           reltol=.Machine$double.eps^.05)

              p_hat_fh_ci[1] <- integral3(det_func_fh, r_min,r_max,phi_min,phi_max,theta_min,theta_max,
                                          mu_fh=resample[i,1+nfac], sigma_fh=resample[i,1+(2*nfac)],
                                          reltol=.Machine$double.eps^.05)

              for (j in 2:nfac){
                p_hat_cov_ci[j]<- integral3(full_det_func, r_min,r_max,phi_min,phi_max,theta_min,theta_max,
                                            sigma_det = resample[i,1]+resample[i,j], mu_fh=resample[i,1+nfac]+resample[i,j+nfac], sigma_fh=resample[i,1+(2*nfac)]+resample[i,j+2*nfac],
                                            reltol=.Machine$double.eps^.05)

                p_hat_fh_ci[j]<- integral3(det_func_fh, r_min,r_max,phi_min,phi_max,theta_min,theta_max,
                                           mu_fh=resample[i,1+nfac]+resample[i,j+nfac], sigma_fh=resample[i,1+(2*nfac)]+resample[i,j+2*nfac],
                                           reltol=.Machine$double.eps^.05)
              }
            }


          } else if (p_meth=="hr"){
            if(fhsigcov==FALSE){
              p_hat_cov_ci[1] <- integral3(full_det_func, r_min,r_max,phi_min,phi_max,theta_min,theta_max,
                                           sigma_det = resample[i,1], beta_det = resample[i,1+nfac], mu_fh=resample[i,2+nfac], sigma_fh=exp(resample[i,2+(2*nfac)]),
                                           reltol=.Machine$double.eps^.05)

              p_hat_fh_ci[1] <- integral3(det_func_fh, r_min,r_max,phi_min,phi_max,theta_min,theta_max,
                                          mu_fh=resample[i,2+nfac], sigma_fh=exp(resample[i,2+(2*nfac)]),
                                          reltol=.Machine$double.eps^.05)

              for (j in 2:nfac){
                p_hat_cov_ci[j]<- integral3(full_det_func, r_min,r_max,phi_min,phi_max,theta_min,theta_max,
                                            sigma_det = resample[i,1]+resample[i,j],beta_det=resample[i,1+nfac], mu_fh=resample[i,2+nfac]+resample[i,1+j+nfac], sigma_fh=exp(resample[i,2+2*nfac]),
                                            reltol=.Machine$double.eps^.05)

                p_hat_fh_ci[j]<- integral3(det_func_fh, r_min,r_max,phi_min,phi_max,theta_min,theta_max,
                                           mu_fh=resample[i,2+nfac]+resample[i,1+j+nfac], sigma_fh=exp(resample[i,2+2*nfac]),
                                           reltol=.Machine$double.eps^.05)
              }
            } else {
              p_hat_cov_ci[1] <- integral3(full_det_func, r_min,r_max,phi_min,phi_max,theta_min,theta_max,
                                           sigma_det = resample[i,1], beta_det = resample[i,1+nfac], mu_fh=resample[i,2+nfac], sigma_fh=resample[i,2+(2*nfac)],
                                           reltol=.Machine$double.eps^.05)

              p_hat_fh_ci[1] <- integral3(det_func_fh, r_min,r_max,phi_min,phi_max,theta_min,theta_max,
                                          mu_fh=resample[i,2+nfac], sigma_fh=resample[i,2+(2*nfac)],
                                          reltol=.Machine$double.eps^.05)

              for (j in 2:nfac){
                p_hat_cov_ci[j]<- integral3(full_det_func, r_min,r_max,phi_min,phi_max,theta_min,theta_max,
                                            sigma_det = resample[i,1]+resample[i,j],beta_det=resample[i,1+nfac], mu_fh=resample[i,2+nfac]+resample[i,1+j+nfac], sigma_fh=resample[i,2+(2*nfac)]+resample[i,2+j+2*nfac],
                                            reltol=.Machine$double.eps^.05)

                p_hat_fh_ci[j]<- integral3(det_func_fh, r_min,r_max,phi_min,phi_max,theta_min,theta_max,
                                           mu_fh=resample[i,2+nfac]+resample[i,1+j+nfac], sigma_fh=resample[i,2+(2*nfac)]+resample[i,2+j+2*nfac],
                                           reltol=.Machine$double.eps^.05)
              }
            }


          }


          vol_sph <- function(r,phi,theta) (r^2*cos(phi))
          volume_obs=integral3(vol_sph,r_min,r_max,phi_min,phi_max,theta_min,theta_max)

          count_fov <- vector()
          for (j in 1:nfac){
            count_fov[j] = sum(data$var_detcov==fac_choose_det[j])
          }


          n_hrs = output[["n_hours"]]

          for (j in 1:nfac){
            dens_ci[[j]][i] = ((count_fov[j]*p_hat_fh_ci[j]*12)/(volume_obs*p_hat_cov_ci[j]*n_hrs[j]))*10000000
          }
        }

        for (j in 1:nfac){
          ci_dens[[j]] <- c(quantile(dens_ci[[j]], prob=0.025,na.rm=T), quantile(dens_ci[[j]], prob=0.975,na.rm=T))
        }
      } else if (detcov==FALSE & fhcov==FALSE) {
        dens_ci <- vector()
        N_ci <- vector()

        for (i in 1:dens_boot){
          r_min <- output[["integral"]][["r_min"]]
          r_max <- output[["integral"]][["r_max"]]
          phi_min <- output[["integral"]][["phi_min"]]
          phi_max <- output[["integral"]][["phi_max"]]
          theta_min <- output[["integral"]][["theta_min"]]
          theta_max <- output[["integral"]][["theta_max"]]


          if (p_meth=="hn"){
            p_hat <- integral3(full_det_func, r_min,r_max,phi_min,phi_max,theta_min,theta_max,
                               sigma_det = exp(resample[i,1]), mu_fh=exp(resample[i,2]), sigma_fh=exp(resample[i,3]),
                               reltol=.Machine$double.eps^.05)



            p_hat_fh <- integral3(det_func_fh, r_min,r_max,phi_min,phi_max,theta_min,theta_max,
                                  mu_fh=exp(resample[i,2]), sigma_fh=exp(resample[i,3]),
                                  reltol=.Machine$double.eps^.05)
          } else if (p_meth=="hr"){
            p_hat <- integral3(full_det_func, r_min,r_max,phi_min,phi_max,theta_min,theta_max,
                               sigma_det = exp(resample[i,1]), beta_det=exp(resample[i,2]), mu_fh=exp(resample[i,3]), sigma_fh=exp(resample[i,4]),
                               reltol=.Machine$double.eps^.05)



            p_hat_fh <- integral3(det_func_fh, r_min,r_max,phi_min,phi_max,theta_min,theta_max,
                                  mu_fh=exp(resample[i,3]), sigma_fh=exp(resample[i,4]),
                                  reltol=.Machine$double.eps^.05)
          }


          vol_sph <- function(r,phi,theta) (r^2*cos(phi))
          volume_obs=integral3(vol_sph,r_min,r_max,phi_min,phi_max,theta_min,theta_max)

          count = nrow(data)
          n_days = output[["n_days"]]
          n_hours = output[["n_hours"]]

          dens_ci[i] = (count*p_hat_fh*12)/(volume_obs*p_hat*n_hours)
          dens_ci[i] = dens_ci[i]*10000000
        }

        ci_low <- quantile(dens_ci, prob=0.025,na.rm=T)
        ci_high <- quantile(dens_ci, prob=0.975,na.rm=T)

        ci_dens <- c(ci_low,ci_high)

      }else if (detcov!=FALSE & fhcov==FALSE){
        nfac <- length(output[["fac_cat"]])

        dens_ci <- vector(mode='list', length=nfac)
        ci_dens <- vector(mode='list', length=nfac)

        for (i in 1:dens_boot){
          r_min <- output[["integral"]][["r_min"]]
          r_max <- output[["integral"]][["r_max"]]
          phi_min <- output[["integral"]][["phi_min"]]
          phi_max <- output[["integral"]][["phi_max"]]
          theta_min <- output[["integral"]][["theta_min"]]
          theta_max <- output[["integral"]][["theta_max"]]


          p_hat_cov_ci <- vector()
          p_hat_fh_ci <- vector()

          if (p_meth=="hn"){
            p_hat_cov_ci[1] <- integral3(full_det_func, r_min,r_max,phi_min,phi_max,theta_min,theta_max,
                                         sigma_det = resample[i,1], mu_fh=resample[i,1+nfac], sigma_fh=exp(resample[i,2+nfac]),
                                         reltol=.Machine$double.eps^.05)

            p_hat_fh_ci <- integral3(det_func_fh, r_min,r_max,phi_min,phi_max,theta_min,theta_max,
                                     mu_fh=resample[i,1+nfac], sigma_fh=exp(resample[i,2+nfac]),
                                     reltol=.Machine$double.eps^.05)

            for (j in 2:nfac){
              p_hat_cov_ci[j]<- integral3(full_det_func, r_min,r_max,phi_min,phi_max,theta_min,theta_max,
                                          sigma_det = resample[i,1]+resample[i,j], mu_fh=resample[i,1+nfac], sigma_fh=exp(resample[i,2+nfac]),
                                          reltol=.Machine$double.eps^.05)

            }
          } else if (p_meth=="hr"){
            p_hat_cov_ci[1] <- integral3(full_det_func, r_min,r_max,phi_min,phi_max,theta_min,theta_max,
                                         sigma_det = resample[i,1], beta_det=resample[i,1+nfac], mu_fh=resample[i,2+nfac], sigma_fh=exp(resample[i,3+nfac]),
                                         reltol=.Machine$double.eps^.05)

            p_hat_fh_ci <- integral3(det_func_fh, r_min,r_max,phi_min,phi_max,theta_min,theta_max,
                                     mu_fh=resample[i,2+nfac], sigma_fh=exp(resample[i,3+nfac]),
                                     reltol=.Machine$double.eps^.05)

            for (j in 2:nfac){
              p_hat_cov_ci[j]<- integral3(full_det_func, r_min,r_max,phi_min,phi_max,theta_min,theta_max,
                                          sigma_det = resample[i,1]+resample[i,j], beta_det=resample[i,1+nfac], mu_fh=resample[i,2+nfac], sigma_fh=exp(resample[i,3+nfac]),
                                          reltol=.Machine$double.eps^.05)

            }
          }



          vol_sph <- function(r,phi,theta) (r^2*cos(phi))
          volume_obs=integral3(vol_sph,r_min,r_max,phi_min,phi_max,theta_min,theta_max)

          count_fov <- vector()
          for (j in 1:nfac){
            count_fov[j] = sum(data$var_detcov==fac_choose_det[j])
          }


          n_hrs = output[["n_hours"]]

          for (j in 1:nfac){
            dens_ci[[j]][i] = ((count_fov[j]*p_hat_fh_ci*12)/(volume_obs*p_hat_cov_ci[j]*n_hrs[j]))*10000000
          }
        }

        for (j in 1:nfac){
          ci_dens[[j]] <- c(quantile(dens_ci[[j]], prob=0.025,na.rm=T), quantile(dens_ci[[j]], prob=0.975,na.rm=T))
        }
      }


      output$ci_dens = ci_dens

    }

  return(output)
}



