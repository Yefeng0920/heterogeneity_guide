#-----------------------heterogeneity calculation-----------------------#
#' @title i2_ml
#' @description I2 (I-squared) for mulilevel meta-analytic models, based on Nakagawa & Santos (2012). Under multilevel models, we can have multiple I2 (see also Senior et al. 2016). Alternatively, the method proposed by Wolfgang Viechtbauer can also be used.
#' @param model Model object of class \code{rma.mv} or \code{rma}. Currently only model objects using the \code{mods} argument work (e.g., \code{mod = ~1}).
#' @param method Method used to calculate I2. Two options exist: a ratio-based calculation proposed by Nakagawa & Santos (\code{"ratio"}), or Wolfgang Viechtbauer's matrix method (\code{"matrix"}).
#' @param boot Number of simulations to run to produce 95 percent confidence intervals for I2. Default is \code{NULL}, where only the point estimate is provided.
#' @return A data frame containing all the model results including mean effect size estimate, confidence, and prediction intervals
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @examples
#' \dontrun{
#' # IMPORTANT NOTE ** boot = 10 is set LOW deliberately to make the models run fast. You should always run for at least boot = 1000
#' # English example
#' data(english)
#' english <- escalc(measure = "SMD", n1i = NStartControl,
#' sd1i = SD_C, m1i = MeanC, n2i = NStartExpt, sd2i = SD_E,
#' m2i = MeanE, var.names=c("SMD","vSMD"),data = english)
#' english_MA <- rma.mv(yi = SMD, V = vSMD,
#' random = list( ~ 1 | StudyNo, ~ 1 | EffectID), data = english)
#' I2_eng_1 <- i2_ml(english_MA, data = english, boot = 10)
#' I2_eng_2 <- i2_ml(english_MA, data = english, method = "ratio")
#' I2_eng_3 <- i2_ml(english_MA, data = english, method = "matrix")
#'
#' ## Fish example
#' data(fish)
#' warm_dat <- fish
#' model <- metafor::rma.mv(yi = lnrr, V = lnrr_vi,
#' random = list(~1 | group_ID, ~1 | es_ID),
#' mods = ~ experimental_design + trait.type + deg_dif + treat_end_days,
#' method = "REML", test = "t", data = warm_dat,
#' control=list(optimizer="optim", optmethod="Nelder-Mead"))
#' I2_fish_1 <- i2_ml(model, data = warm_dat, boot = 10)
#' I2_fish_2 <- i2_ml(model, method = c("matrix"),data = warm_dat)
#' I2_fish_2 <- i2_ml(model, method = c("ratio"),data = warm_dat)
#'
#' # Lim example
#' data(lim)
#' # Add in the sampling variance
#' lim$vi<-(1/sqrt(lim$N - 3))^2
#' # Lets fit a meta-regression - I will do Article non-independence.
#' The phylogenetic model found phylogenetic effects, however, instead we could fit Phylum as a fixed effect and explore them with an Orchard Plot
#' lim_MR<-metafor::rma.mv(yi=yi, V=vi, mods=~Phylum-1, random=list(~1|Article, ~1|Datapoint), data=lim)
#' I2_lim_1 <- i2_ml(lim_MR, data=lim, boot = 10)
#' I2_lim_2 <- i2_ml(lim_MR, data=lim)
#' }
#' @references Senior, A. M., Grueber, C. E., Kamiya, T., Lagisz, M., O’Dwyer, K., Santos, E. S. A. & Nakagawa S. 2016. Heterogeneity in ecological and evolutionary meta-analyses: its magnitudes and implications. Ecology 97(12): 3293-3299.
#'  Nakagawa, S, and Santos, E.S.A. 2012. Methodological issues and advances in biological meta-analysis.Evolutionary Ecology 26(5): 1253-1274.
#' @export

i2_ml <- function(model,
                  method = c("ratio", "matrix"),
                  boot = NULL) {
  
  if(all(class(model) %in% c("robust.rma", "rma.mv", "rma", "rma.uni")) == FALSE) {stop("Sorry, you need to fit a metafor model of class robust.rma, rma.mv, rma, rma.uni")}
  
  if(any(model$tau2 > 0)) { stop("Sorry. At the moment i2_ml cannot take models with heterogeneous variance.")}
  
  ## evaluate choices
  method <- match.arg(method)
  
  if (method == "matrix") {
    # Wolfgang Viechtbauer's method
    I2s <- matrix_i2(model)
  } else {
    # Nakagawa & Santos (2012)
    I2s <- ratio_i2(model)
  }
  
  # Extract the data from the model object
  data <- model$data
  
  # Check if missing values exist and use complete case data
  if(any(model$not.na == FALSE)){
    data <- data[model$not.na,]
  }
  
  if(!is.null(boot)){
    # Simulate the vector of effect sizes
    sim <- metafor::simulate.rma(model, nsim=boot) # Add try catch here? DN
    
    # Get formula from model object.
    random_formula <- model$random
    mods_formula <- metafor::formula.rma(model, type = "mods") #in case moderators
    vi <- model$vi
    
    # Parametric bootstrap
    pb <- progress::progress_bar$new(total = boot,
                                     format = "Bootstrapping [:bar] :percent ETA: :eta",
                                     show_after = 0)
    if(is.null(mods_formula)){
      I2_each <- sapply(sim, function(ysim) {
        tmp <- tryCatch(metafor::rma.mv(ysim, vi,
                                        random = random_formula, data = data))
        pb$tick()
        Sys.sleep(1/boot)
        if (method == "matrix") {
          I2 <- matrix_i2(tmp)
        }
        else {
          I2 <- ratio_i2(tmp)
        }
        return(I2)
      })
    } else{
      I2_each <- sapply(sim, function(ysim) {
        
        # The model
        tmp <- tryCatch(metafor::rma.mv( ysim, vi,
                                         mods = mods_formula,
                                         random = random_formula,
                                         data = data))
        pb$tick()
        Sys.sleep(1 / boot)
        
        if(method == "matrix"){
          I2 <- matrix_i2(tmp)
        } else {
          I2 <- ratio_i2(tmp)
        }
        
        return(I2) })
    }
    # Summarise the bootstrapped distribution.
    I2s_each_95 <- data.frame(t(apply(I2_each, 1, stats::quantile, probs=c(0.5, .025, .975))))
    I2s <-  round(I2s_each_95, digits = 3)
    colnames(I2s) = c("Est.", "2.5%", "97.5%")
  }
  
  return(I2s)
}

#' @title matrix_i2
#' @description Calculated I2 (I-squared) for mulilevel meta-analytic models, based on a matrix method proposed by Wolfgang Viechtbauer.
#' @param model Model object of class \code{rma.mv} or \code{rma}.
#' @examples
#' \dontrun{
#' # English example
#' data(english)
#' english <- escalc(measure = "SMD", n1i = NStartControl, sd1i = SD_C, m1i = MeanC, n2i = NStartExpt, sd2i = SD_E, m2i = MeanE, var.names=c("SMD","vSMD"),data = english)
#' english_MA <- rma.mv(yi = SMD, V = vSMD, random = list( ~ 1 | StudyNo, ~ 1 | EffectID), data = english)
#' I2_eng <- i2_ml(english_MA, data = english, method = "matrix")
#'
#' # Lim example
#' data(lim)
#' # Add in the sampling variance
#' lim$vi<-(1/sqrt(lim$N - 3))^2
#' # Lets fit a meta-regression - I will do Article non-independence. The phylogenetic model found phylogenetic effects, however, instead we could fit Phylum as a fixed effect and explore them with an Orchard Plot
#' lim_MR<-metafor::rma.mv(yi=yi, V=vi, mods=~Phylum-1,
#' random=list(~1|Article, ~1|Datapoint), data=lim)
#' I2_lim <- i2_ml(lim_MR, data=lim, method = "matrix")
#' }
#' @export
matrix_i2 <- function(model){
  
  if(all(class(model) %in% c("robust.rma", "rma.mv", "rma", "rma.uni")) == FALSE) {stop("Sorry, you need to fit a metafor model of class robust.rma, rma.mv, rma, rma.uni")}
  
  W <- solve(model$V)
  X <- model.matrix(model)
  P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
  I2_total <- 100* (sum(model$sigma2) / (sum(model$sigma2) + (model$k - model$p) / sum(diag(P))))
  I2_each <- 100* (model$sigma2 / (sum(model$sigma2) + (model$k - model$p) / sum(diag(P))))
  names(I2_each) <- paste0("I2_", model$s.names)
  names(I2_total) <- "I2_Total"
  I2s <- c(I2_total, I2_each)
  return(I2s)
}


#' @title ratio_i2
#' @description I2 (I-squared) for mulilevel meta-analytic models based on Nakagawa & Santos (2012). Under multilevel models, we can have a multiple I2 (see also Senior et al. 2016).
#' @param model Model object of class \code{rma.mv} or \code{rma}.
#' @examples
#' \dontrun{
#' # English example
#' data(english)
#' english <- escalc(measure = "SMD", n1i = NStartControl,
#' sd1i = SD_C, m1i = MeanC, n2i = NStartExpt,
#' sd2i = SD_E, m2i = MeanE, var.names=c("SMD","vSMD"),data = english)
#' english_MA <- rma.mv(yi = SMD, V = vSMD,
#' random = list( ~ 1 | StudyNo, ~ 1 | EffectID), data = english)
#' I2_eng_1 <- i2_ml(english_MA, data = english, boot = 1000)
#' I2_eng_2 <- i2_ml(english_MA, data = english, method = "ratio")
#' }
#' @export
ratio_i2 <- function(model){
  
  if(all(class(model) %in% c("robust.rma", "rma.mv", "rma", "rma.uni")) == FALSE) {stop("Sorry, you need to fit a metafor model of class robust.rma, rma.mv, rma, rma.uni")}
  
  # sigma2_v = typical sampling error variance
  sigma2_v <- sum(1 / model$vi) * (model$k - 1) /
    (sum(1 / model$vi)^2 - sum((1 / model$vi)^2))
  
  # s^2_t = total variance
  I2_total <- 100 * (sum(model$sigma2) / (sum(model$sigma2) + sigma2_v))
  I2_each <- 100 * (model$sigma2 / (sum(model$sigma2) + sigma2_v))
  names(I2_each) <- paste0("I2_", model$s.names)
  names(I2_total) <- "I2_Total"
  
  I2s <- c(I2_total, I2_each)
  return(I2s)
}


#' @title cvh1_ml
#' @description CVH1 for mulilevel meta-analytic models, based on Yang et al. (2023). Under multilevel models, we can have multiple CVH1. TODO - we need to cite original CVH1 paper
#' @param model Model object of class \code{rma.mv} or \code{rma}. Currently only model objects using the \code{mods} argument work (e.g., \code{mod = ~1}).
#' @param boot Number of simulations to run to produce 95 percent confidence intervals for I2. Default is \code{NULL}, where only the point estimate is provided.
#' @return A data frame containing all the model results including mean effect size estimate, confidence, and prediction intervals
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @examples
#' \dontrun{
#' # IMPORTANT NOTE ** boot = 10 is set LOW deliberately to make the models run fast. You should always run for at least boot = 1000
#' # English example
#' data(english)
#' english <- escalc(measure = "SMD", n1i = NStartControl,
#' sd1i = SD_C, m1i = MeanC, n2i = NStartExpt, sd2i = SD_E,
#' m2i = MeanE, var.names=c("SMD","vSMD"),data = english)
#' english_MA <- rma.mv(yi = SMD, V = vSMD,
#' random = list( ~ 1 | StudyNo, ~ 1 | EffectID), data = english)
#' CVH1_eng_1 <- cvh1_ml(english_MA, boot = 10)
#' CVH1_eng_2 <- cvh1_ml(english_MA)
#'
#' ## Fish example
#' data(fish)
#' warm_dat <- fish
#' model <- metafor::rma.mv(yi = lnrr, V = lnrr_vi,
#' random = list(~1 | group_ID, ~1 | es_ID),
#' mods = ~ experimental_design + trait.type + deg_dif + treat_end_days,
#' method = "REML", test = "t", data = warm_dat,
#' control=list(optimizer="optim", optmethod="Nelder-Mead"))
#' CVH1_fish_1 <- cvh1_ml(model, boot = 10)
#' CVH1_fish_2 <- cvh1_ml(model)
#'
#' # Lim example
#' data(lim)
#' # Add in the sampling variance
#' lim$vi<-(1/sqrt(lim$N - 3))^2
#' # Lets fit a meta-regression - I will do Article non-independence.
#' The phylogenetic model found phylogenetic effects, however, instead we could fit Phylum as a fixed effect and explore them with an Orchard Plot
#' lim_MR<-metafor::rma.mv(yi=yi, V=vi, mods=~Phylum-1, random=list(~1|Article, ~1|Datapoint), data=lim)
#' CVH1_lim_1 <- cvh1_ml(lim_MR, boot = 10)
#' CVH1_lim_2 <- cvh1_ml(lim_MR)
#' }
#' @references TODO
#' @export

cvh1_ml <- function(model,
                    boot = NULL) {
  
  if(all(class(model) %in% c("robust.rma", "rma.mv", "rma", "rma.uni")) == FALSE) {stop("Sorry, you need to fit a metafor model of class robust.rma, rma.mv, rma, rma.uni")}
  
  if(any(model$tau2 > 0)) { stop("Sorry. At the moment cvh1_ml cannot take models with heterogeneous variance.")}
  
  CVH1s <- ml_cvh1(model)
  
  # Extract the data from the model object
  data <- model$data
  
  # Check if missing values exist and use complete case data
  if(any(model$not.na == FALSE)){
    data <- data[model$not.na,]
  }
  
  if(!is.null(boot)){
    # Simulate the vector of effect sizes
    sim <- metafor::simulate.rma(model, nsim=boot) # Add try catch here? DN
    
    # Get formula from model object.
    random_formula <- model$random
    mods_formula <- metafor::formula.rma(model, type = "mods") #in case moderators
    vi <- model$vi
    
    # Parametric bootstrap
    pb <- progress::progress_bar$new(total = boot,
                                     format = "Bootstrapping [:bar] :percent ETA: :eta",
                                     show_after = 0)
    if(is.null(mods_formula)){
      CVH1_each <- sapply(sim, function(ysim) {
        tmp <- tryCatch(metafor::rma.mv(ysim, vi,
                                        random = random_formula, data = data))
        pb$tick()
        Sys.sleep(1/boot)
        CVH1 <- ml_cvh1(tmp)
      })
    } else{
      CVH1_each <- sapply(sim, function(ysim) {
        
        # The model
        tmp <- tryCatch(metafor::rma.mv( ysim, vi,
                                         mods = mods_formula,
                                         random = random_formula,
                                         data = data))
        pb$tick()
        Sys.sleep(1 / boot)
        CVH1 <- ml_cvh1(tmp)
        return(CVH1) })
    }
    # Summarise the bootstrapped distribution.
    CVH1s_each_95 <- data.frame(t(apply(CVH1_each, 1, stats::quantile, probs=c(0.5, .025, .975))))
    CVH1s <-  round(CVH1s_each_95, digits = 3)
    colnames(CVH1s) = c("Est.", "2.5%", "97.5%")
  }
  
  return(CVH1s)
}


#' @title ml_cvh1
#' @description Calculated CVH1 for mulilevel meta-analytic models
#' @param model Model object of class \code{rma.mv} or \code{rma}.
#' @export

ml_cvh1 <- function(model){
  
  # total cvh1
  CVH1_total <- sqrt(sum(model$sigma2)) / abs(model$beta[[1]])
  # cvh1 at different levels
  CVH1_each <-  sqrt(model$sigma2) / abs(model$beta[[1]])
  names(CVH1_each) <- paste0("CVH1_", model$s.names)
  names(CVH1_total) <- "CVH1_Total"
  
  CVH1s <- c(CVH1_total, CVH1_each)
  
  return(CVH1s)
}


#' @title cvh2_ml
#' @description CVH2 for mulilevel meta-analytic models, based on Yang et al. (2023). Under multilevel models, we can have multiple CVH2. TODO - we need to cite original CVH2 paper
#' @param model Model object of class \code{rma.mv} or \code{rma}. Currently only model objects using the \code{mods} argument work (e.g., \code{mod = ~1}).
#' @param boot Number of simulations to run to produce 95 percent confidence intervals for I2. Default is \code{NULL}, where only the point estimate is provided.
#' @return A data frame containing all the model results including mean effect size estimate, confidence, and prediction intervals
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @examples
#' \dontrun{
#' # IMPORTANT NOTE ** boot = 10 is set LOW deliberately to make the models run fast. You should always run for at least boot = 1000
#' # English example
#' data(english)
#' english <- escalc(measure = "SMD", n1i = NStartControl,
#' sd1i = SD_C, m1i = MeanC, n2i = NStartExpt, sd2i = SD_E,
#' m2i = MeanE, var.names=c("SMD","vSMD"),data = english)
#' english_MA <- rma.mv(yi = SMD, V = vSMD,
#' random = list( ~ 1 | StudyNo, ~ 1 | EffectID), data = english)
#' CVH2_eng_1 <- cvh2_ml(english_MA, boot = 10)
#' CVH2_eng_2 <- cvh2_ml(english_MA)
#'
#' ## Fish example
#' data(fish)
#' warm_dat <- fish
#' model <- metafor::rma.mv(yi = lnrr, V = lnrr_vi,
#' random = list(~1 | group_ID, ~1 | es_ID),
#' mods = ~ experimental_design + trait.type + deg_dif + treat_end_days,
#' method = "REML", test = "t", data = warm_dat,
#' control=list(optimizer="optim", optmethod="Nelder-Mead"))
#' CVH2_fish_1 <- cvh2_ml(model, boot = 10)
#' CVH2_fish_2 <- cvh2_ml(model)
#'
#' # Lim example
#' data(lim)
#' # Add in the sampling variance
#' lim$vi<-(1/sqrt(lim$N - 3))^2
#' # Lets fit a meta-regression - I will do Article non-independence.
#' The phylogenetic model found phylogenetic effects, however, instead we could fit Phylum as a fixed effect and explore them with an Orchard Plot
#' lim_MR<-metafor::rma.mv(yi=yi, V=vi, mods=~Phylum-1, random=list(~1|Article, ~1|Datapoint), data=lim)
#' CVH2_lim_1 <- cvh2_ml(lim_MR, boot = 10)
#' CVH2_lim_2 <- cvh2_ml(lim_MR)
#' }
#' @references TODO
#' @export

cvh2_ml <- function(model,
                    boot = NULL) {
  
  if(all(class(model) %in% c("robust.rma", "rma.mv", "rma", "rma.uni")) == FALSE) {stop("Sorry, you need to fit a metafor model of class robust.rma, rma.mv, rma, rma.uni")}
  
  if(any(model$tau2 > 0)) { stop("Sorry. At the moment cvh2_ml cannot take models with heterogeneous variance.")}
  
  CVH2s <- ml_cvh2(model)
  
  # Extract the data from the model object
  data <- model$data
  
  # Check if missing values exist and use complete case data
  if(any(model$not.na == FALSE)){
    data <- data[model$not.na,]
  }
  
  if(!is.null(boot)){
    # Simulate the vector of effect sizes
    sim <- metafor::simulate.rma(model, nsim=boot) # Add try catch here? DN
    
    # Get formula from model object.
    random_formula <- model$random
    mods_formula <- metafor::formula.rma(model, type = "mods") #in case moderators
    vi <- model$vi
    
    # Parametric bootstrap
    pb <- progress::progress_bar$new(total = boot,
                                     format = "Bootstrapping [:bar] :percent ETA: :eta",
                                     show_after = 0)
    if(is.null(mods_formula)){
      CVH2_each <- sapply(sim, function(ysim) {
        tmp <- tryCatch(metafor::rma.mv(ysim, vi,
                                        random = random_formula, data = data))
        pb$tick()
        Sys.sleep(1/boot)
        CVH2 <- ml_cvh2(tmp)
      })
    } else{
      CVH2_each <- sapply(sim, function(ysim) {
        
        # The model
        tmp <- tryCatch(metafor::rma.mv( ysim, vi,
                                         mods = mods_formula,
                                         random = random_formula,
                                         data = data))
        pb$tick()
        Sys.sleep(1 / boot)
        CVH2 <- ml_cvh2(tmp)
        return(CVH2) })
    }
    # Summarise the bootstrapped distribution.
    CVH2s_each_95 <- data.frame(t(apply(CVH2_each, 1, stats::quantile, probs=c(0.5, .025, .975))))
    CVH2s <-  round(CVH2s_each_95, digits = 3)
    colnames(CVH2s) = c("Est.", "2.5%", "97.5%")
  }
  
  return(CVH2s)
}


#' @title ml_cvh2
#' @description Calculated CVH2 for mulilevel meta-analytic models
#' @param model Model object of class \code{rma.mv} or \code{rma}.
#' @export

ml_cvh2 <- function(model){
  
  # total cvh2
  CVH2_total <- sum(model$sigma2) / (model$beta[[1]])^2
  # cvh2 at different levels
  CVH2_each <-  model$sigma2 / (model$beta[[1]])^2
  names(CVH2_each) <- paste0("CVH2_", model$s.names)
  names(CVH2_total) <- "CVH2_Total"
  
  CVH2s <- c(CVH2_total, CVH2_each)
  
  return(CVH2s)
}


#' @title m1_ml
#' @description M1 for mulilevel meta-analytic models, based on Yang et al. (2023). Under multilevel models, we can have multiple M1 - TODO - we need to cite original M1 paper
#' @param model Model object of class \code{rma.mv} or \code{rma}. Currently only model objects using the \code{mods} argument work (e.g., \code{mod = ~1}).
#' @param boot Number of simulations to run to produce 95 percent confidence intervals for M1. Default is \code{NULL}, where only the point estimate is provided.
#' @return A data frame containing all the model results including mean effect size estimate, confidence, and prediction intervals
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @examples
#' \dontrun{
#' # IMPORTANT NOTE ** boot = 10 is set LOW deliberately to make the models run fast. You should always run for at least boot = 1000
#' # English example
#' data(english)
#' english <- escalc(measure = "SMD", n1i = NStartControl,
#' sd1i = SD_C, m1i = MeanC, n2i = NStartExpt, sd2i = SD_E,
#' m2i = MeanE, var.names=c("SMD","vSMD"),data = english)
#' english_MA <- rma.mv(yi = SMD, V = vSMD,
#' random = list( ~ 1 | StudyNo, ~ 1 | EffectID), data = english)
#' M1_eng_1 <- m1_ml(english_MA, boot = 10)
#' M1_eng_2 <- m1_ml(english_MA)
#'
#' ## Fish example
#' data(fish)
#' warm_dat <- fish
#' model <- metafor::rma.mv(yi = lnrr, V = lnrr_vi,
#' random = list(~1 | group_ID, ~1 | es_ID),
#' mods = ~ experimental_design + trait.type + deg_dif + treat_end_days,
#' method = "REML", test = "t", data = warm_dat,
#' control=list(optimizer="optim", optmethod="Nelder-Mead"))
#' M1_fish_1 <- m1_ml(model, boot = 10)
#' M1_fish_2 <- m1_ml(model)
#'
#' # Lim example
#' data(lim)
#' # Add in the sampling variance
#' lim$vi<-(1/sqrt(lim$N - 3))^2
#' # Lets fit a meta-regression - I will do Article non-independence.
#' The phylogenetic model found phylogenetic effects, however, instead we could fit Phylum as a fixed effect and explore them with an Orchard Plot
#' lim_MR<-metafor::rma.mv(yi=yi, V=vi, mods=~Phylum-1, random=list(~1|Article, ~1|Datapoint), data=lim)
#' M1_lim_1 <- m1_ml(lim_MR, boot = 10)
#' M1_lim_2 <- m1_ml(lim_MR)
#' }
#' @references TODO
#' @export

m1_ml <- function(model,
                  boot = NULL) {
  
  if(all(class(model) %in% c("robust.rma", "rma.mv", "rma", "rma.uni")) == FALSE) {stop("Sorry, you need to fit a metafor model of class robust.rma, rma.mv, rma, rma.uni")}
  
  if(any(model$tau2 > 0)) { stop("Sorry. At the moment m1_ml cannot take models with heterogeneous variance.")}
  
  M1s <- ml_m1(model)
  
  # Extract the data from the model object
  data <- model$data
  
  # Check if missing values exist and use complete case data
  if(any(model$not.na == FALSE)){
    data <- data[model$not.na,]
  }
  
  if(!is.null(boot)){
    # Simulate the vector of effect sizes
    sim <- metafor::simulate.rma(model, nsim=boot) # Add try catch here? DN
    
    # Get formula from model object.
    random_formula <- model$random
    mods_formula <- metafor::formula.rma(model, type = "mods") #in case moderators
    vi <- model$vi
    
    # Parametric bootstrap
    pb <- progress::progress_bar$new(total = boot,
                                     format = "Bootstrapping [:bar] :percent ETA: :eta",
                                     show_after = 0)
    if(is.null(mods_formula)){
      M1_each <- sapply(sim, function(ysim) {
        tmp <- tryCatch(metafor::rma.mv(ysim, vi,
                                        random = random_formula, data = data))
        pb$tick()
        Sys.sleep(1/boot)
        M1 <- ml_m1(tmp)
      })
    } else{
      M1_each <- sapply(sim, function(ysim) {
        
        # The model
        tmp <- tryCatch(metafor::rma.mv( ysim, vi,
                                         mods = mods_formula,
                                         random = random_formula,
                                         data = data))
        pb$tick()
        Sys.sleep(1 / boot)
        M1 <- ml_m1(tmp)
        return(M1) })
    }
    # Summarise the bootstrapped distribution.
    M1s_each_95 <- data.frame(t(apply(M1_each, 1, stats::quantile, probs=c(0.5, .025, .975))))
    M1s <-  round(M1s_each_95, digits = 3)
    colnames(M1s) = c("Est.", "2.5%", "97.5%")
  }
  
  return(M1s)
}


#' @title ml_m1
#' @description Calculated CV for mulilevel meta-analytic models
#' @param model Model object of class \code{rma.mv} or \code{rma}.
#' @export

ml_m1 <- function(model){
  
  # total m1
  M1_total <- sqrt(sum(model$sigma2)) / (abs(model$beta[[1]]) + sqrt(sum(model$sigma2)))
  # m1 at different levels
  M1_each <-  sqrt(model$sigma2) / (abs(model$beta[[1]]) + sqrt(sum(model$sigma2)))
  names(M1_each) <- paste0("M1_", model$s.names)
  names(M1_total) <- "M1_Total"
  
  M1s <- c(M1_total, M1_each)
  
  return(M1s)
}





#' @title ml_cv
#' @description Calculated CV for mulilevel meta-analytic models
#' @param model Model object of class \code{rma.mv} or \code{rma}.
#' @export

ml_cv <- function(model){
  
  # total cv
  CV_total <- sqrt(sum(model$sigma2)) / abs(model$beta[[1]])
  # cv at different levels
  CV_each <-  sqrt(model$sigma2) / abs(model$beta[[1]])
  names(CV_each) <- paste0("CV_", model$s.names)
  names(CV_total) <- "CV_Total"
  
  CVs <- c(CV_total, CV_each)
  
  return(CVs)
}


#' @title m1_ml
#' @description M1 for mulilevel meta-analytic models, based on Yang et al. (2023). Under multilevel models, we can have multiple M1 - TODO - we need to cite original M1 paper
#' @param model Model object of class \code{rma.mv} or \code{rma}. Currently only model objects using the \code{mods} argument work (e.g., \code{mod = ~1}).
#' @param boot Number of simulations to run to produce 95 percent confidence intervals for M1. Default is \code{NULL}, where only the point estimate is provided.
#' @return A data frame containing all the model results including mean effect size estimate, confidence, and prediction intervals
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @examples
#' \dontrun{
#' # IMPORTANT NOTE ** boot = 10 is set LOW deliberately to make the models run fast. You should always run for at least boot = 1000
#' # English example
#' data(english)
#' english <- escalc(measure = "SMD", n1i = NStartControl,
#' sd1i = SD_C, m1i = MeanC, n2i = NStartExpt, sd2i = SD_E,
#' m2i = MeanE, var.names=c("SMD","vSMD"),data = english)
#' english_MA <- rma.mv(yi = SMD, V = vSMD,
#' random = list( ~ 1 | StudyNo, ~ 1 | EffectID), data = english)
#' M1_eng_1 <- m1_ml(english_MA, boot = 10)
#' M1_eng_2 <- m1_ml(english_MA)
#'
#' ## Fish example
#' data(fish)
#' warm_dat <- fish
#' model <- metafor::rma.mv(yi = lnrr, V = lnrr_vi,
#' random = list(~1 | group_ID, ~1 | es_ID),
#' mods = ~ experimental_design + trait.type + deg_dif + treat_end_days,
#' method = "REML", test = "t", data = warm_dat,
#' control=list(optimizer="optim", optmethod="Nelder-Mead"))
#' M1_fish_1 <- m1_ml(model, boot = 10)
#' M1_fish_2 <- m1_ml(model)
#'
#' # Lim example
#' data(lim)
#' # Add in the sampling variance
#' lim$vi<-(1/sqrt(lim$N - 3))^2
#' # Lets fit a meta-regression - I will do Article non-independence.
#' The phylogenetic model found phylogenetic effects, however, instead we could fit Phylum as a fixed effect and explore them with an Orchard Plot
#' lim_MR<-metafor::rma.mv(yi=yi, V=vi, mods=~Phylum-1, random=list(~1|Article, ~1|Datapoint), data=lim)
#' M1_lim_1 <- m1_ml(lim_MR, boot = 10)
#' M1_lim_2 <- m1_ml(lim_MR)
#' }
#' @references TODO
#' @export

m1_ml <- function(model,
                  boot = NULL) {
  
  if(all(class(model) %in% c("robust.rma", "rma.mv", "rma", "rma.uni")) == FALSE) {stop("Sorry, you need to fit a metafor model of class robust.rma, rma.mv, rma, rma.uni")}
  
  if(any(model$tau2 > 0)) { stop("Sorry. At the moment m1_ml cannot take models with heterogeneous variance.")}
  
  M1s <- ml_m1(model)
  
  # Extract the data from the model object
  data <- model$data
  
  # Check if missing values exist and use complete case data
  if(any(model$not.na == FALSE)){
    data <- data[model$not.na,]
  }
  
  if(!is.null(boot)){
    # Simulate the vector of effect sizes
    sim <- metafor::simulate.rma(model, nsim=boot) # Add try catch here? DN
    
    # Get formula from model object.
    random_formula <- model$random
    mods_formula <- metafor::formula.rma(model, type = "mods") #in case moderators
    vi <- model$vi
    
    # Parametric bootstrap
    pb <- progress::progress_bar$new(total = boot,
                                     format = "Bootstrapping [:bar] :percent ETA: :eta",
                                     show_after = 0)
    if(is.null(mods_formula)){
      M1_each <- sapply(sim, function(ysim) {
        tmp <- tryCatch(metafor::rma.mv(ysim, vi,
                                        random = random_formula, data = data))
        pb$tick()
        Sys.sleep(1/boot)
        M1 <- ml_m1(tmp)
      })
    } else{
      M1_each <- sapply(sim, function(ysim) {
        
        # The model
        tmp <- tryCatch(metafor::rma.mv( ysim, vi,
                                         mods = mods_formula,
                                         random = random_formula,
                                         data = data))
        pb$tick()
        Sys.sleep(1 / boot)
        M1 <- ml_m1(tmp)
        return(M1) })
    }
    # Summarise the bootstrapped distribution.
    M1s_each_95 <- data.frame(t(apply(M1_each, 1, stats::quantile, probs=c(0.5, .025, .975))))
    M1s <-  round(M1s_each_95, digits = 3)
    colnames(M1s) = c("Est.", "2.5%", "97.5%")
  }
  
  return(M1s)
}


#' @title ml_m1
#' @description Calculated m1_m1 for mulilevel meta-analytic models
#' @param model Model object of class \code{rma.mv} or \code{rma}.
#' @export

ml_m1 <- function(model){
  
  # total m1
  M1_total <- sqrt(sum(model$sigma2)) / (abs(model$beta[[1]]) + sum(unlist(lapply(model$sigma2, sqrt))))
  # m1 at different levels
  M1_each <-  sqrt(model$sigma2) / (abs(model$beta[[1]]) + sum(unlist(lapply(model$sigma2, sqrt))))
  names(M1_each) <- paste0("M1_", model$s.names)
  names(M1_total) <- "M1_Total"
  
  M1s <- c(M1_total, M1_each)
  
  return(M1s)
}


#' @title m2_ml
#' @description M2 for mulilevel meta-analytic models, based on Yang et al. (2023). Under multilevel models, we can have multiple M2 - TODO - we need to cite original M2 paper
#' @param model Model object of class \code{rma.mv} or \code{rma}. Currently only model objects using the \code{mods} argument work (e.g., \code{mod = ~1}).
#' @param boot Number of simulations to run to produce 95 percent confidence intervals for M2. Default is \code{NULL}, where only the point estimate is provided.
#' @return A data frame containing all the model results including mean effect size estimate, confidence, and prediction intervals
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @examples
#' \dontrun{
#' # IMPORTANT NOTE ** boot = 10 is set LOW deliberately to make the models run fast. You should always run for at least boot = 1000
#' # English example
#' data(english)
#' english <- escalc(measure = "SMD", n1i = NStartControl,
#' sd1i = SD_C, m1i = MeanC, n2i = NStartExpt, sd2i = SD_E,
#' m2i = MeanE, var.names=c("SMD","vSMD"),data = english)
#' english_MA <- rma.mv(yi = SMD, V = vSMD,
#' random = list( ~ 1 | StudyNo, ~ 1 | EffectID), data = english)
#' M2_eng_1 <- m2_ml(english_MA, boot = 10)
#' M2_eng_2 <- m2_ml(english_MA)
#'
#' ## Fish example
#' data(fish)
#' warm_dat <- fish
#' model <- metafor::rma.mv(yi = lnrr, V = lnrr_vi,
#' random = list(~1 | group_ID, ~1 | es_ID),
#' mods = ~ experimental_design + trait.type + deg_dif + treat_end_days,
#' method = "REML", test = "t", data = warm_dat,
#' control=list(optimizer="optim", optmethod="Nelder-Mead"))
#' M2_fish_1 <- m2_ml(model, boot = 10)
#' M2_fish_2 <- m2_ml(model)
#'
#' # Lim example
#' data(lim)
#' # Add in the sampling variance
#' lim$vi<-(1/sqrt(lim$N - 3))^2
#' # Lets fit a meta-regression - I will do Article non-independence.
#' The phylogenetic model found phylogenetic effects, however, instead we could fit Phylum as a fixed effect and explore them with an Orchard Plot
#' lim_MR<-metafor::rma.mv(yi=yi, V=vi, mods=~Phylum-1, random=list(~1|Article, ~1|Datapoint), data=lim)
#' M2_lim_1 <- m2_ml(lim_MR, boot = 10)
#' M2_lim_2 <- m2_ml(lim_MR)
#' }
#' @references TODO
#' @export

m2_ml <- function(model,
                  boot = NULL) {
  
  if(all(class(model) %in% c("robust.rma", "rma.mv", "rma", "rma.uni")) == FALSE) {stop("Sorry, you need to fit a metafor model of class robust.rma, rma.mv, rma, rma.uni")}
  
  if(any(model$tau2 > 0)) { stop("Sorry. At the moment m2_ml cannot take models with heterogeneous variance.")}
  
  M2s <- ml_m2(model)
  
  # Extract the data from the model object
  data <- model$data
  
  # Check if missing values exist and use complete case data
  if(any(model$not.na == FALSE)){
    data <- data[model$not.na,]
  }
  
  if(!is.null(boot)){
    # Simulate the vector of effect sizes
    sim <- metafor::simulate.rma(model, nsim=boot) # Add try catch here? DN
    
    # Get formula from model object.
    random_formula <- model$random
    mods_formula <- metafor::formula.rma(model, type = "mods") #in case moderators
    vi <- model$vi
    
    # Parametric bootstrap
    pb <- progress::progress_bar$new(total = boot,
                                     format = "Bootstrapping [:bar] :percent ETA: :eta",
                                     show_after = 0)
    if(is.null(mods_formula)){
      M2_each <- sapply(sim, function(ysim) {
        tmp <- tryCatch(metafor::rma.mv(ysim, vi,
                                        random = random_formula, data = data))
        pb$tick()
        Sys.sleep(1/boot)
        M2 <- ml_m2(tmp)
      })
    } else{
      M2_each <- sapply(sim, function(ysim) {
        
        # The model
        tmp <- tryCatch(metafor::rma.mv( ysim, vi,
                                         mods = mods_formula,
                                         random = random_formula,
                                         data = data))
        pb$tick()
        Sys.sleep(1 / boot)
        M2 <- ml_m2(tmp)
        return(M2) })
    }
    # Summarise the bootstrapped distribution.
    M2s_each_95 <- data.frame(t(apply(M2_each, 1, stats::quantile, probs=c(0.5, .025, .975))))
    M2s <-  round(M2s_each_95, digits = 3)
    colnames(M2s) = c("Est.", "2.5%", "97.5%")
  }
  
  return(M2s)
}


#' @title ml_m2
#' @description Calculated m1_m2 for mulilevel meta-analytic models
#' @param model Model object of class \code{rma.mv} or \code{rma}.
#' @export

ml_m2 <- function(model){
  
  # total m2
  M2_total <- sum(model$sigma2) / ( (model$beta[[1]])^2 + sum(model$sigma2) )
  # m2 at different levels
  M2_each <-  model$sigma2 / ( (model$beta[[1]])^2 + sum(model$sigma2) )
  names(M2_each) <- paste0("M2_", model$s.names)
  names(M2_total) <- "M2_Total"
  
  M2s <- c(M2_total, M2_each)
  
  return(M2s)
}


#-----------------------heterogeneity interpretation-----------------------#
het_interpret <- function(observed_value, het_type, es_type, data) {
  
  if (!het_type %in% c("I2", "CVH", "M", "sigma2", "V_bar")) {
    stop("Invalid heterogeneity type. Choose from 'I2', 'CVH', 'M', 'sigma2', or 'V_bar'.")
  }
  
  if (!es_type %in% unique(data$es.type)) {
    stop("Invalid effect size type. Check the levels of 'es.type'.")
  }
  
  filtered_data <- data %>% filter(es.type == es_type)
  
  if (nrow(filtered_data) == 0) {
    stop("No data available for the selected effect size type.")
  }
  
  het_values <- filtered_data[[het_type]]
  
  percentiles <- quantile(het_values, probs = seq(0, 1, by = 0.05), na.rm = TRUE)
  
  lower_bound <- max(percentiles[percentiles <= observed_value], na.rm = TRUE)
  upper_bound <- min(percentiles[percentiles >= observed_value], na.rm = TRUE)
  lower_percentile <- names(percentiles)[percentiles == lower_bound]
  upper_percentile <- names(percentiles)[percentiles == upper_bound]
  
  percentile_range <- paste0(lower_percentile, "-", upper_percentile, "th percentile")
  
  # return the results
  return(tibble::tibble(
    observed_value = observed_value,
    het_type = het_type,
    es_type = es_type,
    percentile_range = percentile_range
  ))
}


# function to calculate sampling variance
sigma2_v <- function(mod){
  sigma2_v <- sum(1 / mod$vi) * (mod$k - 1) /
    (sum(1 / mod$vi)^2 - sum((1 / mod$vi)^2))
  return(sigma2_v)
}


