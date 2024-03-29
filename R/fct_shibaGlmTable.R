#' shibaGlmTable
#'
#' @description A fct function
#'
#' @return The return value, if any, from executing the function.
#' @importFrom data.table rbindlist
#' @noRd



shibaGlmTable <- function(fit,IC=95, type_glm, seuilTwoIt = NULL, ...) {
  
  nomsModel <- fit$coefficients %>% names()
  noms<-ifelse_perso (type_glm %in% c("poiss", "binom"), c(nomsModel), c(nomsModel,"sigma"))
  
  res <- lapply(noms,
                function(i)  quantile(as.array(fit)[, , i], c(0.5,(1-IC/100)/2,1-(1-IC/100)/2))%>%t%>%as.data.frame(   check.names= F)%>% mutate(var=i, .before=everything()))%>%rbindlist
  
  
  if (is.null(seuilTwoIt)) {
    
    
    if (type_glm %in% c("poiss", "binom")) {
      res[, 2:4] <- exp(res[, 2:4])
      names(res)[1] <- ifelse(type_glm == "poiss", "RR", "OR")
    }
    
    return(res)
  }
  
  
  
  
  
  if (seuilTwoIt$type == "seuil") {
    if (type_glm %in% c("poiss", "binom")) {
      seuils <- log(seuilTwoIt$val)
      seuil_init <- seuilTwoIt$val
    } else {
      seuil_init <- seuils <- (seuilTwoIt$val)
    }
    
    
    if (seuilTwoIt$plusieur_seuils) {
      seuils_df <- data.frame(
        seuils = c(NA, seuils),
        parameter = nomsModel
      )
      
      
      res_s <- sapply(
        nomsModel[-1],
        function(i) length(which(as.array(fit)[, , i] > seuils_df$seuils[seuils_df$parameter == i])) / length(as.array(fit)[, , i])
      )
    } else {
      seuils <- rep(seuils, length(nomsModel))
      seuil_init <- rep(seuil_init, length(nomsModel))
      seuils[1] <- NA
      seuil_init <- seuil_init[-1]
      seuils_df <- data.frame(
        seuils = seuils,
        parameter = nomsModel
      )
      res_s <- sapply(
        nomsModel[-1],
        function(i) length(which(as.array(fit)[, , i] > seuils_df$seuils[seuils_df$parameter == i])) / length(as.array(fit)[, , i])
      )
    }
    
    res_seuil <- data.frame(Seuil = c("", seuil_init), Prob = c("", res_s))
    res_seuil[(nrow(res_seuil) + 1):(nrow(res_seuil) + (nrow(res) - nrow(res_seuil))), ] <- ""
    res <- cbind(res, res_seuil)
    
    
    if (type_glm %in% c("poiss", "binom")) {
      res[, 1:3] <- exp(res[, 1:3])
      names(res)[1] <- ifelse(type_glm == "poiss", "RR", "OR")
    }
    
    return(res)
  } else if (!is.null(seuilTwoIt$val$var)) {
    if (!seuilTwoIt$val$var %in% nomsModel) {
      if (type_glm %in% c("poiss", "binom")) {
        res[, 2:4] <- exp(res[, 2:4])
        names(res)[1] <- ifelse(type_glm == "poiss", "RR", "OR")
      }
      
      return(res)
    }
    
    
    
    
    list_param <- seuilTwoIt$val
    if (type_glm %in% c("poiss", "binom")) {
      twoIt <- twoItStanGlm(fit, list_param$var,
                            HA_diff_l = log(list_param$theta_A_min), HA_diff_u = log(list_param$theta_A_max),
                            HP_diff_l = log(list_param$theta_P_min), HP_diff_u = log(list_param$theta_P_max)
      )
      
      twoIt$names[1:2] <- c(
        paste0("PR(", list_param$theta_A_min, " < diff < ", list_param$theta_A_max, ")"),
        paste0("PR(", list_param$theta_P_min, " < diff < ", list_param$theta_P_max, ")")
      )
    } else {
      twoIt <- twoItStanGlm(fit, list_param$var,
                            HA_diff_l = list_param$theta_A_min, HA_diff_u = list_param$theta_A_max,
                            HP_diff_l = list_param$theta_P_min, HP_diff_u = list_param$theta_P_max
      )
    }
    res<-as.data.frame(res)
    
    res[nrow(res)+1,] <- NA
    
    Pr <- data.frame(c("H|Prior", "H|Données"), twoIt$values[1:2], twoIt$values[3:4])
    names(Pr) <- c("Two It", twoIt$names[1:2])
    
    
    
    
    ligne_para<- which(res$var == list_param$var)
    n_ligne<- nrow(res)
    
    matrice_vide <- matrix("",nrow = n_ligne,ncol = 3)
    
    matrice_vide[ligne_para:(ligne_para+1),]<-as.matrix(Pr)
    matrice_vide<-as.data.frame(matrice_vide)
    names(matrice_vide)<- names(Pr)
    
    Pr<-matrice_vide
    
    
    
    
    
    ligne_vide <- as.data.frame(matrix(NA,nrow = 1,ncol = ncol(res)))
    rownames(ligne_vide)<-""
    names(ligne_vide) <- names(res)
    
    
    
    
    
    res<- rbind(res[1:ligne_para,],ligne_vide, res[(ligne_para+1):n_ligne,])
    res<-res[-nrow(res),]
    
    
    
    
    res <- cbind(res, Pr)
    # res<-res[-nrow(res),]
    
    
    if (type_glm %in% c("poiss", "binom")) {
      res[, 2:4] <- exp(res[, 2:4])
      names(res)[1] <- ifelse(type_glm == "poiss", "RR", "OR")
    }
    
    return(res)
  }
}
