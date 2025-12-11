#' infe_moy 
#'
#' @description A fct function
#'
#' @return The return value, if any, from executing the function.
#' @import tidyr
#' @noRd


estim_moy_gibbs<-  function(x,mu0=NULL,s_m0 =NULL, s0=NULL, n_s0=NULL, n_sample =1000 ) {
  
  moy_data = mean(x,na.rm=T)
  var_data = var(x,na.rm = T)
  n = length(x)
  PHI<-matrix (nrow=n_sample, ncol=2)
  PHI[1,] <- phi <- c(moy_data,1/var_data)
  ###
  ###  Gibbs   sampling
  set.seed (1)
  for(s in 2:n_sample)   {
    # generate a new theta value from its full conditional
    mun <-(mu0/s_m0 + n*moy_data*phi[2])/(1/s_m0  + n*phi[2])
    t2n <- 1/(1/s_m0 + n*phi[2])
    phi[1]<-rnorm (1,mun,sqrt(t2n))
    # generate a new 1/sigmaˆ2 value   from its  full conditional
    nun <- n_s0+n
    s2n <- (n_s0*s0 + (n-1)*var_data  + n*(moy_data - phi[1])^2)/nun
    phi[2] <- rgamma(1,nun/2,nun*s2n/2)
    PHI[s,] <- phi
  }
  return(PHI)
}


compare_moy_gibbs <- function(X, Y,mu0=NULL,s_m0 =NULL, s0=NULL, n_s0=NULL,
                              type =  NULL, plusieurs = NULL, seuild=NULL, twit=NULL,seuil_global = NULL,
                              n_sample =1000,IC=0.95, arr=3) {
  Y = as.factor(Y)
  noms = levels(Y)
  Ngroup=nlevels(Y)
  pourcent_IC2 = (1-IC)/2
  long<-factorial(Ngroup)/(factorial(Ngroup-2)*2)
  
  if(is.null(type)){
    
  }else{
    if((!is.null(type)) &  "seuil" %in%type){
      
      if(!is.null(plusieurs) & !is.null(seuil_global)){
        
        seuild = rep(seuil_global["seuil"],long)
        
      }else{
        
        seuild = unlist(seuild%>%unlist%>%na.omit()%>%as.vector())
        
        
      }
    }
    
  }
  
  if(length(mu0) != Ngroup) stop('Prior mu0 n\' as pas la même nombre de groupe que la variable Y')
  if(length(s_m0) != Ngroup) stop('Prior s_m0 n\' as pas la même nombre de groupe que la variable Y')
  if(length(s0) != Ngroup) stop('Prior s0 n\' as pas la même nombre de groupe que la variable Y')
  if(length(n_s0) != Ngroup) stop('Prior n_s0 n\' as pas la même nombre de groupe que la variable Y')
  
  descriptif<-function(x, probs=c(pourcent_IC2,0.5, 1-pourcent_IC2 )){
    return(c(list(mean(x)),quantile(x,probs=probs)))
  }
  res_infe<-lapply(1:Ngroup, function(x) estim_moy_gibbs(X[which(Y==levels(Y)[x])],
                                                         mu0[x],s_m0[x],s0[x],n_s0[x], n_sample =n_sample
  ))
  res_prior<- lapply(1:Ngroup, function(x) rnorm(n_sample,
                                                 mu0[x],s0[x]
  ))
  
  diffprior<-diffpost<-diff<-list()
  nomsdiff<- NULL
  for(i in 1:Ngroup){
    for(j in 2:Ngroup){
      if(j<=i)next
      # dput(diffprior[[i]])
      # dput(descriptif(res_prior[[i]]-res_prior[[j]]))
      diffprior<-append(diffprior,list(descriptif(res_prior[[i]]-res_prior[[j]])))
      nomsdiff<-c(nomsdiff,paste0("Diff Moy(groupe=",noms[i],")-Moy(Groupe=",noms[j],")"))
      diff<-append( diff,list(res_infe[[i]][,1]-res_infe[[j]][,1]))
      
      diffpost<- append(diffpost,list(descriptif(res_infe[[i]][,1]-res_infe[[j]][,1])))
    }
  }
  
  
  
  nomsx<-paste0("Moy groupe =",noms)
  noms.tab<-c(nomsx,nomsdiff)
  
  
  resprior<-c(lapply(res_prior,descriptif),diffprior)
  dput(diffprior)
  dput(lapply(res_prior,descriptif))
  # resprior<-unlist(resprior,recursive=F)
  # resprior<-lapply(resprior,function(i){round(i,4)})
  resprior<-data.table::rbindlist(resprior)
  resprior<-resprior%>%mutate_all(function(i){round(i,arr)})


  colnames(resprior)<-c("moy", paste0(c(pourcent_IC2,0.5, 1-pourcent_IC2 )*100, "%") )
  rownames(resprior)<-c(noms.tab)
  # print(resprior)
  
  res<-c(lapply(res_infe,descriptif),diffpost)
  
  
  
  res<-rbindlist(res)
  
  res<-res%>%mutate_all(function(i){round(i,arr)})
  
  colnames(res)<-c("moy", paste0(c(pourcent_IC2,0.5, 1-pourcent_IC2 )*100, "%") )
  
  rownames(res)<-c(noms.tab)
  
  
  
  
  
  if(is.null(type)){
  }else if(type =="seuil"){
    seuils<-c(seuild)
    tests<-NULL
    print("a regarder")
    print(length(diff))
    for(i in 1:length(diff)){
      tests<-c(tests,round(mean(diff[[i]]>seuils[i]),arr))
    }
    
    res$seuils<-c(rep(NA,Ngroup),seuils)
    res[,"Pr(X>seuil|D)"]<-c(rep(NA,Ngroup),tests)
    
    
  }else{
    
    sel<-which(apply(utils::combn(noms, 2),2, function(x) c(paste(x, collapse = "vs"),paste(rev(x), collapse = "vs")))%>%as.vector() ==twit$var )
    
    lequel <- ceiling(sel/2)
    signe = 1-sel%%2
    data=twit$data
    print(diff[[lequel]]* ifelse_perso(signe==1, 1, -1))
 
    data$`Pr Ha` = c(round(mean(between(diff[[lequel]]* ifelse_perso(signe==1, 1, -1),data["min_effet_absent"]%>%unlist,data["max_effet_absent"]%>%unlist)),arr))
    
    
    data$`Pr Hr` =  c(round(mean(between(diff[[lequel]]* ifelse_perso(signe==1, 1, -1),data["min_effet_present"]%>%unlist,data["max_effet_present"]%>%unlist)),arr))              
    names(data)<-c("Type","min Abs","max Abs","min pres","max pres",
                   "Pr(effet absent)","Pr(effet présent)")
    
    
    
  }
  
  
  
  # prior<-  data.frame(Loi = "Beta", "Paramètre alpha" =priors[1,] , "Paramètre beta" = priors[2,],row.names = paste("Groupe",noms), check.names = F)
  
  
  
  laliste<-list(
    # prior
    resprior,res)
  names(laliste) <- c(#"prior",
    "Valeurs a priori","Estimations a posteriori" )
  if(is.null(type)){
    
  }else if(type =="twit" & !is.null(type)){
    res_twit= list(TwoIt = data)
    names(res_twit) = twit$var
    laliste<-append(laliste,res_twit)
  }
  
  
  
  return(laliste)
  
  res_infe
  
  
  
  
  
}



Infe_moy2IT<-function(Y,priors,seuil = NULL, twit=NULL, 
                      arr=3,M=1000, IC=0.95, type =  NULL, color_1 ="#DE3163" , color_2="#40E0D0",color_3 = "#EBC341" ,color_4="#9242A6"){
  
  pourcent_IC2 = (1-IC)/2
  
  posteriors<-estim_moy_gibbs(Y,priors[1],priors[2],priors[3],priors[4],n_sample=M)
  posterior_mu <- posteriors[,1]
  posterior_sd <- posteriors[,2]
  
  
  descriptif<-function(x, probs=c(pourcent_IC2,0.5, 1-pourcent_IC2 )){
    return(c(mean(x),quantile(x,probs=probs)))
  }
  
  
  
  # a priori 
  
  xxprior_mu<-rnorm(M,priors[1],priors[2])
  xxprior_sd <- rgamma(M,priors[1],priors[2])
  descrxxprior_mu<-descriptif(xxprior_mu, probs = c(pourcent_IC2,0.5, 1-pourcent_IC2 ))%>%round(arr)
  descrxxprior_sd<-descriptif(xxprior_sd, probs = c(pourcent_IC2,0.5, 1-pourcent_IC2 ))%>%round(arr)
  
  descrxxpost_mu<-descriptif(posterior_mu, probs = c(pourcent_IC2,0.5, 1-pourcent_IC2 ))%>%round(arr)
  descrxxpost_sd<-descriptif(posterior_sd, probs = c(pourcent_IC2,0.5, 1-pourcent_IC2 ))%>%round(arr)
  
  
  res_df_mu = data.frame(prior = xxprior_mu,posterior=posterior_mu)%>%pivot_longer(cols = everything(),names_to = "type",values_to = "x")
  res_df_sd=data.frame(prior = xxprior_sd,posterior=posterior_sd  )%>%pivot_longer(cols = everything(),names_to = "type",values_to = "x")
  
  
  p<- list(
    res_df_mu%>%ggplot(aes(x,fill=type))+geom_histogram(alpha=0.7)+theme_light() +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank()
      )+ ylab("") +
      scale_fill_manual(name = "Distribution",
                        breaks = c("prior", "posterior"),
                        values = c("prior" = color_3, "posterior" = color_4) )
    
    ,
    res_df_sd%>%ggplot(aes(x,fill=type))+geom_histogram(alpha=0.7)+theme_light() +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank()
      )+ ylab("") +
      scale_fill_manual(name = "Distribution",
                        breaks = c("prior", "posterior"),
                        values = c("prior" = color_3, "posterior" = color_4) )
  )
  
  
  
  
  
  
  resprior<-rbind(descrxxprior_mu,descrxxprior_sd)
  
  colnames(resprior)<-c("moy", paste0(c(pourcent_IC2,0.5, 1-pourcent_IC2 )*100, "%") )
  rownames(resprior)<-c("mu","sd")
  
  
  respost<-rbind(descrxxpost_mu,descrxxpost_sd)%>%as.data.frame()
  
  colnames(respost)<-c("moy", paste0(c(pourcent_IC2,0.5, 1-pourcent_IC2 )*100, "%") )
  rownames(respost)<-c("mu","sd")
  
  
  
  
  
  
  
  if(is.null(type)){
  }else if(type =="seuil"){
    
    tests <- round(mean(posterior_mu>seuil),arr)
    
    
    print(seuil)
    print(tests)
    respost$seuils<-c(seuil,NA)
    respost$`Pr(X>seuil|D)`<-c(tests,NA)
    print(respost)
    
    p[[1]]<-  p[[1]]+ geom_segment( aes(x=seuil, xend = seuil,y=-Inf,yend=Inf))
  }else{
    
    
    data=twit$data
    data$`Pr Ha` = round(mean(between(posterior_mu,data[,"min_effet_absent"],data[,"max_effet_absent"])),arr)
    
    
    data$`Pr Hr` =  round(mean(between(posterior_mu,data[,"min_effet_present"],data[,"max_effet_present"])),arr)
    
    
    p[[1]]<-  p[[1]]+geom_rect( mapping= aes( xmin = data[,"min_effet_absent"],xmax = data[,"max_effet_absent"], ymin = -Inf,
                                              ymax = Inf, x=0, y=0,fill ="Acceptée" ), alpha = 0.2 )+
      geom_rect(aes( xmin = data[,"min_effet_present"],
                     xmax = data[,"max_effet_present"], ymin = -Inf,
                     ymax = Inf,x=0,y=0 ,fill ="Rejetée" ), alpha = 0.2 )+
      scale_fill_manual(name = "Hypothèse",
                        labels =c("P : effet présent ","A : Effet absent"),
                        breaks = c("Acceptée", "Rejetée"),
                        values = c("Acceptée" = color_2, "Rejetée" = color_1))
    
    names(data)<-c("min Abs","max Abs","min pres","max pres",
                   "Pr(effet absent)","Pr(effet présent)")
  }
  
  prior<-  data.frame( "Moyenne (mu)" =priors[1], "Moyenne (sd)"=priors[2] ,
                      "Sd (alpha)" =priors[3],"Sd (Beta)"=priors[4],row.names = ("Prior"), check.names = F)
  
  
  
  df<-list(prior,resprior,respost)
  names(df) <- c("prior","Valeurs a priori","Estimations a posteriori" )
  if(is.null(type)){
    
  }else if(type =="twit" & !is.null(type)){
    res_twit= list(TwoIt = data)
    names(res_twit) = twit$var
    df<-append(df,res_twit)
  }
  
  
  
  return(list(df=df, graph = p))
}



