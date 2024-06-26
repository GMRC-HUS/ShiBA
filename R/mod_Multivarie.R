#' Multivarie UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#' @param r variable in the evironnement
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
#' @import broom.mixed
#' @import sortable
#' @import rstanarm
#' @import bayesplot
#' @import shinyWidgets
#' @import shinyalert
#' @import shinycssloaders
#' @import kableExtra
#' @import waiter
#' @import colourpicker

mod_Multivarie_ui <- function(id) {
  ns <- NS(id)
  tagList(
    fluidPage(
      waiter::use_waiter(),
      titlePanel("Multivarié"),
      sidebarLayout(
        sidebarPanel(
          width = 4,
          radioButtons(
            ns("type_glm"), HTML(paste(h2("Type de regression")#, text_aide("Choix type de régression multivarié !")
                                       )),
            c(
              "Linéaire" = "lin",
              Binomial = "binom",
              # Beta = "beta",
              Poisson = "poiss"
            ), "lin"
          ),
          sliderInput(ns("IC"),label = "Intervalle de Crédibilité en %",min = 80,max = 100,step = 1,animate = F,post = " %",value = 95),
          uiOutput(ns("choix_y")),
          uiOutput(ns("propositions_multi")),
          uiOutput(ns("refactorisation")),
          h3("Choix des priors"),
          actionButton(ns("ellicitation"), "Ellicitation"),
          #text_aide("Texte Aide ellicitation multivarié "),
          h3("Seuils/Two IT ?"), #text_aide("Texte Aide Two IT multivarié "),
          # shinyWidgets::materialSwitch(ns("twit"), "", value =FALSE, status = "success", right = T),
          uiOutput(ns("twit_ui")),
          br(),
          actionButton(ns("go"), "Go :")
        ),
        mainPanel(
          # tags$head(tags$style(".butt{background-color:#E9967A;} .butt{color: black;}")),
          uiOutput(ns("result_multi"))
        ) # fin MainPanel
      ) # fin sidebarlayout
    ) # fin fluidpage
  )
}

#' Multivarie Server Functions
#'
#' @noRd
mod_Multivarie_server <- function(id, r) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    output$choix_y <- renderUI({
      
      if( input$type_glm=="poiss"){
      
      choix_possible = apply(r$BDD%>%select_if(is.numeric),2, function(x) ifelse_perso(sum(!x%%1==0, na.rm=T)==0 & min(x,na.rm = T)>=0, return(T), return()))%>%unlist%>%names
      }else if( input$type_glm=="binom"){
        choix_possible = apply(r$BDD,2, function(x) ifelse_perso(nlevels(as.factor(unique(x)))==2, return(T), return()))%>%unlist%>%names
      }else{
        choix_possible = r$noms
      }
      
      
      selectInput(ns("variable"), h3("Variable d'interêt :"),
        choices =
          choix_possible
      )
    })



    # Choix seuil two_it
    output$twit_ui <- renderUI({
      ui_twit("seuil_2it", ns)
     
    })
    #
    seuil_twoit <- reactiveVal(value = NULL)
    
    seuil_twoit_val <- reactiveVal(value = NULL)
    #
    #   # var_input2
    observeEvent(c(input$list_quali, input$list_quanti, input$variable, input$type_glm), {
      output$twit_ui <- renderUI({
        ui_twit("seuil_2it", ns)
       
      })

      model_2(NULL)

      var_quali <- isolate(input$list_quali)
      if (length(var_quali) > 0) {
        nom_var_quali <- lapply(var_quali, function(x) paste(x, levels(factor(r$BDD[, x]))[-1], sep = "")) %>% unlist()
      } else {
        nom_var_quali <- NULL
      }

      var <- c(isolate(input$list_quanti), nom_var_quali)
      seuil_twoit(twitServer("id_i", var, type_glm = input$type_glm))
    })

    observeEvent(input$seuil_2it, {
      model_2(NULL)
      output$twit_ui <- renderUI({
        ui_twit("seuil_2it", ns, color="success", icon = icon("check"))
  
      })

      twitUi(ns("id_i"))
      # seuil_twoit(twitServer("id_i",var))
    })









    # Resultat table of model




    var_input <- reactiveValues(
      choix_base = NULL,
      var_quali = NULL,
      var_quanti = NULL
    )


    var_quali_sel <- reactiveValues(var = NULL)

    var_quanti_sel <- reactiveValues(var = NULL)
    
    observeEvent(c(input$variable,input$list_quali,input$list_quanti), {
      liste_choix <- r$noms
      liste_choix <- liste_choix[-c(which(liste_choix == input$variable),which(liste_choix %in% input$list_quali),which(liste_choix %in% input$list_quanti))]# on met in car dans list_quali on peut avoir plus d'une variable
      var_input$choix_base <- liste_choix
      var_input$var_quali<- var_input$var_quali[var_input$var_quali!= input$variable]
      var_input$var_quanti<-var_input$var_quanti[var_input$var_quanti!= input$variable]  
      
      }
      
      )
    
    observeEvent(input$type_glm,{
      var_input$var_quali <- NULL
      var_input$var_quanti <- NULL
  })

    observeEvent(c(input$choix_base,input$list_quant, input$list_quali), {
      quantis <- input$list_quanti
      qualis <- input$list_quali
      base_encours <- input$choix_base

      base_avant <- var_input$choix_base
      quali_avant <- var_input$var_quali
      quanti_avant <- var_input$var_quanti

      nvx_quantis <- length(quantis) > length(quanti_avant)
      nvx_base <- length(base_encours) > length(base_avant)
      nvx_qualis <- length(qualis) > length(quali_avant)

      if (nvx_base) {
        var_input$choix_base <- base_encours
        var_input$var_quali <- qualis
        var_input$var_quanti <- quantis
      } else if (nvx_quantis) {
        var_sel <- setdiff(quantis, quanti_avant)
        if (sum(is.na(as.numeric(as.character(r$BDD[, var_sel])))) > sum(is.na(r$BDD[, var_sel]))) {
          showNotification(HTML("<b>", var_sel, "</b> ne semble pas être une variable quantitatives.<br>Au moins une donnée non numérique."), type = "warning")
          var_input$choix_base <- base_avant
          var_input$var_quanti <- quanti_avant
          var_input$var_quali <- quali_avant
        } else {
          var_input$choix_base <- base_encours
          var_input$var_quanti <- quantis
          var_input$var_quali <- qualis
        }
      } else if (nvx_qualis) {
        var_sel <- setdiff(qualis, quali_avant)
        if (length(unique(r$BDD[, var_sel])) > 20) {
          showNotification(HTML("<b>", var_sel, "</b> ne semble pas être une variable qualitative.<br>Plus de 20 modalités différentes"), type = "warning")
          var_input$choix_base <- base_avant
          var_input$var_quanti <- quanti_avant
          var_input$var_quali <- quali_avant
        } else {
          var_input$choix_base <- base_encours
          var_input$var_quanti <- quantis
          var_input$var_quali <- qualis
        }
      }
    })
    

      output$propositions_multi <- renderUI({
      bucket_list(
        header = HTML("<h3>Variable explicatives :</h3>"),
        group_name = "bucket_list_group",
        orientation = "vertical",
        add_rank_list(
          text = "Noms des variables",
          labels = var_input$choix_base,
          input_id = ns("choix_base")
        ),
        add_rank_list(
          text = HTML("<strong>Variables quantitatives</strong>"),
          labels = (var_input$var_quanti),
          input_id = ns("list_quanti")
        ),
        add_rank_list(
          text = HTML("<strong>Variables qualitatives</strong>"),
          labels = (var_input$var_quali),
          input_id = ns("list_quali")
        )
      )
    })

    #

    
    output$refactorisation <- renderUI({
      if (length(input$list_quali) == 0) {
        return()
      } else {
        actionButton(ns("refact_button"), "Sélection des modalités de référence")
      }
    })


    observeEvent(input$refact_button, {
      nom_var_quali <- isolate(input$list_quali)

      showModal(
        modalDialog(
          title = "Sélection des références",
          tagList(lapply(1:length(nom_var_quali), function(i) {
            x <- nom_var_quali[i]
            r$BDD[, x] <- as.factor(r$BDD[, x])
            var <- r$BDD[, x]
            var <- as.factor(var)
            noms_levels <- levels(var)


            list(
              radioButtons(inputId = ns(paste("fact_", x, sep = "")), label = x, choices = noms_levels, selected = noms_levels[1])
            )
          })),
          footer = tagList(
            actionButton(ns("ok_fact"), "OK")
          )
        )
      )
    })

    observeEvent(input$ok_fact, {
      nom_var_quali <- isolate(input$list_quali)
      for (x in nom_var_quali) {
        r$BDD[, x] <- relevel(r$BDD[, x], ref = input[[paste("fact_", x, sep = "")]])
      }
      removeModal()
    })

    randomVals <- eventReactive(input$go, {
      runif(n = 1)
    })



    output$model_text <- renderText({
      randomVals()


      formule_default(isolate(input$variable), isolate(input$list_quanti), isolate(input$list_quali))
    })



    model_2 <- reactiveVal(value = NULL)

    output$res_multi <- function() {
      if (is.null(model_2())) {
        return()
      }


      print(     dput(seuil_twoit()$ls))
      
      res <- shibaGlmTable(model_2(),IC = isolate(input$IC), 
                           input$type_glm, 
                           seuilTwoIt = isolate(seuil_twoit()$ls))
      loo<-loo(model_2()) 
      loo$estimates[3,1] # Estimate
      loo$estimates[3,2] # SE
      
    
      opts <- options(knitr.kable.NA = "")
      
      if( input$type_glm == "lin"){
      res %>%
        kbl(digits = 3, caption = paste0("loo : ",round(loo$estimates[3,1],2),
                                         "; se: ",round(loo$estimates[3,2],2),"")) %>%
        kable_styling(full_width = F, bootstrap_options = c("hover"), fixed_thead = T) %>%
        row_spec((nrow(res) ), extra_css = "border-top: thick double #D8D8D8")
      }else{
        res %>%
          kbl(digits = 3, caption = paste0("loo : ",round(loo$estimates[3,1],2),
                                           "; se: ",round(loo$estimates[3,2],2),"")) %>%
          kable_styling(full_width = F, bootstrap_options = c("hover"), fixed_thead = T) 
        
      }
      
      # column_spec(4, extra_css = "border-right: thick double #D8D8D8;")
    }


    prior_glm <- reactiveValues(
      prior_intercept = NULL,
      prior_beta_scale = NULL,
      prior_beta_location = NULL
    )
    observeEvent(input$choix_base, {
      prior_glm$prior_intercept <- NULL
      prior_glm$prior_beta_scale <- NULL
      prior_glm$prior_beta_location <- NULL
    })


    
  
    
    observeEvent(input$go, {
     
      seuil_twoit_val(seuil_twoit())
      showModal(modalDialog(
        title = "Modèle en cours",
        fluidRow(align="center",
        img(src="www/wait.gif", align = "center",height='300px',width='400px')
        )
       
      ))

 
     
      list_quanti <- isolate(input$list_quanti)
      list_quali <- isolate(input$list_quali)
      y <- isolate(input$variable)

      
      data <- r$BDD %>% select(y, list_quanti, list_quali)

      
      if( input$type_glm=="poiss"){
        
        data[,y]=as.numeric(data[,y])
      }else if( input$type_glm=="binom"){
        data[,y]=as.factor(data[,y])
      }else{
        data[,y]=as.numeric(data[,y])}
    
      if (length(list_quanti) > 0) data <- data %>% mutate_at(list_quanti, ~ as.numeric(as.character(.)))

      if (length(list_quali) > 0) data <- data %>% mutate_at(list_quali, as.factor)
      formule <- formule_default(y, list_quanti, list_quali)
     
  
        
        fit <- glm_Shiba(formule,
        data = data,
        prior_intercept = prior_glm$prior_intercept,
        prior = list(scale = prior_glm$prior_beta_scale, location = prior_glm$prior_beta_location),
        type = input$type_glm
      )
       print("ou encore la")
        
      model_2(NULL)
      

      removeModal()
    
      
      
      model_2(fit)
      
    })


    output$model_prior <- renderTable({
      if (is.null(model_2())) {
        return()
      }
      
      df_prior(model_2(), input$type_glm)
      
    })


    # Afficher les diags de convergence :
    output$diag <- renderUI({
      if (is.null(model_2())) {
        return()
      }


      if (diag_convergence(model_2())) {
        actionBttn(
          inputId = ns("convergence"),
          label = "Analyse de convergence",
          style = "gradient",
          color = "success",
          icon = icon("check")
        )
      } else {
        actionBttn(
          inputId = ns("convergence"),
          label = "Analyse de convergence",
          style = "gradient",
          color = "warning",
          icon = icon("exclamation")
        )
      }
    })
    output$graph_model <- renderPlot({
      if (is.null(model_2())) {
        return()
      }

      shibaGlmPlot(model_2(), type_glm = input$type_glm, pars = input$Variable_graph,
                   seuilTwoIt = isolate(seuil_twoit()$ls), prob = isolate(input$IC)/100,
                   hist = ifelse(input$hist_graph =="Histogramme",T,F),color_1 = input$col1,color_2= input$col2)
      # bayesplot::mcmc_areas(model_2() %>% as.matrix(), pars = input$Variable_graph)+theme_light()
    })


    output$result_multi <- renderUI({
      if (is.null(model_2())) {
        return()
      }
    
      var_model <- colnames(model_2()$covmat)
      type_model <- names(which(c(
        "Linéaire" = "lin",
        Binomial = "binom",
        Beta = "beta",
        Poisson = "poiss"
      ) == input$type_glm))


      fluidPage(
        div(
          id = "div_model",
    
          br(),
         fluidRow(
            box(title = paste("Modèle :",type_model ),align = "center",status = "primary", solidHeader = TRUE,
            
            h3(textOutput(ns("model_text")) %>% withSpinner())),
        
          box(width = 12,title = ("Prior :"),status = "primary", solidHeader = TRUE,
              fluidRow(align = "center", tableOutput(ns("model_prior")) %>% withSpinner()))
          ),
         br(),
          fluidRow(align = "center", uiOutput(ns("diag")) %>% withSpinner()),
          br(),
          box(title = 
          ("Résultats :"),width  = 12,status = "primary", solidHeader = TRUE,
          fluidRow(align = "center",
         tableOutput(ns("res_multi")) %>% withSpinner())),
          
         box(title = 
          ("Graphiques :"),width  = 12,status = "primary", solidHeader = TRUE,fluidRow(align = "center",
         
       
            dropdownButton(up=T,
            
            
              fluidPage(awesomeCheckboxGroup(
                inputId = ns("Variable_graph"),
                label = "Variables à afficher",
                choices = var_model,
                selected = ifelse_perso(length(var_model) > 1, var_model[-1], var_model)
              ),
              
              awesomeRadio(
                inputId = ns("hist_graph"),
                label = "Type de graphique :", 
                choices = c("Histogramme", "Densité"),
                selected = "Histogramme",
                inline = TRUE
              ),
              
              
              colourInput(
                ns("col1"), "Couleur 1", "#DE3163",
                showColour = "background"),
              
              colourInput(
                ns("col2"),  "Couleur 2", "#40E0D0",
                showColour = "background")),
              
              
              
              circle = TRUE,
              status = "primary",
              icon = icon("gear"), width = "300px",
              tooltip = tooltipOptions(title = "Modifier les pamètres graphiques")
            )
            ,
            plotOutput(width = 600, height = 400, ns("graph_model")) %>% withSpinner()
          ),
          br()
        )
      ))
    })
    # Afficher les prior


    # diagnostique de convergence
    observeEvent(input$convergence, {
      
      
      
      nom_var_quali <- lapply(input$list_quali, function(x) paste(x, levels(factor(r$BDD[, x]))[-1], sep = "")) %>% unlist()
      variables_conv<- c("(Intercept)", input$list_quanti, nom_var_quali, "sigma")
      if( input$type_glm != "lin")  variables_conv<- c("(Intercept)", input$list_quanti, nom_var_quali)
      showModal(
        modalDialog(
          tagList(
            div(
              align = "center",
              lapply(variables_conv, function(x) {
                list(
                  h2(x),
                  (plotOutput(width = "100%", height = 400, ns(paste(x, "_courbe_diag", sep = "")))) %>% withSpinner()
                )
              })
            )
          ),
          footer = modalButton("", icon = icon("xmark")),
          easyClose = T,
          size = "l",
        )
      )



      lapply(c("(Intercept)", isolate(input$list_quanti), nom_var_quali, "sigma"), function(i) {
        output[[paste(i, "_courbe_diag", sep = "")]] <- renderPlot({
          plot_diag(model_2(), i)
        })
      })
    })


    # Action of Ellicitaion button
    observeEvent(input$ellicitation, ignoreInit = T, {
      if (is.null(prior_glm$prior_intercept) | is.null(prior_glm$prior_beta_scale)) {
        prior_quali_sd <- sd_quali(input$list_quali, r$BDD)

        prior_quanti_sd <- sapply(input$list_quanti, function(x) sd(r$BDD[, x], na.rm = T))
        
        if (length(prior_quanti_sd) == 0) prior_quanti_sd <- NULL

        if (input$type_glm == "lin") {
          prior_glm$prior_intercept <- c(
            round(mean(r$BDD[, input$variable], na.rm = T), 2),
            round(2.5 * sd(r$BDD[, input$variable], na.rm = T), 2)
          )
        } else if (input$type_glm %in% c("binom", "poiss")) {
          prior_glm$prior_intercept <- c(0, 2.5)
        }
        prior_glm$prior_beta_scale <- round(2.5 / c(prior_quanti_sd, prior_quali_sd) * sd(r$BDD[, input$variable], na.rm = T), 2)
        prior_glm$prior_beta_location <- rep(0, length(c(prior_quanti_sd, prior_quali_sd)))
      }



      prior_beta_scale_def <- c(prior_glm$prior_intercept[2], prior_glm$prior_beta_scale) # ), default_prior_beta_scale_def, prior_glm$prior_beta_scale)
      prior_beta_location_def <- c(prior_glm$prior_intercept[1], prior_glm$prior_beta_location) # ), default_prior_beta_location_def, prior_glm$prior_beta_location)

      if (length(input$list_quali) > 0) {
        nom_var_quali <- lapply(input$list_quali, function(x) paste(x, levels(factor(r$BDD[, x]))[-1], sep = "_")) %>% unlist()
      } else {
        nom_var_quali <- NULL
      }


      variables <-c("intercept",input$list_quanti,nom_var_quali)
      if (input$type_glm == "lin") {
      showModal(
        modalDialog(size = "l",
          tagList(lapply(1:length(variables), function(i) ui_choix_prior_norm(i,variables,ns,prior_beta_location_def,prior_beta_scale_def )
            )),
          footer = tagList(
            actionButton(ns("ok"), "OK"),
            actionButton(ns("defaut"), "Défaut")
          )
        )
      )
        
        noms <-variables
        positions <- c(prior_glm$prior_intercept[1], prior_glm$prior_beta_location)
        dispersions <- c(prior_glm$prior_intercept[2], prior_glm$prior_beta_scale)
        
        lapply(1:length(noms), function(i) {
          output[[paste(noms[i], "_courbe", sep = "")]] <- ui_ggplot_prior_norm(i,input,variables)
        })
      }else if (input$type_glm %in% c("binom", "poiss")) {
        
        showModal(
          modalDialog(size = "l",
            tagList(lapply(2:length(variables), 
                           function(i) ui_choix_prior_exp(i,variables,ns,prior_beta_location_def,prior_beta_scale_def ))),
            footer = tagList(
              actionButton(ns("ok"), "OK"),
              actionButton(ns("defaut"), "Défaut")
            )
          )
        )
        
        noms <-variables
        positions <- c(prior_glm$prior_intercept[1], prior_glm$prior_beta_location)
        dispersions <- c(prior_glm$prior_intercept[2], prior_glm$prior_beta_scale)
        
        lapply(1:length(noms), function(i) {
          output[[paste(noms[i], "_courbe", sep = "")]] <- ui_ggplot_prior_exp(i,input,variables)
        })  
      }

   
    })

    observeEvent(input$defaut, {
      
      prior_quanti_sd <- c()
      if (!length(c(input$list_quanti, input$list_quali)) == 0) {
        prior_quali_sd <- sd_quali(input$list_quali, r$BDD)

        prior_quanti_sd <- sapply(input$list_quanti, function(x) sd(r$BDD[, x], na.rm = T))
        
        default_prior_beta_scale_def <- round(2.5 / c(prior_quanti_sd, prior_quali_sd) * sd(r$BDD[, input$variable], na.rm = T), 2)
        default_prior_beta_location_def <- 0
       
      }
      if (length(prior_quanti_sd) == 0) prior_quanti_sd <- NULL
      default_prior_intercept_def <- c(
        round(mean(r$BDD[, input$variable], na.rm = T), 2),
        round(2.5 * sd(r$BDD[, input$variable], na.rm = T)), 2
      )



      if (length(input$list_quali) > 0) {
        nom_var_quali <- lapply(input$list_quali, function(x) paste(x, levels(factor(r$BDD[, x]))[-1], sep = "_")) %>% unlist()
      } else {
        nom_var_quali <- NULL
      }
      
      
      variables <-c("intercept",input$list_quanti,nom_var_quali)
      
      
      for (x in 1:length(variables)) {
        updateNumericInput(session, paste(variables[x], "_mu_0", sep = ""), value = ifelse(x == 1, default_prior_intercept_def[1], default_prior_beta_location_def))
        updateNumericInput(session, paste(variables[x], "_sigma_0", sep = ""), value = ifelse(x == 1, default_prior_intercept_def[2], default_prior_beta_scale_def[x - 1]))
     
        if (!input$type_glm == "lin") {
          transfo_exp = norm_tomin_max_exp(ifelse(x == 1, default_prior_intercept_def[1], default_prior_beta_location_def) , ifelse(x == 1, default_prior_intercept_def[2], default_prior_beta_scale_def[x - 1])   )
          transfo_exp<-as.vector(transfo_exp)
          
          updateNumericInput(session, paste(variables[x], "_min_exp", sep = ""), value = round(transfo_exp[1],3) )
          updateNumericInput(session, paste(variables[x], "_max_exp", sep = ""), value = round(transfo_exp[2],3))
          
        }
        
         }
    })

    observeEvent(input$ok, {
      if (length(input$list_quali) > 0) {
        nom_var_quali <- lapply(input$list_quali, function(x) paste(x, levels(factor(r$BDD[, x]))[-1], sep = "_")) %>% unlist()
      } else {
        nom_var_quali <- NULL
      }
      
      
      variables <-c("intercept",input$list_quanti,nom_var_quali)

      if(input$type_glm == "lin"){
      prior_mu <- unlist(lapply(variables, function(i) {
        input[[paste(i, "_mu_0", sep = "")]]
      }))
      prior_sd <- unlist(lapply(variables, function(i) {
        input[[paste(i, "_sigma_0", sep = "")]]
      }))

      
      }else{
        res_trans <- lapply(variables[-1], function(i){
          
          min = input[[paste(i, "_min_exp", sep = "")]]
          max=  input[[paste(i, "_max_exp", sep = "")]]
          
          trans = min_max_exp_to_norm(min, max)
          return(data.frame(prior_mu = trans[1], prior_sd = trans[2]))
        }
        )%>%rbindlist
     
        
        prior_mu<- c(0,res_trans$prior_mu)
        prior_sd<- c(2.5,res_trans$prior_sd)
      }
      removeModal()
      prior_glm$prior_intercept <- c(prior_mu[1], prior_sd[1])

      prior_glm$prior_beta_scale <- prior_sd[-1]
      prior_glm$prior_beta_location <- prior_mu[-1]
      if (length(prior_sd) == 1) {
        prior_glm$prior_beta_scale <- NULL
        prior_glm$prior_beta_location <- NULL
      }
    })
  })
}



## To be copied in the UI
# mod_Multivarie_ui("Multivarie_1")

## To be copied in the server
# mod_Multivarie_server("Multivarie_1")
