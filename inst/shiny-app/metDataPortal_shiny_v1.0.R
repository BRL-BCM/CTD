require(shiny)
require(shinycssloaders)
require(shinyWidgets)
require(DT)
require(shinydashboard)
require(ggplot2)
require(grid)
require(Hmisc)
require(plotly)
require(grDevices)
require(igraph)
require(networkD3)
require(visNetwork)
require(R.utils)
require(CTD)
require(cowplot)
require(gridExtra)
require(gtools)
data("Thistlethwaite2020")

disMod <<- loadToEnv(system.file(sprintf("shiny-app/disMod.RData"), package = "CTD"))[["disMod"]]
modelChoices <<- tolower(unique(sapply(list.files(system.file("ranks/ind_ranks",package = "CTD")),function(x) sub("[0-9]+-ranks.RData","",x))))

source(system.file("shiny-app/metDataPortal_appFns.r",package = "CTD"))
cohorts_coded <<- lapply(cohorts_coded, mixedsort, decreasing=TRUE)

pwy_choices = c("Choose", "Arginine Metabolism", "Ascorbate Metabolism", "Asp-Glu Metabolism", "BCAA Metabolism",
                "Benzoate Metabolism", "Beta-Oxidation", "Bile-Acid Metabolism", "Carnitine Biosynthesis",
                "Cholesterol Synthesis", "Creatine Metabolism", "Dicarboxylic Acid Metabolism", "Eicosanoids",
                "Endocannabinoid Synthesis", "Fatty Acid Metabolism", "Fibrinogen Cleavage Peptides", "GABA Shunt",
                "Galactose Metabolism", "Glutathione Metabolism", "Gly-Ser-Thr Metabolism", "Glycogen Metabolism",
                "Glycolysis", "Glycosylation", "Hemoglobin-Porphyrin Metabolism", "Histidine Metabolism", "Inositol Metabolism",
                "Ketone Bodies", "Lysine Catabolism", "Met-Cys Metabolism", "Mevalonate Metabolism", "Nicotinate-Nicotinamide Metabolism",
                "Pantothenate Metabolism", "Pentose-Phosphate Metabolism", "Phe-Tyr Metabolism", "Phospholipid Metabolism",
                "Polyamine Metabolism", "Proline Metabolism", "Protein Degradation", "Purine Metabolism", "Pyridoxal Metabolism",
                "Pyrimidine Metabolism", "Riboflavin Metabolism", "Secondary-Bile-Acids", "Sorbitol-Glycerol Metabolism",
                "Sphingolipid-Metabolism","Steroid-Hormone Biosynthesis", "TCA Cycle", "Thyroid Hormone Synthesis",
                "Tryptophan Metabolism", "All")

ui = dashboardPage(
  dashboardHeader(title = "Metabolomics Data Portal"),
  dashboardSidebar(sidebarMenu(id = "tab",style = "position:fixed;",
                               menuItem("View Patient Report", tabName = "ptReport", icon = icon("user-circle-o")),
                               menuItem("Network-Assisted Diagnostics", tabName = "ctd", icon=icon("project-diagram")),
                               menuItem("Inspect Reference Population", tabName = "refPop", icon = icon("bar-chart")),
                               menuItem("Download Data", tabName = "download", icon = icon("download"))
  )),
  dashboardBody(
    tags$script(HTML("$('body').addClass('fixed');")),
    tabItems(
      tabItem(tabName="ptReport",
              fluidRow(h2("Patient Report", align="center"),
                       box(title="Select Patient(s)", status="warning", solidHeader = TRUE, width = 12,
                           fluidRow(column(width = 6,pickerInput(inputId = "diagClass",label = "Select diagnosis.",choices = names(cohorts_coded),selected = names(cohorts_coded)[1], inline = FALSE,width = "100%")),
                                    column(width = 6,pickerInput(inputId = "ptIDs",label = "Select patient.",choices = "",inline = FALSE,width = "100%",multiple=FALSE)))),
                       box(title = "Patient Report", status="info", solidHeader = TRUE, align="left", width=12, collapsible=TRUE,
                           tabsetPanel(type="tabs",tabPanel("Metabolomic Profile", dataTableOutput("patientReport") %>% withSpinner(color="#0dc5c1")),
                                                   tabPanel("Metabolites not Detected", dataTableOutput("missingMets") %>% withSpinner(color="#0dc5c1")))),
                       box(title="Pathway Map", status="primary", solidHeader = TRUE,width = 12,
                           fluidRow(style="padding:20px; height:100px", splitLayout(cellWidths=c("40%", "60%"), align = "left",
                                    selectInput(inputId = "pathwayMapId", label = "Pathway Map",choices = pwy_choices, selected="choose", selectize=FALSE),
                                    sliderInput(inputId = "scalingFactor", label="Node Scaling Factor", min=1, max=5, step=1, value=2))),
                           fluidRow(style="padding:15px; height:110px; ",plotOutput("pmapleg", width = "100%", height = "100px")),
                           fluidRow(width=12, style="padding:10px; height:820px;",visNetworkOutput("pathwayMap")) %>% withSpinner(color="#0dc5c1"), collapsible=TRUE),
                       box(title="Top Perturbed Pathways", status="info", solidHeader=TRUE,  width = 12, collapsible=TRUE,
                           tabsetPanel(type="tabs",tabPanel("Over-representation Analysis", dataTableOutput("oraEnrichment") %>% withSpinner(color="#0dc5c1")),
                                                   tabPanel("Metabolite Set Enrichment Analysis", dataTableOutput("mseaEnrichment")))))), # tabItem ptReport
      tabItem(tabName="ctd", width=12,
              h2("Network-Assisted Diagnostics", align="center"),
              fluidRow(
                box(width=12, title="Select Patient", status="warning", solidHeader = TRUE, splitLayout(cellWidths=c("33%", "33%", "33%"),
                                                                                              selectInput(inputId = "diag_nw_Class", label = "Select diagnosis.",choices = names(cohorts_coded), selected = names(cohorts_coded)[1], selectize=FALSE),
                                                                                              selectInput(inputId = "pt_nw_ID", label = "Select patient.", choices = cohorts_coded[[1]], selected=cohorts_coded[[1]][1], selectize=FALSE, multiple=FALSE),
                                                                                              selectInput(inputId="pvalueType", label="Select method", choices=c("CTD", "CTDdm", "Combined"), selected="Combined", selectize=FALSE)),
                    h4('Click on the cells below to select disease model to interpret patient profile.'),
                    DTOutput('Cohort_pvalRank') %>% withSpinner(color="#0dc5c1")), # box Select Patient
                box(width = 12 ,title = "Network Display", status="info", solidHeader = TRUE, align="left",height="930px", collapsible=FALSE,
                    selectInput(inputId = "bgModel", label = "Select Disease-Specific Background Network.",choices = .GlobalEnv$modelChoices,selected = names(cohorts_coded)[1],selectize = TRUE),
                    div(prettyRadioButtons(inputId = "RangeChoice",label = "Choose range of nodes:",choices = c("Top K perturbed metabolites only", "Abnormal metabolites only", "All Metabolites"),selected = "Top K perturbed metabolites only"),style="display:center-align"),
                    h4(htmlOutput(outputId="selectedPtModel", container = div)),
                    forceNetworkOutput(outputId = "ptNetwork",height = "600px") %>% withSpinner(color="#0dc5c1")) # box Network Display
                #box(title = "genotype", width=NULL, status = "info", solidHeader = TRUE, height = 200)
              ) # fluidRow
      ), # tabItem ctd
      tabItem(tabName="refPop",
              h2("Inspect Reference Population", align="center"),
              fluidRow(box(title = "Inspect the Distribution", status="primary", solidHeader = TRUE, align="left", width=12, collapsible=FALSE,
                           splitLayout(cellWidths=c("33%", "33%", "33%"),
                                       selectInput(inputId="anticoagulant", label="EDTA or Heparin reference population?", choices=c("EDTA", "Heparin"), selected="EDTA", selectize=FALSE),
                                       selectInput(inputId = "metClass", label = "Which metabolite class do you want to select from?",choices = unique(Miller2015$SUPER_PATHWAY), selected="Amino Acid", selectize=FALSE),
                                       selectInput(inputId = "metSelect", label = "Select a metabolite from the chosen class to inspect.", choices = "", selectize=FALSE)),
                           textOutput("estimates"),
                           splitLayout(cellWidths=c("50%", "50%"), plotOutput("referenceReport"), plotOutput("qqplot")),
                           splitLayout(cellWidths=c("50%", "50%"), plotOutput("howRare"), dataTableOutput("refOutliers"))))), # tabItem refPop
      tabItem(tabName="download",
              h2("Download Data", align="center"),
              fluidRow(box(title="Download Data", status="info", solidHeader=TRUE, align="left", width=12, collapsible=FALSE,
                           pickerInput(inputId = "showThese", label = "Diagnoses",choices = names(cohorts_coded)[-which(names(cohorts_coded) %in% c("hep_refs", "edta_refs"))],
                                       selected = names(cohorts_coded)[1],options = list(`actions-box` = TRUE),inline=FALSE,multiple = TRUE),
                           h4(textOutput("st")), downloadButton("downloadButton", "Download"), dataTableOutput("selectedData")))) # tabItem download
    ) # tabItems
  ) # dashboardBody
) # dashboardPage

server = function(input, output, session) {
  observeEvent(input$tab, {
    print(sprintf("%s tab is selected.", input$tab))
    if (input$tab == "ptReport") {
      observeEvent(input$diagClass, priority=1, {
        updatePickerInput(session, "ptIDs", choices = cohorts_coded[[input$diagClass]], selected=cohorts_coded[[input$diagClass]][1])
        report = eventReactive({
          input$diagClass
          input$ptIDs
        },getPatientReport(input))
        output$patientReport = DT::renderDataTable({
          d = report()$patientReport
          DT::datatable(d, rownames=FALSE, options=list(scrollX=TRUE))
        })
        output$missingMets = DT::renderDataTable(report()$missingMets, rownames = FALSE)
        oraDf=eventReactive({
          input$diagClass
          input$ptIDs
        },shiny.getORA_Metabolon(input))
        output$oraEnrichment = renderDataTable({
          datatable(oraDf(),rownames=FALSE, options=list(scrollX=TRUE)) %>%
            formatStyle(c('FDR',"Pvalue"),color = styleInterval(c(0.05),c("red","black")))
        })
        mseaDf=eventReactive({
          input$diagClass
          input$ptIDs
        },getMSEA(input, cohorts_coded))
        output$mseaEnrichment = renderDataTable({
          datatable(mseaDf(),rownames=FALSE, options=list(scrollX=TRUE)) %>%
            formatStyle(colnames(mseaDf)[grepl("val",colnames(mseaDf))],color = styleInterval(c(0.05),c("red","black")))
        })
      })
      updateSelectInput(session,"pathwayMapId", selected="Arginine Metabolism")
      observeEvent(input$pathwayMapId, priority=1, {
        observeEvent(input$scalingFactor, priority=2, {
          pmap = eventReactive({
            input$diagClass
            input$ptIDs
          }, getPathwayMap(input))
          output$pathwayMap = renderVisNetwork({pmap()$pmap})
          output$pmapleg = renderPlot({
            grid.newpage()
            grid.arrange(pmap()$colorbar,pmap()$shapeleg,ncol=2)
          }
          )
        })
      })
    } else if (input$tab == "ctd") {
      observeEvent(input$diag_nw_Class, {
        print(sprintf("nw_Class selected Diagnosis: %s",input$diag_nw_Class))
        print(sprintf("nw_Class selected patient: %s", input$pt_nw_ID))
        print(sprintf("nw_Class selected background graph: %s", input$bgModel))
        updateSelectInput(session, "pt_nw_ID", choices = cohorts_coded[[input$diag_nw_Class]], selected=cohorts_coded[[input$diag_nw_Class]][1])
        updateSelectInput(session, "bgModel", selected=input$diag_nw_Class)
      })
      
      PrankDf=eventReactive({input$diag_nw_Class
        input$pvalueType},
        getPrankDf(input))
      
      print(sprintf(" PrankDf selected patient: %s", input$pt_nw_ID))
      print(sprintf(" Df selected background graph: %s", input$bgModel))
      
      # get ranking table and update cell selection input
      brks = quantile(PrankDf()$df.pranks[PrankDf()$df.pranks<0.05], probs = seq(.05, .95, .05), na.rm = TRUE)
      clrs = round(seq(255, 40, length.out = length(brks) + 1), 0) %>% {paste0("rgb(255,", ., ",", ., ")")}
      clrs = rev(clrs)
      
      output$Cohort_pvalRank = renderDT(datatable(PrankDf()$df.pranks,
                                                  extensions = 'FixedColumns',
                                                  options = list(scrollX = TRUE,fixedColumns = TRUE, #dom = 't',
                                                                 pageLength = 50,
                                                                 lengthMenu = c(10, 25, 50)),
                                                  selection=list(target = 'cell',
                                                                 selected = matrix(c(PrankDf()$model.ind,1),ncol=2),
                                                                 mode = 'single')) %>% 
                                          formatStyle(colnames(PrankDf()$df.pranks),
                                                      backgroundColor = styleInterval(brks,clrs)))
      
      # get patient disease ranking table
      observeEvent(input$Cohort_pvalRank_cells_selected,{
        updateSelectInput(session, "pt_nw_ID", choices = cohorts_coded[[input$diag_nw_Class]], selected = cohorts_coded[[input$diag_nw_Class]][input$Cohort_pvalRank_cells_selected[,2]])
        updateSelectInput(session, "bgModel", choices =.GlobalEnv$modelChoices, selected= rownames(PrankDf()$df.pranks)[input$Cohort_pvalRank_cells_selected[,1]])
        
        print(sprintf("Tb seleted Diagnosis: %s",input$diag_nw_Class))
        print(sprintf("Tb selected patient: %s", input$pt_nw_ID))
        print(sprintf("Tb selected background graph: %s", input$bgModel))
      })
      
      # draw network
      PtResult = reactive({
        shiny::validate(need(try(getVlength(input) != 0), "There are not enough nodes to build network."))
        getPtResult(input)
      })
      output$selectedPtModel=renderText({ paste("Currently viewing", "<font color=\"#FF0000\"><b>",input$diag_nw_Class,"</b></font>",
                                                "patient ","<font color=\"#FF0000\"><b>", input$pt_nw_ID, "</b></font>", 
                                                "in disease model ","<font color=\"#FF0000\"><b>",input$bgModel, "</b></font>",".") })
      output$ptNetwork=renderForceNetwork(PtResult()$ptNetwork)
      
      print(sprintf("selected patient: %s", input$pt_nw_ID))
      print(sprintf("seleted background graph: %s", input$bgModel))
    } else if (input$tab == "refPop") {
      observeEvent(input$metClass, priority=1, {
        observeEvent(input$anticoagulant, priority=0, {
          ch = getMetList(input)
          updateSelectInput(session, "metSelect", choices = ch, selected=ch[1])
          print("metSelect dropdown should be updated now.")
        })
      })
      ref = reactive(getRefPop(input))
      output$estimates = renderText({sprintf("Mean Estimate = %.2f\nStandard Deviation Estimate = %.2f", ref()$ests$mean, ref()$ests$std)})
      output$referenceReport = renderPlot(ref()$hst)
      output$qqplot = renderPlot(ref()$qq)
      output$howRare = renderPlot(ref()$rare)
      output$refOutliers = renderDataTable(ref()$outliers)
    } else if (input$tab == "download") {
      observeEvent(input$showThese, priority = 0, {
        dd = getData(input)
        output$downloadButton = downloadHandler(
          filename = function() { paste(paste(input$showThese, collapse="_"), ".txt", sep="") },
          content = function(file) { write.table(dd, file, sep="\t", col.names = TRUE, row.names = FALSE) }
        )
        output$selectedData = DT::renderDataTable({DT::datatable(dd, rownames=FALSE, options=list(scrollX=TRUE))})
      })
      output$st = renderText({sprintf("Selected Cohort: %s",paste(input$showThese,collapse = ", "))})
    } else {
      print("No tab selected")
    }
  })
}

shinyApp(ui, server)

