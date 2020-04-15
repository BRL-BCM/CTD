require(shiny)
require(DT)
require(shinydashboard)
require(ggplot2)
require(grid)
require(Hmisc)
require(plotly)
require(grDevices)
require(igraph)
require(networkD3)
require(R.utils)
require(CTD)
setwd("/Users/lillian.rosa/Downloads/CTD/inst/shiny-app")
source("/Users/lillian.rosa/Downloads/CTD/inst/shiny-app/metDataPortal_appFns.r")

ui = dashboardPage(
  dashboardHeader(title = "Metabolomics Data Portal"),
  dashboardSidebar(sidebarMenu(id = "tab",
                               menuItem("View Patient Report", tabName = "ptReport", icon = icon("user-circle-o")),
                               menuItem("Inspect Reference Population", tabName = "refPop", icon = icon("bar-chart")),
                               menuItem("Network-Assisted Diagnostics", tabName = "ctd", icon=icon("project-diagram")))),
  dashboardBody(
                tabItems(
                  tabItem(tabName="ptReport",
                          fluidRow(h2("Patient Report", align="center"),
                          box(title="Select Patient(s)", status="warning", solidHeader = TRUE,
                              splitLayout(cellWidths=c("25%", "75%"),
                                       selectInput(inputId = "diagClass", label = "Select diagnosis.", 
                                                   choices = names(cohorts_coded), selected = names(cohorts_coded)[1], selectize=FALSE),
                                       selectInput(inputId = "ptIDs", label = "Select patients.", choices = "", selectize=TRUE, multiple=TRUE)), width=4),
                          box(title="Top Perturbed Pathways", status="info", solidHeader=TRUE, 
                              tabsetPanel(type="tabs", 
                                          tabPanel("Over-representation Analysis", dataTableOutput("oraEnrichment")),
                                          tabPanel("Metabolite Set Enrichment Analysis", dataTableOutput("mseaEnrichment"))), width=8, collapsible=TRUE),
                          box(title = "Patient Report", status="info", solidHeader = TRUE,
                              #downloadButton("downloadPatientReport", "Download Patient Report"),
                              splitLayout(cellWidths=c("60%", "40%"), dataTableOutput("patientReport"), dataTableOutput("missingMets")),
                              align="left", width=12, collapsible=TRUE),
                          box(title="Pathway Map", status="primary", solidHeader = TRUE,
                              fluidRow(style="padding:10px; height:80px;", 
                                       splitLayout(cellWidths=c("20%", "40%", "40%"),
                                                   selectInput(inputId = "pathwayMapId", label = "Pathway Map", 
                                                               choices = c("Choose", "Arginine Metabolism", "Ascorbate Metabolism", "Asp-Glu Metabolism", "BCAA Metabolism", 
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
                                                                           "Tryptophan Metabolism", "All"), selectize=FALSE),
                                                   sliderInput(inputId = "scalingFactor", label="Node Scaling Factor", min=1, max=5, step=1, value=3),
                                                   plotOutput("colorbar"))),
                              fluidRow(splitLayout(cellWidths=c("15%", "85%"),
                                                   textOutput("dummy"),
                                                   imageOutput("pathwayMap", height="100%", width="100%"))), width=12, collapsible=TRUE)
                          )),
                  tabItem(tabName="refPop",
                          h2("Inspect Reference Population", align="center"),
                          fluidRow(box(title = "Inspect the Distribution", status="primary", solidHeader = TRUE,
                                       splitLayout(cellWidths=c("33%", "33%", "33%"),
                                          selectInput(inputId="anticoagulant", label="EDTA or Heparin reference population?", choices=c("EDTA", "Heparin"), selected="EDTA", selectize=FALSE),
                                          selectInput(inputId = "metClass", label = "Which metabolite class do you want to select from?",
                                                      choices = unique(Miller2015$SUPER_PATHWAY), selected="Amino Acid", selectize=FALSE),
                                          selectInput(inputId = "metSelect", label = "Select a metabolite from the chosen class to inspect.", choices = "", selectize=FALSE)),
                                       textOutput("estimates"),
                                       splitLayout(cellWidths=c("50%", "50%"), plotOutput("referenceReport"), plotOutput("qqplot")),
                                       splitLayout(cellWidths=c("50%", "50%"), plotOutput("howRare"), dataTableOutput("refOutliers")),
                                       align="left", width=12, collapsible=TRUE)
                                   ),
                          fluidRow(box(title="Download Data", status="info", solidHeader=TRUE,
                                       splitLayout(cellWidths=c("85%", "15%"),
                                                   checkboxGroupInput(inputId = "showThese", label = "Diagnoses", 
                                                                      choices = names(cohorts_coded)[-which(names(cohorts_coded) %in% c("hep_refs", "edta_refs"))], 
                                                                      selected = names(cohorts_coded)[1], inline=TRUE),
                                                   selectInput(inputId = "raworZscore", label = "Data Processing Level", choices = list("Raw", "Normalized", "Zscored"), 
                                                               selected = "Zscored", selectize = FALSE)),
                                       textOutput("st"), textOutput("rz"),
                                       downloadButton("downloadButton", "Download"), dataTableOutput("selectedData"),
                                       align="left", width=12, collapsible=TRUE)
                                   )
                          ),
                  tabItem(tabName="ctd",
                          h2("Network-Assisted Diagnostics"),
                          box(title="Select Patient", status="warning", solidHeader = TRUE,
                              splitLayout(cellWidths=c("25%", "50%"),
                                          selectInput(inputId = "diag_nw_Class", label = "Select diagnosis.", 
                                                      choices = names(cohorts_coded), selected = names(cohorts_coded)[1], selectize=FALSE),
                                          selectInput(inputId = "pt_nw_ID", label = "Select patient.", choices = cohorts_coded[[1]], selected=cohorts_coded[[1]][1], selectize=FALSE, multiple=FALSE),
                                          selectInput(inputId="kmx", label="Top K Metabolites", choices=seq(5,50), selected=30)
                              ), width=12),
                          fluidRow(
                            column(width = 4,
                                   #box(title = "genotype", width=NULL, status = "info", solidHeader = TRUE, height = 200),
                                   box(title = "P-value Rankings",width=NULL,status = "info", solidHeader = TRUE, height = 700,
                                       align = "left",
                                       tableOutput("pvalRank"))
                            ),
                            column(width = 8 ,
                                   box(title = "Network", status="info", solidHeader = TRUE, width = NULL,
                                       selectInput(inputId = "bgModel", label = "Select Disease-Specific Background Network.",
                                                   choices = names(cohorts_coded),
                                                   selected = names(cohorts_coded)[1],selectize = TRUE),
                                       h4(htmlOutput(outputId="selectedPtModel", container = div)),
                                       forceNetworkOutput(outputId = "ptNetwork"),
                                       #h5(htmlOutput("pctd")),
                                       align="left", 
                                       height=700, collapsible=FALSE))
                          ))
                        )
  )
)

server = function(input, output, session) {
  observe({
    print(sprintf("%s tab is selected.", input$tab))
  })

  observeEvent(input$tab, {
    if (input$tab == "ptReport") {
      observeEvent(input$diagClass, priority=1, {
        print(input$diagClass)
        updateSelectInput(session, "ptIDs", choices = cohorts_coded[[input$diagClass]], selected=cohorts_coded[[input$diagClass]][1])

        report = reactive(getPatientReport(input))
        output$patientReport = DT::renderDataTable({
          d = report()$patientReport
          DT::datatable(d, rownames=FALSE, options=list(scrollX=TRUE))
        })
        output$missingMets = DT::renderDataTable(report()$missingMets, rownames = FALSE)
        #output$downloadPatientReport <- downloadHandler(
        #  filename = function() { paste(input$biofluid, "-", input$patientID, ".txt", sep="") },
        #  content = function(file) { write.table(report()$patientReport, file, sep="\t", col.names = TRUE, row.names = FALSE) }
        #)
        output$oraEnrichment = renderDataTable(shiny.getORA_Metabolon(input), rownames=FALSE)
        output$mseaEnrichment = renderDataTable(shiny.getMSEA_Metabolon(input, cohorts_coded), rownames=FALSE)
        #output$geneticVars = DT::renderDataTable(getGeneticVariants(input), rownames = FALSE)

        observeEvent(input$pathwayMapId, priority=0, {
          print(input$pathwayMapId)
          observeEvent(input$scalingFactor, priority=-1, {
            pmap = reactive(isolate(getPathwayMap(input)))
            output$pathwayMap = renderImage({pmap()$pmap})
            output$colorbar = renderPlot({
              grid.newpage()
              grid.layout(nrow = 1, ncol = 1, just = c("right", "top"))
              grid.draw(pmap()$colorbar)
            }, height = 20)
          })
        })
      })
    } else if (input$tab == "refPop") {
      observeEvent(input$metClass, priority = 1, {
        ch = getMetList(input)
        updateSelectInput(session, "metSelect", choices = ch, selected=ch[1])
        print("metSelect dropdown should be updated now.")
      })
      ref = reactive(getRefPop(input))
      output$estimates = renderText({sprintf("Mean Estimate = %.2f\nStandard Deviation Estimate = %.2f", ref()$ests$mean, ref()$ests$std)})
      output$referenceReport = renderPlot(ref()$hst)
      output$qqplot = renderPlot(ref()$qq)
      output$howRare = renderPlot(ref()$rare)
      output$refOutliers = renderDataTable(ref()$outliers)

      observeEvent(input$showThese, priority = 0, { 
        observeEvent(input$raworZscore, priority = -1, {
          dd = getData(input)
          output$downloadButton = downloadHandler(
            filename = function() { paste(paste(input$showThese, collapse="_"), "-", input$raworZscore, ".txt", sep="") },
            content = function(file) { write.table(dd, file, sep="\t", col.names = TRUE, row.names = FALSE) }
          )
          output$selectedData = DT::renderDataTable({DT::datatable(dd, rownames=FALSE, options=list(scrollX=TRUE))})
        })
      })
      output$st = renderText({input$showThese})
      output$rz = renderText({input$raworZscore})
    } else if (input$tab == "ctd") {
        observeEvent(input$diag_nw_Class, priority=1, {
          print(sprintf("selected Diagnosis: %s",input$diag_nw_Class))
          updateSelectInput(session, "pt_nw_ID", choices = cohorts_coded[[input$diag_nw_Class]], selected=cohorts_coded[[input$diag_nw_Class]][1])
          updateSelectInput(session, "bgModel", choices = names(cohorts_coded), selected=input$diag_nw_Class)
          print(sprintf("selected patient: %s", input$pt_nw_ID))
          print(sprintf("selected background graph: %s", input$bgModel))
          output$selectedPtModel=renderText({ paste("Currently viewing", "<font color=\"#FF0000\"><b>",input$diag_nw_Class,"</b></font>",
                                                    "patient ","<font color=\"#FF0000\"><b>", input$pt_nw_ID, "</b></font>", 
                                                    "in disease model ","<font color=\"#FF0000\"><b>",input$bgModel, "</b></font>",".") })
          PtResult = reactive({getPtResult(input)})
          PrankDf=reactive({getPrankDf(input)})
          output$pvalRank = renderTable(PrankDf())
          output$ptNetwork=renderForceNetwork(PtResult()[["ptNetwork"]])
        })
    } else {
      print("No tab selected")
    }
  })



}

shinyApp(ui, server)
