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
kmx=15
setwd("/Users/lillian.rosa/Downloads/CTD/vignette/")
source("/Users/lillian.rosa/Downloads/CTD/vignette/shiny-app/metDataPortal_appFns.r")
data("Miller2015")
cohorts = list()
cohorts$mcc = diagnoses$id[which(diagnoses$diagnosis=="3-methylcrotonyl CoA carboxylase")]
cohorts$arg = diagnoses$id[which(diagnoses$diagnosis=="Argininemia")]
cohorts$cit = diagnoses$id[which(diagnoses$diagnosis=="Citrullinemia")]
cohorts$cob = diagnoses$id[which(diagnoses$diagnosis=="Cobalamin biosynthesis")]
cohorts$ga = diagnoses$id[which(diagnoses$diagnosis=="Glutaric Aciduria")]
cohorts$gamt = diagnoses$id[which(diagnoses$diagnosis=="Guanidinoacetate methyltransferase")]
cohorts$msud = diagnoses$id[which(diagnoses$diagnosis=="Maple syrup urine disease")]
cohorts$mma = diagnoses$id[which(diagnoses$diagnosis=="Methylmalonic aciduria")]
cohorts$otc = diagnoses$id[which(diagnoses$diagnosis=="Ornithine transcarbamoylase")]
cohorts$pa = diagnoses$id[which(diagnoses$diagnosis=="Propionic aciduria")]
cohorts$pku = diagnoses$id[which(diagnoses$diagnosis=="Phenylketonuria")]
cohorts$tmhle = diagnoses$id[which(diagnoses$diagnosis=="Trimethyllysine hydroxylase epsilon")]

ui = dashboardPage(
  dashboardHeader(title = "Metabolomics Data Portal"),
  dashboardSidebar(sidebarMenu(id = "tab",
                               menuItem("View Patient Report", tabName = "ptReport", icon = icon("user-circle-o")),
                               menuItem("Inspect Reference Population", tabName = "refPop", icon = icon("bar-chart")),
                               menuItem("Examine Patient Similarity", tabName="similarity", icon= icon("search")))),
  dashboardBody(height="100%",
                tabItems(
                  tabItem(tabName="ptReport",
                          h2("Patient Report"),
                          fluidRow(box(title="Select Patient(s)", status="warning", solidHeader = TRUE,
                                       selectInput(inputId = "diagClass", label = "Select diagnosis.", choices = names(cohorts), selected = names(cohorts)[1]),
                                       checkboxGroupInput(inputId = "ptIDs", label = "Select patients.", choices = ""),
                                       align="left", width=2),
                                   box(title="Pathway Map", status="primary", solidHeader = TRUE,
                                       splitLayout(cellWidths=c("33%", "33%", "33%"),
                                                   selectInput(inputId = "pathwayMapId", label = "Pathway Map", choices = ""),
                                                   sliderInput(inputId = "scalingFactor", label="Node Scaling Factor", min=1, max=5, step=1, value=1),
                                                   plotOutput("colorbar")),
                                       imageOutput("pathwayMap"),
                                       align="left", width=10, collapsible=TRUE)
                                   ),
                          box(title = "Patient Report", status="info", solidHeader = TRUE,
                              downloadButton("downloadPatientReport", "Download Patient Report"),
                              splitLayout(cellWidths=c("60%", "40%"), dataTableOutput("patientReport"), dataTableOutput("missingMets")),
                              align="left", width=12, collapsible=TRUE),
                          fluidRow(box(title="Top Perturbed Pathways", status="info", solidHeader=TRUE, width=4, collapsible=TRUE), #tableOutput("pathwayEnrichment")
                                   box(title="Inspect Genetic Variants", status="info", solidHeader=TRUE, width=8, collapsible=TRUE)) #dataTableOutput("geneticVars")
                          ),
                  tabItem(tabName="refPop",
                          h2("Inspect Reference Population"),
                          fluidRow(box(title = "Inspect the Distribution", status="primary", solidHeader = TRUE,
                                       splitLayout(cellWidths=c("50%", "50%"),
                                          selectInput(inputId = "metClass", label = "Which metabolite class do you want to select from?",
                                                      choices = sort(unique(.GlobalEnv$metClass)), selected="Amino Acid"),
                                          selectInput(inputId = "metSelect", label = "Select a metabolite from the chosen class to inspect.", choices = "")),
                                       textOutput("estimates"),
                                       splitLayout(cellWidths=c("50%", "50%"), plotOutput("referenceReport"), plotOutput("qqplot")),
                                       splitLayout(cellWidths=c("50%", "50%"), plotOutput("howRare"), dataTableOutput("refOutliers")),
                                       align="left", width=12, collapsible=TRUE)
                                   ),
                          fluidRow(box(title="Download Data", status="info", solidHeader=TRUE,
                                       splitLayout(cellWidths=c("85%", "15%"),
                                                   checkboxGroupInput(inputId = "showThese", label = "Diagnoses", choices = names(cohorts), selected = names(cohorts)[1], inline=TRUE),
                                                   selectInput(inputId = "raworZscore", label = "Data Processing Level", choices = list("Raw", "Normalized", "Zscored"), selected = "Zscored")),
                                       textOutput("st"), textOutput("rz"),
                                       downloadButton("downloadButton", "Download"), dataTableOutput("selectedData"),
                                       align="left", width=12, collapsible=TRUE)
                                   )
                          ),
                  tabItem(tabName="similarity",
                          h2("Examine Patient Similarity"),
                          fluidRow(tabBox(title="", id="diagOrExtract",
                                          tabPanel("Multi-dimensional Scaling",
                                                   splitLayout(cellWidths=c("33%", "33%", "33%"),
                                                               selectInput(inputId = "diagnosis", label="Select Disease-Specific Knowledge Graph", choices = names(cohorts), selected = names(cohorts)[1]),
                                                               selectInput(inputId = "dim", label="#Dimensions", choices = c(2,3), selected=3),
                                                               textOutput("print_diagnosis")),
                                                   splitLayout(cellWidths=c("60%", "40%"), dataTableOutput("algSig_dscore"), plotlyOutput("mds"))),
                                          tabPanel("Modular Feature Extraction",
                                                   splitLayout(cellWidths=c("25%", "25%", "25%", "25%"),
                                                               selectInput(inputId = "diagnosis", label="Select Disease-Specific Knowledge Graph", choices = names(cohorts), selected = names(cohorts)[1]),
                                                               selectInput(inputId = "ptID", label="Select Patient 1", choices=NULL),
                                                               selectInput(inputId = "ptID2", label="Select Patient 2", choices=NULL),
                                                               selectInput(inputId = "kmax", label="Select k", choices = c(2:kmx))),
                                                   splitLayout(cellWidths=c("60%", "40%"), dataTableOutput("sim"), forceNetworkOutput("pt_cmp"))), width=12))
                          )
                )
              )
)

server = function(input, output, session) {
  observe({
    print(sprintf("%s tab is selected.", input$tab))
  })

  observeEvent(input$tab, {
    if (input$tab == "ptReport") {
      # Pathway Analysis code
      observeEvent(input$diagClass, priority=1, {
        print(input$diagClass)
        updateCheckboxGroupInput(session, "ptIDs", choices = cohorts[[input$diagClass]], selected = cohorts[[input$diagClass]])

        report = reactive(getPatientReport(input))
        output$patientReport = DT::renderDataTable({
          d = report()$patientReport
          DT::datatable(d, rownames=FALSE, options=list(scrollX=TRUE))
        })
        output$missingMets = DT::renderDataTable(report()$missingMets, rownames = FALSE)
        output$downloadPatientReport <- downloadHandler(
          filename = function() { paste(input$biofluid, "-", input$patientID, ".txt", sep="") },
          content = function(file) { write.table(report()$patientReport, file, sep="\t", col.names = TRUE, row.names = FALSE) }
        )
        #output$pathwayEnrichment = renderTable(getPathwayEnrichment2(input))
        #output$geneticVars = DT::renderDataTable(getGeneticVariants(input), rownames = FALSE)

        observeEvent(input$pathwayMapId, priority=0, {
          print(input$pathwayMapId)
          updateSelectInput(session, "pathwayMapId", choices = c("All", "Arginine Metabolism", "Ascorbate Metabolism", "Asp-Glu Metabolism",
                                        "BCAA Metabolism", "Benzoate Metabolism", "Beta-Oxidation", "Bile-Acid Metabolism",
                                        "Carnitine Biosynthesis", "Cholesterol Synthesis", "Creatine Metabolism", "Dicoarboxylic Acid Metabolism",
                                        "Eicosanoids", "Endocannabinoid Synthesis", "Fatty Acid Metabolism", "Fibrinogen Cleavage Peptides",
                                        "GABA Shunt", "Galactose Metabolism", "Glutathione Metabolism", "Gly-Ser-Thr Metabaolism", "Glycogen Metabolism",
                                        "Glycolysis", "Glycosylation", "Hemoglobin-Porphyrin Metabolism", "Histidine Metabolism", "Inositol Metabolism",
                                        "Ketone Bodies", "Lysine Catabolism", "Met-Cys Metabolism", "Mevalonate Metabolism", "Nicotinate-Nicotinamide Metabolism",
                                        "Pantothenate Metabolism", "Pentose-Phosphate Metabolism", "Phe-Tyr Metabolism", "Phospholipid Metabolism", "Polyamine Metabolism",
                                        "Proline Metabolism", "Protein Degradation", "Purine Metabolism", "Pyridoxal Metabolism", "Pyrimidine Metabolism",
                                        "Riboflavin Metabolism", "Secondary-Bile-Acids", "Sorbitol-Glycerol Metabolism", "Sphingolipid-Metabolism",
                                        "Steroid-Hormone Biosynthesis", "TCA Cycle", "Thyroid Hormone Synthesis", "Tryptophan Metabolism"),
                            selected="Arginine Metabolism")
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
      ref = reactive(getRefPop(input, .GlobalEnv$all_norm_data))
      output$estimates = renderText({sprintf("Mean Estimate = %.2f\nStandard Deviation Estimate = %.2f", ref()$ests$mean, ref()$ests$std)})
      output$referenceReport = renderPlot(ref()$hst)
      output$qqplot = renderPlot(ref()$qq)
      output$howRare = renderPlot(ref()$rare)
      output$refOutliers = renderDataTable(ref()$outliers)

      observeEvent(c(input$raworZscore, input$showThese), priority = -1, {
        data = Miller2015[,grep("IEM_", colnames(Miller2015))]
        output$downloadButton = downloadHandler(
          filename = function() { paste(paste(input$showThese, collapse="_"), "-", input$raworZscore, ".txt", sep="") },
          content = function(file) { write.table(data, file, sep="\t", col.names = TRUE, row.names = FALSE) }
        )
        output$selectedData = DT::renderDataTable({DT::datatable(data, rownames=FALSE, options=list(scrollX=TRUE))})
      })
      output$st = renderText({input$showThese})
      output$rz = renderText({input$raworZscore})
    } else if (input$tab == "similarity") {
      data = eval(parse(text=sprintf("graphs$%s$data", input$diagnosis)))
      updateSelectInput(session, "ptID", label="Select a patient to diagnose", choices = colnames(data), selected=colnames(data)[1])
      updateSelectInput(session, "ptID2", label="Select a patient to diagnose", choices = colnames(data), selected=colnames(data)[2])

      output$print_diagnosis = renderText(input$diagnosis)
      output$algSig_dscore = renderDataTable(comparePatientModPerts(input), rownames=FALSE)
      output$mds = renderPlotly(getMDS(input))
      res = reactive(extractModPerts(input))
      output$sim = renderDataTable(res()$sim, rownames=FALSE)
      output$pt_cmp = renderForceNetwork(res()$pt_ig)
    } else {
      print("No tab selected")
    }
  })



}

shinyApp(ui, server)
