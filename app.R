#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
# R script to build Shiny APP showing ensemble model procedure and output for Sole GSA17
# @author Francesco Masnadi
# @email francesco.masnadi@outlook.it
#
# This app has been developed as the Ph.D. thesis of Francesco Masnadi 
# Ph.D. Program "Innovative Technologies and Sustainable Use of Mediterranean Sea Fishery and Biological Resources" (www.FishMed-PhD.org). 
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>

# Main source: https://rstudio-education.github.io/shiny-course/

library(readr)
library(shiny)
library(readxl)
library(dplyr)

# Data pre-processing ----
run_list <- read_excel("single_run/run_list.xlsx", sheet = "list") %>% select(-Weighting)
table_diags <- read.csv2("ensemble/Diags_table.csv")
# ensemble output to be download
ensF <- read.csv2("ensemble/Ensemble_F.csv")
ensSSB <- read.csv2("ensemble/Ensemble_SSB.csv")
ensRec <- read.csv2("ensemble/Ensemble_Recr.csv")
ensF_Ftrg <- read.csv2("ensemble/Ensemble_F_Ftrg.csv")
ensSSB_SSBtrg <- read.csv2("ensemble/Ensemble_SSB_SSBtrg.csv")
tableICES <-   read.csv2("forecast/Forecast_table_ICES.csv")
  
  
# Define UI  ----
ui <- fluidPage(titlePanel("Topic 2 - Stock assessment of common sole in GSA17: an ensemble approach using data-rich model"),
                
                tabsetPanel(type = "tabs",
                            tabPanel("Workflow",
                                     plotOutput("workf",  inline = T) 
                            ),
                            tabPanel("Single runs", 
                                     
                                     # App title ----
                                     titlePanel("Single Runs Outputs and Diags"),
                                     
                                     # Sidebar layout with input and output definitions ----
                                     sidebarLayout(
                                       
                                       # Sidebar panel for inputs ----
                                       sidebarPanel(
                                         
                                         # Input: Selector for run to plot  ----
                                         selectInput("run", "Select run from ensemble grid:", choices =  run_list$Name, multiple = F),
                                       ),   #end sidebarPanel
                                       
                                       # Main panel for displaying outputs ----
                                       mainPanel(
                                         
                                         tabsetPanel(type = "tabs",
                                                     tabPanel("SA Summary", 
                                                              # Output: Plot of the SA summary ----
                                                              # Output: Formatted text for caption ----
                                                              h3(textOutput("caption")),
                                                              plotOutput("runsPlotsSSB",  inline = T),
                                                              plotOutput("runsPlotsF",  inline = T),
                                                              plotOutput("runsPlotsRec",  inline = T),
                                                              plotOutput("runsPlotsCurve",  inline = T),
                                                              plotOutput("runsPlotsCATCH",  inline = T)
                                                     ),  #end tabPanel SA summary
                                                     tabPanel("Diagnostics", 
                                                              # Output: Plot of the diags ----
                                                              # Output: Formatted text for caption ----
                                                              h3("Joint Residulas"),
                                                              plotOutput("diagsJres",  inline = T),
                                                              h3("Run Test"),
                                                              plotOutput("diagsRun",  inline = T),
                                                              h3("Retrospective"),
                                                              plotOutput("diagsRetro",  inline = T),
                                                              h3("Hindcasting CPUE"),
                                                              plotOutput("diagsHCxCPUE",  inline = T),
                                                              h3("Hindcasting Length"),
                                                              plotOutput("diagsHCxLEN",  inline = T)
                                                     )  #end tabPanel Diags
                                         )  # tabsetPanel 2
                                       )   #end mainPanel
                                     ),    #end sidebarLayout
                                     # add table on the page
                                     fluidRow(column(12,tableOutput('table'))
                                     )  #end fluidRow table
                            ), #end of tabPanel Single run output
                            
                            tabPanel("Model weighting",
                                     h4("- Summary table of the diagnostics used in the weighting procedure:"),
                                     tableOutput('tablediags'),
                                     downloadButton("downloadDataDT", "Download diagnostic table")
                            ), #end model weighting
                            
                            tabPanel("Ensemble",
                                     titlePanel("Ensemble Output"),
                                     # Output: Plot of the ensemble----         
                                     fluidRow(column(12,offset=1, h3("Single run Trajectories Comparison 
                          ---------------------------------------------------> 
                                                            Final Ensemble Trajectories")),       
                                              column(6, plotOutput("ensetrj_sing",  inline = T)),
                                              # column(6, h3("Final Ensemble Trajectories",  inline = T)),       
                                              column(6, plotOutput("ensetrj_all",  inline = T)),
                                              # Button
                                              fluidRow(column(12 , offset = 2,downloadButton("downloadDataF", "Download F timeseries"),
                                                              downloadButton("downloadDataSSB", "Download SSB timeseries"),
                                                              downloadButton("downloadDataRec", "Download Rec timeseries"),
                                                              downloadButton("downloadDataF_Ftrg", "Download F_Ftrg timeseries"),
                                                              downloadButton("downloadDataSSB_SSBtrg", "Download SSB_SSBtrg timeseries"))),
                                              column(12,offset = 5, h3("Ensemble kobe")),
                                              column(12,offset = 3, plotOutput("kobe",  inline = T))
                                     ), 
                                     
                            ),#end tabPanel Ensemble
                            
                            tabPanel("Forecast",
                                     mainPanel(tabsetPanel(type = "tabs",
                                                           tabPanel("Short-term projection",      
                                                           fluidRow(column(6, plotOutput("MLVN_Ensemble_forecast",  inline = T)),
                                                                    column(3, offset= 3,plotOutput("Kobe_Ensemble_forecast",  inline = T)))),
                                                           tabPanel("Probabilities",
                                                                    plotOutput("forecast_percent_FminorFtarget",  inline = T),
                                                                    plotOutput("forecast_percent_BmajorBtarget",  inline = T),
                                                                    plotOutput("forecast_percent_BmajorBlim",  inline = T),  
                                                           downloadButton("downloadICES", "Download ICES"))
                                     )),
                            ),     #end tabPanel forecast                  
                ), #end tabsetPanel 1
) #end fluidPage




# Define server logic to plot various runs ----
server <- function(input, output,session) {
  
  #############################
  # Introduction              #
  #############################
  # Plot F
  output$workf <- renderImage({
    filename <-  file.path(paste("./ensemble/", sep=''), paste("wf.jpg" , sep=''))
    # Return a list containing the filename and alt text
    list(src = filename,
         width = 1300,
         height = 630,
         alt = "z"
    )
  }, deleteFile = F)
  
  
  
  #############################
  # Single run output         #
  #############################
  # show run list
  output$table <- renderTable(run_list)
  
  # Return the requested dataset ----
  datasetInput <- reactive({
    (input$run_list)
  })
  
  # Compute the formula text ----
  # This is in a reactive expression since it is shared by the
  # output$caption and output$runsPlots functions
  formulaText <- reactive({
    paste("Summary output of ", input$run)
  })
  
  # Return the formula text for printing as a caption ----
  output$caption <- renderText({
    formulaText()
  })
  
  # Generate a plot for SA summary----
  # Plot SSB
  output$runsPlotsSSB<- renderImage({
    filename <-  file.path(paste("./single_run/", input$run,  sep=''),'ts7_Spawning_biomass_(mt)_with_95_asymptotic_intervals_intervals.png' )
    # Return a list containing the filename and alt text
    list(src = filename,
         width = 400,
         height = 300,
         alt = "tot SBB"
    )
  }, deleteFile = F)
  
  # Plot F
  output$runsPlotsF <- renderImage({
    filename <-  file.path(paste("./single_run/", input$run,  sep=''),'ts_summaryF.png' )
    # Return a list containing the filename and alt text
    list(src = filename,
         width = 400,
         height = 300,
         alt = "Total F"
    )
  }, deleteFile = F)
  
  #plot Recr
  output$runsPlotsRec<- renderImage({
    filename <-  file.path(paste("./single_run/", input$run,  sep=''),'ts11_Age-0_recruits_(1000s)_with_95_asymptotic_intervals.png' )
    # Return a list containing the filename and alt text
    list(src = filename,
         width = 400,
         height = 300,
         alt = "tot SBB"
    )
  }, deleteFile = F)
  
  
  #plot Curve
  output$runsPlotsCurve <- renderImage({
    filename <-  file.path(paste("./single_run/", input$run,  sep=''),'yield2_yield_curve_with_refpoints.png' )
    # Return a list containing the filename and alt text
    list(src = filename,
         width = 400,
         height = 300,
         alt = "refp"
    )
  }, deleteFile = F)
  
  #plot Catch
  output$runsPlotsCATCH <- renderImage({
    filename <-  file.path(paste("./single_run/", input$run,  sep=''),'catch1 landings.png' )
    # Return a list containing the filename and alt text
    list(src = filename,
         width = 400,
         height = 300,
         alt = "Landings by fleet"
    )
  }, deleteFile = F)
  
  #
  # Generate a plot for Diags----
  # Plot Hind CPUE
  output$diagsHCxCPUE<- renderImage({
    filename <-  file.path(paste("./single_run/", input$run, "/Plotdiags", sep=''), paste("HCxval_CPUE_", input$run,"_update.jpg" , sep=''))
    list(src = filename,
         width = 500,
         height = 500,
         alt = "aa"
    )
  }, deleteFile = F)
  # plot Hind Length
  output$diagsHCxLEN<- renderImage({
    filename <-  file.path(paste("./single_run/", input$run, "/Plotdiags", sep=''), paste("HCxval_length_", input$run,"_update.jpg" , sep=''))
    list(src = filename,
         width = 800,
         height = 800,
         alt = "aa"
    )
  }, deleteFile = F)
  # plot joint Residual
  output$diagsJres<- renderImage({
    filename <-  file.path(paste("./single_run/", input$run, "/Plotdiags", sep=''), paste("JointResiduals_", input$run,"_update.jpg" , sep=''))
    list(src = filename,
         width = 800,
         height = 400,
         alt = "aa"
    )
  }, deleteFile = F)
  # plot Run test
  output$diagsRun<- renderImage({
    filename <-  file.path(paste("./single_run/", input$run, "/Plotdiags", sep=''), paste("RunTestResidual_", input$run,"_update.jpg" , sep=''))
    list(src = filename,
         width = 800,
         height = 800,
         alt = "aa"
    )
  }, deleteFile = F)
  
  # plot Retro
  output$diagsRetro<- renderImage({
    filename <-  file.path(paste("./single_run/", input$run, "/Plotdiags", sep=''), paste("Retro_", input$run,"_update.jpg" , sep=''))
    list(src = filename,
         width = 800,
         height = 400,
         alt = "aa"
    )
  }, deleteFile = F)
  
  ####################
  # Model weighting         #
  ####################
  # table diags
  output$tablediags <- renderTable(table_diags, width = "100%")
  output$downloadDataDT <- downloadHandler(
    filename = function() { 
      paste("Diagnostic_table.csv", sep="")
    },
    content = function(file) {
      write.csv(table_diags, file)
    })
  
  ####################
  # Ensemble         #
  ####################
  # plot Stock status trajectories single run
  output$ensetrj_sing <- renderImage({
    filename <-  file.path(paste("./ensemble/", sep=''), paste("MLVN_Compare.jpg" , sep=''))
    list(src = filename,
         width = 700,
         height = 550,
         alt = "-"
    )
  }, deleteFile = F)
  
  # plot ALL trajectories 
  output$ensetrj_all <- renderImage({
    filename <-  file.path(paste("./ensemble/", sep=''), paste("MLVN_All_wSole.jpg" , sep=''))
    list(src = filename,
         width = 700,
         height = 550,
         alt = "-"
    )
  }, deleteFile = F)
  
  # plot kobe
  output$kobe <- renderImage({
    filename <-  file.path(paste("./ensemble/", sep=''), paste("Kobe_final_wSolea.jpg" , sep=''))
    list(src = filename,
         width = 700,
         height = 550,
         alt = "-"
    )
  }, deleteFile = F)
  
  
  # Downloadable csv of selected dataset ----
  output$downloadDataF <- downloadHandler(
    filename = function() { 
      paste("Ensemble_F.csv", sep="")
    },
    content = function(file) {
      write.csv(ensF, file)
    })
  output$downloadDataSSB <- downloadHandler(
    filename = function() { 
      paste("Ensemble_SSB.csv", sep="")
    },
    content = function(file) {
      write.csv(ensSSB, file)
    })
  output$downloadDataRec <- downloadHandler(
    filename = function() { 
      paste("Ensemble_Rec.csv", sep="")
    },
    content = function(file) {
      write.csv(ensRec, file)
    })
  output$downloadDataF_Ftrg <- downloadHandler(
    filename = function() { 
      paste("Ensemble_F_Ftrg.csv", sep="")
    },
    content = function(file) {
      write.csv(ensF_Ftrg, file)
    })
  output$downloadDataSSB_SSBtrg <- downloadHandler(
    filename = function() { 
      paste("Ensemble_SSB_SSBtrg.csv", sep="")
    },
    content = function(file) {
      write.csv(ensSSB_SSBtrg, file)
    })
  
  ####################
  # Forecast         #
  ####################
  # plot trj
  output$MLVN_Ensemble_forecast <- renderImage({
    filename <-  file.path(paste("./forecast/", sep=''), paste("MLVN_Ensemble_forecast.jpg" , sep=''))
    list(src = filename,
         width = 750,
         height = 580,
         alt = "-"
    )
  }, deleteFile = F)
  # plot kobe
  output$Kobe_Ensemble_forecast <- renderImage({
    filename <-  file.path(paste("./forecast/", sep=''), paste("Kobe_Ensemble_forecast.jpg" , sep=''))
    list(src = filename,
         width = 580,
         height = 580,
         alt = "-"
    )
  }, deleteFile = F)
  # plot prob f
  output$forecast_percent_FminorFtarget <- renderImage({
    filename <-  file.path(paste("./forecast/", sep=''), paste("forecast_percent_FminorFtarget__.png" , sep=''))
    list(src = filename,
         width = 700,
         height = 300,
         alt = "-"
    )
  }, deleteFile = F)
  # plot prob ssb
  output$forecast_percent_BmajorBtarget <- renderImage({
    filename <-  file.path(paste("./forecast/", sep=''), paste("forecast_percent_BmajorBtarget__.png" , sep=''))
    list(src = filename,
         width = 700,
         height = 300,
         alt = "-"
    )
  }, deleteFile = F)
  # plot prob blim
  output$forecast_percent_BmajorBlim <- renderImage({
    filename <-  file.path(paste("./forecast/", sep=''), paste("forecast_percent_BmajorBlim__.png" , sep=''))
    list(src = filename,
         width = 700,
         height = 300,
         alt = "-"
    )
  }, deleteFile = F)
  #ices table
  output$downloadICES <- downloadHandler(
    filename = function() { 
      paste("tableICES.csv", sep="")
    },
    content = function(file) {
      write.csv(tableICES, file)
    })
  
}

shinyApp(ui, server)



