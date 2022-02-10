#########################################################################################################################
#                                                                                                                       #
#   SS3 ENSEMBLE MODEL script                                                                                           #
#   Run all scenarios at once and plot results, do diagnostic, do Ensemble on final weighted grid                       #
#                                                                                                                       #
#   September 2021                                                                                                      #
#                                                                                                                       #
#   Authors:                                                                                                            #
#   Francesco Masnadi (CNR-IRBIM & UNIBO, Ancona)                                                                       #
#   Massimiliano Cardinale (SLU Aqua, Lysekil)                                                                          #
#   mainly based on ss3diags functions developed by Hennig Winker (JRC-EC) and Felipe Carvalho (NOAA)                   #
#                                                                                                                       #
#  This script is now tailored on ensemble grid of N runs with 5 survey,15 Length data slot and 10 Age data slot.       #
#                                                                                                                       #
#  - ATTENTION: Age slot have to be activate in the Diagnostic loop before running the script (now only Length active)  #
#  - Also "sspar(mfrow=c(X,X)" in Diagnostic loop must be set accordingly                                               #
#  - Retro analysis is now set to -5 year                                                                               #
#  - Mohn Rho value bound set now for long-live species, please change for short-live species                           #  
#  - In Ensemble part make sure you set correct Fref depending on your choice (es. MSY,Btgt,..)                         #
#  - ATTENTION: Weighting vector automatically produce is a pure mean of all diagnostic scores as they are,             #
#    if the user want to use different Weighting methods (ex. merge multiple similar diagnostics into a single score)   #
#    have to produce the vector externally and save it in the main dir as "weight_vector.csv"                           #
#                                                                                                                       #
#                                                                                                                       #
#########################################################################################################################
library(r4ss)
library(ss3diags)
library(dplyr)
library(readr)

main.dir <- "C:/Users/f.masnadi/Desktop/Stock Assessment/SS3/Ensemble_scripts/Es" # set working directory in which you crated the runs subfolders
runs = paste0("run_",1:4) # 1:N Number of run

###############################################################################
# Run of all grid of models together and make SS plots (Skip if already done!)
###############################################################################

 dir <- main.dir
 setwd(dir)
 for(i in 1:length(runs)){
    dir.runN <- paste0(dir,"/",runs[i])
    run_SS_models(dirvec =paste0(dir.runN), extras = "ss", skipfinished = F)
    tmp = SS_output(dir=file.path(runs[i]))
    SS_plots(tmp, uncertainty=T, datplot=F, forecast=T, maxyr=2025, minbthresh=0.25,fitrange = F)
}

###########################################################
# build the empty diags table for weighting the models later
############################################################
results <- data.frame(Run=runs, Convergence=NA, Total_LL=NA, N_Params=NA, Runs_test_cpue1=factor(NA, levels = c("Passed","Failed")), Runs_test_cpue2=factor(NA, levels = c("Passed","Failed")),Runs_test_cpue3=factor(NA, levels = c("Passed","Failed")),Runs_test_cpue4=factor(NA, levels = c("Passed","Failed")), Runs_test_cpue5=factor(NA, levels = c("Passed","Failed")), Runs_test_len1=factor(NA, levels = c("Passed","Failed")), Runs_test_len2=factor(NA, levels = c("Passed","Failed")), Runs_test_len3=factor(NA, levels = c("Passed","Failed")), Runs_test_len4=factor(NA, levels = c("Passed","Failed")), Runs_test_len5=factor(NA, levels = c("Passed","Failed")), Runs_test_len6=factor(NA, levels = c("Passed","Failed")), Runs_test_len7=factor(NA, levels = c("Passed","Failed")), Runs_test_len8=factor(NA, levels = c("Passed","Failed")), Runs_test_len9=factor(NA, levels = c("Passed","Failed")), Runs_test_len10=factor(NA, levels = c("Passed","Failed")), Runs_test_len11=factor(NA, levels = c("Passed","Failed")), Runs_test_len12=factor(NA, levels = c("Passed","Failed")), Runs_test_len13=factor(NA, levels = c("Passed","Failed")), Runs_test_len14=factor(NA, levels = c("Passed","Failed")), Runs_test_len15=factor(NA, levels = c("Passed","Failed")), Runs_test_age1=factor(NA, levels = c("Passed","Failed")), Runs_test_age2=factor(NA, levels = c("Passed","Failed")),Runs_test_age3=factor(NA, levels = c("Passed","Failed")),Runs_test_age4=factor(NA, levels = c("Passed","Failed")),Runs_test_age5=factor(NA, levels = c("Passed","Failed")),Runs_test_age6=factor(NA, levels = c("Passed","Failed")),Runs_test_age7=factor(NA, levels = c("Passed","Failed")),Runs_test_age8=factor(NA, levels = c("Passed","Failed")),Runs_test_age9=factor(NA, levels = c("Passed","Failed")), Runs_test_age10=factor(NA, levels = c("Passed","Failed")), RMSE_Perc=NA, RMSE_Perc_1=NA,RMSE_Perc_2=NA, Retro_Rho_SSB= NA, Forecast_Rho_SSB=NA, Retro_Rho_F= NA, Forecast_Rho_F=NA, MASE_cpue1=NA, MASE_cpue2=NA, MASE_cpue3=NA, MASE_cpue4=NA, MASE_cpue5=NA, MASE_len1=NA, MASE_len2=NA, MASE_len3=NA, MASE_len4=NA, MASE_len5=NA, MASE_len6=NA, MASE_len7=NA,MASE_len8=NA,MASE_len9=NA,MASE_len10=NA,MASE_len11=NA,MASE_len12=NA,MASE_len13=NA,MASE_len14=NA,MASE_len15=NA,MASE_len16=NA,MASE_len17=NA,MASE_len18=NA,MASE_len19=NA,MASE_len20=NA, MASE_age1=NA,MASE_age2=NA,MASE_age3=NA,MASE_age4=NA,MASE_age5=NA,MASE_age6=NA,MASE_age7=NA,MASE_age8=NA,MASE_age9=NA,MASE_age10=NA )
     


##########################################################       
#*********************************************************
# Diagnostic loop
#********************************************************
########################################################## 

for(i in runs){
  # load all scenarios as list  
  dir <- main.dir
  setwd(dir)
  readme<-paste0("/",i)
  tmp <- SS_output(dir=paste0(dir,readme),covar=T,ncols=1000)
  
  ##############################
  # check for Convergence 
  ##############################
  # put in final diags table
  results$Convergence[results$Run== i] <- tmp$maximum_gradient_component
  results$Total_LL[results$Run== i] <-   tmp$likelihoods_used$values[1]
  results$N_Params[results$Run== i] <-   tmp$N_estimated_parameters
  # Par_nearBound <- as.data.frame(ifelse(tmp$parameters$Status  == "OK" | tmp$parameters$Status  == "act", 1, 0) %>% na.omit())
  
  #*********************************************************
  # Basic Residual Diagnostics (joint-residual and run test)
  #********************************************************
 # dir.create("Plotdiags",showWarnings = F) 
  dir.runN <- paste0(dir,"/",i)
  dir.diag <- paste0(dir.runN,'/Plotdiags')
  dir.create(path=dir.diag, showWarnings = T, recursive = T)
  
  # Joint-Residual (JABBAres)
  # Check conflict between indices and mean length or age
  sspar(mfrow=c(1,2),plot.cex = 0.8)   # change to 3 if also Age present
  jr.cpue <- SSplotJABBAres(tmp,subplots="cpue",add=T,col=sscol(3)[c(1,3,2)])
  jr.len <- SSplotJABBAres(tmp,subplots="len",add=T,col=sscol(3)[c(1,3,2)])
# jr.age <-SSplotJABBAres(tmp,subplots="age",add=T,col=sscol(3)[c(1,3,2)])  # activate if age data is present!
  dev.print(jpeg,paste0(dir = dir.diag,"/JointResiduals_",i,".jpg"), width = 8, height = 3.5, res = 300, units = "in")
  # put in final diags table
  results$RMSE_Perc[results$Run== i] <-   jr.cpue$RMSE.perc[jr.cpue$indices=="Combined"]
  results$RMSE_Perc_1[results$Run== i] <-   jr.len$RMSE.perc[jr.len$indices=="Combined"]
 # results$RMSE_Perc_2[results$Run== i] <-   jr.age$RMSE.perc[age.len$indices=="Combined"]
  
  # Check Runs Test 
  sspar(mfrow=c(3,2),plot.cex = 0.8)    # change based on number of fleet and length/age data
  rt.cpue <- SSplotRunstest(tmp,subplots="cpue",add=T,mixing="two.sided")
  rt.len <-SSplotRunstest(tmp,subplots="len",add=T,mixing="two.sided")
# rt.age <- SSplotRunstest(tmp,subplots="age",add=T)   # activate if age data is present!
  dev.print(jpeg,paste0(dir = dir.diag,"/RunTestResidual_",i,".jpg"), width = 8, height = 9, res = 300, units = "in")
  # put in final diags table
  #CPUE
  results$Runs_test_cpue1[results$Run== i] <-   rt.cpue$test[1]
  results$Runs_test_cpue2[results$Run== i] <-   rt.cpue$test[2]
  results$Runs_test_cpue3[results$Run== i] <-   rt.cpue$test[3]
  results$Runs_test_cpue4[results$Run== i] <-   rt.cpue$test[4]
  results$Runs_test_cpue5[results$Run== i] <-   rt.cpue$test[5]
  # length
  results$Runs_test_len1[results$Run== i] <-   rt.len$test[1]
  results$Runs_test_len2[results$Run== i] <-   rt.len$test[2]
  results$Runs_test_len3[results$Run== i] <-   rt.len$test[3]
  results$Runs_test_len4[results$Run== i] <-   rt.len$test[4]
  results$Runs_test_len5[results$Run== i] <-   rt.len$test[5]
  results$Runs_test_len6[results$Run== i] <-   rt.len$test[6]
  results$Runs_test_len7[results$Run== i] <-   rt.len$test[7]
  results$Runs_test_len8[results$Run== i] <-   rt.len$test[8]
  results$Runs_test_len9[results$Run== i] <-   rt.len$test[9]
  results$Runs_test_len10[results$Run== i] <-   rt.len$test[10]
  results$Runs_test_len11[results$Run== i] <-   rt.len$test[11]
  results$Runs_test_len12[results$Run== i] <-   rt.len$test[12]
  results$Runs_test_len13[results$Run== i] <-   rt.len$test[13]
  results$Runs_test_len14[results$Run== i] <-   rt.len$test[14]
  results$Runs_test_len15[results$Run== i] <-   rt.len$test[15]
  # Age
#  results$Runs_test_age1[results$Run== i] <-   rt.age$test[1]
#  results$Runs_test_age2[results$Run== i] <-   rt.age$test[2]
#  results$Runs_test_age3[results$Run== i] <-   rt.age$test[3]
#  results$Runs_test_age4[results$Run== i] <-   rt.age$test[4]
#  results$Runs_test_age5[results$Run== i] <-   rt.age$test[5]
#  results$Runs_test_age6[results$Run== i] <-   rt.age$test[6]
#  results$Runs_test_age7[results$Run== i] <-   rt.age$test[7]
#  results$Runs_test_age8[results$Run== i] <-   rt.age$test[8]
#  results$Runs_test_age9[results$Run== i] <-   rt.age$test[9]
#  results$Runs_test_age10[results$Run== i] <-   rt.age$test[10]


  ##############################
  # Retrospective analyses
  ##############################
  ## Automatically running the retrospective analyses
  start.retro <- 0    #end year of model
  end.retro   <- 5    #number of years for retrospective 
  # Identify the directory where a completed model run is located
  dirname.completed.model.run <- file.path(dir.runN)
  # Create a subdirectory for the Retrospectives
  dirname.Retrospective <- paste0(dir.runN,'/Retrospective')
  dir.create(path=dirname.Retrospective, showWarnings = TRUE, recursive = TRUE)
  setwd(dirname.Retrospective)
  #----------------- copy model run files ----------------------------------------
  file.copy(paste(dirname.completed.model.run,       "starter.ss_new", sep="/"),
            paste(dirname.Retrospective, "starter.ss", sep="/"))
  file.copy(paste(dirname.completed.model.run,       "control.ss_new", sep="/"),
            paste(dirname.Retrospective, "CONTROL.SS", sep="/"))
  file.copy(paste(dirname.completed.model.run,       "data.ss_new", sep="/"),
            paste(dirname.Retrospective, "DATA.SS", sep="/"))	
  file.copy(paste(dirname.completed.model.run,       "forecast.ss", sep="/"),
            paste(dirname.Retrospective, "forecast.ss", sep="/"))
  file.copy(paste(dirname.completed.model.run,       "SS.exe", sep="/"),
            paste(dirname.Retrospective, "SS.exe", sep="/"))
  file.copy(paste(dirname.completed.model.run,       "wtatage.ss", sep="/"),
            paste(dirname.Retrospective, "wtatage.ss", sep="/"))
  #------------Make Changes to the Starter.ss file (DC Example) ------------------------------- 
  starter <- readLines(paste(dirname.Retrospective, "/starter.ss", sep=""))
  # Starter File changes to speed up model runs
  # Run Display Detail
  #[8] "2 # run display detail (0,1,2)" 
  linen <- grep("# run display detail", starter)
  starter[linen] <- paste0( 1 , " # run display detail (0,1,2)" )
  write(starter, paste(dirname.Retrospective, "starter.ss", sep="/"))
  # Run the retrospective analyses with r4SS function "SS_doRetro"
  # Switch off Hessian "-nohess" (much faster)
  SS_doRetro(masterdir=dirname.Retrospective, oldsubdir="", newsubdir="", years=start.retro:-end.retro,extras="-nohess")
  # Read "SS_doRetro" output
  retroModels <- SSgetoutput(dirvec=file.path(paste0(dirname.Retrospective),paste("retro",start.retro:-end.retro,sep="")))
  # Summarize output
  retroSummary <- r4ss::SSsummarize(retroModels)
  endyrvec <- retroSummary$endyrs + start.retro:-end.retro
  sspar(mfrow=c(1,2),plot.cex = 0.8)   
  retro.ssb<- SSplotRetro(retroSummary,add=T,legendcex=0.8,tickEndYr=F,xylabs=F,legendloc = "bottomleft",uncertainty = T,showrho = T,forecast = T,labels="SSB (t)", endyrvec = c(endyrvec), subplots = c("SSB"))
    retro.f<- SSplotRetro(retroSummary,add=T,legendcex=0.8,tickEndYr=F,xylabs=F,legendloc = "bottomleft",uncertainty = T,showrho = T,forecast = T,labels="SSB (t)", endyrvec = c(endyrvec), subplots = c("F"))
  dev.print(jpeg,paste0(dir = dir.diag,"/Retro_",i,".jpg"), width = 12, height = 5, res = 300, units = "in")
   # put in final diags table
  results$Retro_Rho_SSB[results$Run== i] <-   retro.ssb$Rho[retro.ssb$peel=="Combined"]
  results$Forecast_Rho_SSB[results$Run== i] <-   retro.ssb$ForecastRho[retro.ssb$peel=="Combined"]
  results$Retro_Rho_F[results$Run== i] <-   retro.f$Rho[retro.f$peel=="Combined"]
  results$Forecast_Rho_F[results$Run== i] <-   retro.f$ForecastRho[retro.f$peel=="Combined"]

  ##############################
  # Hindcasting analyses
  ##############################
  # Do Hindcast with Cross-Validation of CPUE observations
  sspar(mfrow=c(1,1),plot.cex = 0.9)     # change based on number of cpue or survey
  mase_cpue <- SSplotHCxval(retroSummary,add=T)
  dev.print(jpeg,paste0(dir = dir.diag,"/HCxval_CPUE_",i,".jpg"), width = 8, height = 9, res = 300, units = "in")
  results$MASE_cpue1[results$Run== i] <-   mase_cpue$MASE[1]
  results$MASE_cpue2[results$Run== i] <-   mase_cpue$MASE[2]
  results$MASE_cpue3[results$Run== i] <-   mase_cpue$MASE[3]
  results$MASE_cpue4[results$Run== i] <-   mase_cpue$MASE[4]
  results$MASE_cpue5[results$Run== i] <-   mase_cpue$MASE[5]
  # Also Hindcast with Cross-Validation for mean length or age
   # Use new converter fuction SSretroComps()
  hccomps = SSretroComps(retroModels)
  # Specify subplots = "age" or "len" in SSplotHCxval
  sspar(mfrow=c(3,2),plot.cex = 0.7) # change based on number of fleet or survey per length or age data
  mase_len.plot <- SSplotHCxval(hccomps,add=T,subplots = "len",legendloc="topright",legend = FALSE, indexUncertainty = T,legendcex = 1)
#  mase_age <- SSplotHCxval(hccomps,add=T,subplots = "age",legendloc="topright",legend = FALSE, indexUncertainty = T,legendcex = 1) # activate if age data is present!
  dev.print(jpeg,paste0(dir = dir.diag,"/HCxval_length_",i,".jpg"), width = 8, height = 9, res = 300, units = "in")
  # put in final diags table
  mase_len <- SSmase(hccomps,quants = "len")
  results$MASE_len1[results$Run== i] <-   mase_len$MASE.adj[1]
  results$MASE_len2[results$Run== i] <-   mase_len$MASE.adj[2]
  results$MASE_len3[results$Run== i] <-   mase_len$MASE.adj[3]
  results$MASE_len4[results$Run== i] <-   mase_len$MASE.adj[4]
  results$MASE_len5[results$Run== i] <-   mase_len$MASE.adj[5]
  results$MASE_len6[results$Run== i] <-   mase_len$MASE.adj[6]
  results$MASE_len7[results$Run== i] <-   mase_len$MASE.adj[7]
  results$MASE_len8[results$Run== i] <-   mase_len$MASE.adj[8]
  results$MASE_len9[results$Run== i] <-   mase_len$MASE.adj[9]
  results$MASE_len10[results$Run== i] <-   mase_len$MASE.adj[10]
  results$MASE_len11[results$Run== i] <-   mase_len$MASE.adj[11]
  results$MASE_len12[results$Run== i] <-   mase_len$MASE.adj[12]
  results$MASE_len13[results$Run== i] <-   mase_len$MASE.adj[13]
  results$MASE_len14[results$Run== i] <-   mase_len$MASE.adj[14]
  results$MASE_len15[results$Run== i] <-   mase_len$MASE.adj[15]
  results$MASE_len16[results$Run== i] <-   mase_len$MASE.adj[16]
  results$MASE_len17[results$Run== i] <-   mase_len$MASE.adj[17]
  results$MASE_len18[results$Run== i] <-   mase_len$MASE.adj[18]
  results$MASE_len19[results$Run== i] <-   mase_len$MASE.adj[19]
  results$MASE_len20[results$Run== i] <-   mase_len$MASE.adj[20]
  # age
#  results$MASE_age1[results$Run== i] <-   mase_age$MASE.adj[1]
#  results$MASE_age2[results$Run== i] <-   mase_age$MASE.adj[2]
#  results$MASE_age3[results$Run== i] <-   mase_age$MASE.adj[3]
#  results$MASE_age4[results$Run== i] <-   mase_age$MASE.adj[4]
#  results$MASE_age5[results$Run== i] <-   mase_age$MASE.adj[5]
#  results$MASE_age6[results$Run== i] <-   mase_age$MASE.adj[6]
#  results$MASE_age7[results$Run== i] <-   mase_age$MASE.adj[7]
#  results$MASE_age8[results$Run== i] <-   mase_age$MASE.adj[8]
#  results$MASE_age9[results$Run== i] <-   mase_age$MASE.adj[9]
#  results$MASE_age10[results$Run== i] <-   mase_age$MASE.adj[10]
}


#*********************************************************
# save the table of diagnostic with real value
#********************************************************
diags_table <- as.data.frame(t(results) %>% na.omit())
write.csv(diags_table, file = paste0(main.dir,"/Diags_table.csv"))


#*********************************************************
# save the table of diagnostic for weigthing purpose (value 0 to 1)
#********************************************************
results2 <- results %>% select(-Convergence,-Total_LL,-N_Params,-Run)

results2$Runs_test_cpue1 <-   ifelse(results2$Runs_test_cpue1 == 'Passed', 1, 0)
results2$Runs_test_cpue2 <-   ifelse(results2$Runs_test_cpue2 == 'Passed', 1, 0)
results2$Runs_test_cpue3 <-   ifelse(results2$Runs_test_cpue3 == 'Passed', 1, 0)
results2$Runs_test_cpue4 <-   ifelse(results2$Runs_test_cpue4 == 'Passed', 1, 0)
results2$Runs_test_cpue5 <-   ifelse(results2$Runs_test_cpue5 == 'Passed', 1, 0)
# length
results2$Runs_test_len1 <-   ifelse(results2$Runs_test_len1 == 'Passed', 1, 0)
results2$Runs_test_len2 <-   ifelse(results2$Runs_test_len2 == 'Passed', 1, 0)
results2$Runs_test_len3 <-   ifelse(results2$Runs_test_len3 == 'Passed', 1, 0)
results2$Runs_test_len4 <-   ifelse(results2$Runs_test_len4 == 'Passed', 1, 0)
results2$Runs_test_len5 <-   ifelse(results2$Runs_test_len5 == 'Passed', 1, 0)
results2$Runs_test_len6 <-   ifelse(results2$Runs_test_len6 == 'Passed', 1, 0)
results2$Runs_test_len7 <-   ifelse(results2$Runs_test_len7 == 'Passed', 1, 0)
results2$Runs_test_len8 <-   ifelse(results2$Runs_test_len8 == 'Passed', 1, 0)
results2$Runs_test_len9 <-   ifelse(results2$Runs_test_len9 == 'Passed', 1, 0)
results2$Runs_test_len10 <-   ifelse(results2$Runs_test_len10 == 'Passed', 1, 0)
results2$Runs_test_len11 <-   ifelse(results2$Runs_test_len11 == 'Passed', 1, 0)
results2$Runs_test_len12 <-   ifelse(results2$Runs_test_len12 == 'Passed', 1, 0)
results2$Runs_test_len13 <-   ifelse(results2$Runs_test_len13 == 'Passed', 1, 0)
results2$Runs_test_len14 <-   ifelse(results2$Runs_test_len14 == 'Passed', 1, 0)
results2$Runs_test_len15 <-   ifelse(results2$Runs_test_len15 == 'Passed', 1, 0)
# Age
#  results2$Runs_test_age1 <-   ifelse(results2$Runs_test_age1 == 'Passed', 1, 0)
#  results2$Runs_test_age2 <-   ifelse(results2$Runs_test_age2 == 'Passed', 1, 0)
#  results2$Runs_test_age3 <-   ifelse(results2$Runs_test_age3 == 'Passed', 1, 0)
#  results2$Runs_test_age4 <-   ifelse(results2$Runs_test_age4 == 'Passed', 1, 0)
#  results2$Runs_test_age5 <-   ifelse(results2$Runs_test_age5 == 'Passed', 1, 0)
#  results2$Runs_test_age6 <-   ifelse(results2$Runs_test_age6 == 'Passed', 1, 0)
#  results2$Runs_test_age7 <-   ifelse(results2$Runs_test_age7 == 'Passed', 1, 0)
#  results2$Runs_test_age8 <-   ifelse(results2$Runs_test_age8 == 'Passed', 1, 0)
#  results2$Runs_test_age9 <-   ifelse(results2$Runs_test_age9 == 'Passed', 1, 0)
#  results2$Runs_test_age10 <-   ifelse(results2$Runs_test_age10 == 'Passed', 1, 0)
results2$RMSE_Perc <- ifelse(results2$RMSE_Perc < 30, 1, 0)
results2$RMSE_Perc_1 <- ifelse(results2$RMSE_Perc_1 < 30, 1, 0)
# Age
#  results2$RMSE_Perc_2 <- ifelse(results2$RMSE_Perc_2 < 30, 1, 0)
results2$Retro_Rho_SSB <- ifelse(results2$Retro_Rho_SSB > -0.15 & results2$Retro_Rho_SSB < 0.2, 1, 0) # change to -0.22-0.30 for shorter-lived species
results2$Forecast_Rho_SSB <- ifelse(results2$Forecast_Rho_SSB > -0.15 & results2$Forecast_Rho_SSB < 0.2, 1, 0) # change to -0.22-0.30 for shorter-lived species
results2$Retro_Rho_F <- ifelse(results2$Retro_Rho_F > -0.15 & results2$Retro_Rho_F < 0.2, 1, 0) # change to -0.22-0.30 for shorter-lived species
results2$Forecast_Rho_F <- ifelse(results2$Forecast_Rho_F > -0.15 & results2$Forecast_Rho_F < 0.2, 1, 0) # change to -0.22-0.30 for shorter-lived species
results2$MASE_cpue1 <- ifelse(results2$MASE_cpue1 < 1, 1, 0)
results2$MASE_cpue2 <- ifelse(results2$MASE_cpue2 < 1, 1, 0)
results2$MASE_cpue3 <- ifelse(results2$MASE_cpue3 < 1, 1, 0)
results2$MASE_cpue4 <- ifelse(results2$MASE_cpue4 < 1, 1, 0)
results2$MASE_cpue5 <- ifelse(results2$MASE_cpue5 < 1, 1, 0)
results2$MASE_len1 <- ifelse(results2$MASE_len1 < 1, 1, 0)
results2$MASE_len2 <- ifelse(results2$MASE_len2 < 1, 1, 0)
results2$MASE_len3 <- ifelse(results2$MASE_len3 < 1, 1, 0)
results2$MASE_len4 <- ifelse(results2$MASE_len4 < 1, 1, 0)
results2$MASE_len5 <- ifelse(results2$MASE_len5 < 1, 1, 0)
results2$MASE_len6 <- ifelse(results2$MASE_len6 < 1, 1, 0)
results2$MASE_len7 <- ifelse(results2$MASE_len7 < 1, 1, 0)
results2$MASE_len8 <- ifelse(results2$MASE_len8 < 1, 1, 0)
results2$MASE_len9 <- ifelse(results2$MASE_len9 < 1, 1, 0)
results2$MASE_len10 <- ifelse(results2$MASE_len10 < 1, 1, 0)
results2$MASE_len11 <- ifelse(results2$MASE_len11 < 1, 1, 0)
results2$MASE_len12 <- ifelse(results2$MASE_len12 < 1, 1, 0)
results2$MASE_len13 <- ifelse(results2$MASE_len13 < 1, 1, 0)
results2$MASE_len14 <- ifelse(results2$MASE_len14 < 1, 1, 0)
results2$MASE_len15 <- ifelse(results2$MASE_len15 < 1, 1, 0)
results2$MASE_len16 <- ifelse(results2$MASE_len16 < 1, 1, 0)
results2$MASE_len17 <- ifelse(results2$MASE_len17 < 1, 1, 0)
results2$MASE_len18 <- ifelse(results2$MASE_len18 < 1, 1, 0)
results2$MASE_len19 <- ifelse(results2$MASE_len19 < 1, 1, 0)
results2$MASE_len20 <- ifelse(results2$MASE_len20 < 1, 1, 0)
# Age
# results2$MASE_age1 <- ifelse(results2$MASE_age1 < 1, 1, 0)
# results2$MASE_age2 <- ifelse(results2$MASE_age2 < 1, 1, 0)
# results2$MASE_age3 <- ifelse(results2$MASE_age3 < 1, 1, 0)
# results2$MASE_age4 <- ifelse(results2$MASE_age4 < 1, 1, 0)
# results2$MASE_age5 <- ifelse(results2$MASE_age5 < 1, 1, 0)
# results2$MASE_age6 <- ifelse(results2$MASE_age6 < 1, 1, 0)
# results2$MASE_age7 <- ifelse(results2$MASE_age7 < 1, 1, 0)
# results2$MASE_age8 <- ifelse(results2$MASE_age8 < 1, 1, 0)
# results2$MASE_age9 <- ifelse(results2$MASE_age9 < 1, 1, 0)
# results2$MASE_age10 <- ifelse(results2$MASE_age10 < 1, 1, 0)

weight_table <-  as.data.frame(t(results2) %>% na.omit())
write.csv(weight_table, file = paste0(main.dir,"/Weight_table.csv"))


####################################################
#### weight_vector to be used in Ensemble
####################################################
weight_vector <- rowMeans(as.data.frame(lapply(as.data.frame(t(weight_table)), as.numeric)))
write.csv(t(weight_vector), file = paste0(main.dir,"/weight_vector.csv"),row.names=FALSE)



#****************************************************************
# Approximate uncertainty with MVLN (hessian)
#****************************************************************
setwd(main.dir)  # set dir to the initial one 
#### weight_vector from diagnostic scores coming from "Ensemble_diags_runs.R" script
weight_vector <- unname(unlist(read_csv(paste0(main.dir,"/weight_vector.csv"))[1,]))

kbproj = NULL
# Compile MVLN posteriors by scenario run
for(i in 1:length(runs)){
  # load all scenarios as list  
  run = SS_output(dir=file.path(runs[i])) 
  # get MVLN mvn.temp for each scenario
  # Make sure you set correct Fref depending on your choice (in this case F at Btarget)
  mvn.temp = SSdeltaMVLN(run,run = runs[i],years = run$startyr:run$endyr,addprj = T,mc=5000,weight = weight_vector[i],Fref = "Btgt",plot=F)
  kbproj = rbind(kbproj,data.frame(mvn.temp$kb,model=runs[i]))
  # save labels once
  if(i==1 ) labels = mvn.temp$labels
} 
#### save Ensemble r.data
save(kbproj,file="Ensemble_model.rdata")

##########################
# FINAL Kobe_plot ensemble
##########################
sspar(mfrow=c(1,1),plot.cex = 0.9)
SSplotKobe(kbproj,fill=T,joint=F,posterior="kernel",ylab=labels[2],xlab=labels[1])
dev.print(jpeg,paste0("Kobe_final_wProva.jpg"), width = 12, height = 8, res = 300, units = "in") 

# Show trajectories one by one
sspar(mfrow=c(3,2),plot.cex = 0.9)
# zoom <- c(2010,2020)
quants = c("stock", "harvest", "SSB", "F", "Recr", "Catch")
for(i in 1:1){
  SSplotEnsemble(kbproj,quantiles = c(0.05, 0.95),add=T,subplots = quants[i],ylabs=labels[i],
                 legendcex = 0.6,legendloc = "topleft",legend=ifelse(i==1,T,T))
  abline(h= 0.5, lty = 2)   # in this case h= 0.5 because Blim is half of Btarget
}
for(i in 2:6){
  SSplotEnsemble(kbproj,quantiles = c(0.05, 0.95),add=T,subplots = quants[i],ylabs=labels[i],
                 legendcex = 0.8,legendloc = "topleft",legend=ifelse(i==1,F,F))
}
dev.print(jpeg,paste0("MLVN_Compare_wProva.jpg"), width = 12, height = 8, res = 300, units = "in")

# Show trajectories All together 
kbproj2  <- kbproj 
kbproj2$run <- "All_runs"
sspar(mfrow=c(3,2),plot.cex = 0.9)
# zoom <- c(2010,2020)
quants = c("stock", "harvest", "SSB", "F", "Recr", "Catch")
for(i in 1:1){
  SSplotEnsemble(kbproj2,quantiles = c(0.05, 0.95),add=T,subplots = quants[i],ylabs=labels[i],
                 legendcex = 0.6,legendloc = "topleft",legend=ifelse(i==1,T,T))
  abline(h= 0.5, lty = 2) # in this case h= 0.5 because Blim is half of Btarget
}
for(i in 2:6){
  SSplotEnsemble(kbproj2,quantiles = c(0.05, 0.95),add=T,subplots = quants[i],ylabs=labels[i],
                 legendcex = 0.8,legendloc = "topleft",legend=ifelse(i==1,F,F))
}
dev.print(jpeg,paste0("MLVN_All_wProva.jpg"), width = 12, height = 8, res = 300, units = "in")



# get median + 90% CIs for SSB/SSBmsy across all models
prj_SBB_SBBtrg = aggregate(stock~year,kbproj,quantile,c(0.5,0.05,0.95))
prj_F_Btrg = aggregate(harvest~year,kbproj,quantile,c(0.5,0.05,0.95))
prj_SSB = aggregate(SSB~year,kbproj,quantile,c(0.5,0.05,0.95))
prj_F = aggregate(F~year,kbproj,quantile,c(0.5,0.05,0.95))
prj_Recr = aggregate(Recr~year,kbproj,quantile,c(0.5,0.05,0.95))
write.csv(prj_SBB_SBBtrg , "Ensemble_SBB_SBBtrg.csv")
write.csv(prj_F_Btrg , "Ensemble_F_Btrg.csv")
write.csv(prj_SSB , "Ensemble_SBB.csv")
write.csv(prj_F , "Ensemble_F.csv")
write.csv(prj_Recr , "Ensemble_Recr.csv")




