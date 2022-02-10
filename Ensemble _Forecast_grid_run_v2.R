#+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + 
# 
#   SS3 ENSEMBLE MODEL script                                                                                           
#   Copy and paste forecast scenarios SS files,                                                
#   run all forecast scenarios models at once and plot results                                                                                 
#   December 2021                                                                                                      
#                                                                                                                       
#   Authors:                                                                                                            
#   Francesco Masnadi (CNR-IRBIM & UNIBO, Ancona)                                                                       
#   Massimiliano Cardinale (SLU Aqua, Lysekil)                                                                          
#   mainly based on ss3diags functions developed by Hennig Winker (JRC-EC) and Felipe Carvalho (NOAA)  
#
#+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + 


library(r4ss)
library(ss3diags)
library(readr)
library(kobe)
library(plyr)
library(reshape)

# Defining Blim e Btrigger for later (line 170 and 172 of the script)
###################
RATIOLIM <- 2     # ratio of Btarget on Blim: e.g. if target=SSB40% and lim=SSB20% --> 40/20=2  
RATIOTRIG <- 1.25 # ratio of 1 on trigger percentage: e.g. if trigger=80% of target --> 1/0.8=1.25  
###################

#set working directorty in wich you want to crate the runs subfolder:
dir = "C:/Users/f.masnadi/Desktop/Stock Assessment/SS3/Ensemble_scripts/Es" 
setwd(dir)

# Subfolder for the runs from 1 to N:
runs = paste0("run_",1:2)  # here 2 subfolder
# Subfolder within TAC levels scenarios for forecast:
tacs = paste0("TAC",seq(80,120,20))

# link to forecast TAC subfolder where you manually created the forecast files with different catch levels simulating the TACs; name the different files "forecast_TAC..."
dir.forecastTAC <-  paste0(dir,"/forecast_TAC")

# create and run forecast folder and subfolders (if the first time)
for(i in 1:length(runs)){
  dir.runN <- paste0(dir,"/",runs[i])
  dir.runN.new <- paste0(dir,"/Forecast/",runs[i])
  dir.create(path=dir.runN.new, showWarnings = T, recursive = T)
  for(j in 1:length(tacs)){
    dir.tacN <- paste0(dir.runN.new,"/",tacs[j])
    dir.create(path=dir.tacN, showWarnings = T, recursive = T)
    # copy the SS base files in every TAC subfolder 
    file.copy(paste(dir.runN,       "starter.ss_new", sep="/"),
              paste(dir.tacN, "starter.ss", sep="/"))
    file.copy(paste(dir.runN,       "control.ss_new", sep="/"),
              paste(dir.tacN, "control.ss", sep="/"))
    file.copy(paste(dir.runN,       "data.ss_new", sep="/"),
              paste(dir.tacN, "data.ss", sep="/"))	
    file.copy(paste(dir.runN,       "SS.exe", sep="/"),
              paste(dir.tacN, "SS.exe", sep="/"))
    file.copy(paste(dir.runN,       "wtatage.ss", sep="/"),
              paste(dir.tacN, "wtatage.ss", sep="/"))
    # copy the right forecast file from the "forecast_TAC" folder
    file.copy(paste(dir.forecastTAC,  paste0("forecast_",tacs[j], ".ss") , sep="/"),
              paste(dir.tacN, "forecast.ss", sep="/"))
  }
}

############################################################
# Run of all grid of models togheter
############################################################
for(i in 1:length(runs)){
  for(j in 1:length(tacs)){
    dir.runN.new <- paste0(dir,"/Forecast/",runs[i])
    dir.tacN <- paste0(dir.runN.new,"/",tacs[j])
  run_SS_models(dirvec =paste0(dir.tacN), extras = "ss", skipfinished = F)
 }
}

##################################
# ENSEMBLE FORECAST 
#################################
dir = "C:/Users/f.masnadi/Desktop/Stock Assessment/SS3/Ensemble_scripts/Es/Forecast" 
setwd(dir)

#### weight_vector from diagnostic scores coming from "Ensemble_diags_runs.R" script
weight_vector <- unname(unlist(read_csv("C:/Users/f.masnadi/Desktop/Stock Assessment/SS3/Ensemble_scripts/Es/weight_vector.csv")[1,]))

kbproj.fr = NULL
# Compile MVLN posteriors by scenario run
for(i in 1:length(runs)){
  # load all scenarios as list  
  run = SSgetoutput(dirvec=file.path(runs[i],tacs))  
  # get MVLN mvn.temp for each scenario
  for(j in 1:length(tacs)){  
    # this is a super long projection so I decrease mc to 1000 to speed up 
    # Make sure you set correct Fref depending on your choice
    # Make sure you name run = according to TACs
    mvn.temp = SSdeltaMVLN(run[[j]],run=tacs[j],years = run[[j]]$startyr:run[[j]]$endyr+run[[j]]$nforecastyears,addprj = T,mc=5000,weight = 1,Fref = "Btgt",plot=F)
    # build kbproj.fr and add model name [j] for each TAC i 
    kbproj.fr = rbind(kbproj.fr,data.frame(mvn.temp$kb,model=runs[i]))
    # save labels once
    if(i==1 & j ==1) labels = mvn.temp$labels
  } # end j loop
} # end k loop
#### save Ensemble r.data
save(kbproj.fr,file="Ensemble_Forecast_model.rdata")


############################################
# Plotting forecast results
############################################
# Name labels correctly if Btgt
zoom <- c(2010,2024) # put the range you want to visualize in the forecast plot (zoom) 
labels <- expression("SSB/SSB"[40], "F/F"[SB ~ 40], "SSB", 
                   "F", "Recruits", "Catch")
sspar(mfrow=c(3,2),plot.cex = 0.9)
quants = c("stock", "harvest", "SSB", "F", "Recr", "Catch")
for(i in 1:1){
  SSplotEnsemble(kbproj.fr,xlim=zoom,add=T,quantiles = c(0.05, 0.95),subplots = quants[i],ylabs=labels[i],
                 legendcex = 0.8,legendloc = "bottomleft",legend=ifelse(i==1,T,F))
  abline(h= 0.5, lty = 2)   # h = % of Blim respect to Btarget
}
for(i in 2:6){
  SSplotEnsemble(kbproj.fr,xlim=zoom,add=T,quantiles = c(0.05, 0.95),subplots = quants[i],ylabs=labels[i],
                 legendcex = 0.8,legendloc = "topleft",legend=ifelse(i==2,F,F))
}  
dev.print(jpeg,paste0("MLVN_Ensemble_forecast.jpg"), width = 12, height = 8, res = 300, units = "in") 

# check Kobe for final year of forecast (2024; need to fix)
sspar(mfrow=c(1,1),plot.cex = 1)
yrs = c(2024)
for(i in 1:1){
  SSplotKobe(kbproj.fr,joint=F,year=yrs[i],legendruns=T,fill=F,ylab=labels[1],xlab=labels[2])
  legend("top",paste(yrs[i]),bty="n",cex=1.2,y.intersp = -0.5, x.intersp = -0.5)
}
dev.print(jpeg,paste0("Kobe_Ensemble_forecast.jpg"), width = 8, height = 8, res = 300, units = "in")



####################################################################
# Probabilistic Forecast table and results (Ftarget, Btarget, Blim)
####################################################################

prob.dir <- paste0(dir,"/prob_forecast")
dir.create(path=prob.dir, showWarnings = T, recursive = T)

# Load forecast grid rdata file if not already loaded
# load("C:/Users/f.masnadi/Desktop/Stock Assessment/SS3/Ensemble_scripts/Es/Forecast/Ensemble_Forecast_model.rdata", verbose=T)

nTAC = length(tacs)
yrs = 2021:2024 # Projection horizon   XXX DG: add 2024?

### data.frame saved in the loaded .rdata-file miss the following columns:
###
### $Btriggerratio
### $Fmsyratio
### moreover added Blim and Btrigger for ICES forecast table
#################################################################################
kbproj.fr$Ftrgtratio <- kbproj.fr$harvest
kbproj.fr$Btrftratio <- kbproj.fr$stock
# adding Blim 
kbproj.fr$Blimratio <- kbproj.fr$Btrftratio * RATIOLIM   
# adding Btrigger 
kbproj.fr$Btriggerratio <- kbproj.fr$Btrftratio * RATIOTRIG 
kbproj.fr$tac <- kbproj.fr$run       # rename to fit the script
#################################################################################



########################################
# Table for SSB target
# Prepare matrices using library(kobe) 
t. = plyr::ddply(kbproj.fr, .(year,tac), with, kobe:::smry(Btrftratio, Ftrgtratio )) 
# Compile projection matrices in tables
k2smTab=list()
k2smTab[[1]]=cast(t., tac~year, value="underFishing")
k2smTab[[2]]=cast(t., tac~year, value="underFished")
k2smTab[[3]]=cast(t., tac~year, value="green")

########################################
# Table for SSB LIM
# Prepare matrices using library(kobe) 
t.lim = plyr::ddply(kbproj.fr, .(year,tac), with, kobe:::smry(Blimratio, Ftrgtratio )) 
# Compile projection matrices in tables
k2smTabt.lim=list()
k2smTabt.lim[[1]]=cast(t.lim, tac~year, value="underFished")

########################################
# Table for SSB Btrigger
# Prepare matrices using library(kobe) 
t.trig = plyr::ddply(kbproj.fr, .(year,tac), with, kobe:::smry(Btriggerratio, Ftrgtratio )) 


# Probabilistic plot Target
# text size
tx = 0.75
titles <- list("under Ftarget (F<Ftrgt)","above Btarget (SSB>SSBtrgt)","in the green zone (Kobeplot: SSB>SSBtrgt & F<Ftrgt)")
#Loop through B > Btarget, F < Ftarget and Kobe K2SM tables
for(k in 1:3){
  # Define Projection matrix
  pjm = as.matrix(k2smTab[[k]])
  pjm = round(pjm[, which(colnames(pjm) %in% paste(yrs))]*100,0)
  ypr = as.numeric(colnames(pjm))
  npr = length(ypr)
  mat.names = c("/forecast_percent_FminorFtarget_","/forecast_percent_BmajorBtarget_","/forecast_percent_KobeGreen_Btarget")
  
  # Write table
  pjm.save = k2smTab[[k]]
  pjm.save[,2:ncol(pjm.save)] = round(pjm.save[,2:ncol(pjm.save)]*100,1)
  write.csv(pjm.save,paste0(prob.dir,paste(mat.names[k]),".csv"),row.names = FALSE)
  
  op <- par(mfrow = c(1,1),mar = rep(1, 4),mai = rep(1, 4),omi = rep(0, 4))
  png(file = paste0(prob.dir,"/",paste(mat.names[k]),"_",".png"), width = 6.5, height = 6.5*nTAC/npr*0.4 ,
      res = 200, units = "in")
  par(op)
  plot(1, 1, axes=FALSE, frame.plot=FALSE, xlab="", ylab="", type="n", ylim=c(0,nTAC+1), xlim=c(-1,npr+1),main=paste0("Probability of being " ,titles[k]) ,cex.main=0.75)
  # first line
  rect(-1,1:nTAC+1,1,0:nTAC);
  text(rep(0,nTAC+1),seq(0.5,nTAC+0.5,1),c(rev(paste(tacs)),"TAC | Year"),cex=tx)
  
  # Set grey shading
  # Set grey shading
  mat = pjm/100
  mat[mat<0.5]=-1
  mat=(1-mat)*2
  mat[mat>1]=1
  
  rect(1:npr,rep(nTAC+1,npr),1:(npr+1),rep(nTAC,npr))
  text(c(seq(1.5,npr+0.5,1)),rep(nTAC+0.5,npr),c(paste(ypr)),cex=tx)
  for(t in 1:nTAC) rect(1:npr,rep(nTAC+1-t,npr),2:(npr+1),rep(nTAC-t,npr),col=grey(mat[tacs[t],],0.5))
  for(t in 1:nTAC) text(c(seq(1.5,npr+0.5,1)),rep(nTAC+0.5-t,npr),paste(pjm[tacs[t],]),cex=tx)
  dev.off()
}



# Probabilistic plot Limit
# text size
tx = 0.75
titles <- list("above Blimit (SSB>SSBlim)")
#Loop through B > Btarget, F < Ftarget and Kobe K2SM tables
for(k in 1:1){
  # Define Projection matrix
  pjm = as.matrix(k2smTabt.lim[[k]])
  pjm = round(pjm[, which(colnames(pjm) %in% paste(yrs))]*100,0)
  ypr = as.numeric(colnames(pjm))
  npr = length(ypr)
  mat.names = c("/forecast_percent_BmajorBlim_")
  
  # Write table
  pjm.save = k2smTabt.lim[[k]]
  pjm.save[,2:ncol(pjm.save)] = round(pjm.save[,2:ncol(pjm.save)]*100,1)
  write.csv(pjm.save,paste0(prob.dir,paste(mat.names[k]),".csv"),row.names = FALSE)
  
  op <- par(mfrow = c(1,1),mar = rep(1, 4),mai = rep(1, 4),omi = rep(0, 4))
  png(file = paste0(prob.dir,"/",paste(mat.names[k]),"_",".png"), width = 6.5, height = 6.5*nTAC/npr*0.4 ,
      res = 200, units = "in")
  par(op)
  plot(1, 1, axes=FALSE, frame.plot=FALSE, xlab="", ylab="", type="n", ylim=c(0,nTAC+1), xlim=c(-1,npr+1),main=paste0("Probability of being " ,titles[k]) ,cex.main=0.75)
  # first line
  rect(-1,1:nTAC+1,1,0:nTAC);
  text(rep(0,nTAC+1),seq(0.5,nTAC+0.5,1),c(rev(paste(tacs)),"TAC | Year"),cex=tx)
  
  # Set grey shading
  # Set grey shading
  mat = pjm/100
  mat[mat<0.5]=-1
  mat=(1-mat)*2
  mat[mat>1]=1
  
  rect(1:npr,rep(nTAC+1,npr),1:(npr+1),rep(nTAC,npr))
  text(c(seq(1.5,npr+0.5,1)),rep(nTAC+0.5,npr),c(paste(ypr)),cex=tx)
  for(t in 1:nTAC) rect(1:npr,rep(nTAC+1-t,npr),2:(npr+1),rep(nTAC-t,npr),col=grey(mat[tacs[t],],0.5))
  for(t in 1:nTAC) text(c(seq(1.5,npr+0.5,1)),rep(nTAC+0.5-t,npr),paste(pjm[tacs[t],]),cex=tx)
  dev.off()
}


#####################################
# csv forecast table for ICES
#####################################

prj_SBB_SBBtrg = aggregate(stock~year+tac,kbproj.fr,median)%>% dplyr::filter(year %in% yrs)
prj_F_Btrg = aggregate(harvest~year+tac,kbproj.fr,median) %>% dplyr::filter(year %in% yrs)
prj_SSB = aggregate(SSB~year+tac,kbproj.fr,median) %>% dplyr::filter(year %in% yrs)
prj_F = aggregate(F~year+tac,kbproj.fr,median) %>% dplyr::filter(year %in% yrs)
prj_Recr = aggregate(Recr~year+tac,kbproj.fr,median) %>% dplyr::filter(year %in% yrs)
prj_Catch = aggregate(Catch~year+tac,kbproj.fr,median) %>% dplyr::filter(year %in% yrs)
# join the quantities for the forecast yrs
prj_join <- join_all(list(prj_SBB_SBBtrg,prj_F_Btrg,prj_SSB,prj_F,prj_Recr,prj_Catch), type='left') 
# recall and modify the target quantities table
target <- t. %>% dplyr::filter(year %in% yrs) %>% rename("overFished(SSBtrg)" = overFished ) %>% rename("underFished(SSBtrg)" = underFished )%>% rename("overFishing(Ftrg)" = overFishing  ) %>% rename("underFishing(Ftrg)" = underFishing )
# recall and modify the limit quantities table
lim <- t.lim %>% dplyr::filter(year %in% yrs) %>% rename("overFished(Blim)" = overFished ) %>% rename("underFished(Blim)" = underFished ) %>% rename(Flim_ratio = harvest)%>% rename(SBBlim_ratio = stock) %>% select(year,tac,SBBlim_ratio,"overFished(Blim)","underFished(Blim)")
trig <- t.trig %>% dplyr::filter(year %in% yrs) %>% rename("overFished(Btrig)" = overFished ) %>% rename("underFished(Btrig)" = underFished ) %>% rename(Ftrig_ratio = harvest)%>% rename(SBBtrig_ratio = stock) %>% select(year,tac,SBBtrig_ratio,"overFished(Btrig)","underFished(Btrig)")
# final table
tab_join <- join_all(list(prj_join,target,lim,trig), type='left')%>% rename(SBBtarget_ratio = stock) %>% rename(Ftarget_ratio = harvest)
tab_join <- tab_join[c(1,2,8,5,6,7,3,4,16,19,12,14,13,15,17,18,20,21)]
write.csv(tab_join , "Forecast_table_ICES.csv")





