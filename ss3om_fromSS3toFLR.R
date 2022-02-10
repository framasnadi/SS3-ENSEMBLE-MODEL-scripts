#https://github.com/iagomosqueira/ss3om/
# Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="TRUE") 
# options(devtools.install.args = "--no-multiarch")
#devtools::install_github("flr/ss3om", INSTALL_opts=c("--no-multiarch"), force=TRUE)
library(readr)
library(ss3om)
library(rlist)
library(dplyr)
library(tidyverse)
library(gtools)
###Working directory

##Adriatic sole
dir <- "C:/Users/f.masnadi/Desktop/Stock Assessment/SS3/SOLEA_SS3/SOLE_2021/BENCHMARK_update/"
setwd(dir)
FLRdir <- "FLRfolder"
readme = paste0("run",1:18, "_update") # 1:N Number of run

# LOOP: transfer all runs in FLR
for(i in 1:length(readme)){
  out = SS_output(dir=file.path(readme[i]))
  rstk <- readFLSss3(dir=file.path(readme[i]))
#Checks and settings
rstk@name <- readme[i]
range(rstk)["minfbar"] <-1
range(rstk)["maxfbar"] <-4
R0 <- out$parameters$Value[out$parameters$Label=="SR_LN(R0)"]
s <- out$parameters$Value[out$parameters$Label=="SR_BH_steep"]
sigmaR <- out$parameters$Value[out$parameters$Label=="SR_sigmaR"]
rho <- out$parameters$Value[out$parameters$Label=="SR_autocorr"]
B0 <- out$derived_quants$Value[out$derived_quants$Label=="SSB_unfished"]
##Add attributes to the FLstock
attr(rstk,"benchmark") = c(rho=rho,sigmaR=sigmaR,s=s,sd.logit.s=0.3, R0=R0, B0=B0)
save(rstk,file=paste0(FLRdir,"/FLR_Sole_",readme[i],".RData"))
}

# file list from FLR run folder defined above
fileL = split(mixedsort(list.files(FLRdir,pattern ="FLR_Sole_")),seq(length(readme)))
# Load runs (Big File)
runs = lapply(fileL,function(x){load(file.path(FLRdir,x))
  rstk})
# names(runs) = str_remove(str_remove(str_remove(fileL , "FLR_") , ".RData"),"_update")


# save final ensemble grid list with all runs
save(runs,file=file.path(FLRdir,"list_sol_WGSAD2021B.Rdata"))


#####
load(file=file.path(FLRdir,"list_sol_WGSAD2021B.Rdata"))

