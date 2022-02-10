#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
# R script to build BRP short-mse for Sole GSA17 ensemble model
# @author Henning Winker & Massimiliano Cardinale & Francesco Masnadi
# @email henning.winker@ec.europa.eu
# Licence: EUPL
# mse version: https://github.com/flr/mse/tree/FLombf
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
#devtools::install_github("henning-winker/FLRef", INSTALL_opts=c("--no-multiarch"), force=TRUE)
library(FLCore)
library(ggplotFL)
library(mse)
library(FLBRP)
library(gridExtra)
library(FLRef)
library(FLSRTMB)
library(mseviz)

library(readr)
library(ss3om)
library(rlist)
library(dplyr)
library(gtools)

##VENDACE gulf of bothnia
dir <- "C:/Users/f.masnadi/Desktop/BIOLOGIA della PESCA/DOTTORATO/overseas_SLU/VENDACE/"
setwd(dir)

load(file=file.path(dir,"/list_vendace.Rdata"))
stks <- FLStocks(runs); rm(runs)
names(stks) <- do.call(rbind,lapply(stks,function(x) x@name))
do.call(rbind,lapply(stks,function(x) x@name)) ==stks@names 

dir.create("MSE")
dir.create("MSE/plots")
setwd(paste0(dir,"MSE/plots"))
#-----------------------------------------------------------------
# 1:  Get the SRR
#-----------------------------------------------------------------
# Bevholt
sr <- FLSRs(lapply(stks, function(x) { 
  # Get Pars
  s=   x@benchmark[["s"]]
  R0 =  exp(x@benchmark[["R0"]])
  v=  x@benchmark[["B0"]]
  srr = as.FLSR(x,model=bevholtSV)
  return(ab(fmle(as.FLSR(x,model=bevholtSV),fixed=list(s=s,v=v,spr0=v/R0))))
}))  
srs = list(sr=sr)

# Blim for plot
lim <- (lapply(stks, function(x) { 
  # Get Pars
  v=  x@benchmark[["B0"]]
  BLIM01 <- v*0.1   # set Blim 10%B0
  BLIM015 <- v*0.15 # set Blim 15%B0
  BLIM020 <- v*0.20 # set Blim 20%B0
  return(list(BLIM01,BLIM015,BLIM020))
}))
dbBLim <- as.data.frame(t(as.data.frame(matrix(unlist(lim), ncol=length(lim) ))))
# plot S-R
plot(FLSRs(sr))+ theme(legend.position = "right")+ geom_vline(xintercept = median(dbBLim$V1), linetype = "longdash", colour = "red")+geom_vline(xintercept = median(dbBLim$V2), linetype = "longdash", colour = "blue")+geom_vline(xintercept = median(dbBLim$V3), linetype = "longdash", colour = "black")+
  xlim(0, 10000)
dev.print(jpeg,paste0("SR_allrunszomm.jpg"), width = 12, height = 8, res = 300, units = "in")


#-----------------------------------------------------------------
# 2:  Set up refpoints
#-----------------------------------------------------------------
# OM: "deterministic" MSY
omplot = Map(function(x,y){
  computeFbrp(x,y,proxy="msy",blim = 0.1)  # put Blim in agreement with the best percentage from the plot above
},  stks,sr)
#plot Fbrp for every runs
pdf(paste0("determMSYBRP.pdf"))
for(i in 1:length(stks)){
  p <-ploteq(omplot[[i]]) + ggtitle(paste0("BRP ",names(stks)[i]))
  print(p)
  #ggsave(paste0("determMSYBRP_",names(stks)[i] ,".jpg") , p)
}
dev.off()

#fb20
b20plot = Map(function(x,y){
  computeFbrp(x,y,proxy = "bx",x=20,blim = 0.1)  # put Blim in agreement with the best percentage from the plot above
},  stks,sr)
#fb25
b25plot = Map(function(x,y){
  computeFbrp(x,y,proxy = "bx",x=25,blim = 0.1)  # put Blim in agreement with the best percentage from the plot above
},  stks,sr)
#fb30
b30plot = Map(function(x,y){
  computeFbrp(x,y,proxy = "bx",x=30,blim = 0.1)  # put Blim in agreement with the best percentage from the plot above
},  stks,sr)
#fb35
b35plot = Map(function(x,y){
  computeFbrp(x,y,proxy = "bx",x=35,blim = 0.1)  # put Blim in agreement with the best percentage from the plot above
},  stks,sr)
#fb40
b40plot = Map(function(x,y){
  computeFbrp(x,y,proxy = "bx",x=40,blim = 0.1)  # put Blim in agreement with the best percentage from the plot above
},  stks,sr)

brps = list()
for(i in 1:length(omplot)){
  brps[[i]] = FLBRPs(list(  om=omplot[[i]],fb20=b20plot[[i]], fb25=b25plot[[i]], fb30=b30plot[[i]],fb35=b35plot[[i]],fb40=b40plot[[i]]  ))
}
names(brps) = stks@names

# compute "deterministic" MSY values
om.ref = Map(function(x,y){
  Fbrp(computeFbrp(x,y,proxy="msy",blim = 0.1))
},  stks,sr)
# to dataset
om.ref.db= as.data.frame(t(as.data.frame(om.ref,col.names = c(names(stks)))))
rownames(om.ref.db) <- names(stks)
om.ref.db$stock <- names(stks)
om.ref.db$type <- "detMSY"

# compute fb20 value
fb20.ref = Map(function(x,y){
  Fbrp(computeFbrp(x,y,proxy = "bx",x=20,blim = 0.1))
},  stks,sr)
# to dataset
om.fb20.db= as.data.frame(t(as.data.frame(fb20.ref,col.names = c(names(stks)))))
rownames(om.fb20.db) <- names(stks)
om.fb20.db$stock <- names(stks)
om.fb20.db$type <- "fb20"

# compute fb25 value
fb25.ref = Map(function(x,y){
  Fbrp(computeFbrp(x,y,proxy = "bx",x=25,blim = 0.1))
},  stks,sr)
# to dataset
om.fb25.db= as.data.frame(t(as.data.frame(fb25.ref,col.names = c(names(stks)))))
rownames(om.fb25.db) <- names(stks)
om.fb25.db$stock <- names(stks)
om.fb25.db$type <- "fb25"

# compute fb30 value
fb30.ref = Map(function(x,y){
  Fbrp(computeFbrp(x,y,proxy = "bx",x=30,blim = 0.1))
},  stks,sr)
# to dataset
om.fb30.db= as.data.frame(t(as.data.frame(fb30.ref,col.names = c(names(stks)))))
rownames(om.fb30.db) <- names(stks)
om.fb30.db$stock <- names(stks)
om.fb30.db$type <- "fb30"

# compute fb35 value
fb35.ref = Map(function(x,y){
  Fbrp(computeFbrp(x,y,proxy = "bx",x=35,blim = 0.1))
},  stks,sr)
# to dataset
om.fb35.db= as.data.frame(t(as.data.frame(fb35.ref,col.names = c(names(stks)))))
rownames(om.fb35.db) <- names(stks)
om.fb35.db$stock <- names(stks)
om.fb35.db$type <- "fb35"

# compute fb40 value
fb40.ref = Map(function(x,y){
  Fbrp(computeFbrp(x,y,proxy = "bx",x=40,blim = 0.1))
},  stks,sr)
# to dataset
om.fb40.db= as.data.frame(t(as.data.frame(fb40.ref,col.names = c(names(stks)))))
rownames(om.fb40.db) <- names(stks)
om.fb40.db$stock <- names(stks)
om.fb40.db$type <- "fb40"

mp.refs = rbind(om.ref.db,om.fb20.db,om.fb25.db,om.fb30.db,om.fb35.db,om.fb40.db)
#mp.refs = rbind(om.ref.db,om.fb20.db,om.fb25.db,om.fb30.db, .....)
refs = mp.refs[mp.refs$type=="detMSY",]

#-------------------------------------------------------------------------------
# 3 : setting dimensions
#-------------------------------------------------------------------------------
set.seed(1234)

it <- 250   # number of iterations #250                      
ny <- 30    # number of years to project #60 
nsqy <- 3 # numbers of year to average for internal forecasts

#------------------------------------------------------------------
# 4: Set up OM structure
#------------------------------------------------------------------

# Project...takes a second
foms <- FLStocks(lapply(stks,function(x){
  stf(x, ny, nsqy, nsqy)
})) 

# propagate
foms <- FLStocks(lapply(foms,function(x){
  propagate(x,it)
})) 


# Assign attributes of mpargs to all OMs
foms = FLStocks(Map(function(x,y){
  mpargs = list(
    y0=an(range(x)["minyear"]),
    dy=an(range(x)["maxyear"])-1,
    ay=an(range(x)["maxyear"])-1,
    iy=an(range(x)["maxyear"]),
    fy=an(range(x)["maxyear"])+ny,
    it=it,ny=ny,nsqy=nsqy)
  attr(y,"args") = mpargs
  y
},stks,foms)    
)

plot(foms[[1]]) #example plot run1
foms[[1]]@args

#-----------------------------------------------------------
#  5: Build OM recruitment residuals 
#----------------------------------------------------------
# Beverton Holt
set.seed(123)
om = Map(function(x,y,z){
  rho <-  z@benchmark[["rho"]] # or fixed to 0.341 (fishbase prior)
  sigmaR <- z@benchmark[["sigmaR"]]
  res =window(rec(x),end=x@args$fy+1)   
  res =ar1rlnorm(rho=rho,years=ac(range(x)["minyear"]:range(x)["maxyear"]),it=x@args$it,sdlog = sigmaR) 
  res_om = propagate(y,x@args$it)
  residuals(res_om) <- res
  res_om
},foms,sr,stks)

plot(om[[1]]@residuals) #example plot run1

#------------------------------------------------------------------
# 6: Set projection method
#------------------------------------------------------------------

  proj <- mseCtrl(method=fwd.om, args=list(maxF=3.))
  
  om =  Map(function(x,y,z){
    FLom(stock=x, sr=y,refpts=z,projection=proj)  }
    ,foms,om,om.ref)
  
 # Check that args are still there
  om[[1]]@stock@args
  om[[1]]@stock
  om[[1]]@sr
  om[[1]]@refpts

#-------------------------------------------------------
# 7: harvest control rule function 
#------------------------------------------------------
  
  # F-HCR will be converted to TACy+1 by mse::tac.is
  fhcr = function (stk, ftrg, blim, btrigger, fmin = 0.001,as.spr=FALSE, args, tracking) 
  {
    ay <- args$ay
    data_lag <- args$data_lag
    man_lag <- args$management_lag
    ssb <- ssb(stk)[, ac(ay - data_lag)]
    if(as.spr){ # Translate biomass into SB/R, R recruitment is geomean of last 3 years
      ssb <-  ssb(stk)[, ac(ay - data_lag)]/exp(mean(log(stock.n(stk)[1,ac((ay - data_lag-2):(ay - data_lag))])))
    }
    fout <- FLQuant(fmin, dimnames = list(iter = dimnames(ssb)$iter))
    fout[ssb >= btrigger] <- ftrg # Fish at Ftrg if > Btrigger
    inbetween <- (ssb < btrigger) & (ssb > blim) # here Blim at origin
    gradient <- (ftrg - fmin)/(btrigger - blim)
    fout[inbetween] <- (ssb[inbetween] - blim) * gradient + fmin
    ctrl <- fwdControl(year = ay + man_lag, quant = "fbar", 
                       value = c(fout))
    list(ctrl = ctrl, tracking = tracking)
  }
  
  # Example for Fb40% for Sole (actual ref point!)
  fb = Fbrp(computeFbrp(stks[[1]],sr[[1]],proxy="bx",x=40,blim = 0.1))
  fb
  plotICES(fmsy=fb[[1]],blim=fb[[3]],btrigger = fb[[2]]*1,kobe=T,obs=stks[[1]]) # trigger very close to the target (solo per specie mediterranee tipo sardina)
  dev.print(jpeg,paste0("ActualBRP_HCR_run1.jpg"), width = 6, height = 6, res = 300, units = "in")
# example of FMSY
  fb = Fbrp(computeFbrp(stks[[10]],sr[[10]],proxy="msy",blim = 0.1))
  plotICES(fmsy=fb[[1]],blim=fb[[3]],btrigger = fb[[2]]*1,kobe=T,obs=stks[[10]])
  dev.print(jpeg,paste0("MSYBRP_HCR_runCS.jpg"), width = 6, height = 6, res = 300, units = "in")
#-------------------------------------------------------
# 7: Parameterize Harvest Control Rule for MPs
#------------------------------------------------------
  #-------------------------------------------------------
  # 7.A: Parameterize Harvest Control Rule for MPs
  #------------------------------------------------------
  
  # OM: BRP deterministic MSY
  ref = mp.refs[mp.refs$type=="detMSY",]
  y= setNames(split(ref, seq(nrow(ref))), rownames(ref))
  det.MSY.hcr = lapply(y,function(x){
    mseCtrl(method=fhcr,args=list(fmin=0.001,ftrg=x$Fbrp,blim=0.001,btrigger=1*x$Btrg))
  }) 
  
  # fb20
  ref = mp.refs[mp.refs$type=="fb20",]
  y= setNames(split(ref, seq(nrow(ref))), rownames(ref))
  # trigger 0.6
  fb20.06.hcr = lapply(y,function(x){
    mseCtrl(method=fhcr,args=list(fmin=0.001,ftrg=x$Fbrp,blim=0.001,btrigger=0.6*x$Btrg))
  }) 
  # trigger 0.8
  fb20.08.hcr = lapply(y,function(x){
    mseCtrl(method=fhcr,args=list(fmin=0.001,ftrg=x$Fbrp,blim=0.001,btrigger=0.8*x$Btrg))
  }) 
  # trigger 1
  fb20.1.hcr = lapply(y,function(x){
    mseCtrl(method=fhcr,args=list(fmin=0.001,ftrg=x$Fbrp,blim=0.001,btrigger=1*x$Btrg))
  })
  
  # fb25
  ref = mp.refs[mp.refs$type=="fb25",]
  y= setNames(split(ref, seq(nrow(ref))), rownames(ref))
  # trigger 0.6
  fb25.06.hcr = lapply(y,function(x){
    mseCtrl(method=fhcr,args=list(fmin=0.001,ftrg=x$Fbrp,blim=0.001,btrigger=0.6*x$Btrg))
  }) 
  # trigger 0.8
  fb25.08.hcr = lapply(y,function(x){
    mseCtrl(method=fhcr,args=list(fmin=0.001,ftrg=x$Fbrp,blim=0.001,btrigger=0.8*x$Btrg))
  }) 
  # trigger 1
  fb25.1.hcr = lapply(y,function(x){
    mseCtrl(method=fhcr,args=list(fmin=0.001,ftrg=x$Fbrp,blim=0.001,btrigger=1*x$Btrg))
  }) 
  
  # fb30
  ref = mp.refs[mp.refs$type=="fb30",]
  y= setNames(split(ref, seq(nrow(ref))), rownames(ref))
  # trigger 0.6
  fb30.06.hcr = lapply(y,function(x){
    mseCtrl(method=fhcr,args=list(fmin=0.001,ftrg=x$Fbrp,blim=0.001,btrigger=0.6*x$Btrg))
  }) 
  # trigger 0.8
  fb30.08.hcr = lapply(y,function(x){
    mseCtrl(method=fhcr,args=list(fmin=0.001,ftrg=x$Fbrp,blim=0.001,btrigger=0.8*x$Btrg))
  }) 
  # trigger 1
  fb30.1.hcr = lapply(y,function(x){
    mseCtrl(method=fhcr,args=list(fmin=0.001,ftrg=x$Fbrp,blim=0.001,btrigger=1*x$Btrg))
  }) 
  
  # fb35
  ref = mp.refs[mp.refs$type=="fb35",]
  y= setNames(split(ref, seq(nrow(ref))), rownames(ref))
  # trigger 0.6
  fb35.06.hcr = lapply(y,function(x){
    mseCtrl(method=fhcr,args=list(fmin=0.001,ftrg=x$Fbrp,blim=0.001,btrigger=0.6*x$Btrg))
  }) 
  # trigger 0.8
  fb35.08.hcr = lapply(y,function(x){
    mseCtrl(method=fhcr,args=list(fmin=0.001,ftrg=x$Fbrp,blim=0.001,btrigger=0.8*x$Btrg))
  }) 
  # trigger 1
  fb35.1.hcr = lapply(y,function(x){
    mseCtrl(method=fhcr,args=list(fmin=0.001,ftrg=x$Fbrp,blim=0.001,btrigger=1*x$Btrg))
  }) 
  
  # fb40
  ref = mp.refs[mp.refs$type=="fb40",]
  y= setNames(split(ref, seq(nrow(ref))), rownames(ref))
  # trigger 0.6
  fb40.06.hcr = lapply(y,function(x){
    mseCtrl(method=fhcr,args=list(fmin=0.001,ftrg=x$Fbrp,blim=0.001,btrigger=0.6*x$Btrg))
  }) 
  # trigger 0.8
  fb40.08.hcr = lapply(y,function(x){
    mseCtrl(method=fhcr,args=list(fmin=0.001,ftrg=x$Fbrp,blim=0.001,btrigger=0.8*x$Btrg))
  }) 
  # trigger 1
  fb40.1.hcr = lapply(y,function(x){
    mseCtrl(method=fhcr,args=list(fmin=0.001,ftrg=x$Fbrp,blim=0.001,btrigger=1*x$Btrg))
  }) 
  
  
  #-------------------------------------------------------
  # 7.B ## HCRs
  #------------------------------------------------------
  # om mp (self-test)
  sc.ctrl <- list(
    detMSY = lapply(det.MSY.hcr,function(x){ # deterministic MSY
      mpCtrl(list(
        est = mseCtrl(method=perfect.sa),
        hcr  = x,
        iysy = mseCtrl(method=tac.is)
      ))}),
    # fb20
      fb20.bt06 = lapply(fb20.06.hcr,function(x){ # fb20 trigger 0.6
      mpCtrl(list(
        est = mseCtrl(method=perfect.sa),
        hcr  = x,
        iysy = mseCtrl(method=tac.is)
      ))}),
    fb20.bt08 = lapply(fb20.08.hcr,function(x){ # fb20 trigger 0.8
      mpCtrl(list(
        est = mseCtrl(method=perfect.sa),
        hcr  = x,
        iysy = mseCtrl(method=tac.is)
      ))}),
    fb20.bt1 = lapply(fb20.1.hcr,function(x){ # fb20 trigger 1
      mpCtrl(list(
        est = mseCtrl(method=perfect.sa),
        hcr  = x,
        iysy = mseCtrl(method=tac.is)
      ))}),
    # fb25
    fb25.bt06 = lapply(fb25.06.hcr,function(x){ # fb25 trigger 0.6
      mpCtrl(list(
        est = mseCtrl(method=perfect.sa),
        hcr  = x,
        iysy = mseCtrl(method=tac.is)
      ))}),
    fb25.bt08 = lapply(fb25.08.hcr,function(x){ # fb25 trigger 0.8
      mpCtrl(list(
        est = mseCtrl(method=perfect.sa),
        hcr  = x,
        iysy = mseCtrl(method=tac.is)
      ))}),
    fb25.bt1 = lapply(fb25.1.hcr,function(x){ # fb25 trigger 1
      mpCtrl(list(
        est = mseCtrl(method=perfect.sa),
        hcr  = x,
        iysy = mseCtrl(method=tac.is)
      ))}),
    # fb30
    fb30.bt06 = lapply(fb30.06.hcr,function(x){ # fb30 trigger 0.6
      mpCtrl(list(
        est = mseCtrl(method=perfect.sa),
        hcr  = x,
        iysy = mseCtrl(method=tac.is)
      ))}),
    fb30.bt08 = lapply(fb30.08.hcr,function(x){ # fb30 trigger 0.8
      mpCtrl(list(
        est = mseCtrl(method=perfect.sa),
        hcr  = x,
        iysy = mseCtrl(method=tac.is)
      ))}),
    fb30.bt1 = lapply(fb30.1.hcr,function(x){ # fb30 trigger 1
      mpCtrl(list(
        est = mseCtrl(method=perfect.sa),
        hcr  = x,
        iysy = mseCtrl(method=tac.is)
      ))}),
    # fb35
    fb35.bt06 = lapply(fb35.06.hcr,function(x){ # fb35 trigger 0.6
      mpCtrl(list(
        est = mseCtrl(method=perfect.sa),
        hcr  = x,
        iysy = mseCtrl(method=tac.is)
      ))}),
    fb35.bt08 = lapply(fb35.08.hcr,function(x){ # fb35 trigger 0.8
      mpCtrl(list(
        est = mseCtrl(method=perfect.sa),
        hcr  = x,
        iysy = mseCtrl(method=tac.is)
      ))}),
    fb35.bt1 = lapply(fb35.1.hcr,function(x){ # fb35 trigger 1
      mpCtrl(list(
        est = mseCtrl(method=perfect.sa),
        hcr  = x,
        iysy = mseCtrl(method=tac.is)
      ))}),
    # fb40
    fb40.bt06 = lapply(fb40.06.hcr,function(x){ # fb40 trigger 0.6
      mpCtrl(list(
        est = mseCtrl(method=perfect.sa),
        hcr  = x,
        iysy = mseCtrl(method=tac.is)
      ))}),
    fb40.bt08 = lapply(fb40.08.hcr,function(x){ # fb40 trigger 0.8
      mpCtrl(list(
        est = mseCtrl(method=perfect.sa),
        hcr  = x,
        iysy = mseCtrl(method=tac.is)
      ))}),
    fb40.bt1 = lapply(fb40.1.hcr,function(x){ # fb40 trigger 1
      mpCtrl(list(
        est = mseCtrl(method=perfect.sa),
        hcr  = x,
        iysy = mseCtrl(method=tac.is)
      ))})
  ) # end of hcr list
  
  
  #------------------------------------------------------------
  # Save OM files
  setwd(paste0(dir,"MSE"))
  save(om,stks,sc.ctrl,mp.refs,srs,brps,omplot,om.ref, file="omVENDACE2022.rdata")
  

  #====================================================================
  # RUN 
  #====================================================================
  dir <- "C:/Users/f.masnadi/Desktop/BIOLOGIA della PESCA/DOTTORATO/overseas_SLU/VENDACE/"
  setwd(paste0(dir,"MSE"))
  output = "runs" 
  dir.create(output,showWarnings = F)
  load(file="omVENDACE2022.rdata",verbose=T)  # load the file if already produced
  
  stks@names==names(om)
  names(sc.ctrl$detMSY) == names(om)
  om[[1]]@stock@args 
  #-------------------------------------------------------
  # 8: Set up implementation error
  #------------------------------------------------------
  iem <- FLiem(method=noise.iem, args=list(fun="rlnorm", mean=0.1, sd=0.15, multiplicative=TRUE))  #error between expected TAC and actual landing, if you do not have true TAC put plausible mean and sd base on prior experience   
  
  #---!!!---!!!---!!!---!!!---!!!---!!!---!!!---!!!---!!!---!!!
  #
  #   WARNING BEFORE RUNNING THIS LOOP:                     !!!
  #
  # depending on the number of runs and HCR,                !!!
  # this step can be time consuming (hours to days)
  #---!!!---!!!---!!!---!!!---!!!---!!!---!!!---!!!---!!!---!!!
  for(s in 1:length(stks)) {
    system.time({
    # Only run if not exist already
    if(file.exists(paste0(output,"/mps.",om[[s]]@stock@name,".rdata"))==F){ # F to re-write
      cat("Running",om[[s]]@stock@name, "__")
      mps = FLStocks(lapply(sc.ctrl,function(x){ 
        out = mp(om[[s]], ctrl=x[[s]], args=om[[s]]@stock@args,iem=iem,parallel=F)
        stock(out@om)
      }))
      save(mps,file=paste0(output,"/mps.",om[[s]]@stock@name,".rdata"))
    }})
  } # end of loop 
  

#====================================================================
# 9.Eval 
#====================================================================
  setwd(dir <-"C:/Users/f.masnadi/Desktop/BIOLOGIA della PESCA/DOTTORATO/overseas_SLU/VENDACE/MSE")
# load stocks
load("omVENDACE2022.rdata",verbose=T) # load the file if already produced

# Set run folder
runs.dir = "runs"
yr.eval = 2040:2050 # 10 years

# file list from MSE run folder defined above
run.names = mixedsort(list.files(runs.dir))
# To use only ceryain model use [1:3] after list.files(runs.dir) #><> use model 1-3

# list of file names
fileL = split(run.names,seq(length(run.names)))

# Load runs (Big File)
runs = lapply(fileL,function(x){load(
  file.path(runs.dir,x))
  mps})

# Name stks uniquely
names(runs) = str_remove(str_remove(fileL , "mps.") , ".rdata")

#><>><>><>><>><>><>><>><>><>><>
#> Performance statistics
#><>><>><>><>><>><>><>><>><>><>
library(mseviz)

# Define Metrics for Evaluation
metrics <- list(SB = ssb, F = fbar, C = landings,TC=catch,Rec=rec)

Catch = list(C = list(~yearMeans(C),name="mean(C)",desc="Average Catch"))

# compute TRUE values in a list of FLPars
# Loop through stks 
rpts <- FLPars(Map(function(x,y){
  
  # Get Pars
  s= x@benchmark[["s"]]
  R0 = exp(x@benchmark[["R0"]])
  v=  x@benchmark[["B0"]]
  sr = ab(fmle(as.FLSR(x,model=bevholtSV),fixed=list(s=s,v=v,spr0=v/R0)))
  
  #><> set Blim
  rp= Fbrp(computeFbrp(x,sr,proxy="msy",blim=0.1))
  
  # define MSY as median yield under true Fmsy for each OM
  msy = performance(y[1], statistics=Catch,metrics=metrics, years=list(yr.eval))$data
  FLPar(
    SBmsy = rp["Btrg"], 
    Fmsy =  rp["Fbrp"], 
    Blim =  rp["Blim"],
    MSY =  msy,
    SB0 = rp["B0"])
},x=stks,y=runs))

# rpts@names = paste0("Run",1:length(run.names))

# Performance Statistics
inds=list(
  a = list(~yearMeans(F/Fmsy),name="F/Fmsy",desc="Median annual F/Fmsy"),
  b = list(~yearMeans(SB/SBmsy),name="B/Bmsy",desc="Median annual B/Bmsy"),
  c = list(~yearMeans(C/MSY),name="Catch/MSY",desc="Mean Catch/MSY over years"),
  d = list(~yearMeans(iav(C)),name="AAV",desc="Average annual variation in catches"),
  e = list(~apply(iterMeans((SB/Blim) < 1),1,max),name="P3(B<Blim)",desc="Probability that SSB < Blim Hockey-Stick"),
  x = list(~yearMeans(SB/SBmsy),name="SSB/SSB[MSY]",desc="Average annual SSB/SSBmsy"),
  y = list(~yearMeans(F/Fmsy),name="F/F[MSY]",desc="Average annual F/Fmsy")
)

perfrun = Map(function(x,y,z){
  pb = performance(x, refpts = y, statistics=inds, # Now statistics
                   metrics = metrics, years=list(yr.eval))
  mp = pb$mp
  pb = pb[,!"mp"]
  # Add things
  pb$run = paste0("Run",z) # allow to identify models
  pb$mp = mp # Add mp back in last column
  return(pb)
},x=runs,y=rpts,z=seq(length(rpts@names)))

# Combine in data table
perf = rbindlist(c(perfrun),idcol = "stock")

#Median AAV
muAAV=median(perf[perf$name=="AAV",]$data)


###########################
## Save performance plots 
###########################
png(paste0("plots/mseperf_joint.png"), width = 8, height =9., res = 250, units = "in")
pbp = plotBPs(perf,c("a","b","c","e","d"),
              target = c(b=1,c=1,d=muAAV),
              limit= c( e=0.05,c=0.95),
              yminmax = c(0.05, 0.95))+theme_bw()+
  facet_wrap(~name,scales = "free_y",ncol=2)+
  ggtitle(paste0("Performance: All runs"))+
  ylab("Performance statistics")+
  theme(axis.text.x=element_blank())+xlab("Candidates")
pbp
dev.off()

##########################
## Save single trj plots 
##########################
nm = unique(perf$mp)

library(scales) 
gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}
cols = c("black",gg_color_hue(11))

pdf("plots/shortcut_ouput.pdf")
# Loop trj plots 
for(i in 1:length(names(runs))){ 
run = FLStocks(c(FLStocks(om=stks[[i]]),runs[[i]]))
run@names = c("stock",paste(nm))
p1 = plot(run,metrics=list(Landings=landings,F=fbar,SSB=ssb,Recruitment=rec))+theme_bw()+ 
  geom_flpar(data=FLPars(SSB  = FLPar(Bmsy=median(an(rpts[[i]]["SBmsy"])),B0=median(an(rpts[[i]]["SB0"])),Blim=median(an(rpts[[i]]["Blim"]))),
                         F = FLPar(Fmsy=median(an(rpts[[i]]["Fmsy"]))),Landings=FLPar(MSY=median(an(rpts[[i]]["MSY"])))),x=c(1995,rep(1976,4)),colour=c("darkgreen","blue","red","darkgreen","darkgreen"))+
  facet_wrap(~qname, scales="free")+ggtitle(paste0("  Run",i))+
  scale_colour_manual(values=cols) # Make sure colors stay consistent by using stock = "black"
print(p1)
}
dev.off()

##################
# Save kobe plots 
##################
png(paste0("plots/msekobe.png"), width = 8, height =7, res = 250, units = "in")
kbcex =function(){theme(plot.title = element_text(size=10),
                        legend.key.size = unit(0.3, 'cm'), #change legend key size
                        legend.key.height = unit(0.4, 'cm'), #change legend key height
                        legend.key.width = unit(0.4, 'cm'), #change legend key width
                        legend.text = element_text(size=10)) #change legend text font size
}
kobeMPs(perf,x="x", y="y")+
  ylab(expression(F/F[MSY]))+xlab(expression(B/B[MSY]))+ylim(0,2.5)+kbcex()+theme()
dev.off()


#############################################


