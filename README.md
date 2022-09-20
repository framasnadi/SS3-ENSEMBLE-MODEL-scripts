# SS3 ENSEMBLE MODEL scripts

This repo contains all the script that has been developed during the work done for the **Ph.D. thesis of Francesco Masnadi** to support the methods described in Chapter 4:_"Stock assessment of common sole in GSA17: an ensemble approach using data-rich model"_

PhD Project in _"Evaluation and management of demersal stocks in the Adriatic: the case of common sole (Solea solea, L.) in the northern and central Adriatic Sea"_ conceived under the International Ph.D. Program _“Innovative Technologies and Sustainable Use of Mediterranean Sea Fishery and Biological Resources”_ (www.FishMed-PhD.org)

- Supervisor: Giuseppe Scarcella (CNR-IRBIM)

- Co-supervisors: Stefano Goffredo (UNIBO); Massimiliano Cardinale (SLU)


File description:

1) ```Ensemble_grid.R``` : Run all scenarios at once and plot results, do diagnostic, do Ensemble on final weighted grid

2) ```Ensemble _Forecast_grid_run_v2.R``` : run all forecast scenarios models at once and plot results

3) ```ss3om_fromSS3toFLR.R``` : convert multiple SS3 models to flr (FLRStock)

4) ```MSE_ensemble_loop.R``` : R script to build BioRefPoint short-mse for ensemble model (based on ICES WKREF1&2 scripts)

5) ```app.R``` :  R script to build Shiny APP showing ensemble model procedure and output for Sole in GSA17 (https://framasnadi.shinyapps.io/AppSOL)

* Authors:                                                                                          
Francesco Masnadi (CNR-IRBIM & UNIBO, Ancona)                                                     
Massimiliano Cardinale (SLU Aqua, Lysekil)                                                        
mainly based on functions developed by Hennig Winker (JRC-EC), Felipe Carvalho (NOAA) and Iago Mosqueira (WUR-IMARES)

 
