##### A L L    P A C K A G E S    N E E D E D  #########



#### UN -COMMENT BELOW IF NEED TO INSTALL ANY/ALL OF THESE PACKAGES USED IN CODE ###
#install.packages(pkgs = c("ggplot2","dplyr","tidyr","stringr","gplots","plotrix"),dependencies = T)
#install.packages("gridExtra")
#install.packages("viridis")
#install.packages('multipanelfigure')
#install.packages('magrittr')
#install.packages('MASS')
#install.packages('ggpubr', dependencies = T)
#install.packages('corpcor')
#install.packages('matrixStats')
#install.packages('readr')
#install.packages('vctrs')
#install.packages('colorRamps')
#install.packages('fitdistrplus')

library('MASS')#
library(viridis)#
library(ggplot2) # ggplot() for plotting
library(dplyr) # data reformatting
library(tidyr) # data reformatting
library(stringr) # string manipulation
library(gridExtra)#
library(multipanelfigure)#
library(magrittr)#
library(ggpubr)#
library("corpcor")
library('matrixStats')
library('readr')
library('vctrs')
library('colorRamps')
library('fitdistrplus')




#### FUNCTIONS NEEDED FOR WHOLE CODE ####
### make sure to run these before running any other section of code

## function for running the evolution sim
# inputs: svalALlTheoMuts (col 1 is the DFE1 col2 is the DFE2), numGens, num iterations, popsize, U, yplotmin, yplotmax, wantPLot (a T/F for if want fit traj plot)
# outputs: graph of fit traj in home/non-home (if wantPLot == TRUE), list of (fit traj values, variance over time, covariance over time, r2, d22, d12)

simEvo <- function(sValallTheoMuts, numGens, numIterations, popsize, U, yplotmin, yplotmax, wantPlot)
{
  
  sValAllMuts = sValallTheoMuts
  numDrugs = 2
  binSize = 0.1
  minS = min(sValAllMuts)
  maxS = max(sValAllMuts)
  numBins = length(seq(floor(minS),ceiling(maxS*10)/10, by = binSize)) - 1
  binSeq = seq(floor(minS),ceiling(maxS*10)/10, by = binSize) # length = 1+numBins
  lastNegBin = which.min(abs(binSeq)) - 1
  
  
  allHomeSValLists= list()
  allHomeGaussianParameters = list()
  for(home in 1:numDrugs)
  {
    sVallistALLRespDrugs = list()
    GaussianParametersALLRespDrugs = list()
    for(resp in 1:numDrugs)
    {
      
      singleDrugList = list() # 1 entry per bin
      singleDrugGaussianParameters = list()
      for(bin in 1:numBins)
      {
        binMIN = binSeq[bin]
        binMAX = binSeq[bin+1]
        
        
        sValsinBin = (sValAllMuts[,home]>=binMIN & sValAllMuts[,home]<=binMAX)
        
        respSvals = sValAllMuts[,resp][sValsinBin!=0]
        
        singleDrugList[[bin]] = respSvals
        
        
        #####  NOW: SAVE THE Gaussian mean AND sd PARAMETERS FOR ALL FITS YOU CAN 
        
        if(length(respSvals) >= 2 )
        {
          dat <- data.frame(x=c(respSvals)) 
          fitW <- fitdistr(dat$x, densfun = "normal")
          
          singleDrugGaussianParameters[[bin]] = fitW$estimate
          
        }else{
          singleDrugGaussianParameters[[bin]] = 0 
        }
        
        
      }
      
      sVallistALLRespDrugs[[resp]] = singleDrugList
      GaussianParametersALLRespDrugs[[resp]] = singleDrugGaussianParameters
    }
    allHomeSValLists[[home]] = sVallistALLRespDrugs
    allHomeGaussianParameters[[home]] = GaussianParametersALLRespDrugs
  }
  
  
  
  #==============================================================================================================================================#
  
  beneficialCutoffVec = numeric(numDrugs)  ## in theo JDFEs, no meaurement error/no need for cuttoff other than 0
  deleteriousCutoffVec = numeric(numDrugs)
  
  
  
  
  allGaussianParameters = list()
  for( i in 1:numDrugs)
  {
    fit = fitdistr(sValAllMuts[,i], "normal")
    mean = fit$estimate[1]
    sd = fit$estimate[2]
    
    
    allGaussianParameters[[i]] = c(mean,sd)
  }
  
  
  
  
  
  ### calc R and D values ###
  r1Easy= array(0, c(numDrugs,numDrugs))
  r2Easy = array(0, c(numDrugs,numDrugs))
  D11Easy = array(0, c(numDrugs,numDrugs))
  D22Easy = array(0, c(numDrugs,numDrugs))
  D12Easy = array(0, c(numDrugs,numDrugs))
  
  for(home in 1:numDrugs)
  {
    ## make all negative svals in home which are negative 0 so dont dont contribute to calc!
    homeSvals = sValAllMuts[,home]*(sValAllMuts[,home]>=0)
    
    for(resp in 1:numDrugs)
    {
      r1Easy[home,resp] = mean((homeSvals)^2)
      r2Easy[home,resp] = mean((homeSvals)*(sValAllMuts[,resp]))
      D11Easy[home,resp] = mean((homeSvals)^3)
      D12Easy[home,resp] = mean((homeSvals)^2*(sValAllMuts[,resp]))
      D22Easy[home,resp] = mean((homeSvals)*(sValAllMuts[,resp])^2)
    }
  }
  
  
  
  
  #=====================================================================================================================#
  #========================================================================================================#
  #========================================================================================================#
  
  #### RUN SIMULATION ####
  
  slopeArray = array(data = 0, dim = c(numDrugs, numDrugs))
  covarArray = array(data = 0, dim = c(numDrugs, numDrugs))
  varArray = array(data = 0, dim = c(numDrugs, numDrugs))
  
  
  numBenMutsAllDrugs = numeric(numDrugs)
  benMutEffects = numeric(numDrugs)
  allUbs = array() 
  
  
  allCovars = array(data = 0, dim = c(numGens, numDrugs, numDrugs))
  allCors = array(data = 0, dim = c(numGens, numDrugs, numDrugs))
  for( yay in 1)
  {
    enviroEvolveID = yay
    
    
    allWbarsAllIterations = array(data = 0, dim = c(numGens, (numDrugs), numIterations))
    fullDrugVariDataSet = array(data = 0, dim = c(numGens, numDrugs, numIterations))
    
    for(w in 1:numIterations)
    { 
      
      print('iteration number')
      print(w)
      
      startIDarray = array(data = c(1, popsize, 0.01), dim = c(1, 3))
      mutFitnesses = array(data = 0.01, dim = c(1, numDrugs))
      allWbars = array(data = 0, dim = c(numGens, (numDrugs)))
      itVar = array(data = 0, dim = c(numGens, numDrugs))
      for( x in 1:numGens)
      {
        numMuts = dim(startIDarray)[1]
        
        #### Selection and Drift ####
        
        Wbar = sum(startIDarray[ , 3] * (startIDarray[ ,2] / popsize))
        
        # generate probability vector
        probvec = numeric(numMuts)
        for( i in 1:numMuts)
        {
          
          if((startIDarray[i,2] / popsize) + (startIDarray[i,2] / popsize)* (startIDarray[i,3] - Wbar) >0)
          {
            
            probvec[i] = (startIDarray[i,2] / popsize) + (startIDarray[i,2] / popsize)* (startIDarray[i,3] - Wbar)
          }else{
            probvec[i] = 0
          }
        }
        
        # generate IDarray (with only sel/drift) for next gen
        nextGenIDarray = array(data = 0, dim = c(dim(startIDarray)[1], dim(startIDarray)[2]))
        nextGenOffspring = rmultinom(1, size = popsize, prob = probvec)
        
        
        for( i in 1:numMuts)
        {
          nextGenIDarray[i, 1] = startIDarray[i, 1]
          nextGenIDarray[i, 2] = nextGenOffspring[i]
          nextGenIDarray[i, 3] = startIDarray[i, 3]
        }
        
        
        
        ## Mutations ##
        M = rpois(1, popsize*U)
        parentIDs = sample(nextGenIDarray[ , 1], M, prob = (nextGenIDarray[ ,2]/popsize), replace = TRUE)
        
        for( i in parentIDs)
        {
          # pick s in home
          homeS = rnorm(1, mean = allGaussianParameters[[enviroEvolveID]][1], sd = allGaussianParameters[[enviroEvolveID]][2]) 
          
          if(homeS > beneficialCutoffVec[enviroEvolveID])
          {
            numBenMutsAllDrugs[enviroEvolveID] = numBenMutsAllDrugs[enviroEvolveID] + 1
            benMutEffects[enviroEvolveID] = benMutEffects[enviroEvolveID] + homeS
          }
          
          
          killChance = sample(c(0,1), 1, prob = c(0, 1)) # dont allow chance of death
          if(killChance == 0)
          {
            svec = numeric(numDrugs) - 1
          }else{
            
            # find bin of homeS
            if(homeS >= -1 & homeS <= -0.9)
            {
              binHome = 1
            }else if(homeS >= -0.9 & homeS < -0.8){
              binHome = 2
            }else if(homeS >= -0.8 & homeS < -0.7){
              binHome = 3
            }else if(homeS >= -0.7 & homeS < -0.6){
              binHome = 4
            }else if(homeS >= -0.6 & homeS < -0.5){
              binHome = 5
            }else if(homeS >= -0.5 & homeS < -0.4){
              binHome = 6
            }else if(homeS >= -0.4 & homeS < -0.3){
              binHome = 7
            }else if(homeS >= -0.3 & homeS < -0.2){
              binHome = 8
            }else if(homeS >= -0.2 & homeS < -0.1){
              binHome = 9
            }else if(homeS >= -0.1 & homeS < 0){
              binHome = 10
            }else if(homeS >= 0 & homeS < 0.1){
              binHome = 11
            }else if(homeS >= 0.1 & homeS < 0.2){
              binHome = 12
            }else if(homeS >= 0.2){
              binHome = 13
            }
            
            #### make svec ####
            
            svec = numeric(numDrugs)
            svec[enviroEvolveID] = homeS
            for(respit in 1:numDrugs)
            {
              if(respit != enviroEvolveID) #only need to sample new value from NON homes
              {
                if(length(allHomeGaussianParameters[[enviroEvolveID]][[respit]][[binHome]]) == 2) # only fit dist to sets of svalues with more than 10 vals
                {
                  
                  svec[respit] = rnorm(1, mean = allHomeGaussianParameters[[enviroEvolveID]][[respit]][[binHome]][1], sd = allHomeGaussianParameters[[enviroEvolveID]][[respit]][[binHome]][2]) 
                  
                  
                }else{
                  
                  
                  svec[respit] = 0
                  
                  
                  
                  
                }
                
              }
            }
            
            
            
            
            
          }
          
          
          
          
          
          #remove indiv from its parent class if can be removed
          if(nextGenIDarray[i,2] != 0)
          {
            nextGenIDarray[i, 2] = nextGenIDarray[i, 2]- 1
            
            #make new row for this new mutant; mutliplicative fitness change (mult parent fit by (1+s), so dec if s neg or inc if s pos)
            newRow = c(nextGenIDarray[dim(nextGenIDarray)[1], 1]+ 1, 1, startIDarray[i,3] + svec[enviroEvolveID])
            
            #bind newrow with nextgenIDarray
            nextGenIDarray = rbind(nextGenIDarray, newRow)
            
            #add row with fitness of new mutant in all environments to mutFitnesses
            newFit = mutFitnesses[i , ]+ svec
            mutFitnesses = rbind(mutFitnesses, newFit)
          }
          
          
        }
        
        output = list(nextGenIDarray, mutFitnesses)
        
        
        #====================================================================================================================#
        
        startIDarray = output[[1]]
        mutFitnesses = output[[2]]
        
        
        # find and save Wbars in all enviros 
        
        for( n in 1:(numDrugs))
        {
          allWbars[x, n] = sum(mutFitnesses[ , n] * (startIDarray[ ,2] / popsize))
        }
        
        # remove all muttypes with frequency of 0 from idarray and mutfitarray
        toRemove = c()
        if(length(startIDarray) > 3) # only have chance to remove if have more than 1 line in the startIDarray
        {
          for(q in 1:(dim(startIDarray)[1]))
          {
            if(startIDarray[q,2] == 0)
            {
              toRemove = c(toRemove, q)
              
            }
          }
        }
        
        
        if(length(toRemove) > 0)
        {
          startIDarray = startIDarray[-toRemove, ]
          mutFitnesses = mutFitnesses[ -toRemove, ]
          
          
          
          if(length(startIDarray) == 3) 
          {
            startIDarray = array(data = c(startIDarray), dim = c(1,3))
            mutFitnesses = array( data = c(mutFitnesses), dim= c(1, numDrugs))
          }
          
          count = 0
          
          for(q in 1:(dim(startIDarray)[1]))
          {
            count = count + 1
            startIDarray[q, 1] = count
          }
          
        }
        
        
      }
      
      allWbarsAllIterations[ , , w] = allWbars
    }
    
    # variance across all iterations for this evolve enviro in this data set
    avgDrugVariDataSet = array(data = 0, dim = c(numGens, numDrugs))
    for( i in 1:numGens)
    {
      for(y in 1:(numDrugs))
      {
        avgDrugVariDataSet[i,y] = var((allWbarsAllIterations[i,y, ])) #
      }
    }
    
    
    #COVARIANCE
    avgDrugCoVariDataSet = array(data = 0, dim = c(numGens, numDrugs))
    for( i in 1:numGens)
    {
      for(y in 1:(numDrugs))
      {
        avgDrugCoVariDataSet[i,y] = cov((allWbarsAllIterations[i,y, ]), (allWbarsAllIterations[i,enviroEvolveID, ])) 
      }
    }
    
    allCovars[,,yay] = avgDrugCoVariDataSet
    
    # CORREALTIONS 
    avgDrugCoriDataSet = array(data = 0, dim = c(numGens, numDrugs))
    for( i in 1:numGens)
    {
      for(y in 1:(numDrugs))
      {
        avgDrugCoriDataSet[i,y] = cor((allWbarsAllIterations[i,y, ]), (allWbarsAllIterations[i,enviroEvolveID, ]))
      }
    }
    
    
    
    
    # average Wbar across all iterations for this evolve enviro in this data set 
    averageAllWbar = array(data = 0, dim = c(numGens, (numDrugs)))
    for( i in 1:numGens)
    {
      for(y in 1:(numDrugs))
      {
        averageAllWbar[i,y] = mean(allWbarsAllIterations[i,y, ])
      }
    }
    
    
    
    
    
    
    ##### PLot results for this EnviroEvolveID ######
    WbarData = averageAllWbar
    varData = avgDrugVariDataSet
    covarData = avgDrugCoVariDataSet
    corData = avgDrugCoriDataSet
    evo = enviroEvolveID
    
    
    
    
    genVec = c()
    nameVec = c()
    fitnessVec = c()
    sdVec_df = c()
    
    for( resp in 1:numDrugs)
    {
      
      # make drug column have drugs actual names
      
      if(resp == 1){
        name = "Home"
      }else if(resp == 2){
        name = "Non-Home"
      }
      
      
      
      genVec = c(genVec, 1:numGens)
      nameVec = c(nameVec,rep(name, times = numGens) )
      fitnessVec = c(fitnessVec, (WbarData[,resp]))
      sdVec_df = c(sdVec_df, sqrt(varData[,resp]))
      
      
      
      
    }
    
    minVec = fitnessVec - sdVec_df
    maxVec= fitnessVec + sdVec_df
    
    
    
    
    # create data frame from the matrix
    WbarDf = data.frame(genVec, nameVec, fitnessVec, minVec, maxVec)
    colnames(WbarDf) = c("time", "Env", "fitness", 'min', 'max')
    
    
    
  
    
    # for  fit trajectory plot
    if(wantPlot == TRUE)
    {
      print( ggplot(WbarDf, aes(x=as.numeric(time),y= as.numeric(fitness), colour = Env )) + geom_line()  +
               geom_ribbon(aes(x= time, ymin=min, ymax=max,fill = Env,group = Env), alpha = 0.3, linetype = 0) + 
               scale_color_manual(values=c('chocolate1','dodgerblue2'))+  geom_hline(yintercept =  0, linetype = 'dashed')+
               theme_bw()+theme(panel.border = element_rect(colour = 'black', fill = NA, size = 2), 
                                panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
                                axis.text =element_text(size = 20,family = 'Helvetica'), axis.title = element_blank())+ 
               xlab("Time")+ ylab("Fitness")+ylim(yplotmin,yplotmax))
    }
    
    ## save average slope
    for(slope in 1:numDrugs)
    {                           
      #slopeArray[evo, slope] = ((averageAllWbar[numGens/2, slope]) - (averageAllWbar[numGens, slope]))/(numGens/2 - numGens)
      slopeArray[evo, slope] = lm(averageAllWbar[100:numGens,2]~c(100:numGens))$coefficients[2]
    }
    
    
    
    ##save variance
    for(iEnv in 1:numDrugs)
    {
      varArray[evo, iEnv] = avgDrugVariDataSet[numGens, iEnv]
    }
    
    
    ##savecovariance
    for(iEnv in 1:numDrugs)
    {
      covarArray[evo, iEnv] = avgDrugCoVariDataSet[numGens, iEnv]
    }
    
    ## SAVE SLOPE OF VAR AND COVAR TRAJ'S 
    
    mVarSlopeResp = lm(avgDrugVariDataSet[100:numGens,2]~c(100:numGens))$coefficients[2]
    mCoVarSlopeResp = lm(avgDrugCoVariDataSet[100:numGens,2]~c(100:numGens))$coefficients[2]
    
    
  }
  
  
  output = list(slopeArray[1,], varArray[1,], covarArray[1,], r2Easy[1,2], D22Easy[1,2], D12Easy[1,2], mVarSlopeResp, mCoVarSlopeResp)
  return(output)
  
  
  
}


## NOTE: only fits the var traj for if have > 100 generations, else there is error, this is to correct for the 'lag time' in WF

scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}



get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y,...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}




######################################### F I G U R E     1  #################################################


# (optional) set working directory to where you'd like to save the outputs of this section, example below 
# setwd("~/Documents/JDFE Project/Figures/Figure example JDFES")

########### Simple 2D Gaussian (main text) ###############

# vars that are changed for text figures
corr = 0
respMean = 0


# vars that can be changed, but are constant for text figures 
homeMean = -0.1
homeSd = 0.1
respSd = 0.1
xplotmin = -0.35
xplotmax = 0.15
yplotmin = -0.38
yplotmax = 0.38


# Based on above parameters, make the JDFE
#-------------#
homeVar = homeSd^2
respVar = respSd^2
covar = corr*homeSd*respSd

meanVec <- c(homeMean, respMean)
sigma <- matrix(c(homeVar,covar,covar,respVar), nrow=2)
data.grid <- expand.grid(vals1 = seq(xplotmin,xplotmax, length.out=200), vals2 = seq(yplotmin, yplotmax, length.out=200))
df <- cbind(data.grid, prob = mvtnorm::dmvnorm(data.grid, mean = meanVec, sigma = sigma))


## plot the JDFE
print(ggplot(df, aes(x=vals1, y=vals2, z=prob)) + 
        geom_contour(binwidth = 3, col = 'black') +
        xlim(xplotmin, xplotmax)+ ylim(yplotmin, yplotmax) +
        theme_bw()+theme(panel.border = element_rect(colour = 'black', fill = NA, size = 2), 
                         panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
                         axis.text =element_text(size = 20,family = 'Helvetica'), axis.title = element_blank())+ 
        geom_hline(yintercept =  0, linetype = 'dashed') + 
        geom_vline(xintercept =  0, linetype = 'dashed')+
        geom_point(aes(weighted.mean(vals1, prob), weighted.mean(vals2, prob)), pch = 4))







# MAKE THE DFES FOR USE IN THE SIM
DFE1 = sort(rnorm(10000, mean = homeMean,sd = homeSd)) # constant for all JDFEs, based on no drug
DFE2 = sort(rnorm(10000, mean = respMean, sd = respSd)) # based off iteration

cor_mat = matrix(corr, ncol = 2, nrow = 2)
diag(cor_mat) = 1

# make data sets with perfect desired corr
mvdat = mvrnorm(10000, mu= c(0,0), Sigma = cor_mat, empirical = TRUE)

# compute ranks of random data 
rx <- rank(mvdat[ , 1], ties.method = "first")
ry <- rank(mvdat[ , 2], ties.method = "first")

# cor corrected DFEs 
DFE1 = DFE1[rx]
DFE2 = DFE2[ry]

sValallTheoMuts = array(0, dim= c(10000, 2))
sValallTheoMuts[,1] = DFE1
sValallTheoMuts[,2] = DFE2



# EXAMPLE OF EVOLUTION FUNCTION CALL #
# outputs graph of mean fitness trajectory in home and non-home; ribbons show +/- 1sd
# also saves results in list (in order):
# list(slopeArray, varArray, covarArray, r2Easy, D22Easy, D12Easy, mrespvarslope, mrespcovarslope)

yplotmin = -8
yplotmax = 8
output = simEvo(sValallTheoMuts, numGens=1000, numIterations = 10, popsize = 10^6, U = 10^-4, yplotmin, yplotmax, wantPlot = TRUE)




########### MORE COMPLEX 'DOUBLE GAUSSIAN' FOR SUPPlement ###############


## vals that are changed
resp1Mean = 0.5
resp2Mean = -0.17
corr = 0
x1plotmin = -1  ### need to not use overlapping y limits else there will be mutiple of the same x/y pair with different density values in the compiled data fram
x1plotmax = 1  ### so, just stick to resp1mean > resp2mean and youll be golden 
y1plotmin = 0
y1plotmax = 1
x2plotmin = -1
x2plotmax = 1
y2plotmin = -1
y2plotmax = 0
xTOTplotmin = -0.1
xTOTplotmax = 0.75
yTOTplotmin = -0.75
yTOTplotmax = 0.75

# vars that can be changed, but are constant for supp figure 
homeMean = 0.4
homeSd = 0.1
respSd = 0.1


# Make the JDFE
#-------------#
homeVar = homeSd^2
respVar = respSd^2
covar = corr*homeSd*respSd

## resp 1
meanVec <- c(homeMean, resp1Mean)
sigma <- matrix(c(homeVar,covar,covar,respVar), nrow=2)
data.grid <- expand.grid(vals1 = seq(x1plotmin,x1plotmax, length.out=200), vals2 = seq(y1plotmin, y1plotmax, length.out=200))
df1 <- cbind(data.grid, prob = mvtnorm::dmvnorm(data.grid, mean = meanVec, sigma = sigma))

## resp 2
meanVec <- c(homeMean, resp2Mean)
sigma <- matrix(c(homeVar,covar,covar,respVar), nrow=2)
data.grid <- expand.grid(vals1 = seq(x2plotmin,x2plotmax, length.out=200), vals2 = seq(y2plotmin, y2plotmax, length.out=200))
df2 <- cbind(data.grid, prob = mvtnorm::dmvnorm(data.grid, mean = meanVec, sigma = sigma))

## together
df = data.frame(c(df1$vals1, df2$vals1), c(df1$vals2, df2$vals2),c(df1$prob, df2$prob))
colnames(df) = c('vals1', 'vals2', 'prob')

## plot the JDFE
print(ggplot(df, aes(x=vals1, y=vals2, z=prob)) + 
        geom_contour(bins = 5 ,col = 'black') +
        xlim(xTOTplotmin, xTOTplotmax)+ ylim(yTOTplotmin, yTOTplotmax) +
        theme_bw()+theme(panel.border = element_rect(colour = 'black', fill = NA, size = 2), 
                         panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
                         axis.text =element_text(size = 20,family = 'Helvetica'), axis.title = element_blank())+ 
        geom_hline(yintercept =  0, linetype = 'dashed') + 
        geom_vline(xintercept =  0, linetype = 'dashed')+
        geom_point(aes(weighted.mean(vals1, prob), weighted.mean(vals2, prob)), pch = 4))




# MAKE THE DFES FOR USE IN THE SIM

## 1st env ##
DFE1_1 = sort(rnorm(10000, mean = homeMean,sd = homeSd)) # constant for all JDFEs, based on no drug
DFE2_1 = sort(rnorm(10000, mean = resp1Mean, sd = respSd)) # based off iteration

cor_mat = matrix(corr, ncol = 2, nrow = 2)
diag(cor_mat) = 1

# make data sets with perfect desired corr
mvdat = mvrnorm(10000, mu= c(0,0), Sigma = cor_mat, empirical = TRUE)

# compute ranks of random data 
rx <- rank(mvdat[ , 1], ties.method = "first")
ry <- rank(mvdat[ , 2], ties.method = "first")

# cor corrected DFEs 
DFE1_1 = DFE1_1[rx]
DFE2_1 = DFE2_1[ry]


## 2nd env ##
DFE1_2 = sort(rnorm(10000, mean = homeMean,sd = homeSd)) # constant for all JDFEs, based on no drug
DFE2_2 = sort(rnorm(10000, mean = resp2Mean, sd = respSd)) # based off iteration

cor_mat = matrix(corr, ncol = 2, nrow = 2)
diag(cor_mat) = 1

# make data sets with perfect desired corr
mvdat = mvrnorm(10000, mu= c(0,0), Sigma = cor_mat, empirical = TRUE)

# compute ranks of random data 
rx <- rank(mvdat[ , 1], ties.method = "first")
ry <- rank(mvdat[ , 2], ties.method = "first")

# cor corrected DFEs 
DFE1_2 = DFE1_2[rx]
DFE2_2 = DFE2_2[ry]


sValallTheoMuts = array(0, dim= c(20000, 2))
sValallTheoMuts[,1] = c(DFE1_1, DFE1_2)
sValallTheoMuts[,2] = c(DFE2_1, DFE2_2)






#### EXAMPLE OF EVOLUTION FUNCTION CALL ###
# outputs graph of mean fitness trajectory in home and non-home; ribbons show +/- 1sd
# also saves results in list (in order):
# list(slopeArray, varArray, covarArray, r2Easy, D22Easy, D12Easy)

yplotmin = -50
yplotmax = 50
output = simEvo(sValallTheoMuts, numGens=1000, numIterations = 5, 10^6, 10^-4, yplotmin, yplotmax, wantPlot = TRUE)





#### END FIGURE 1 ####

############################################ F I G U R E    2  ################################################################

# (optional) set working directory to where you'd like to save the outputs of this section, example below 
# setwd("~/Documents/JDFE Project/Figures/Figure Theo JDFES")

#### NON SSWM ####


# parameters to cycle through for making JDFES #
m2Vec = seq(-0.15,0.05, by = 0.05)
sdVec = seq(0.06, 0.1, by = 0.03)
corVec = seq(-0.9,0.9, by = 0.3)

# parameters that can be changed

popsize = 10^4
N = popsize
numIterations = 100
numDrugs = 2
numGens = 1000
homeMean = -0.1
homeSd = 0.1

mu2Vec_sim  = c()
c0Vec_sim = c()
sdVec_sim = c()

for(U in c(10^-4,10^-3,10^-2))
{
  
  
  print(U)
  
  # establish vectors for saving results
  r2Vals = c()
  D2Vals= c()
  D12Vals = c()
  meanRespSlopes = c()
  respVars = c()
  covars = c()
  meanVarSlopeResp = c()
  meanCoVarSlopeResp = c()
  
  
  numJDFEs = 0
  for(mu2 in m2Vec)
  {
    for(sd in sdVec)
    {
      for(c0 in corVec)
      {
        
        mu2Vec_sim  = c(mu2Vec_sim, mu2)
        c0Vec_sim = c(c0Vec_sim, c0)
        sdVec_sim = c(sdVec_sim, sd)
        
        
        numJDFEs = numJDFEs+1
        print('JDFE')
        print(numJDFEs)
        
        
        respMean = mu2
        respSd = sd
        corr = c0
        
        #-------------#
        homeVar = homeSd^2
        respVar = respSd^2
        covar = corr*homeSd*respSd
        
        
        
        # MAKE THE DFES FOR USE IN THE SIM
        DFE1 = sort(rnorm(10000, mean = homeMean,sd = homeSd)) # constant for all JDFEs, based on no drug
        DFE2 = sort(rnorm(10000, mean = respMean, sd = respSd)) # based off iteration
        
        cor_mat = matrix(corr, ncol = 2, nrow = 2)
        diag(cor_mat) = 1
        
        # make data sets with perfect desired corr
        mvdat = mvrnorm(10000, mu= c(0,0), Sigma = cor_mat, empirical = TRUE)
        
        # compute ranks of random data 
        rx <- rank(mvdat[ , 1], ties.method = "first")
        ry <- rank(mvdat[ , 2], ties.method = "first")
        
        # cor corrected DFEs 
        DFE1 = DFE1[rx]
        DFE2 = DFE2[ry]
        
        sValallTheoMuts = array(0, dim= c(10000, 2))
        sValallTheoMuts[,1] = DFE1
        sValallTheoMuts[,2] = DFE2
        
        
        output = simEvo(sValallTheoMuts, numGens, numIterations, popsize, U, yplotmin = -20, yplotmax = 20, wantPlot = FALSE)
        
        
        #list(slopeArray, varArray, covarArray, r2Easy, D22Easy, D12Easy, mVarSlope, mCovarslope)
        r2Vals = c(r2Vals, output[[4]])
        D2Vals= c(D2Vals, output[[5]])
        D12Vals = c(D12Vals, output[[6]])
        meanRespSlopes = c(meanRespSlopes, output[[1]][2])
        respVars = c(respVars, output[[2]][2])
        covars = c(covars, output[[3]][2])
        meanVarSlopeResp = c(meanVarSlopeResp, output[[7]] )
        meanCoVarSlopeResp = c(meanCoVarSlopeResp, output[[8]])
        
        
      }
    }
  }
  
  assign(paste("meanRespSlopes", U*popsize, sep = ""), meanRespSlopes)
  assign(paste("respVars", U*popsize, sep = ""), respVars)
  assign(paste("covars", U*popsize, sep = ""), covars)
  assign(paste("meanVarSlopeResp", U*popsize, sep = ""), meanVarSlopeResp)
  assign(paste("meanCoVarSlopeResp", U*popsize, sep = ""), meanCoVarSlopeResp)
  
  #Nu 1 palegreen4
  #Nu 10 darkgoldenrod
  #Nu 100 hotpink4
  
  meanRespSlopes = meanRespSlopes
  respVars = respVars
  covars = covars
  meanVarSlopeResp = meanVarSlopeResp
  meanCoVarSlopeResp = meanCoVarSlopeResp
  
  
  N = popsize
  Ub = U*(1- pnorm(0, -0.1,0.1))
  
  
  
  
  if(U*popsize == 1){
    
    color = 'palegreen4'
  }else if(U*popsize == 10)
  {
    color = 'darkgoldenrod4'
  }else if(U*popsize == 100)
  {
    color = 'hotpink4' 
  }
  
  
  # ADJUST ALL MEASURES BY 2NUb (as discussed in text)
  adjR2Vals = r2Vals*2*N*Ub
  adjD22Vals = D2Vals*2*N*Ub 
  adjD12Vals = D12Vals*2*N*Ub
  adjD22endTimeVals= D2Vals*2*N*Ub*numGens
  adjD12endTimeVals = D12Vals*2*N*Ub*numGens
  
  
  data_thisNu = data.frame(adjR2Vals,adjD22Vals,adjD12Vals,adjD22endTimeVals,adjD12endTimeVals,meanRespSlopes, respVars, covars, meanVarSlopeResp, meanCoVarSlopeResp)
  
  
  scientific_10 <- function(x) {
    parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
  }
  
  ################### make and save plots
  ## MEAN SLOPE
  assign(paste('r2vSlope_Nu', N*U,sep = ""), ggplot(data_thisNu, aes(x=adjR2Vals,y=meanRespSlopes)) + geom_point(col = color, size = 3) +
           geom_smooth(method='lm', formula= y~x, col = color, lwd = 0.5) + geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey')+
           theme_bw()+theme(panel.border = element_rect(colour = 'black', fill = NA, size = 1), 
                            panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
                            axis.text.x =element_text(size = 35,family = 'Helvetica', color = c('black', 'transparent', 'black', 'transparent', 'transparent')),axis.text.y =element_text(size = 35,family = 'Helvetica', color = c('black')), axis.title = element_blank()) +
           scale_y_continuous(label=scientific_10) + scale_x_continuous(label=scientific_10)  )
  
  
  
  ## SLOPE OF VARIANCE TRAJECTORY
  assign(paste('D22vVarSlope_Nu', N*U,sep = ""), ggplot(data_thisNu, aes(x=adjD22Vals,y=meanVarSlopeResp)) + geom_point(col = color, size = 3) +
           geom_smooth(method='lm', formula= y~x, col = color, lwd = 0.5) + geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey')+
           theme_bw()+theme(panel.border = element_rect(colour = 'black', fill = NA, size = 1), 
                            panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
                            axis.text.x =element_text(size = 35,family = 'Helvetica', color = c('black', 'transparent', 'black', 'transparent', 'transparent')),axis.text.y =element_text(size = 35,family = 'Helvetica', color = c('black')), axis.title = element_blank()) +
           scale_y_continuous(label=scientific_10) + scale_x_continuous(label=scientific_10)  )
  
  
  
  ## SLOPE OF COVARIANCE TRAJECTORY
  assign(paste('D12vCoVarSlope_Nu', N*U,sep = ""),  ggplot(data_thisNu, aes(x=adjD12Vals,y=meanCoVarSlopeResp)) + geom_point(col = color, size = 3) +
           geom_smooth(method='lm', formula= y~x, col = color, lwd = 0.5) + geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey')+
           theme_bw()+theme(panel.border = element_rect(colour = 'black', fill = NA, size = 1), 
                            panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
                            axis.text.x =element_text(size = 35,family = 'Helvetica', color = c('black', 'transparent', 'black', 'transparent', 'transparent')),axis.text.y =element_text(size = 35,family = 'Helvetica', color = c('black')), axis.title = element_blank()) +
           scale_y_continuous(label=scientific_10) + scale_x_continuous(label=scientific_10)  )
  
  
  
  
}


### SAVE ALL DATA
allData_nonSSWM_TheoJDFEs = data.frame(c(rep(1, times = 70), rep(10, times = 70), rep(100, times = 70)), c(meanRespSlopes1, meanRespSlopes10, meanRespSlopes100), 
                                        c(respVars1, respVars10,respVars100),c(covars1,covars10,covars100), c(meanVarSlopeResp1,meanVarSlopeResp10,meanVarSlopeResp100), c(meanCoVarSlopeResp1, meanCoVarSlopeResp10, meanCoVarSlopeResp100))


colnames(allData_nonSSWM_TheoJDFEs) = c('Nu', 'meanRespSlope','respVar', 'covar', 'meanrespVarSlope', 'meancovarSlope' )

write.csv(allData_nonSSWM_TheoJDFEs, "allData_nonSSWM_TheoJDFEs.csv")











#### SSWM THEO JDFE ####
m2Vec = seq(-0.15,0.05, by = 0.05)
sdVec = seq(0.06, 0.1, by = 0.03)
corVec = seq(-0.9,0.9, by = 0.3)
NuVec = c(10^-7, 10^-6,10^-5,10^-4)
numIterations = 500 
numDrugs = 2
numGens = 1500

meanHomeSlopes = c()
meanRespSlopes = c()
homeVars = c()
respVars = c()
covars = c()
cors = c()

r1Vals = c()
r2Vals = c()
varOfr2Dist = c()
D1Vals = c()
D2Vals= c()
D12Vals = c()

allM2 = c()
allsd = c()
allCor = c()
allNu = c()
allUbs = c()
allFixRate = c()
meanVarSlopeResp = c()
meanCoVarSlopeResp = c()



allEndFits = c()

popsize = 10^6
N = popsize
U = 10^-7
Nu = N*U
sizeOfBin = 0.02
Ub = U*(1- pnorm(0, -0.1,0.1)) # prob of beneficial mut (as definined in theory)

count = 0
for(mu2 in m2Vec)
{
  for(sd in sdVec)
  {
    for(c0 in corVec)
    {
      count = count+1
      print('JDFE num')
      print(count)
      
      
      # save vars for this it
      allM2 = c(allM2, mu2)
      allsd = c(allsd, sd)
      allCor = c(allCor, c0)
      allNu = c(allNu, Nu*popsize)
      
      
      ####  generate the sorted two DFES --> must be sorted to apply ranks to ####
      
      DFE1 = sort(rnorm(10000, mean = -0.1,sd = 0.1)) # constant for all JDFEs, based on no drug
      DFE2 = sort(rnorm(10000, mean = mu2, sd = sd)) # based off iteration
      
      corr = c0
      cor_mat = matrix(corr, ncol = 2, nrow = 2)
      diag(cor_mat) = 1
      
      # make data sets with perfect desired corr
      mvdat = mvrnorm(10000, mu= c(0,0), Sigma = cor_mat, empirical = TRUE)
      
      # compute ranks of random data 
      rx <- rank(mvdat[ , 1], ties.method = "first")
      ry <- rank(mvdat[ , 2], ties.method = "first")
      
      # cor corrected DFEs --> of course, this isnt 100% perfect but corr always w/in +/- 0.001 of desired
      DFE1 = DFE1[rx]
      DFE2 = DFE2[ry]
      
      
      ### Find R and D-values ####
      
      
      r1 = mean((DFE1>0)*DFE1^2)
      r2 = mean((DFE1>0)*DFE1*DFE2)
      D1 = mean((DFE1>0)*DFE1^3)
      D2 = mean((DFE1>0)*DFE1*DFE2^2)
      D12 = mean((DFE1>0)*DFE1^2*DFE2)
      
      
      r1Vals = c(r1Vals, r1)
      r2Vals = c(r2Vals, r2)
      D1Vals = c(D1Vals, D1)
      D2Vals = c(D2Vals, D2)
      D12Vals = c(D12Vals, D12)
      
      
      
      WbarAllIterationsHome = array(0, dim= c(numGens, numIterations))
      WbarAllIterationsResp = array(0, dim = c(numGens, numIterations))
      allWbarsAllIterations = array(data = 0, dim = c(numGens, (numDrugs), numIterations))
      
      allNumben = c()
      allNumFix = c()
      allIntraItCorrs = c()
      
      homeFitMean = mean(DFE1)
      homeFitSD = sd(DFE1)
      for(it in 1:numIterations)
      {
        print(it)
        wBarHome = 1
        wBarResp = 1
        numBen = 0
        numFix = 0
        
        ### make vec of times in which get mutations
        end = 0
        current = 0
        mutTimesVec = c()
        Ub = U*(1- pnorm(0, -0.1,0.1)) 
        while(end < numGens)
        {
          current = end
          nextMutTime = ceiling(rexp(1, rate = popsize*Ub))
          end = current + nextMutTime
          
          mutTimesVec = c(mutTimesVec, end)
          
        }
        
        
        for(gen in 1:numGens)
        {
          if(sum(gen==mutTimesVec)!=0)
          {
            
            
            # home S-val
            homeS = rnorm(1, mean = homeFitMean, sd = homeFitSD)
            
            while(homeS > max(DFE1) | homeS < min(DFE1))
            {
              homeS = rnorm(1, mean = homeFitMean, sd = homeFitSD)
            }
            
            if(homeS > 0)
            {
              ## 'flip coin' to see if fixes##
              prob = 2*homeS
              if(prob > 1)
              {
                prob = 1
              }
              flip = sample(c(0,1), 1, prob = c(1-prob,prob))
              numBen = numBen + 1
              
              if(flip == 1) ## if mut does fix
              {
                
                
                # find fit of S-vals in bin in env 2
                minofBin = homeS - sizeOfBin
                maxofBin = homeS + sizeOfBin
                sValsinBin = DFE2[DFE1>=minofBin & DFE1<=maxofBin]
                
                # if only have 1 sval in bin, respFit = the direct translation from that one fit to env 2
                if(length(sValsinBin)==0)
                {
                  # just propogate last Wbar in both envs to this timepoint
                  WbarAllIterationsHome[gen,it] = wBarHome
                  WbarAllIterationsResp[gen,it] = wBarResp
                  
                  allWbarsAllIterations[gen,1, it] = wBarHome
                  allWbarsAllIterations[gen,2, it] = wBarResp
                }else if(length(sValsinBin) == 1)
                {
                  numFix = numFix + 1
                  location = which.min(abs(DFE1 - sValsinBin))
                  respS = DFE2[location]
                  wBarHome = (homeS)+wBarHome
                  wBarResp = (respS)+wBarResp
                  
                  WbarAllIterationsHome[gen, it] = wBarHome
                  WbarAllIterationsResp[gen, it] = wBarResp
                  
                  allWbarsAllIterations[gen,1, it] = wBarHome
                  allWbarsAllIterations[gen,2, it] = wBarResp
                }else{
                  
                  numFix = numFix + 1
                  
                  
                  # fit dist to sValsinBin and sample resp Sval
                  respDistFit = fitdistr(sValsinBin, 'normal')
                  respS = rnorm(1, respDistFit$estimate[1], respDistFit$estimate[2])
                  wBarHome = (homeS)+wBarHome
                  wBarResp = (respS)+wBarResp
                  
                  WbarAllIterationsHome[gen, it] = wBarHome
                  WbarAllIterationsResp[gen, it] = wBarResp
                  
                  allWbarsAllIterations[gen,1, it] = wBarHome
                  allWbarsAllIterations[gen,2, it] = wBarResp
                }
                
                
              }else{
                WbarAllIterationsHome[gen,it] = wBarHome
                WbarAllIterationsResp[gen,it] = wBarResp
                
                allWbarsAllIterations[gen,1, it] = wBarHome
                allWbarsAllIterations[gen,2, it] = wBarResp
              }
            }else{
              
              # homeS < 0 
              
              # just propogate last Wbar in both envs to this timepoint
              WbarAllIterationsHome[gen,it] = wBarHome
              WbarAllIterationsResp[gen,it] = wBarResp
              
              allWbarsAllIterations[gen,1, it] = wBarHome
              allWbarsAllIterations[gen,2, it] = wBarResp
            }
            
          }else{
            # dont get mut
            
            # just propogate last Wbar in both envs to this timepoint
            WbarAllIterationsHome[gen,it] = wBarHome
            WbarAllIterationsResp[gen,it] = wBarResp
            
            allWbarsAllIterations[gen,1, it] = wBarHome
            allWbarsAllIterations[gen,2, it] = wBarResp
          }
        }# end gen loop
        allNumben = c(allNumben, numBen)
        allNumFix = c(allNumFix, numFix)
      }#end iteration loop 
      
      averageAllWbar = array(data = 0, dim = c(numGens, (numDrugs)))
      for( gen2 in 1:numGens)
      {
        for(d2 in 1:(numDrugs))
        {
          averageAllWbar[gen2,d2] = mean(allWbarsAllIterations[gen2,d2, ])
        }
      }
      
      
      DrugVariDataSet = array(data = 0, dim = c(numGens, numDrugs))
      for( gen2 in 1:numGens)
      {
        for(d2 in 1:(numDrugs))
        {
          DrugVariDataSet[gen2,d2] = var((allWbarsAllIterations[gen2,d2, ]))
        }
      }
      
      DrugCoVariDataSet = array(data = 0, dim = c(numGens, numDrugs))
      for( gen2 in 1:numGens)
      {
        for(d2 in 1:(numDrugs))
        {
          DrugCoVariDataSet[gen2,d2] = cov((allWbarsAllIterations[gen2,d2, ]),(allWbarsAllIterations[gen2,1, ]) )
        }
      }
      
      
      
      
      mHomeSlope = ((averageAllWbar[numGens, 1]) - (averageAllWbar[numGens/2, 1]))/(numGens-numGens/2)
      mRespSlope = ((averageAllWbar[numGens, 2]) - (averageAllWbar[numGens/2, 2]))/(numGens-numGens/2)
      
      ## SAVE SLOPE OF VAR AND COVAR TRAJ'S 
      mVarSlopeResp = ((DrugVariDataSet[numGens, 2]) - (DrugVariDataSet[numGens/2, 2]))/(numGens-numGens/2)
      mCoVarSlopeResp = ((DrugCoVariDataSet[numGens, 2]) - (DrugCoVariDataSet[numGens/2, 2]))/(numGens-numGens/2)
      
      
      
      # save results
      meanHomeSlopes = c(meanHomeSlopes, mHomeSlope)
      meanRespSlopes = c(meanRespSlopes, mRespSlope)
      homeVars = c(homeVars, DrugVariDataSet[numGens, 1])
      respVars = c(respVars, DrugVariDataSet[numGens, 2])
      covars = c(covars, DrugCoVariDataSet[numGens,2])
      
      meanVarSlopeResp = c(meanVarSlopeResp, mVarSlopeResp)
      meanCoVarSlopeResp = c(meanCoVarSlopeResp, mCoVarSlopeResp)
      
      
      fixRate = mean(allNumFix)
      Ub = mean(allNumben)/numGens
      allUbs = c(allUbs, Ub) 
      allFixRate = c(allFixRate, fixRate)
      
      endFit = mean((averageAllWbar[numGens,2]))
      allEndFits = c(allEndFits, endFit)
      
      
      
      
    }
  }
}

N = popsize
Ub = U*(1- pnorm(0, -0.1,0.1)) # prob of beneficial mut (as definined in theory)


assign(paste("meanRespSlopes", U*popsize, sep = ""), meanRespSlopes)
assign(paste("respVars", U*popsize, sep = ""), respVars)
assign(paste("covars", U*popsize, sep = ""), covars)
assign(paste("meanVarSlopeResp", U*popsize, sep = ""), meanVarSlopeResp)
assign(paste("meanCoVarSlopeResp", U*popsize, sep = ""), meanCoVarSlopeResp)



N = popsize
Ub = U*(1- pnorm(0, -0.1,0.1))



color = 'dodgerblue4'

adjR2Vals = r2Vals*2*N*Ub
adjD22Vals = D2Vals*2*N*Ub 
adjD12Vals = D12Vals*2*N*Ub
adjD22endTimeVals= D2Vals*2*N*Ub*numGens
adjD12endTimeVals = D12Vals*2*N*Ub*numGens


data_thisNu = data.frame(adjR2Vals,adjD22Vals,adjD12Vals,adjD22endTimeVals,adjD12endTimeVals,meanRespSlopes, respVars, covars, meanVarSlopeResp, meanCoVarSlopeResp)


scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}


# make and save plots
## MEAN SLOPE

assign(paste('r2vSlope_Nu', N*U,sep = ""), ggplot(data_thisNu, aes(x=adjR2Vals,y=meanRespSlopes)) + geom_point(col = color, size = 4) +
         geom_smooth(method='lm', formula= y~x, col = color, lwd = 0.5) + geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey')+
         theme_bw()+theme(panel.border = element_rect(colour = 'black', fill = NA, size = 1), 
                          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
                          axis.text.x =element_text(size = 35,family = 'Helvetica', color = c('black', 'transparent', 'black', 'transparent', 'transparent')),axis.text.y =element_text(size = 35,family = 'Helvetica', color = c('black')), axis.title = element_blank()) +
         scale_y_continuous(label=scientific_10) + scale_x_continuous(label=scientific_10)  )



## SLOPE OF VARIANCE TRAJECTORY


assign(paste('D22vVarSlope_Nu', N*U,sep = ""), ggplot(data_thisNu, aes(x=adjD22Vals,y=meanVarSlopeResp)) + geom_point(col = color, size = 4) +
         geom_smooth(method='lm', formula= y~x, col = color, lwd = 0.5) + geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey')+
         theme_bw()+theme(panel.border = element_rect(colour = 'black', fill = NA, size = 1), 
                          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
                          axis.text.x =element_text(size = 35,family = 'Helvetica', color = c('black', 'transparent', 'black', 'transparent', 'transparent')),axis.text.y =element_text(size = 35,family = 'Helvetica', color = c('black')), axis.title = element_blank()) +
         scale_y_continuous(label=scientific_10) + scale_x_continuous(label=scientific_10)  )







## SLOPE OF COVARIANCE TRAJECTORY
assign(paste('D12vCoVarSlope_Nu', N*U,sep = ""),  ggplot(data_thisNu, aes(x=adjD12Vals,y=meanCoVarSlopeResp)) + geom_point(col = color, size = 4) +
         geom_smooth(method='lm', formula= y~x, col = color, lwd = 0.5) + geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey')+
         theme_bw()+theme(panel.border = element_rect(colour = 'black', fill = NA, size = 1), 
                          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
                          axis.text.x =element_text(size = 35,family = 'Helvetica', color = c('black', 'transparent', 'black', 'transparent', 'transparent')),axis.text.y =element_text(size = 35,family = 'Helvetica', color = c('black')), axis.title = element_blank()) +
         scale_y_continuous(label=scientific_10) + scale_x_continuous(label=scientific_10)  )





#### SAVE AND PRINT ALL PLOTS FROM ALL NUs ####

## save all individual plots

#SSWM
ggsave(r2vSlope_Nu0.1, file = 'r2vSlope_Nu0.1.pdf', width = 8.99, height = 7.25 )
ggsave(D22vVarSlope_Nu0.1, file = 'D22vVarSlope_Nu0.1.pdf' , width = 8.99, height = 7.25 )
ggsave(D12vCoVarSlope_Nu0.1, file = 'D12vCoVarSlope_Nu0.1.pdf', width = 8.99, height = 7.25 )


#Nu 1
ggsave(r2vSlope_Nu1, file = 'r2vSlope_Nu1.pdf' , width = 8.99, height = 7.25 )
ggsave(D22vVarSlope_Nu1, file = 'D22vVarSlope_Nu1.pdf' , width = 8.99, height = 7.25 )
ggsave(D12vCoVarSlope_Nu1, file = 'D12vCoVarSlope_Nu1.pdf' ,width = 8.99, height = 7.25 )


#Nu 10
ggsave(r2vSlope_Nu10, file = 'r2vSlope_Nu10.pdf', width = 8.99, height = 7.25 )
ggsave(D22vVarSlope_Nu10, file = 'D22vVarSlope_Nu10.pdf' ,width = 8.99, height = 7.25 )
ggsave(D12vCoVarSlope_Nu10, file = 'D12vCoVarSlope_Nu10.pdf' , width = 8.99, height = 7.25 )

#Nu 100
ggsave(r2vSlope_Nu100, file = 'r2vSlope_Nu100.pdf' , width = 8.99, height = 7.25 )
ggsave(D22vVarSlope_Nu100, file = 'D22vVarSlope_Nu100.pdf' , width = 8.99, height = 7.25 )
ggsave(D12vCoVarSlope_Nu100, file = 'D12vCoVarSlope_Nu100.pdf' , width = 8.99, height = 7.25 )


## print all plots

print(r2vSlope_Nu0.1)
print(D22vVarSlope_Nu0.1)
print(D12vCoVarSlope_Nu0.1)

print(r2vSlope_Nu1)
print(D22vVarSlope_Nu1)
print(D12vCoVarSlope_Nu1)


print(r2vSlope_Nu10)
print(D22vVarSlope_Nu10)
print(D12vCoVarSlope_Nu10)


print(r2vSlope_Nu100)
print(D22vVarSlope_Nu100)
print(D12vCoVarSlope_Nu100)





##############################################  F I G U R E    3   #################################################


# (optional) set working directory to where you'd like to save the outputs of this section, example below 
# setwd("~/Documents/JDFE Project/Figures/Figure Epistasis")



########### SSWM WITH EPISTASIS SIM #########

# variables that can be changed 
numIterations = 100
numDrugs = 2
numGens = 1000

r1Vals = c()
r2Vals = c()

allUbs = c()

### MOST IMPORTANT PARAMETERS
U = 10^-7
gamma1 = 0.05
#####

N = 1.6*10^5
popsize = 1.6*10^5
xmax = 2
ymax = 2
gamma2 = 1
sd = 0.01
mu2 = -0.1
x0 = 1
y0 = 1
rhoVec = seq(-0.9,0.9, by = 0.2)
endGens = c()
sizeOfBin = sd/2

allF1theo = c()
allF2theo = c()
allavgWbar1 = c()
allavgWbar2 = c()


simfit1_tEND = c()
simfit2_tEND = c()

theofit1_tEND = c()
theofit2_tEND = c()


r2means = c()

r1vec = c()

numBenallRhos = c()


totalIT= 0


endFit = c()
meanendfit1 = c()
meanendFit2 = c()
meanr2vec = c()

for(rhoDesired in rhoVec)
{
  
  WbarAllIterationsHome = array(0, dim= c(numGens, numIterations))
  WbarAllIterationsResp = array(0, dim = c(numGens, numIterations))
  allWbarsAllIterations = array(data = 0, dim = c(numGens, (numDrugs), numIterations))
  
  allNumben = c()
  
  print('rho')
  print(rhoDesired)
  allrhos = c()
  numBensallIts = c()
  totalIT = totalIT + 1
  r2vec = c()
  for(it in 1:numIterations)
  {
    
    #print(it)
    wBarHome = x0
    wBarResp = y0
    numBen = 0
    
    gen = 0
    
    numBenthisIt= 0 
    
    ##make vec of times in which get a mut by sampling from waiting time distribution
    end = 0
    current = 0
    mutTimesVec = c()
    while(end < numGens)
    {
      current = end
      nextMutTime = ceiling(rexp(1, rate = popsize*U))
      end = current + nextMutTime
      
      mutTimesVec = c(mutTimesVec, end)
      
    }
    
    # remove last entry bc is past end num gens doing in sim
    
    
    for(try in 1:numGens)
    {
      
      gen = gen+ 1
      
      ### only do process of mutation if is time in which get a mutation
      if(sum(try==mutTimesVec)!=0)
      {
        mu1 = gamma1*(xmax - wBarHome)
        
        
        ####  generate the sorted two DFES --> must be sorted to apply ranks to ####
        #### THESE CHANGE WITH FITNESS (hence)
        
        DFE1 = sort(rexp(10000, rate = 1/abs(mu1))) # exponential dist, changes with fitness
        DFE2 = sort(rnorm(10000, mean = mu2, sd = sd)) # stays SAME for all its
        
        corr = rhoDesired
        cor_mat = matrix(corr, ncol = 2, nrow = 2)
        diag(cor_mat) = 1
        
        # make data sets with perfect desired corr
        mvdat = mvrnorm(10000, mu= c(0,0), Sigma = cor_mat, empirical = TRUE)
        
        # compute ranks of random data 
        rx <- rank(mvdat[ , 1], ties.method = "first")
        ry <- rank(mvdat[ , 2], ties.method = "first")
        
        # cor corrected DFEs --> of course, this isnt 100% perfect but corr always w/in +/- 0.001 of desired
        DFE1 = DFE1[rx]
        DFE2 = DFE2[ry]
        
        rho = cor(DFE1, DFE2)
        allrhos = c(allrhos, rho)
        
        
        ### save INITIal r values
        if(gen == 1)
        {
          r1vec = c(r1vec, mean(DFE1^2*(DFE1>0)))
          r2vec = c(r2vec, mean(DFE2*DFE1*(DFE1>0)) )
          
        }
        
        
        
        # home S-val
        
        if(mu1<0)
        {
          homeS = -1*rexp(1, rate = 1/abs(mu1))
        }else{
          homeS = rexp(1, rate = 1/mu1)
        }
        
        
        if(homeS > 0)
        {
          numBen = numBen+1
          ## 'flip coin' to see if fixes##
          
          prob = 2*homeS
          if(prob > 1){flip = 1}else{
            flip = sample(c(0,1), 1, prob = c(1-prob,prob))
          }
          
          
          if(flip == 1) ## if mut does fix
          {
            numBenthisIt = numBenthisIt + 1
            ## find respS
            
            muThishomeS = mu2 + ((rho*sd)/mu1)*(homeS - mu1)
            respS = rnorm(1, muThishomeS, sd)
            
            
            wBarHome = (homeS)+wBarHome
            wBarResp = ( respS)+wBarResp
            
            WbarAllIterationsHome[gen, it] = wBarHome
            WbarAllIterationsResp[gen, it] = wBarResp
            
            allWbarsAllIterations[gen,1, it] = wBarHome
            allWbarsAllIterations[gen,2, it] = wBarResp
          }else{
            # just propogate last Wbar in both envs to this timepoint
            WbarAllIterationsHome[gen,it] = wBarHome
            WbarAllIterationsResp[gen,it] = wBarResp
            
            allWbarsAllIterations[gen,1, it] = wBarHome
            allWbarsAllIterations[gen,2, it] = wBarResp
          }
        }else if(homeS < 0 & wBarHome > xmax){ ## to let fit go back down if goes over
          
          prob = 2*abs(homeS)
          if(prob > 1){flip = 1}else{
            flip = sample(c(0,1), 1, prob = c(1-prob,prob))
          }
          
          
          if(flip == 1) ## if mut does fix
          {
            
            ## find respS   ### asumming neg home with result in same mag for selection of resp s, but will just have opposite value as if homeS was same vsalue but pos
            
            muThishomeS = mu2 + ((rho*sd)/mu1)*(homeS - mu1)
            respS = rnorm(1, muThishomeS, sd)
            
            wBarHome = (homeS)+wBarHome
            wBarResp = ( respS)+wBarResp
            
            WbarAllIterationsHome[gen, it] = wBarHome
            WbarAllIterationsResp[gen, it] = wBarResp
            
            allWbarsAllIterations[gen,1, it] = wBarHome
            allWbarsAllIterations[gen,2, it] = wBarResp
          }else{
            # just propogate last Wbar in both envs to this timepoint
            WbarAllIterationsHome[gen,it] = wBarHome
            WbarAllIterationsResp[gen,it] = wBarResp
            
            allWbarsAllIterations[gen,1, it] = wBarHome
            allWbarsAllIterations[gen,2, it] = wBarResp
          }
          
          
        }else{
          # just propogate last Wbar in both envs to this timepoint
          WbarAllIterationsHome[gen,it] = wBarHome
          WbarAllIterationsResp[gen,it] = wBarResp
          
          allWbarsAllIterations[gen,1, it] = wBarHome
          allWbarsAllIterations[gen,2, it] = wBarResp
        }
        
        
      }else{
        # just propogate last Wbar in both envs to this timepoint
        WbarAllIterationsHome[gen,it] = wBarHome
        WbarAllIterationsResp[gen,it] = wBarResp
        
        allWbarsAllIterations[gen,1, it] = wBarHome
        allWbarsAllIterations[gen,2, it] = wBarResp
      }
      
    }# end generation loop
    allNumben = c(allNumben, numBen)
    endGens = c(endGens, gen)
    
    
    numBensallIts = c(numBensallIts,numBenthisIt)
    
  }#end iteration loop 
  
  numBenallRhos = c(numBenallRhos, mean(numBensallIts))
  
  UbthisRho = mean(numBensallIts)/numGens
  
  
  
  averageAllWbar = array(data = 0, dim = c(numGens, (numDrugs)))
  for( gen2 in 1:numGens)
  {
    for(d2 in 1:(numDrugs))
    {
      
      
      averageAllWbar[gen2,d2] = mean((allWbarsAllIterations[gen2,d2, ]))
    }
  }
  meanendfit1 = c(meanendfit1, averageAllWbar[numGens,1])
  meanendFit2 = c(meanendFit2, averageAllWbar[numGens,2])
  meanr2vec = c(meanr2vec, mean(r2vec))
  
  DrugSDiDataSet = array(data = 0, dim = c(numGens, numDrugs))
  for( gen2 in 1:numGens)
  {
    for(d2 in 1:(numDrugs))
    {
      DrugSDiDataSet[gen2,d2] = sd(((allWbarsAllIterations[gen2,d2, ])))
    }
  }
  
  
  endGens =numGens
  
  f1theo = c((x0))
  f2theo = c((y0))
  min1Vec = c(x0)
  max1Vec = c(x0)
  
  min2Vec = c(y0)
  max2Vec = c(y0)
  
  rho = mean(allrhos)
  
  N = 1
  
  tadj = 2*popsize*U
  for(t in 1:(min(endGens) -1 ))
  {
    f1theo = c(f1theo, xmax - 1/(2*gamma1^2*t*tadj + 1/(xmax-(x0))))
    min1Vec = c(min1Vec,averageAllWbar[t,1]-DrugSDiDataSet[t+1,1] )  # max and min are +/- 1 sd
    max1Vec = c(max1Vec,averageAllWbar[t,1]+DrugSDiDataSet[t+1,1])
    
    
    f2theo = c(f2theo, (y0) + ((mu2 + rho*sd)/(2*gamma1))*(log(1 + 2*gamma1^2*(xmax - (x0))*t*tadj)))
    min2Vec = c(min2Vec,averageAllWbar[t,2]-DrugSDiDataSet[t+1,2]  )
    max2Vec = c(max2Vec,averageAllWbar[t,2]+DrugSDiDataSet[t+1,2]  )
    
  }
  
  
  
  
  
  
  plotDf = data.frame(c(rep('Env1', times = min(endGens)),rep('Env2', times = min(endGens)),rep('Env1', times = min(endGens))
                        ,rep('Env2', times = min(endGens))),
                      c(rep('Sim', times = min(endGens*2)),rep('Theo', times = min(endGens*2))),
                      c((c(x0,averageAllWbar[1:min(endGens-1), 1])), (c(y0,averageAllWbar[1:min(endGens-1), 2])), f1theo, f2theo),
                      rep(1:min(endGens), times= 4), c(min1Vec, min2Vec, min1Vec, min2Vec), c(max1Vec, max2Vec,max1Vec,max2Vec))
  colnames(plotDf) = c('Env', 'Type', "Fitness", 'Time','min', 'max')
  
  plt = ggplot(plotDf, aes(x = Time, y = Fitness)) + geom_line(aes(color = Env, lty = Type))+  ylim(0.5,1.5)
  
  plt = plt + theme_bw()+theme(axis.title = element_text(size = 20), axis.text = element_text(size = 25),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill = NA, color = 'black', size = 1.5), axis.line = element_line(colour = "black"), legend.position = 'none')+
    geom_ribbon(aes(x = Time, y = Fitness,ymin = min, ymax = max, fill = Env,group = Env), alpha = 0.3, linetype = 0)
  print(plt)
  
  name = paste0('SUPP_FitTrajPlot2_tadj2_Nu', popsize*U,'_SSWM_rho_',rhoDesired, '.pdf')
  ggsave(plt,file= name, width = 6.18, height = 5.5)
  
  ### save plotDF
  if(rho > 0)
  {
    assign(paste("plotDf_tadj2_rho", rhoDesired, '_Nu_', 'SSWM', sep = ""), plotDf) 
  }else{
    assign(paste("plotDf_tadj2_rho_neg_",abs(rhoDesired),'_Nu_' ,'SSWM', sep = ""), plotDf) 
  }
  
  allF1theo = c(allF1theo, f1theo)
  allF2theo = c(allF2theo, f2theo)
  allavgWbar1 = c(allavgWbar1, (averageAllWbar[1:min(endGens), 1]) )
  allavgWbar2 = c(allavgWbar2,  (averageAllWbar[1:min(endGens), 2]))
  
  
  simfit1_tEND = c(simfit1_tEND, averageAllWbar[numGens, 1])
  simfit2_tEND = c(simfit2_tEND, averageAllWbar[numGens, 2])
  
  theofit1_tEND = c(theofit1_tEND, f1theo[numGens])
  theofit2_tEND = c(theofit2_tEND, f2theo[numGens])
  
  
  r2means = c(r2means, mean(r2vec))
  lastHomeSArray[totalIT, 1] = rho
  lastHomeSArray[totalIT, 2] = gamma1
  lastHomeSArray[totalIT, 3] = homeS
  
}




U = U
meanr2vec = meanr2vec
meanendFit2 = meanendFit2

df = data_frame(meanr2vec, meanendFit2-y0)
colnames(df) = c('r2', 'MeandeltaY')
newPlt = ggplot(df, aes(x = r2, y = MeandeltaY)) + 
  geom_point() +
  geom_smooth(method='lm',formula = y~x, se = T, color = 'black', lwd = 0.8)+ 
  scale_x_continuous(label=scientific_10,guide = guide_axis(check.overlap = TRUE)) +
  theme(axis.title = element_text(size = 15), axis.text = element_text(size = 35),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black', size = 1.5), 
        axis.line = element_line(colour = "black"))

print(newPlt)  

name = paste0('r2vsDeltay_Nu_', popsize*U, '.pdf')
ggsave(newPlt, file = name, width = 6.18, height = 5.5) 

name = paste0('r2vsDeltay_Nu_', popsize*U, '.csv')
write.csv(df, name)






#####  EPISTASIS SIM, NOT SSWM  #####
allMeanR2s = c()
allMeanEndFits2 = c()
allrhos = c()
allUs = c()
for(U in c( 10^-6, 10^-5, 10^-4))
{
  if(U == 10^-6)
  {
    gamma1 = 0.05
  }else if(U == 10^-5){
    gamma1 = 0.025
    
  }else if(U == 10^-4){
    gamma1 = 0.0125
  }
  
  numIterations = 1
  numDrugs = 2
  numGens = 1000
  
  r1Vals = c()
  r2Vals = c()
  
  allUbs = c()
  
  popsize = 1.6*10^5
  
  
  N = popsize
  
  sizeOfBin = 0.01
  
  xmax = 3
  ymax = 3
  gamma2 = 1
  sd = 0.01
  mu2 = -0.1
  x0 = 1
  y0 = 1
  rhoVec = seq(-0.9,0.9, by = 0.2)
  endGens = c()
  
  
  allF1theo = c()
  allF2theo = c()
  allavgWbar1 = c()
  allavgWbar2 = c()
  
  
  simfit1_t200 = c()
  simfit2_t200 = c()
  meanendfit1 = c()
  meanendFit2 = c()
  meanr2vec = c()
  
  
  wBarHome = x0
  wBarResp = y0
  numBen = 0
  
  gen = 0
  WbarAllIterationsHome = array(0, dim= c(numGens, numIterations))
  WbarAllIterationsResp = array(0, dim = c(numGens, numIterations))
  allWbarsAllIterations = array(data = 0, dim = c(numGens, (numDrugs), numIterations))
  
  allNumben = c()
  
  meanrhovec = c()
  
  
  for(rhoDesired in rhoVec)
  {
    print('rhoDesired')
    print(rhoDesired)
    r1vec = c()
    r2vec = c()
    allrhos = c()
    
    
    
    
    for(it in 1:numIterations)
    {
      enviroEvolveID = 1
      print(it)
      
      
      ######################################################
      ## pre-run from start fit 0 til mean fit 1
      enviroEvolveID = 1
      #print(it)
      
      startIDarray = array(data = c(1, popsize, 0.5), dim = c(1, 3))
      mutFitnesses = array(data = 0.5, dim = c(1, numDrugs))
      allWbars = array(data = 0.5, dim = c(1000000, (numDrugs)))
      numBen = 0 
      
      x = 0
      Wbar = 0
      while(Wbar < 1)
      {
        x = x+1
        if(x == 1)
        {
          r1vec = c(r1vec, mean(DFE1^2*(DFE1>0)))
          r2vec = c(r2vec, mean(DFE2*DFE1*(DFE1>0)) )
          
        }
        
        
        
        
        
        ### Simulate evo on DFE1, trans to traj on DFE2
        numMuts = dim(startIDarray)[1]
        
        #### Selection and Drift ####
        
        Wbar = sum(startIDarray[ , 3] * (startIDarray[ ,2] / popsize))
        
        # generate probability vector
        probvec = numeric(numMuts)
        for( i in 1:numMuts)
        {
          #probvec[i] = (startIDarray[i,3] * (startIDarray[i,2] / popsize)) / Wbar
          
          if((startIDarray[i,2] / popsize) + (startIDarray[i,2] / popsize)* (startIDarray[i,3] - Wbar) >0)
          {
            ## SK way new below
            probvec[i] = (startIDarray[i,2] / popsize) + (startIDarray[i,2] / popsize)* (startIDarray[i,3] - Wbar)
          }else{
            probvec[i] = 0
          }
        }
        
        # generate IDarray (with only sel/drift) for next gen
        nextGenIDarray = array(data = 0, dim = c(dim(startIDarray)[1], dim(startIDarray)[2]))
        nextGenOffspring = rmultinom(1, size = popsize, prob = probvec)
        for( i in 1:numMuts)
        {
          nextGenIDarray[i, 1] = startIDarray[i, 1]
          nextGenIDarray[i, 2] = nextGenOffspring[i]
          nextGenIDarray[i, 3] = startIDarray[i, 3]
        }
        
        
        
        #### Mutations ####
        gen = gen+ 1
        
        mu1 = gamma1*(xmax - Wbar)
        mu2 = mu2
        
        ####  generate the sorted two DFES --> must be sorted to apply ranks to ####
        #### THESE CHANGE WITH FITNESS (hence)
        DFE1 = sort(rexp(10000, rate = 1/abs(mu1))) # exponential dist, changes with fitness
        DFE2 = sort(rnorm(10000, mean = mu2, sd = sd)) #
        
        
        corr = rhoDesired
        cor_mat = matrix(corr, ncol = 2, nrow = 2)
        diag(cor_mat) = 1
        
        # make data sets with perfect desired corr
        mvdat = mvrnorm(10000, mu= c(0,0), Sigma = cor_mat, empirical = TRUE)
        
        # compute ranks of random data 
        rx <- rank(mvdat[ , 1], ties.method = "first")
        ry <- rank(mvdat[ , 2], ties.method = "first")
        
        # cor corrected DFEs --> of course, this isnt 100% perfect but corr always w/in +/- 0.001 of desired
        DFE1 = DFE1[rx]
        DFE2 = DFE2[ry]
        
        
        rho = cor(DFE1,DFE2) ## use the actual, precise rho of the JDFE for all calcs, since the above strategy is an approximation
        allrhos = c(allrhos, rho)
        M = rpois(1, popsize*U)
        parentIDs = sample(nextGenIDarray[ , 1], M, prob = (nextGenIDarray[ ,2]/popsize), replace = TRUE)
        
        for( i in parentIDs)
        {
          ######pick s in home
          
          # home S-val
          
          # home S-val
          
          if(mu1<0)
          {
            homeS = -1*rexp(1, rate = 1/abs(mu1))
          }else{
            homeS = rexp(1, rate = 1/abs(mu1))
          }
          
          
          if(homeS > 0)
          {
            numBen = numBen + 1 
          }
          
          
          
          muThishomeS = mu2 + ((rho*sd)/mu1)*(homeS - mu1)
          respS = rnorm(1, muThishomeS, sd)
          
          svec = c(homeS, respS)
          
          
          ####################
          
          for( q in 1:length(svec))
          {
            if(svec[q] < -1)
            {
              svec[q] = -1
            }
          }
          
          
          #remove indiv from its parent class if can be removed
          if(nextGenIDarray[i,2] != 0)
          {
            nextGenIDarray[i, 2] = nextGenIDarray[i, 2]- 1
            
            #make new row for this new mutant''s home data; 
            newRow = c(nextGenIDarray[dim(nextGenIDarray)[1], 1]+ 1, 1, startIDarray[i,3] + svec[enviroEvolveID])
            
            #bind newrow with nextgenIDarray
            nextGenIDarray = rbind(nextGenIDarray, newRow)
            
            #add row with fitness of new mutant in all environments to mutFitnesses
            newFit = mutFitnesses[i , ]+ svec
            mutFitnesses = rbind(mutFitnesses, newFit)
          }
          
          
        }
        
        output = list(nextGenIDarray, mutFitnesses)
        startIDarray = output[[1]]
        mutFitnesses = output[[2]]
        
        
        # find and save Wbars in all enviros 
        
        for( n in 1:(numDrugs))
        {
          allWbars[x, n] = sum(mutFitnesses[ , n] * (startIDarray[ ,2] / popsize))
        }
        
        # remove all muttypes with frequency of 0 from idarray and mutfitarray
        toRemove = c()
        if(length(startIDarray) > 3) # only have chance to remove if have more than 1 line in the startIDarray
        {
          for(q in 1:(dim(startIDarray)[1]))
          {
            if(startIDarray[q,2] == 0)
            {
              toRemove = c(toRemove, q)
              
            }
          }
        }
        
        
        if(length(toRemove) > 0)
        {
          startIDarray = startIDarray[-toRemove, ]
          mutFitnesses = mutFitnesses[ -toRemove, ]
          
          
          
          if(length(startIDarray) == 3) ## line 266 turns startIDarray to vector if go back down to 1 row only, need it to stay array for appropriate indexing
          {
            startIDarray = array(data = c(startIDarray), dim = c(1,3))
            mutFitnesses = array( data = c(mutFitnesses), dim= c(1, numDrugs))
          }
          
          countq = 0
          
          for(q in 1:(dim(startIDarray)[1]))
          {
            countq = countq + 1
            startIDarray[q, 1] = countq
          }
          
        }
        
        
      }# end pre run generation loop 
      
      
      
      ##############################################
      
      
      
      
      x0 = mean(mutFitnesses[,1])
      y0 = mean(mutFitnesses[,2])
      
      # startIDarray = array(data = c(1, popsize, 1), dim = c(1, 3))
      # mutFitnesses = array(data = 1, dim = c(1, numDrugs))
      allWbars = array(data = 0, dim = c(numGens, (numDrugs)))
      numBen = 0 
      
      for(x in 1:numGens)
      {
        if(x == 1)
        {
          r1vec = c(r1vec, mean(DFE1^2*(DFE1>0)))
          r2vec = c(r2vec, mean(DFE2*DFE1*(DFE1>0)) )
          
        }
        
        
        
        
        
        ### Simulate evo on DFE1, trans to traj on DFE2
        numMuts = dim(startIDarray)[1]
        
        #### Selection and Drift ####
        
        Wbar = sum(startIDarray[ , 3] * (startIDarray[ ,2] / popsize))
        
        # generate probability vector
        probvec = numeric(numMuts)
        for( i in 1:numMuts)
        {
          
          if((startIDarray[i,2] / popsize) + (startIDarray[i,2] / popsize)* (startIDarray[i,3] - Wbar) >0)
          {
            
            probvec[i] = (startIDarray[i,2] / popsize) + (startIDarray[i,2] / popsize)* (startIDarray[i,3] - Wbar)
          }else{
            probvec[i] = 0
          }
        }
        
        # generate IDarray (with only sel/drift) for next gen
        nextGenIDarray = array(data = 0, dim = c(dim(startIDarray)[1], dim(startIDarray)[2]))
        nextGenOffspring = rmultinom(1, size = popsize, prob = probvec)
        for( i in 1:numMuts)
        {
          nextGenIDarray[i, 1] = startIDarray[i, 1]
          nextGenIDarray[i, 2] = nextGenOffspring[i]
          nextGenIDarray[i, 3] = startIDarray[i, 3]
        }
        
        
        
        #### Mutations ####
        gen = gen+ 1
        
        mu1 = gamma1*(xmax - Wbar)
        mu2 = mu2
        
        ####  generate the sorted two DFES --> must be sorted to apply ranks to ####
        #### THESE CHANGE WITH FITNESS (hence)
        DFE1 = sort(rexp(10000, rate = 1/abs(mu1))) # exponential dist, changes with fitness
        DFE2 = sort(rnorm(10000, mean = mu2, sd = sd)) #
        
        
        corr = rhoDesired
        cor_mat = matrix(corr, ncol = 2, nrow = 2)
        diag(cor_mat) = 1
        
        # make data sets with perfect desired corr
        mvdat = mvrnorm(10000, mu= c(0,0), Sigma = cor_mat, empirical = TRUE)
        
        # compute ranks of random data 
        rx <- rank(mvdat[ , 1], ties.method = "first")
        ry <- rank(mvdat[ , 2], ties.method = "first")
        
        # cor corrected DFEs --> of course, this isnt 100% perfect but corr always w/in +/- 0.001 of desired
        DFE1 = DFE1[rx]
        DFE2 = DFE2[ry]
        
        
        rho = cor(DFE1,DFE2) ## use the actual, precise rho of the JDFE for all calcs, since the above strategy is an approximation
        allrhos = c(allrhos, rho)
        M = rpois(1, popsize*U)
        parentIDs = sample(nextGenIDarray[ , 1], M, prob = (nextGenIDarray[ ,2]/popsize), replace = TRUE)
        
        for( i in parentIDs)
        {
          ######pick s in home
          
          if(mu1<0)
          {
            homeS = -1*rexp(1, rate = 1/abs(mu1))
          }else{
            homeS = rexp(1, rate = 1/abs(mu1))
          }
          
          
          if(homeS > 0)
          {
            numBen = numBen + 1 
          }
          
          
          
          muThishomeS = mu2 + ((rho*sd)/mu1)*(homeS - mu1)
          respS = rnorm(1, muThishomeS, sd)
          
          svec = c(homeS, respS)
          
          
          ####################
          
          for( q in 1:length(svec))
          {
            if(svec[q] < -1)
            {
              svec[q] = -1
            }
          }
          
          
          #remove indiv from its parent class if can be removed
          if(nextGenIDarray[i,2] != 0)
          {
            nextGenIDarray[i, 2] = nextGenIDarray[i, 2]- 1
            
            #make new row for this new mutant''s home data; 
            newRow = c(nextGenIDarray[dim(nextGenIDarray)[1], 1]+ 1, 1, startIDarray[i,3] + svec[enviroEvolveID])
            
            #bind newrow with nextgenIDarray
            nextGenIDarray = rbind(nextGenIDarray, newRow)
            
            #add row with fitness of new mutant in all environments to mutFitnesses
            newFit = mutFitnesses[i , ]+ svec
            mutFitnesses = rbind(mutFitnesses, newFit)
          }
          
          
        }
        
        output = list(nextGenIDarray, mutFitnesses)
        startIDarray = output[[1]]
        mutFitnesses = output[[2]]
        
        
        # find and save Wbars in all enviros 
        
        for( n in 1:(numDrugs))
        {
          allWbars[x, n] = sum(mutFitnesses[ , n] * (startIDarray[ ,2] / popsize))
        }
        
        # remove all muttypes with frequency of 0 from idarray and mutfitarray
        toRemove = c()
        if(length(startIDarray) > 3) # only have chance to remove if have more than 1 line in the startIDarray
        {
          for(q in 1:(dim(startIDarray)[1]))
          {
            if(startIDarray[q,2] == 0)
            {
              toRemove = c(toRemove, q)
              
            }
          }
        }
        
        
        if(length(toRemove) > 0)
        {
          startIDarray = startIDarray[-toRemove, ]
          mutFitnesses = mutFitnesses[ -toRemove, ]
          
          
          
          if(length(startIDarray) == 3) 
          {
            startIDarray = array(data = c(startIDarray), dim = c(1,3))
            mutFitnesses = array( data = c(mutFitnesses), dim= c(1, numDrugs))
          }
          
          countq = 0
          
          for(q in 1:(dim(startIDarray)[1]))
          {
            countq = countq + 1
            startIDarray[q, 1] = countq
          }
          
        }
        
        
      }# end generation loop 
      allNumben = c(allNumben, numBen)
      
      
      allWbarsAllIterations[ , , it] = allWbars
    } # end iteration loop
    
    averageAllWbar = array(data = 0, dim = c(numGens, (numDrugs)))
    for( gen2 in 1:numGens)
    {
      for(d2 in 1:(numDrugs))
      {
        averageAllWbar[gen2,d2] = mean(allWbarsAllIterations[gen2,d2, ])
      }
    }
    
    meanendfit1 = c(meanendfit1, averageAllWbar[numGens,1])
    meanendFit2 = c(meanendFit2, averageAllWbar[numGens,2])
    meanr2vec = c(meanr2vec, mean(r2vec))
    meanrhovec = c(meanrhovec, mean(allrhos))
    DrugSDiDataSet = array(data = 0, dim = c(numGens, numDrugs))
    for( gen2 in 1:numGens)
    {
      for(d2 in 1:(numDrugs))
      {
        DrugSDiDataSet[gen2,d2] = sd((allWbarsAllIterations[gen2,d2, ]))
      }
    }
    
    
    
    y0 = averageAllWbar[1,2]
    
    endGens = numGens
    f1theo = c((x0))
    f2theo = c((y0))
    min1Vec = c(x0)
    max1Vec = c(x0)
    
    min2Vec = c(y0)
    max2Vec = c(y0)
    
    
    N = popsize
    Ub = U
    rho = mean(allrhos)
    tadj = 2*N*Ub 
    for(t in 1:(min(endGens) -1 ))
    {
      f1theo = c(f1theo, xmax - 1/(2*gamma1^2*t*tadj + 1/(xmax-(x0))))
      min1Vec = c(min1Vec,averageAllWbar[t,1]-DrugSDiDataSet[t+1,1] )  # max and min are +/- 1 sd
      max1Vec = c(max1Vec,averageAllWbar[t,1]+DrugSDiDataSet[t+1,1])
      
      
      f2theo = c(f2theo, (y0) + ((mu2 + rho*sd)/(2*gamma1))*(log(1 + 2*gamma1^2*(xmax - (x0))*t*tadj)))
      min2Vec = c(min2Vec,averageAllWbar[t,2]-DrugSDiDataSet[t+1,2]  )
      max2Vec = c(max2Vec,averageAllWbar[t,2]+DrugSDiDataSet[t+1,2]  )
      
    }
    
    
    
    
    plotDf = data.frame(c(rep('Env1', times = min(endGens)),rep('Env2', times = min(endGens)),rep('Env1', times = min(endGens))
                          ,rep('Env2', times = min(endGens))),
                        c(rep('Sim', times = min(endGens*2)),rep('Theo', times = min(endGens*2))),
                        c((c(x0,averageAllWbar[1:min(endGens-1), 1])), (c(y0,averageAllWbar[1:min(endGens-1), 2])), f1theo, f2theo),
                        rep(1:min(endGens), times= 4), c(min1Vec, min2Vec, min1Vec, min2Vec), c(max1Vec, max2Vec,max1Vec,max2Vec))
    colnames(plotDf) = c('Env', 'Type', "Fitness", 'Time','min', 'max')
    
    

    
    
    
    plt = ggplot(plotDf, aes(x = Time, y = Fitness)) + geom_line(aes(color = Env, lty = Type)) +ylim(-15,3) 
    
    plt = plt + theme_bw()+theme(axis.title = element_text(size = 15), axis.text = element_text(size = 35),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill = NA, color = 'black', size = 1.5), axis.line = element_line(colour = "black"), legend.position = 'none')+
      geom_ribbon(aes(x = Time, y = Fitness,ymin = min, ymax = max, fill = Env,group = Env), alpha = 0.3, linetype = 0)+ggtitle( as.character(N*U))
    print(plt)
    
    
    name = paste0('FitTrajPlot2_tadj2NUb_Nu_', popsize*U,'_rho_',rhoDesired, '.pdf')
    ggsave(plt,file= name, width = 6.18, height = 5.5)
    
    
    
    
    
    ### save the plotDfs for all 
    
    if(rho > 0)
    {
      assign(paste("plotDf_tadj2Nub__rho", rhoDesired, '_Nu_' ,popsize*U, sep = ""), plotDf) 
    }else{
      assign(paste("plotDf_tadj2Nub_rho_neg_",abs(rhoDesired),'_Nu_' ,popsize*U, sep = ""), plotDf) 
    }
    
    
    y0 = averageAllWbar[1,2] 
    
    
    
  }
  
  
  
  U = U
  meanr2vec = meanr2vec
  meanendFit2 = meanendFit2
  
  df = data_frame(meanr2vec, meanendFit2-y0)
  colnames(df) = c('r2', 'MeandeltaY')
  newPlt = ggplot(df, aes(x = r2, y = MeandeltaY)) + 
    geom_point() +
    geom_smooth(method='lm',formula = y~x, se = T, color = 'black', lwd = 0.8)+ 
    scale_x_continuous(label=scientific_10,guide = guide_axis(check.overlap = TRUE)) +
    theme(axis.title = element_text(size = 15), axis.text = element_text(size = 35),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = NA, color = 'black', size = 1.5), 
          axis.line = element_line(colour = "black"))
  
  print(newPlt)  
  
  name = paste0('r2vsDeltay_Nu_', popsize*U, '.pdf')
  ggsave(newPlt, file = name, width = 6.18, height = 5.5) 
  
  summary(lm(meanendFit2~meanr2vec))
  
  name = paste0('r2vsDeltay_Nu_', popsize*U, '.csv')
  write.csv(df, name)
  

  
  allMeanR2s = c(allMeanR2s, meanr2vec)
  allMeanEndFits2 = c(allMeanEndFits2, meanendFit2)
  allrhos = c(allrhos, rhoVec)
  allUs = c(allUs, rep(U, times = length(allMeanR2s)))
}








####### PANEL A: how JDFES change over time ####

timeVec = c(1, 100, 500, 1000)

# parameters
N = 1.6*10^5
U = 10^-4
Ub = 10^-4
xmax = 3
ymax = 3
gamma1 = 0.05
gamma2 = 1
sd = 0.01
mu2 = -0.1
x0 = 1
y0 = 1
rho = -0.9

n = 0
JDFEpltList = list()
f1Vec = c()
f2Vec = c()
for( t in timeVec)
{
  n = n+1
  #theo Fit 
  tadj = 1 #want the inputted times to be the adjusted times 
  
  f1theo =  xmax - 1/(2*gamma1^2*t*tadj + 1/(xmax-(x0)))
  f2theo = (y0) + ((mu2 + rho*sd)/(2*gamma1))*(log(1 + 2*gamma1^2*(xmax - (x0))*t*tadj))
  
  
  
  
  f1Vec = c(f1Vec, f1theo)
  f2Vec = c(f2Vec, f2theo)
  
  #use to make DFE1 and 2
  mu1 = gamma1*(xmax - f1theo)
  mu2 = mu2
  
  ####  generate the sorted two DFES --> must be sorted to apply ranks to ####
  #### THESE CHANGE WITH FITNESS (hence)
  DFE1 = sort(rexp(10000, rate = 1/abs(mu1))) # exponential dist, changes with fitness
  DFE2 = sort(rnorm(10000, mean = mu2, sd = sd)) #
  
  rhoDesired = rho
  corr = rhoDesired
  cor_mat = matrix(corr, ncol = 2, nrow = 2)
  diag(cor_mat) = 1
  
  # make data sets with perfect desired corr
  mvdat = mvrnorm(10000, mu= c(0,0), Sigma = cor_mat, empirical = TRUE)
  
  # compute ranks of random data 
  rx <- rank(mvdat[ , 1], ties.method = "first")
  ry <- rank(mvdat[ , 2], ties.method = "first")
  
  # cor corrected DFEs --> of course, this isnt 100% perfect but corr always w/in +/- 0.001 of desired
  DFE1 = DFE1[rx]
  DFE2 = DFE2[ry]
  
  
  # ggplot with points colored by denisty  and the built in drawing of contour lines
  # is there better way?well canjust do the contours.. but these will not be smooth 
  df = data.frame(DFE1, DFE2)
  
  
  
  if(n == 1)
  {
    JDFEpltList[[n]] = ggplot(df, aes(x = DFE1, y = DFE2)) + stat_density_2d(n = 200, bins = 5, color = 'dodgerblue4', lwd = 2)+
      theme_bw()+theme(panel.border = element_rect(colour = 'black', fill = NA, size = 1), 
                       panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                       axis.line = element_line(colour = "black"),   
                       axis.text.x =element_text(size = 35,family = 'Helvetica',color = c('black','transparent', 'black','transparent','black')), 
                       axis.text.y =element_text(size = 35,family = 'Helvetica', color = 'black'),
                       axis.title = element_blank())+  ylim(-0.125,-0.075) + xlim(-0.005,0.12)
    
    
    name = paste0('JDFE_TheoTime_',t,'_gamma_',gamma1, '.pdf')
    ggsave(JDFEpltList[[n]], file = name, width = 6.18, height = 5.5)
    
  }else if(n == 2){
    
    JDFEpltList[[n]] = ggplot(df, aes(x = DFE1, y = DFE2)) + stat_density_2d(n = 200, bins = 5, color = 'palegreen4', lwd = 2)+
      theme_bw()+theme(panel.border = element_rect(colour = 'black', fill = NA, size = 1), 
                       panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                       axis.line = element_line(colour = "black"),   
                       axis.text.x =element_text(size = 35,family = 'Helvetica',color = c('black','transparent', 'black','transparent','black')), 
                       axis.text.y =element_text(size = 35,family = 'Helvetica', color = 'black'),
                       axis.title = element_blank())+  ylim(-0.125,-0.075) + xlim(-0.005,0.12)
    
    name = paste0('JDFE_TheoTime_',t,'_gamma_',gamma1, '.pdf')
    ggsave(JDFEpltList[[n]], file = name, width = 6.18, height = 5.5)
    
    
  }else if(n ==3){
    JDFEpltList[[n]] = ggplot(df, aes(x = DFE1, y = DFE2)) + stat_density_2d(n = 200, bins = 5, color = 'palegreen4', lwd = 2)+
      theme_bw()+theme(panel.border = element_rect(colour = 'black', fill = NA, size = 1), 
                       panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                       axis.line = element_line(colour = "black"),   
                       axis.text.x =element_text(size = 35,family = 'Helvetica',color = c('black','transparent', 'black','transparent','black')), 
                       axis.text.y =element_text(size = 35,family = 'Helvetica', color = 'black'),
                       axis.title = element_blank())+  ylim(-0.125,-0.075) + xlim(-0.005,0.12)
    
    ## to make bottom left the coordiante plane of fit in home vs non-home
    # run code below this loop
    
    name = paste0('JDFE_TheoTime_',t,'_gamma_',gamma1, '.pdf')
    ggsave(JDFEpltList[[n]], file = name, width = 6.18, height = 5.5)
    
  }else if(n ==4){
    JDFEpltList[[n]] = ggplot(df, aes(x = DFE1, y = DFE2)) + stat_density_2d(n = 200, bins = 5, color = 'hotpink4', lwd = 2)+
      theme_bw()+theme(panel.border = element_rect(colour = 'black', fill = NA, size = 1), 
                       panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                       axis.line = element_line(colour = "black"),   
                       axis.text.x =element_text(size = 35,family = 'Helvetica',color = c('black','transparent', 'black','transparent','black')), 
                       axis.text.y =element_text(size = 35,family = 'Helvetica', color = 'black'),
                       axis.title = element_blank())+  ylim(-0.125,-0.075) + xlim(-0.005,0.12)
    
    name = paste0('JDFE_TheoTime_',t,'_gamma_',gamma1, '.pdf')
    ggsave(JDFEpltList[[n]], file = name, width = 6.18, height = 5.5)
  }
  
  
  ## add text for num gens in the top right of each panel 
  
}


df = data.frame(f1Vec[c(1,2,4)], f2Vec[c(1,2,4)], timeVec[c(1,2,4)])
colnames(df) = c('f1', 'f2', 'time')

JDFEpltList[[3]] = ggplot(df, aes(x = f1, y = f2)) + geom_point(color = c('dodgerblue4', 'palegreen4', 'hotpink4'), shape = 17, size = 10)+
  theme_bw()+theme(panel.border = element_rect(colour = 'black', fill = NA, size = 1), 
                   panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                   axis.line = element_line(colour = "black"),   
                   axis.text.x =element_text(size = 35,family = 'Helvetica',color ='black'), 
                   axis.text.y =element_text(size = 35,family = 'Helvetica', color = 'black'),
                   axis.title = element_blank())

name = paste0('Fitness1vsFitness2_gamma_',gamma1, '.pdf')
ggsave(JDFEpltList[[3]], file = name, width = 6.18, height = 5.5)



f1Vec = formatC(f1Vec, digits = 2, format = "f")


jdfeOverTimeEpistasisFig = ggarrange(JDFEpltList[[1]], JDFEpltList[[2]], JDFEpltList[[3]], JDFEpltList[[4]], 
                                     ncol = 2, nrow = 2, widths = c(1.25,1), 
                                     labels =f1Vec,
                                     font.label = list(size = 16, color = "black" ),
                                     hjust = c(-7.8,-6,-20,-6),
                                     vjust = 2)


print(jdfeOverTimeEpistasisFig)

ggsave(jdfeOverTimeEpistasisFig, file = 'jdfeOverTimewEpistasis.pdf')



###################################### F I G U R E     4  ################################################################


# (optional) set working directory to where you'd like to save the outputs of this section, example below 
# setwd("~/Documents/JDFE Project/Figures/Figure ABR JDFES")

### NOTE: WHEREVER YOU SET THE WORKING DIRECTORY MUST CONTAIN THE FILES "KOGrDATAPrtoStop.csv' AND 'WTReplicate.csv' FOR THIS CODE TO RUN 


####Import data ####

# Import KO growth rate data  so can have full data (including KO names)
KOGrowthRates = read.csv('KOGrDATAPrtoStop.csv', header = 1, sep = ",")
# remove KOID and FOX and AMP (last 2 drugs), as dont have WT replicate data for them
# remove KOID and FOX and AMP (last 2 drugs), as dont have WT replicate data for them
KOPrtoStop = KOGrowthRates[,11]
KOGrowthRates = as.matrix(KOGrowthRates[,2:8])


##### Import Data #####
WTReplicateData = read.csv('WTReplicate.csv', header = 1, sep = ":")
# remove excess rows; keep only  No Drug and 1st 6 drug enviros and remove header
WTReplicateData = as.matrix(WTReplicateData[,2:8])

# Important Variables #
WTGrowthRateMeans = colMeans(WTReplicateData)
WTGrowthRateSDs   = colSds(WTReplicateData) 
numKnockouts = dim(KOGrowthRates)[1]
numDrugs = dim(KOGrowthRates)[2]


sValAllMuts = array(0, dim(KOGrowthRates))
for( i in 1:dim(KOGrowthRates)[1])
{
  sValAllMuts[i, ] = KOGrowthRates[i,] - WTGrowthRateMeans
}



## find the right cutoff values for calling beneficial muts so have FDR ~ 25% for beneficial muts in  all drugs

testCutoffs = seq(0.01, 0.5, length.out = 100)

nameVec = c()
fdrVec = c()
cutoffVec = c()

for(home in 1:7)
{
  
  fdrsThisHome = c()
  cutoffsthishome = c()
  for(cutoff in testCutoffs)
  {
    
    evo = home
    if(evo == 1){
      title = "No Drug"
    }else if(evo == 2){
      title = "CHL"
    }else if(evo == 3){
      title = "CPR"
    }else if(evo ==4){
      title = "MEC"
    }else if(evo == 5){
      title = "NIT"
    }else if(evo == 6){
      title = "TET"
    }else{
      title = "TMP"
    }
    
    
    benCutoff = qnorm((1-cutoff), mean = 0, sd = WTGrowthRateSDs[home]) 
    deltCutoff = qnorm(cutoff, mean = 0, sd = WTGrowthRateSDs[home])
    
    DFE = sValAllMuts[,home]
    errorDist = rnorm(400, mean = 0, sd = WTGrowthRateSDs[home])
    breaks = seq(min(c(DFE, errorDist)), max(c(DFE, errorDist)), length.out = 50)

    
    
    mutTypeVec = (sValAllMuts[,(home)]>benCutoff)*2 +  (sValAllMuts[,(home)]<deltCutoff)
    
    beneficialMutCount = sum(mutTypeVec == 2)
    
    expectedNumFalsePositive = sum(sValAllMuts[,home]>0)*cutoff*2 ## 
    
    fdr = expectedNumFalsePositive/beneficialMutCount ## want to be < 1 ,, ie want to be calling more than would be expected just by noise 
    
    fdrsThisHome = c(fdrsThisHome, fdr)
    cutoffsthishome = c(cutoffsthishome, cutoff)
  }
  
  cutofftouse = which.min(abs(fdrsThisHome-0.25))
  
  
  nameVec = c(nameVec, title)
  cutoffVec = c(cutoffVec, cutoffsthishome[cutofftouse])
  fdrVec = c(fdrVec, fdrsThisHome[cutofftouse])
}

allFDR_Results = data.frame(nameVec, cutoffVec, fdrVec)







# CREATE BENEFICIAL AND DELTERIOUS FITNESS VALUE CUTTOFF VECS 
numDrugs = 7
beneficialCutoffVec = numeric(numDrugs)
deleteriousCutoffVec = numeric(numDrugs)

cutoffVec = allFDR_Results$cutoffVec#the values of cutoff so have lowest false disocvery rate in all 

for( i in 1:numDrugs)
{
  cutoff = cutoffVec[i]
  ### mean = 0 bc svals are calc as subtract wt gr
  beneficialCutoffVec[i] = qnorm((1-cutoff), mean = 0, sd = WTGrowthRateSDs[i]) 
  deleteriousCutoffVec[i] = qnorm((cutoff), mean = 0, sd = WTGrowthRateSDs[i]) 
  
}




# CREATE DATA FRAME WHICH IDENTIFIES ALL KOS AS 'BENEFICIAL', 'DELETERIOUS', OR 'NEUTRAL' IN ALL ENVIRONMENTS
mutTypeDf = data.frame(KOGrowthRates)
beneficialMutCount = numeric(7)
deleteriousMutCount = numeric(7)
totalPossibleBenMuts = numeric(7)
for(home in 1:7)
{
  
  benCutoff = beneficialCutoffVec[(home)]
  deltCutoff = deleteriousCutoffVec[(home)]
  
  ##home -1 bc in svalallmuts dont have 1st col of mutID
  mutTypeDf[,home] = (sValAllMuts[,(home)]>benCutoff)*2 +  (sValAllMuts[,(home)]<deltCutoff)
  
  
  ## also count Num Muts which are beneficial in each env and deleterious in each env
  beneficialMutCount[(home)] = sum(mutTypeDf[,(home)] == 2)
  deleteriousMutCount[(home)] =sum(mutTypeDf[,(home)] == 1)
  
  totalPossibleBenMuts[(home)] = sum(sValAllMuts[,(home)]>0)
  
}

KOnames = read.csv('KOGrDATAPrtoStop.csv', header = 1, sep = ",")[,1]
mutTypeDf$KOnames = KOnames
mutTypeDf = mutTypeDf[c(8, 1:7)]
mutTypeDf = mutTypeDf[c(1,4,5,6,7)]  # only save data on drugs focus on in text

## save 
write.csv(mutTypeDf, "mutTypeDf.csv")


#
MutTypeCounts = data.frame(c('No Drug', 'CHL', 'CPR','MEC', 'NIT', 'TET', 'TMP'), beneficialMutCount, deleteriousMutCount, totalPossibleBenMuts)
colnames(MutTypeCounts) = c('Drug', "Num Benefical Muts Called", 'Num Deleterious Muts Called', 'NumMutsGreaterThan0')


MutTypeCounts$ExpectedNumFalseDiscovery = (MutTypeCounts$NumMutsGreaterThan0)*cutoffVec*2 ## expect to call cutoff% *2 of the muts beneficial that are actuallt WT with a cutoff at top 10% of WT dist
MutTypeCounts$FalseDiscoveryRate = ( MutTypeCounts$ExpectedNumFalseDiscovery/MutTypeCounts$`Num Benefical Muts Called`)






write.csv(MutTypeCounts, 'MutTypeCounts.csv')







### graph the JDFEs




homeMeanVec = numeric(16)
respMeanVec = numeric(16)
homeVarVec = numeric(16)
respVarVec = numeric(16)
coVarVec = numeric(16)
r1Vec = numeric(16)
r2Vec = numeric(16)
D11Vec = numeric(16)
D12Vec = numeric(16)
D22Vec = numeric(16)
cVec = numeric(16)
timeTo90pSameasR2PredVec = numeric(16)
homeNameVec = numeric(16)
respNameVec = numeric(16)





homeMeanVec_sig = numeric(16)
respMeanVec_sig = numeric(16)
homeVarVec_sig = numeric(16)
respVarVec_sig = numeric(16)
coVarVec_sig = numeric(16)
r1Vec_sig = numeric(16)
r2Vec_sig = numeric(16)
D11Vec_sig = numeric(16)
D12Vec_sig = numeric(16)
D22Vec_sig = numeric(16)
cVec_sig = numeric(16)
timeTo90pSameasR2PredVec_sig = numeric(16)





pltListHM = list()


## drug order for plot
## CHL, CPR, NIT, TET, TMP, MEC, No Drug
## (which is, in numbers of the cols from the original excel sheet c(2,3,5,6,7,4,1))
#3 get rid of 2 and 7 though

no = 0
for(resp in c(3,5,6,4)){
  
  for(home in c(3,5,6,4)) 
  {
    no = no+1
    
    ### save homeMean, respMean, homeVariance, respVariance, covariance, R1, R2, D11, D12, D22, c, and timeto90%sameasr2sign values 
    homeSvals = sValAllMuts[,home]
    
    homeMeanVec[no] = mean(homeSvals)
    respMeanVec[no] = mean(sValAllMuts[,resp])
    homeVarVec[no] = var(homeSvals)
    respVarVec[no] = var(sValAllMuts[,resp])
    coVarVec[no] = cov(sValAllMuts[,home], sValAllMuts[,resp])
    
    homeSvals = sValAllMuts[,home]*(sValAllMuts[,home]>=0) # make all neg values 0 so dont contrib to calc
    r1Vec[no] = mean((homeSvals)^2)
    r2Vec[no] = mean((homeSvals)*(sValAllMuts[,resp]))
    D11Vec[no] = mean((homeSvals)^3)
    D12Vec[no] = mean((homeSvals)^2*(sValAllMuts[,resp]))
    D22Vec[no] = mean((homeSvals)*(sValAllMuts[,resp])^2)
    cVec[no] = sqrt(mean((homeSvals)*(sValAllMuts[,resp])^2))/(mean((homeSvals)*(sValAllMuts[,resp])))
    timeTo90pSameasR2PredVec[no] = (1.28^2*mean((homeSvals)*(sValAllMuts[,resp])^2))/mean((homeSvals)*(sValAllMuts[,resp]))^2
    
    
    
    ## using only signficant values 
    mutsToUse = sValAllMuts[sValAllMuts[,home]>beneficialCutoffVec[home] & (sValAllMuts[,resp]>beneficialCutoffVec[resp] | sValAllMuts[,resp]<deleteriousCutoffVec[resp]), c(home,resp)]
    
    if(length(mutsToUse) < 4)
    {
      
      r1Vec_sig[no] = NA
      r2Vec_sig[no] = NA
      D11Vec_sig[no] = NA
      D12Vec_sig[no] = NA
      D22Vec_sig[no] = NA
      cVec_sig[no] = NA
      timeTo90pSameasR2PredVec_sig[no] = NA
      
    }else{
      
      homeSvals =  mutsToUse[,1]
      respSvals = mutsToUse[,2]
      r1Vec_sig[no] = mean((homeSvals)^2)
      r2Vec_sig[no] = mean((homeSvals)*(respSvals))
      D11Vec_sig[no] = mean((homeSvals)^3)
      D12Vec_sig[no] = mean((homeSvals)^2*(respSvals))
      D22Vec_sig[no] = mean((homeSvals)*(respSvals)^2)
      cVec_sig[no] = sqrt(mean((homeSvals)*(respSvals)^2))/(mean((homeSvals)*(respSvals)))
      timeTo90pSameasR2PredVec_sig[no] = (1.28^2*mean((homeSvals)*(respSvals)^2))/mean((homeSvals)*(respSvals))^2
      
    }
    
    
    
    
    evo = home
    if(evo == 1){
      title = "No Drug"
    }else if(evo == 2){
      title = "CHL"
    }else if(evo == 3){
      title = "CPR"
    }else if(evo ==4){
      title = "MEC"
    }else if(evo == 5){
      title = "NIT"
    }else if(evo == 6){
      title = "TET"
    }else{
      title = "TMP"
    }
    
    evo = resp
    if(evo == 1){
      title2 = "No Drug"
    }else if(evo == 2){
      title2 = "CHL"
    }else if(evo == 3){
      title2 = "CPR"
    }else if(evo ==4){
      title2 = "MEC"
    }else if(evo == 5){
      title2 = "NIT"
    }else if(evo == 6){
      title2 = "TET"
    }else{
      title2 = "TMP"
    }
    
    
    
    homeNameVec[no] = title
    respNameVec[no] = title2
    
    
    homeBenCutoff = beneficialCutoffVec[home]
    homeDeltCutoff = deleteriousCutoffVec[home]
    respBenCutoff = beneficialCutoffVec[resp]
    respDeltCutoff = deleteriousCutoffVec[resp]
    
    
    ## 2 is beneficial, 1 is delterious, 0 is neutral
    mutTypeHome = (sValAllMuts[,home]>homeBenCutoff)*2 +  (sValAllMuts[,home]<homeDeltCutoff)
    mutTypeResp = (sValAllMuts[,resp]>homeBenCutoff)*2 +  (sValAllMuts[,resp]<homeDeltCutoff)
    
    if(home!=resp)
    {
      

      
      df2 = data.frame(sValAllMuts[,home], sValAllMuts[,resp], as.factor(mutTypeHome),as.factor(mutTypeResp))
      colnames(df2) = c("HomeS", "RespS", 'mutTypeHome', 'mutTypeResp')
      
      
      df2$mutTypeBoth = interaction(df2$mutTypeHome, df2$mutTypeResp)
      mutTypeBoth = df2$mutTypeBoth
      
      df2 = df2[order(mutTypeBoth),] 
      
      
      ## if in 1st col AND in last row  save with both x and y axis label and x and y tick labels
      if(title == 'CPR' & title2 == 'MEC'){
        
        # make df2 so all grey points are first 
        
        
        
        pltListHM[[no]] = ggplot(df2, aes(x=HomeS, y=RespS, color = interaction(mutTypeHome,mutTypeResp,sep="-",lex.order=TRUE))) +  xlim(-1,0.3) +ylim(-1, 0.3)+
          geom_point(aes(HomeS, RespS), alpha = 0.5)+ scale_color_manual(values = c('grey', 'grey', 'grey', 'grey', 'grey', 'grey', '#13944a', 'blue','orange'), drop = 'FALSE') +
          theme_bw()+ labs(colour="Home-Resp") + xlab(title) + ylab(title2)+
          theme(panel.border = element_rect(color = 'black', fill = NA, size = 1.5), panel.grid.major = element_blank(),
                legend.position = 'none',axis.text.x = element_text(size = 15, color = c('black', 'transparent', 'black')), axis.text.y = element_text(size = 15, color = 'black'),panel.grid.minor = element_blank(), 
                axis.line = element_line(colour = "black"), axis.title.x = element_text(size = 20, color = 'black'), axis.title.y = element_text( size = 20, color = 'black')) +
          geom_hline(yintercept =  0, linetype = 'dashed') + 
          geom_vline(xintercept =  0, linetype = 'dashed')+
          annotate(geom="text", fontface = 'bold', size = 3,x=-0.95, y=0.3, label=as.character(sum(mutTypeHome==2 & mutTypeResp == 2)),color="orange")+
          annotate(geom="text", fontface = 'bold', size = 3,x=-0.95, y=0.1, label=as.character(sum(mutTypeHome==2 & mutTypeResp == 1)),color="blue")+
          annotate(geom="text", fontface = 'bold', size = 3,x=-0.95, y=0.2, label=as.character(sum(mutTypeHome==2 & mutTypeResp == 0)),color="#13944a")
        
        
        ggsave(pltListHM[[no]], file = paste0('JDFE_', no, '.pdf') )
      }else if(title == 'CPR')## if in 1st col , save with y axis label and y tick labels
      {
        
        pltListHM[[no]] = ggplot(df2, aes(x=HomeS, y=RespS, color = interaction(mutTypeHome,mutTypeResp,sep="-",lex.order=TRUE))) +  xlim(-1,0.3) +ylim(-1, 0.3)+
          geom_point(aes(HomeS, RespS), alpha = 0.5)+ scale_color_manual(values = c('grey', 'grey', 'grey', 'grey', 'grey', 'grey', '#13944a', 'blue','orange'), drop = 'FALSE') +
          theme_bw()+ labs(colour="Home-Resp") + xlab(title) + ylab(title2)+
          theme(panel.border = element_rect(color = 'black', fill = NA, size = 1.5), panel.grid.major = element_blank(),
                legend.position = 'none',axis.text.x = element_blank(), axis.text.y = element_text(size = 15, color = 'black'),panel.grid.minor = element_blank(), 
                axis.line = element_line(colour = "black"), axis.title.x = element_blank(), axis.title.y = element_text( size = 20, color = 'black')) +
          geom_hline(yintercept =  0, linetype = 'dashed') + 
          geom_vline(xintercept =  0, linetype = 'dashed')+
          annotate(geom="text", fontface = 'bold', size = 3,x=-0.95, y=0.3, label=as.character(sum(mutTypeHome==2 & mutTypeResp == 2)),color="orange")+
          annotate(geom="text", fontface = 'bold', size = 3,x=-0.95, y=0.1, label=as.character(sum(mutTypeHome==2 & mutTypeResp == 1)),color="blue")+
          annotate(geom="text", fontface = 'bold', size = 3,x=-0.95, y=0.2, label=as.character(sum(mutTypeHome==2 & mutTypeResp == 0)),color="#13944a")
        
        
        ggsave(pltListHM[[no]], file = paste0('JDFE_', no, '.pdf') )
        
        
      }else if (title2 == 'MEC') { ## or in last row , save with x axis label and x ticks
        
        pltListHM[[no]] = ggplot(df2, aes(x=HomeS, y=RespS, color = interaction(mutTypeHome,mutTypeResp,sep="-",lex.order=TRUE))) +  xlim(-1,0.3) +ylim(-1, 0.3)+
          geom_point(aes(HomeS, RespS), alpha = 0.5)+ scale_color_manual(values = c('grey', 'grey', 'grey', 'grey', 'grey', 'grey', '#13944a', 'blue','orange'), drop = 'FALSE') +
          theme_bw()+ labs(colour="Home-Resp") +xlab(title) + ylab(title2)+
          theme(panel.border = element_rect(color = 'black', fill = NA, size = 1.5), panel.grid.major = element_blank(),
                legend.position = 'none',axis.text.y = element_blank(), axis.text.x = element_text(size = 15, color = c('black', 'transparent', 'black')),panel.grid.minor = element_blank(), 
                axis.line = element_line(colour = "black"), axis.title.y = element_blank(), axis.title.x= element_text( size = 20, color = 'black')) +
          geom_hline(yintercept =  0, linetype = 'dashed') + 
          geom_vline(xintercept =  0, linetype = 'dashed')+
          annotate(geom="text", fontface = 'bold', size = 3,x=-0.95, y=0.3, label=as.character(sum(mutTypeHome==2 & mutTypeResp == 2)),color="orange")+
          annotate(geom="text", fontface = 'bold', size = 3,x=-0.95, y=0.1, label=as.character(sum(mutTypeHome==2 & mutTypeResp == 1)),color="blue")+
          annotate(geom="text", fontface = 'bold', size = 3,x=-0.95, y=0.2, label=as.character(sum(mutTypeHome==2 & mutTypeResp == 0)),color="#13944a")
        
        
        ggsave(pltListHM[[no]], file = paste0('JDFE_', no, '.pdf') )
        
      }else{
        pltListHM[[no]] = ggplot(df2, aes(x=HomeS, y=RespS, color = interaction(mutTypeHome,mutTypeResp,sep="-",lex.order=TRUE))) +  xlim(-1,0.3) +ylim(-1, 0.3)+
          geom_point(aes(HomeS, RespS), alpha = 0.5)+ scale_color_manual(values = c('grey', 'grey', 'grey', 'grey', 'grey', 'grey', '#13944a', 'blue','orange'), drop = 'FALSE') +
          theme_bw()+ labs(colour="Home-Resp") +xlab(title) + ylab(title2)+
          theme(panel.border = element_rect(color = 'black', fill = NA, size = 1.5), panel.grid.major = element_blank(),
                legend.position = 'none',axis.text.x = element_blank(), axis.text.y = element_blank(),panel.grid.minor = element_blank(), 
                axis.line = element_line(colour = "black"), axis.title.x = element_blank(), axis.title.y = element_blank()) +
          geom_hline(yintercept =  0, linetype = 'dashed') + 
          geom_vline(xintercept =  0, linetype = 'dashed')+
          annotate(geom="text", fontface = 'bold', size = 3,x=-0.95, y=0.3, label=as.character(sum(mutTypeHome==2 & mutTypeResp == 2)),color="orange")+
          annotate(geom="text", fontface = 'bold', size = 3,x=-0.95, y=0.1, label=as.character(sum(mutTypeHome==2 & mutTypeResp == 1)),color="blue")+
          annotate(geom="text", fontface = 'bold', size = 3,x=-0.95, y=0.2, label=as.character(sum(mutTypeHome==2 & mutTypeResp == 0)),color="#13944a")
        
        
        ggsave(pltListHM[[no]], file = paste0('JDFE_', no, '.pdf') )
        
      }
      
      
      
      
      
      
      
      
      
      
      
    }else{
      df2 = data.frame(c(sValAllMuts[,home],(WTReplicateData[,home]-WTGrowthRateMeans[home])), c(rep('KO', times = length(sValAllMuts[,home])), rep('WT', times = length(WTReplicateData[,home]))))
      colnames(df2) = c("Fitness", 'Type')
      
      if(title == 'CPR')
      {
        
        
        
        pltListHM[[no]] = ggplot(df2, aes(x = Fitness, fill = Type)) + geom_histogram(bins = 50) +
          theme_bw() + xlab(title) + ylab(title2)+ scale_fill_manual(values = c('grey40', 'red'))+
          theme(panel.border = element_rect(color = 'black', fill = NA, size = 1.5), panel.grid.major = element_blank(),axis.text.y = element_text(size = 15, color = 'black'), axis.text.x = element_blank(),
                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_blank(), axis.title.y = element_text( size = 20, color = 'black'), legend.position = 'none') +
          geom_vline(xintercept =  0, linetype = 'dashed')+ coord_cartesian(ylim=c(0,1000), xlim = c(-1,0.3)) +
          scale_y_continuous(expand = c(0,0)) +
          annotate(geom="text", fontface = 'bold',size = 3,x=-0.9, y=940, label=as.character(sum(mutTypeHome==2 & mutTypeResp == 2)),color="black")
        
        
      }else if (title2 == 'MEC'){
        pltListHM[[no]] =ggplot(df2, aes(x = Fitness, fill = Type)) + geom_histogram(bins = 50) +
          theme_bw() + xlab(title) + ylab(title2)+ scale_fill_manual(values = c('grey40', 'red'))+
          theme(panel.border = element_rect(color = 'black', fill = NA, size = 1.5), panel.grid.major = element_blank(),axis.text.x = element_text(size = 15, color = c('black', 'transparent', 'black')), axis.text.y = element_blank(),
                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_text( size = 20, color = 'black'), axis.title.y = element_blank(), legend.position = 'none') +geom_vline(xintercept =  0, linetype = 'dashed')+ 
          coord_cartesian(ylim=c(0,1000), xlim = c(-1,0.3)) +
          scale_y_continuous(expand = c(0,0)) +
          annotate(geom="text", fontface = 'bold', size = 3,x=-0.9, y=940, label=as.character(sum(mutTypeHome==2 & mutTypeResp == 2)),color="black")
        
        
      }else{
        pltListHM[[no]] = ggplot(df2, aes(x = Fitness, fill = Type)) + geom_histogram(bins = 50) +
          theme_bw() +  scale_fill_manual(values = c('grey40', 'red'))+
          theme(panel.border = element_rect(color = 'black', fill = NA, size = 1.5), panel.grid.major = element_blank(),axis.text.x = element_blank(), axis.text.y = element_blank(),
                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = 'none') +
          geom_vline(xintercept =  0, linetype = 'dashed')+ 
          coord_cartesian(ylim=c(0,1000), xlim = c(-1,0.3)) +
          scale_y_continuous(expand = c(0,0)) +
          annotate(geom="text", fontface = 'bold', size = 3,x=-0.9, y=940, label=as.character(sum(mutTypeHome==2 & mutTypeResp == 2)),color="black")
        
        
      }
      
      ggsave(pltListHM[[no]], file = paste0('JDFE_', no, '.pdf') )
      
      
      
      
    }
    
  }
}


## make all the plots have individual names and save those names to a vector
pltNameVec =c()
for(i in 1:16)
{
  
  assign(paste("PLT_", i, sep = ""), pltListHM[[i]])
  
  pltNameVec = c(pltNameVec, paste("PLT_", i, sep = ""))
  
}



figure4 <- ggarrange(PLT_1, PLT_2,PLT_3, PLT_4,PLT_5, PLT_6,PLT_7, PLT_8,PLT_9, PLT_10,PLT_11, PLT_12,PLT_13, PLT_14,PLT_15, PLT_16, 
                     ncol = 4, nrow = 4, widths = c(1.3,1,1,1), heights = c(1,1,1,1.18))

print(figure4)
ggsave(figure4, file = 'compiledPlot.pdf', height = 11, width = 11)






#### save data frame of all JDFEs stats as CSV for supplemental file ####
allABRJDFEStats = data.frame(homeNameVec, respNameVec, homeMeanVec, respMeanVec,
                             homeVarVec, respVarVec, coVarVec, r1Vec, r2Vec, D11Vec,
                             D22Vec, D12Vec, cVec)

colnames(allABRJDFEStats) = c('Home', 'Response', 'Home_Mean', 'Response_Mean',
                              'Home_Var', 'Resp_Var', 'Covariance', 'r1', 'r2',
                              'D11', 'D22', 'D12', 'c')

write.csv(allABRJDFEStats, "allABRJDFEStats.csv")




allABRJDFEStats_sigVals = data.frame(homeNameVec, respNameVec, homeMeanVec, respMeanVec,
                                     homeVarVec, respVarVec, coVarVec, r1Vec_sig, r2Vec_sig, D11Vec_sig,
                                     D22Vec_sig, D12Vec_sig, cVec_sig)

colnames(allABRJDFEStats_sigVals) = c('Home', 'Response', 'Home_Mean', 'Response_Mean',
                                      'Home_Var', 'Resp_Var', 'Covariance', 'r1', 'r2',
                                      'D11', 'D22', 'D12', 'c')

write.csv(allABRJDFEStats_sigVals, "allABRJDFEStats_sigVals.csv")








##########################################################################################################################################################################################################
##########################################################################################################################################################################################################
##########################################################################################################################################################################################################




######  F I G U R E     5     P L E I O T R O P Y     P A R A M E T E R S #######

# NOTE: THIS FIGURE RELIES ON RESULTS GENERATED FROM FIGURE 4 CODE, YOU MUST RUN FIGURE 4 CODE BEFORE THIS SECTION WILL WORK 



########################################
#### ALL VALUES  ####


# r2 rank ordered
allABRJDFEStats$Response <- with(allABRJDFEStats,factor(Response,levels = rev(sort(unique(Response)))))

# just have colors be 1 for each rank (16 total colors, with 6 blues and 10 oranges)
limit = c(1, max(rank(allABRJDFEStats$r2)))
r2HeatMap = (ggplot(allABRJDFEStats, aes(x = Home,y=Response))+ 
               geom_tile(aes(fill = as.numeric(rank(r2))))+
               xlab("Home") +ylab("Non-home")+ 
               theme_bw()+ 
               geom_text( parse = TRUE, aes(label = gsub("e", " %*% 10^", scales::scientific_format()(r2)) ))+
               theme(panel.border = element_rect(color = 'black', fill = NA, size = 1), legend.position = 'bottom', legend.direction = 'horizontal', legend.text = element_text(size=15), legend.title = element_blank(),axis.text = element_text(size = 20), axis.title.x = element_text(size = 25), axis.title.y = element_text(size = 25))+
               scale_fill_gradientn(colors = c('#14a6ff','#3bb5ff', '#62c4ff','#89d3ff','#b1e1ff','#d8f0ff' ,'#fff3eb','#ffdac4','#ffc29d','#ffb689', '#ffa976','#ff9d62','#ff914e','#ff853b','#ff7827','#ff6c14'), limit = limit)) 

print(r2HeatMap)

ggsave(r2HeatMap, file = 'r2HeatMap_AllValues_rankOrdered.pdf', width = 5.5, height = 6.28)




## c rank ordered
#allABRJDFEStats$Response <- with(allABRJDFEStats,factor(Response,levels = rev(sort(unique(Response)))))

limit <- c(1,16)
cHeatMap = (ggplot(allABRJDFEStats, aes(x = Home,y=Response))+ geom_tile(aes(fill = as.numeric(rank(abs(c)))))+
              xlab("Home") +ylab("Non-home")+ 
              theme_bw()+ geom_text( parse = TRUE, aes(label = gsub("e", " %*% 10^", scales::scientific_format()(r2)) ))+
              theme(panel.border = element_rect(color = 'black', fill = NA, size = 1), legend.position = 'bottom', legend.direction = 'horizontal', legend.text = element_text(size=15),legend.title = element_blank() ,axis.text = element_text(size = 20), axis.title.x = element_text(size = 25), axis.title.y = element_blank())+
              scale_fill_gradientn(colors = c( '#ffffff','#fff3eb','#ffe7d8','#ffdac4','#ffceb1','#ffc29d','#ffb689', '#ffa976','#ff9d62','#ff914e','#ff853b','#ff7827','#ff6c14'), limit = limit))


print(cHeatMap)

ggsave(cHeatMap, file = 'cHeatMap_AllValues_rankOrdered.pdf', width = 5.5, height = 6.28)





## make the absr2vsqrtD22 scatterplot


allABRJDFEStats_sub = subset(allABRJDFEStats, Home!=Response)
## solve linear regression

group1_lm = summary(lm(sqrt(allABRJDFEStats_sub$D22) ~ 0+ abs(allABRJDFEStats_sub$r2)))
group1_slope = group1_lm$coefficients[1]


## plot scatter plot + the regressions

ymax = max(sqrt(allABRJDFEStats_sub$D22))

absr2VsqrtD22 = ggplot(allABRJDFEStats_sub, aes(x = abs(r2), y = sqrt(D22), color = as.factor(sign(r2))), group = group) + 
  geom_abline(intercept = 0, slope = group1_slope, color = '#282828', size = 1.3) +
  geom_point(size = 5) + ylim(0,ymax)+
  scale_x_continuous(label = scientific_10)+
  theme_bw()+ scale_color_manual(values = c('blue', 'orange'))+
  theme(panel.border = element_rect(color = 'black', fill = NA, size = 1.5), panel.grid.major = element_blank(),
        legend.position = 'none',axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) 





print(absr2VsqrtD22)

ggsave(absr2VsqrtD22, file = 'absr2vsqrtD22_AllVals.pdf', width = 5.5, height = 5.33)







##########  EXTRAS - ALL VALS######


## D12 (supp fig)

allABRJDFEStats$Response <- with(allABRJDFEStats,factor(Response,levels = rev(sort(unique(Response)))))

limit <- max(abs(allABRJDFEStats$D12)) * c(-1, 1)
D12HeatMap = (ggplot(allABRJDFEStats, aes(x = Home,y=Response))+ geom_tile(aes(fill = D12))+
                xlab("Home") +ylab("Non-Home")+ 
                theme_bw()+
                theme(panel.border = element_rect(color = 'black', fill = NA, size = 1), legend.position = 'bottom', legend.direction = 'horizontal', legend.text = element_text(size=15), axis.text = element_text(size = 20), axis.title.x = element_text(size = 25), axis.title.y = element_text(size = 25))+
                scale_fill_gradientn(colors = c('#14a6ff','#27a9ff','#3bb5ff', '#4ebdff','#62c4ff','#76cbff','#89d3ff', '#9ddaff','#b1e1ff','#c4e9ff','#d8f0ff' ,'#ebf8ff', '#ffffff','#fff3eb','#ffe7d8','#ffdac4','#ffceb1','#ffc29d','#ffb689', '#ffa976','#ff9d62','#ff914e','#ff853b','#ff7827','#ff6c14'), limit = limit))

print(D12HeatMap)

#ggsave(D12HeatMap, file = 'D12HeatMap.pdf')


## D22 (supp fig)

allABRJDFEStats$Response <- with(allABRJDFEStats,factor(Response,levels = rev(sort(unique(Response)))))

limit <-c(min(abs(allABRJDFEStats$D22)), max(abs(allABRJDFEStats$D22)) )
D22HeatMap = (ggplot(allABRJDFEStats, aes(x = Home,y=Response))+ geom_tile(aes(fill = D22))+
                xlab("Home") +ylab("Non-Home")+ 
                theme_bw()+
                theme(panel.border = element_rect(color = 'black', fill = NA, size = 1), legend.position = 'bottom', legend.direction = 'horizontal', legend.text = element_text(size=15), axis.text = element_text(size = 20), axis.title.x = element_text(size = 25), axis.title.y = element_text(size = 25))+
                scale_fill_gradientn(colors = c('#ffffff','#fff3eb','#ffe7d8','#ffdac4','#ffceb1','#ffc29d','#ffb689', '#ffa976','#ff9d62','#ff914e','#ff853b','#ff7827','#ff6c14'), limit = limit))

print(D22HeatMap)

#ggsave(D22HeatMap, file = 'D22HeatMap.pdf')





####  PLOT R2 VS D22, R2 VS D12, D12 VS D22 ,, sig values####

r2VD22 = ggplot(allABRJDFEStats, aes(x = r2, y = D22)) + geom_point(size = 3) +
  theme_bw()+ 
  theme(panel.border = element_rect(color = 'black', fill = NA, size = 1.5), panel.grid.major = element_blank(),
        legend.position = 'none',axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  geom_vline(xintercept =  0, linetype = 'dashed')



r2VD12 = ggplot(allABRJDFEStats, aes(x = r2, y = D12)) + geom_point(size = 3) +
  theme_bw()+ 
  theme(panel.border = element_rect(color = 'black', fill = NA, size = 1.5), panel.grid.major = element_blank(),
        legend.position = 'none',axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  geom_vline(xintercept =  0, linetype = 'dashed')



D12VD22 = ggplot(allABRJDFEStats, aes(x = D12, y = D22)) + geom_point(size = 3) +
  theme_bw()+ 
  theme(panel.border = element_rect(color = 'black', fill = NA, size = 1.5), panel.grid.major = element_blank(),
        legend.position = 'none',axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  geom_vline(xintercept =  0, linetype = 'dashed')


r2Vc = ggplot(allABRJDFEStats, aes(x = r2, y = c)) + geom_point(size = 3) +
  theme_bw()+ 
  theme(panel.border = element_rect(color = 'black', fill = NA, size = 1.5), panel.grid.major = element_blank(),
        legend.position = 'none',axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), axis.title.x = element_text(size = 25), axis.title.y = element_text(size = 25)) +
  geom_vline(xintercept =  0, linetype = 'dashed') + geom_hline(yintercept = 0, linetype = 'dashed')


csqVTto90 = ggplot(allABRJDFEStats, aes(x = c^2, y = timeTo90)) + geom_point(size = 3) +
  theme_bw()+ 
  theme(panel.border = element_rect(color = 'black', fill = NA, size = 1.5), panel.grid.major = element_blank(),
        legend.position = 'none',axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), axis.title.x = element_text(size = 25), axis.title.y = element_text(size = 25)) +
  geom_smooth(formula = y~x, method = 'lm', lwd = 0.5, color = 'black')



## for extreme time value, reduce 
df = subset(allABRJDFEStats, timeTo90 <=1000)

r2VtTo90 = ggplot(df, aes(x = r2, y = timeTo90)) + geom_point(size = 3) +
  theme_bw()+ 
  theme(panel.border = element_rect(color = 'black', fill = NA, size = 1.5), panel.grid.major = element_blank(),
        legend.position = 'none',axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  geom_vline(xintercept =  0, linetype = 'dashed')




#ggsave(r2Vc, file = 'r2vc.pdf')
#ggsave(r2VD22, file = 'r2VD22.pdf')
#ggsave(r2VD12, file = 'r2VD12.pdf')
#ggsave(D12VD22, file = 'D12VD22.pdf')






#####################################################################
#### MAIN TEXT-  O N L Y   S I G    V A L U E S ####


# r2 rank ordered

allABRJDFEStats_sigVals$Response <- with(allABRJDFEStats,factor(Response,levels = rev(sort(unique(Response)))))


# just have colors be 1 for each rank (16 total colors, with 6 blues and 10 oranges)

## to get geom text to accept the power notation, need to put in the guts of the scientific 10 function and let geom text parse it itself (the function was parsing it for me, which is fine for scale_y_continuous, for example, but not for geom_text)
limit = c(1, max(rank(allABRJDFEStats_sigVals$r2)))
r2HeatMap_sigValues = (ggplot(allABRJDFEStats_sigVals, aes(x = Home,y=Response))+ 
                         geom_tile(aes(fill = as.numeric(rank(r2))))+
                         xlab("Home") +ylab("Non-home")+ 
                         theme_bw()+
                         geom_text( parse = TRUE, aes(label = gsub("e", " %*% 10^", scales::scientific_format()(r2)) ))+
                         theme(panel.border = element_rect(color = 'black', fill = NA, size = 1), legend.position = 'bottom', legend.direction = 'horizontal', legend.text = element_text(size=15), legend.title = element_blank(),axis.text = element_text(size = 20), axis.title.x = element_text(size = 25), axis.title.y = element_text(size = 25))+
                         scale_fill_gradientn(colors = c('#14a6ff','#3bb5ff', '#62c4ff','#89d3ff','#b1e1ff','#d8f0ff' ,'#fff3eb','#ffdac4','#ffc29d','#ffb689', '#ffa976','#ff9d62','#ff914e','#ff853b','#ff7827','#ff6c14'), limit = limit)) 

print(r2HeatMap_sigValues)

ggsave(r2HeatMap_sigValues, file = 'r2HeatMap_sigValues_rankOrdered.pdf', width = 5.5, height = 6.28)





## c rank ordered
#allABRJDFEStats_sigVals$Response <- with(allABRJDFEStats,factor(Response,levels = rev(sort(unique(Response)))))

limit <- c(1,16)

cHeatMap_sigValues = (ggplot(allABRJDFEStats_sigVals, aes(x = Home,y=Response))+ 
                        geom_tile(aes(fill = as.numeric(rank(abs(c)))))+
                        xlab("Home") +ylab("Non-home")+ 
                        theme_bw()+
                        geom_text( parse = TRUE, aes(label = gsub("e", " %*% 10^", scales::scientific_format()(r2)) ))+
                        theme(panel.border = element_rect(color = 'black', fill = NA, size = 1), legend.position = 'bottom', legend.direction = 'horizontal', legend.text = element_text(size=15),legend.title = element_blank() ,axis.text = element_text(size = 20), axis.title.x = element_text(size = 25), axis.title.y = element_blank())+
                        scale_fill_gradientn(colors = c( '#ffffff','#fff3eb','#ffe7d8','#ffdac4','#ffceb1','#ffc29d','#ffb689', '#ffa976','#ff9d62','#ff914e','#ff853b','#ff7827','#ff6c14'), limit = limit))


print(cHeatMap_sigValues)

ggsave(cHeatMap_sigValues, file = 'cHeatMap_sigValues_rankOrdered.pdf', width = 5.5, height = 6.28)










## make the absr2vsqrtD22 scatterplot

allABRJDFEStats_sigVals_sub = subset(allABRJDFEStats_sigVals, Home!=Response)

## solve linear regressions for 3 groups and save parameters
group1 = subset(allABRJDFEStats_sigVals_sub, abs(r2)<0.005)

group1_lm = summary(lm(sqrt(group1$D22) ~ 0+ abs(group1$r2)))
group1_slope = group1_lm$coefficients[1]




group3 = subset(allABRJDFEStats_sigVals_sub, abs(r2)<0.016 & abs(r2) > 0.005)

group3_lm = summary(lm(sqrt(group3$D22) ~ 0+abs(group3$r2)))
group3_slope = group3_lm$coefficients[1]




## plot scatter plot + the regressions

scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}


ymax = max(sqrt(allABRJDFEStats_sigVals_sub$D22))
xmax = max(abs(allABRJDFEStats_sigVals_sub$r2))
absr2VsqrtD22 = ggplot(allABRJDFEStats_sigVals_sub, aes(x = abs(r2), y = sqrt(D22), color = as.factor(sign(r2))), group = group) + 
  geom_abline(intercept = 0, slope = group1_slope, color = '#282828', size = 1.3) +
  geom_abline(intercept = 0, slope = group3_slope, color = '#767676', size = 1.3) +
  geom_point(size = 5) + 
  ylim(0.001, ymax)+ 
  scale_x_continuous(label = scientific_10, limits = c(0.0005,xmax))+
  theme_bw()+ 
  scale_color_manual(values = c('blue', 'orange'))+
  theme(panel.border = element_rect(color = 'black', fill = NA, size = 1.5), panel.grid.major = element_blank(),
        legend.position = 'none',axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) 






print(absr2VsqrtD22)

ggsave(absr2VsqrtD22, file = 'absr2vsqrtD22_sigValues.pdf', width = 5.5, height = 5.33)





######### EXTRAS - SIG VALS ####


## D12 (supp fig)

allABRJDFEStats_sigVals$Response <- with(allABRJDFEStats,factor(Response,levels = rev(sort(unique(Response)))))

limit <- max(abs(allABRJDFEStats_sigVals$D12)) * c(-1, 1)
D12HeatMap_sigValues = (ggplot(allABRJDFEStats_sigVals, aes(x = Home,y=Response))+ geom_tile(aes(fill = D12))+
                          xlab("Home") +ylab("Non-Home")+ 
                          theme_bw()+
                          theme(panel.border = element_rect(color = 'black', fill = NA, size = 1), legend.position = 'bottom', legend.direction = 'horizontal', legend.text = element_text(size=15), axis.text = element_text(size = 20), axis.title.x = element_text(size = 25), axis.title.y = element_text(size = 25))+
                          scale_fill_gradientn(colors = c('#14a6ff','#27a9ff','#3bb5ff', '#4ebdff','#62c4ff','#76cbff','#89d3ff', '#9ddaff','#b1e1ff','#c4e9ff','#d8f0ff' ,'#ebf8ff', '#ffffff','#fff3eb','#ffe7d8','#ffdac4','#ffceb1','#ffc29d','#ffb689', '#ffa976','#ff9d62','#ff914e','#ff853b','#ff7827','#ff6c14'), limit = limit))

print(D12HeatMap_sigValues)

#ggsave(D12HeatMap_sigValues, file = 'D12HeatMap_sigValues.pdf')


## D22 (supp fig)

allABRJDFEStats_sigVals$Response <- with(allABRJDFEStats,factor(Response,levels = rev(sort(unique(Response)))))

limit <-c(min(abs(allABRJDFEStats_sigVals$D22)), max(abs(allABRJDFEStats_sigVals$D22)) )
D22HeatMap_sigValues = (ggplot(allABRJDFEStats_sigVals, aes(x = Home,y=Response))+ geom_tile(aes(fill = D22))+
                          xlab("Home") +ylab("Non-Home")+ 
                          theme_bw()+
                          theme(panel.border = element_rect(color = 'black', fill = NA, size = 1), legend.position = 'bottom', legend.direction = 'horizontal', legend.text = element_text(size=15), axis.text = element_text(size = 20), axis.title.x = element_text(size = 25), axis.title.y = element_text(size = 25))+
                          scale_fill_gradientn(colors = c('#ffffff','#fff3eb','#ffe7d8','#ffdac4','#ffceb1','#ffc29d','#ffb689', '#ffa976','#ff9d62','#ff914e','#ff853b','#ff7827','#ff6c14'), limit = limit))

print(D22HeatMap_sigValues)

#ggsave(D22HeatMap_sigValues, file = 'D22HeatMap_sigValues.pdf')





####  PLOT R2 VS D22, R2 VS D12, D12 VS D22 ,, sig values####

absr2VsqrtD22 = ggplot(allABRJDFEStats_sigVals, aes(x = abs(r2), y = sqrt(D22), color = as.factor(sign(r2)))) + geom_point(size = 3) +
  theme_bw()+ scale_color_manual(values = c('blue', 'orange'))+
  theme(panel.border = element_rect(color = 'black', fill = NA, size = 1.5), panel.grid.major = element_blank(),
        legend.position = 'none',axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  geom_vline(xintercept =  0, linetype = 'dashed')



r2VD12 = ggplot(allABRJDFEStats_sigVals, aes(x = r2, y = D12)) + geom_point(size = 3) +
  theme_bw()+ 
  theme(panel.border = element_rect(color = 'black', fill = NA, size = 1.5), panel.grid.major = element_blank(),
        legend.position = 'none',axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  geom_vline(xintercept =  0, linetype = 'dashed')



D12VD22 = ggplot(allABRJDFEStats_sigVals, aes(x = D12, y = D22)) + geom_point(size = 3) +
  theme_bw()+ 
  theme(panel.border = element_rect(color = 'black', fill = NA, size = 1.5), panel.grid.major = element_blank(),
        legend.position = 'none',axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  geom_vline(xintercept =  0, linetype = 'dashed')


r2Vc = ggplot(allABRJDFEStats_sigVals, aes(x = r2, y = c, col = as.factor(sign(r2)))) + geom_point(size = 3) +
  theme_bw()+ scale_color_manual(values = c('blue', 'orange'))+
  theme(panel.border = element_rect(color = 'black', fill = NA, size = 1.5), panel.grid.major = element_blank(),
        legend.position = 'none',axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), axis.title.x = element_text(size = 25), axis.title.y = element_text(size = 25)) +
  geom_vline(xintercept =  0, linetype = 'dashed') + geom_hline(yintercept = 0, linetype = 'dashed')


csqVTto90 = ggplot(allABRJDFEStats_sigVals, aes(x = c^2, y = timeTo90)) + geom_point(size = 3) +
  theme_bw()+ 
  theme(panel.border = element_rect(color = 'black', fill = NA, size = 1.5), panel.grid.major = element_blank(),
        legend.position = 'none',axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), axis.title.x = element_text(size = 25), axis.title.y = element_text(size = 25)) +
  geom_smooth(formula = y~x, method = 'lm', lwd = 0.5, color = 'black')



## for extreme time value, reduce 
df = subset(allABRJDFEStats_sigVals, timeTo90 <=1000)

r2VtTo90 = ggplot(df, aes(x = r2, y = timeTo90)) + geom_point(size = 3) +
  theme_bw()+ 
  theme(panel.border = element_rect(color = 'black', fill = NA, size = 1.5), panel.grid.major = element_blank(),
        legend.position = 'none',axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  geom_vline(xintercept =  0, linetype = 'dashed')




#ggsave(r2Vc, file = 'r2vc_sigVals.pdf', width = 5.5, height = 5.33)
#ggsave(absr2VsqrtD22, file = 'absr2VsqrtD22_sigVals.pdf', width = 5.5, height = 5.33)
#ggsave(r2VD12, file = 'r2VD12_sigVals.pdf', width = 5.5, height = 5.33)
#ggsave(D12VD22, file = 'D12VD22_sigVals.pdf', width = 5.5, height = 5.33)







################### SUPP - ANTIBIOTIC RESISTANCE DATA SIMULATION ###############
## this section can run independently, so long as working directory has the correct files below to import

#setwd("~/Documents/Antibiotic Resistance JDFE Project/ Figures/SI fig ABR params vs sims")

##### Import Data #####
WTReplicateData = read.csv('WTReplicate.csv', header = 1, sep = ":")
# remove excess rows; keep only  No Drug and 1st 6 drug enviros and remove header
WTReplicateData = as.matrix(WTReplicateData[,2:8])


# Import KO growth rate data 
KOGrowthRates = read.csv('KOGrDATAPrtoStop.csv', header = 1, sep = ",")
# remove KOID and FOX and AMP (last 2 drugs), as dont have WT replicate data for them
KOPrtoStop = KOGrowthRates[,11]
KOGrowthRates = as.matrix(KOGrowthRates[,2:8])



# Important Variables #
WTGrowthRateMeans = colMeans(WTReplicateData)
WTGrowthRateSDs   = colSds(WTReplicateData) 
numKnockouts = dim(KOGrowthRates)[1]
numDrugs = dim(KOGrowthRates)[2]





#### Get Svalue of all Muts from KOGrowthRates ####

# first, eliminate all 0 values from the KOGrowthRate matrix and record probability of
# KILLER mut (should be same for all drug environments)
toRemove = c()
count = 0
for( i in 1:dim(KOGrowthRates)[1])
{
  if(KOGrowthRates[i,1] == 0)
  {
    count = count + 1
    toRemove = c(toRemove, i)
  }
}

probKiller = count/dim(KOGrowthRates)[1]
KOGrowthRates = KOGrowthRates[-toRemove, ]


# turn all values from growth rate to selection coefficients
sValAllMuts = array(0, dim(KOGrowthRates))
for( i in 1:dim(KOGrowthRates)[1])
{
  sValAllMuts[i, ] = KOGrowthRates[i,] - WTGrowthRateMeans
}



## data showing weibull is better fit for majority than normal or exponential
weibullLogLik = c()
gaussLogLik = c()
expLogLik = c()

for(drug in 1:7)
{
  weibullLogLik = c(weibullLogLik, fitdist(sValAllMuts[,drug]+1, distr = 'weibull')$loglik)
  gaussLogLik = c(gaussLogLik, fitdist(sValAllMuts[,drug]+1, distr = 'norm')$loglik)
  expLogLik = c(expLogLik, fitdist(sValAllMuts[,drug]+1, distr = 'exp')$loglik)
}

names = c('No Drug', 'CHL', 'CPR', "MEC", 'NIT', 'TET', 'TMP')
dfeFitdata = data.frame(names, weibullLogLik, gaussLogLik, expLogLik)

write.csv(dfeFitdata, file = 'DFE_DistFit_Data.csv')

##### Create Binned Conditional Probability Matricies for ALL drug pairs ######

#### create list for each resp enviro with the sVals in RESP for each bin in HOME ####

#### this is a 3 layered list with 1 layer being HOME enviro
# 2nd being RESP enviro
# 3rd being BIN
# each entry at this third level is a vector of all the Svalues in a given RESP envuro
# which correspond to a HOME svalue in the given BIN
numDrugs = 7
binSize = 0.1
minS = min(sValAllMuts)
maxS = max(sValAllMuts)
numBins = length(seq(floor(minS),ceiling(maxS*10)/10, by = binSize)) - 1
binSeq = seq(floor(minS),ceiling(maxS*10)/10, by = binSize) # length = 1+numBins
lastNegBin = which.min(abs(binSeq)) - 1


allHomeSValLists= list()
allHomeWeibullParameters = list()
for(home in 1:numDrugs)
{
  sVallistALLRespDrugs = list()
  weibullParametersALLRespDrugs = list()
  for(resp in 1:numDrugs)
  {
    
    singleDrugList = list() # 1 entry per bin
    singleDrugWeibullParameters = list()
    for(bin in 1:numBins)
    {
      binMIN = binSeq[bin]
      binMAX = binSeq[bin+1]
      
      
      sValsinBin = (sValAllMuts[,home]>=binMIN & sValAllMuts[,home]<=binMAX)
      
      respSvals = sValAllMuts[,resp][sValsinBin!=0]
      
      singleDrugList[[bin]] = respSvals
      
      
      #####  NOW: SAVE THE WEIBULL SHAPE AND SCALE PARAMETERS FOR ALL FITS YOU CAN ####
      
      if(length(respSvals) >= 50 & var(respSvals >0))
      {
        dat <- data.frame(x=c(respSvals + 1)) ## +1 bc some dists cant do negatives
        fitW <- fitdistr(dat$x, densfun = "weibull")
        
        singleDrugWeibullParameters[[bin]] = fitW$estimate
        
      }else{
        singleDrugWeibullParameters[[bin]] = 0 
      }
      
      
    }
    
    sVallistALLRespDrugs[[resp]] = singleDrugList
    weibullParametersALLRespDrugs[[resp]] = singleDrugWeibullParameters
  }
  allHomeSValLists[[home]] = sVallistALLRespDrugs
  allHomeWeibullParameters[[home]] = weibullParametersALLRespDrugs
}



#==============================================================================================================================================#
##### Find cutoffs for calling mutation neutral, deleterious or beneficial for each drug #####

## find the right cutoff values for calling beneficial muts so have FDR ~ 25% for beneficial muts in  all drugs
## same code to find be/delt cutoff as ABRJDFE_HEatmap JDFE
testCutoffs = seq(0.01, 0.5, length.out = 100)

nameVec = c()
fdrVec = c()
cutoffVec = c()

for(home in 1:7)
{
  
  fdrsThisHome = c()
  cutoffsthishome = c()
  for(cutoff in testCutoffs)
  {
    
    evo = home
    if(evo == 1){
      title = "No Drug"
    }else if(evo == 2){
      title = "CHL"
    }else if(evo == 3){
      title = "CPR"
    }else if(evo ==4){
      title = "MEC"
    }else if(evo == 5){
      title = "NIT"
    }else if(evo == 6){
      title = "TET"
    }else{
      title = "TMP"
    }
    
    
    benCutoff = qnorm((1-cutoff), mean = 0, sd = WTGrowthRateSDs[home]) 
    deltCutoff = qnorm(cutoff, mean = 0, sd = WTGrowthRateSDs[home])
    
    DFE = sValAllMuts[,home]
    errorDist = rnorm(400, mean = 0, sd = WTGrowthRateSDs[home])
    breaks = seq(min(c(DFE, errorDist)), max(c(DFE, errorDist)), length.out = 50)
    
    
    mutTypeVec = (sValAllMuts[,(home)]>benCutoff)*2 +  (sValAllMuts[,(home)]<deltCutoff)
    
    beneficialMutCount = sum(mutTypeVec == 2)
    
    expectedNumFalsePositive = sum(sValAllMuts[,home]>0)*cutoff*2 ## 
    
    fdr = expectedNumFalsePositive/beneficialMutCount ## want to be < 1 ,, ie want to be calling more than would be expected just by noise 
    
    fdrsThisHome = c(fdrsThisHome, fdr)
    cutoffsthishome = c(cutoffsthishome, cutoff)
  }
  
  cutofftouse = which.min(abs(fdrsThisHome-0.25))
  
  
  nameVec = c(nameVec, title)
  cutoffVec = c(cutoffVec, cutoffsthishome[cutofftouse])
  fdrVec = c(fdrVec, fdrsThisHome[cutofftouse])
}

allFDR_Results = data.frame(nameVec, cutoffVec, fdrVec)








numDrugs = 7
beneficialCutoffVec = numeric(numDrugs)
deleteriousCutoffVec = numeric(numDrugs)

cutoffVec = allFDR_Results$cutoffVec#the values of cutoff so have lowest false disocvery rate in all 

for( i in 1:numDrugs)
{
  cutoff = cutoffVec[i]
  ### mean = 0 bc svals are calc as subtract wt gr
  beneficialCutoffVec[i] = qnorm((1-cutoff), mean = 0, sd = WTGrowthRateSDs[i]) 
  deleteriousCutoffVec[i] = qnorm((cutoff), mean = 0, sd = WTGrowthRateSDs[i]) 
  
}





mutTypeDf = data.frame(KOGrowthRates)
beneficialMutCount = numeric(7)
deleteriousMutCount = numeric(7)
totalPossibleBenMuts = numeric(7)
for(home in 1:7)
{
  
  benCutoff = beneficialCutoffVec[(home)]
  deltCutoff = deleteriousCutoffVec[(home)]
  
  ##home -1 bc in svalallmuts dont have 1st col of mutID
  mutTypeDf[,home] = (sValAllMuts[,(home)]>benCutoff)*2 +  (sValAllMuts[,(home)]<deltCutoff)
  
  
  ## also count Num Muts which are beneficial in each env and deleterious in each env
  beneficialMutCount[(home)] = sum(mutTypeDf[,(home)] == 2)
  deleteriousMutCount[(home)] =sum(mutTypeDf[,(home)] == 1)
  
  totalPossibleBenMuts[(home)] = sum(sValAllMuts[,(home)]>0)
  
}





allWeibullParameters = list()
for( i in 1:numDrugs)
{
  fit = fitdistr(sValAllMuts[,i] + 1,"weibull")
  shape = fit$estimate[1]
  scale = fit$estimate[2]
  
  
  allWeibullParameters[[i]] = c(shape,scale)
}





### calc R and D values ###
r1Easy= array(0, c(7,7))
r2Easy = array(0, c(7,7))
D11Easy = array(0, c(7,7))
D22Easy = array(0, c(7,7))
D12Easy = array(0, c(7,7))
cEasy = array(0, c(7,7))
timeTo90pSameasR2Pred = array(0, c(7,7))

# standard error for all measures
r1EasySE= array(0, c(7,7))
r2EasySE = array(0, c(7,7))
D11EasySE = array(0, c(7,7))
D22EasySE = array(0, c(7,7))
D12EasySE = array(0, c(7,7))

for(home in 1:7)
{
  ## make all negative svals in home which are negative 0 so dont dont contribute to calc!
  homeSvals = sValAllMuts[,home]*(sValAllMuts[,home]>=0)
  
  for(resp in 1:7)
  {
    r1Easy[home,resp] = mean((homeSvals)^2)
    r2Easy[home,resp] = mean((homeSvals)*(sValAllMuts[,resp]))
    D11Easy[home,resp] = mean((homeSvals)^3)
    D12Easy[home,resp] = mean((homeSvals)^2*(sValAllMuts[,resp]))
    D22Easy[home,resp] = mean((homeSvals)*(sValAllMuts[,resp])^2)
    cEasy[home, resp] = sqrt(mean((homeSvals)*(sValAllMuts[,resp])^2))/(mean((homeSvals)*(sValAllMuts[,resp])))
    timeTo90pSameasR2Pred[home, resp] = (1.28^2*mean((homeSvals)*(sValAllMuts[,resp])^2))/mean((homeSvals)*(sValAllMuts[,resp]))^2
    
    r1EasySE[home,resp] = sd((homeSvals)^2)/sqrt(numKnockouts)
    r2EasySE[home,resp] = sd((homeSvals)*(sValAllMuts[,resp]))/sqrt(numKnockouts)
    D11EasySE[home,resp] = sd((homeSvals)^3)/sqrt(numKnockouts)
    D12EasySE[home,resp] = sd((homeSvals)^2*(sValAllMuts[,resp]))/sqrt(numKnockouts)
    D22EasySE[home,resp] = sd((homeSvals)*(sValAllMuts[,resp])^2)/sqrt(numKnockouts)
  }
}




#=====================================================================================================================#
#========================================================================================================#
#========================================================================================================#

#### RUN SIMULATION ####
popsize = 10^4
U = 10^-4     # mutation rate
numGens = 1000 # per evolution simulation
numIterations = 5 # number of times to simulate evo in the given environment 
numDrugs = 7
slopeArray = array(data = 0, dim = c(numDrugs, numDrugs))
errorBarSize = array(data = 0, dim = c(numDrugs, numDrugs))
covarArray = array(data = 0, dim = c(numDrugs, numDrugs))
#corrArray = array(data = 0, dim = c(numDrugs, numDrugs))
varArray_noLog = array(data = 0, dim = c(numDrugs, numDrugs))
varArray = array(data = 0, dim = c(numDrugs, numDrugs))
errorBarsAllDrugs = array(data = 0, dim = c(1, 6))

numBenMutsAllDrugs = numeric(numDrugs)
benMutEffects = numeric(numDrugs)
allUbs = array() 


allCovars = array(data = 0, dim = c(numGens, numDrugs, numDrugs))
allCors = array(data = 0, dim = c(numGens, numDrugs, numDrugs))
for( yay in c(3,5,6,4))
{
  enviroEvolveID = yay
  
  
  allWbarsAllIterations = array(data = 0, dim = c(numGens, (numDrugs), numIterations))
  fullDrugVariDataSet = array(data = 0, dim = c(numGens, numDrugs, numIterations))
  
  for(w in 1:numIterations)
  { 
    print('drug num')
    print(yay)
    print('iteration number')
    print(w)
    
    startIDarray = array(data = c(1, popsize, WTGrowthRateMeans[yay]), dim = c(1, 3))
    mutFitnesses = array(data = WTGrowthRateMeans, dim = c(1, numDrugs))
    allWbars = array(data = 0, dim = c(numGens, (numDrugs)))
    itVar = array(data = 0, dim = c(numGens, numDrugs))
    for( x in 1:numGens)
    {
      numMuts = dim(startIDarray)[1]
      
      #### Selection and Drift ####
      
      Wbar = sum(startIDarray[ , 3] * (startIDarray[ ,2] / popsize))
      
      # generate probability vector
      probvec = numeric(numMuts)
      for( i in 1:numMuts)
      {
        
        if((startIDarray[i,2] / popsize) + (startIDarray[i,2] / popsize)* (startIDarray[i,3] - Wbar) >0)
        {
          ## SK way new below
          probvec[i] = (startIDarray[i,2] / popsize) + (startIDarray[i,2] / popsize)* (startIDarray[i,3] - Wbar)
        }else{
          probvec[i] = 0
        }
      }
      
      # generate IDarray (with only sel/drift) for next gen
      nextGenIDarray = array(data = 0, dim = c(dim(startIDarray)[1], dim(startIDarray)[2]))
      nextGenOffspring = rmultinom(1, size = popsize, prob = probvec)
      
      
      for( i in 1:numMuts)
      {
        nextGenIDarray[i, 1] = startIDarray[i, 1]
        nextGenIDarray[i, 2] = nextGenOffspring[i]
        nextGenIDarray[i, 3] = startIDarray[i, 3]
      }
      
      
      
      ## Mutations ##
      M = rpois(1, popsize*U)
      parentIDs = sample(nextGenIDarray[ , 1], M, prob = (nextGenIDarray[ ,2]/popsize), replace = TRUE)
      
      for( i in parentIDs)
      {
        # pick s in home
        homeS = rweibull(1, shape = allWeibullParameters[[enviroEvolveID]][1], scale = allWeibullParameters[[enviroEvolveID]][2]) -1
        
        if(homeS > beneficialCutoffVec[enviroEvolveID])
        {
          numBenMutsAllDrugs[enviroEvolveID] = numBenMutsAllDrugs[enviroEvolveID] + 1
          benMutEffects[enviroEvolveID] = benMutEffects[enviroEvolveID] + homeS
        }
        
        ### add in chance of killer muatant
        killChance = sample(c(0,1), 1, prob = c(probKiller, (1-probKiller))) # 0 if get killer, 1 otherwise
        if(killChance == 0)
        {
          svec = numeric(numDrugs) - 1
        }else{
          
          # find bin of homeS
          if(homeS >= -1 & homeS <= -0.9)
          {
            binHome = 1
          }else if(homeS >= -0.9 & homeS < -0.8){
            binHome = 2
          }else if(homeS >= -0.8 & homeS < -0.7){
            binHome = 3
          }else if(homeS >= -0.7 & homeS < -0.6){
            binHome = 4
          }else if(homeS >= -0.6 & homeS < -0.5){
            binHome = 5
          }else if(homeS >= -0.5 & homeS < -0.4){
            binHome = 6
          }else if(homeS >= -0.4 & homeS < -0.3){
            binHome = 7
          }else if(homeS >= -0.3 & homeS < -0.2){
            binHome = 8
          }else if(homeS >= -0.2 & homeS < -0.1){
            binHome = 9
          }else if(homeS >= -0.1 & homeS < 0){
            binHome = 10
          }else if(homeS >= 0 & homeS < 0.1){
            binHome = 11
          }else if(homeS >= 0.1 & homeS < 0.2){
            binHome = 12
          }else if(homeS >= 0.2 & homeS < 0.3){
            binHome = 13
          }
          
          #### make svec####
          
          svec = numeric(numDrugs)
          svec[enviroEvolveID] = homeS
          for(respit in 1:numDrugs)
          {
            if(respit != enviroEvolveID) #only need to sample new value from NON homes
            {
              if(length(allHomeWeibullParameters[[enviroEvolveID]][[respit]][[binHome]]) == 2) # only fit dist to sets of svalues with more than 10 vals
              {
                
                svec[respit] = rweibull(1, shape = allHomeWeibullParameters[[enviroEvolveID]][[respit]][[binHome]][1], scale = allHomeWeibullParameters[[enviroEvolveID]][[respit]][[binHome]][2]) -1
                
                
              }else{
                
                ### WHAT TO DO IF CANT FIT THE DISTRIBUTION ####
                
                if(length(allHomeSValLists[[enviroEvolveID]][[respit]][[binHome]]) != 0)
                {
                  svec[respit] = rnorm(1, mean = mean(allHomeSValLists[[enviroEvolveID]][[respit]][[binHome]]), sd = WTGrowthRateSDs[[respit]])
                }else{
                  svec[respit] = 0
                  
                }
                
                
              }
              
            }
          }
          
          
          ### adjust for neutral mutations
          for(sVal in 1:numDrugs)
          {
            if(svec[sVal] < beneficialCutoffVec[sVal] & svec[sVal]>deleteriousCutoffVec[sVal])
            {
              svec[sVal] = 0
            }
          }
          
          
        }
        
        
        
        ####################
        
        #remove indiv from its parent class if can be removed
        if(nextGenIDarray[i,2] != 0)
        {
          nextGenIDarray[i, 2] = nextGenIDarray[i, 2]- 1
          
          #make new row for this new mutant; mutliplicative fitness change (mult parent fit by (1+s), so dec if s neg or inc if s pos)
          newRow = c(nextGenIDarray[dim(nextGenIDarray)[1], 1]+ 1, 1, startIDarray[i,3] + svec[enviroEvolveID])
          
          #bind newrow with nextgenIDarray
          nextGenIDarray = rbind(nextGenIDarray, newRow)
          
          #add row with fitness of new mutant in all environments to mutFitnesses
          newFit = mutFitnesses[i , ]+ svec
          mutFitnesses = rbind(mutFitnesses, newFit)
        }
        
        
      }
      
      output = list(nextGenIDarray, mutFitnesses)
      
      
      #====================================================================================================================#
      
      startIDarray = output[[1]]
      mutFitnesses = output[[2]]
      
      
      # find and save Wbars in all enviros 
      
      for( n in 1:(numDrugs))
      {
        allWbars[x, n] = sum(mutFitnesses[ , n] * (startIDarray[ ,2] / popsize))
      }
      
      # remove all muttypes with frequency of 0 from idarray and mutfitarray
      toRemove = c()
      if(length(startIDarray) > 3) # only have chance to remove if have more than 1 line in the startIDarray
      {
        for(q in 1:(dim(startIDarray)[1]))
        {
          if(startIDarray[q,2] == 0)
          {
            toRemove = c(toRemove, q)
            
          }
        }
      }
      
      
      if(length(toRemove) > 0)
      {
        startIDarray = startIDarray[-toRemove, ]
        mutFitnesses = mutFitnesses[ -toRemove, ]
        
        
        
        if(length(startIDarray) == 3) 
        {
          startIDarray = array(data = c(startIDarray), dim = c(1,3))
          mutFitnesses = array( data = c(mutFitnesses), dim= c(1, numDrugs))
        }
        
        count = 0
        
        for(q in 1:(dim(startIDarray)[1]))
        {
          count = count + 1
          startIDarray[q, 1] = count
        }
        
      }
      
      
    }
    
    allWbarsAllIterations[ , , w] = allWbars
  }
  
  # variance across all iterations for this evolve enviro in this data set
  avgDrugVariDataSet = array(data = 0, dim = c(numGens, numDrugs))
  for( i in 1:numGens)
  {
    for(y in 1:(numDrugs))
    {
      avgDrugVariDataSet[i,y] = var((allWbarsAllIterations[i,y, ])) #
    }
  }
  
  avgDrugVariDataSet_noLog = array(data = 0, dim = c(numGens, numDrugs))
  for( i in 1:numGens)
  {
    for(y in 1:(numDrugs))
    {
      avgDrugVariDataSet_noLog[i,y] = var((allWbarsAllIterations[i,y, ]))
    }
  }
  
  #COVARIANCE
  avgDrugCoVariDataSet = array(data = 0, dim = c(numGens, numDrugs))
  for( i in 1:numGens)
  {
    for(y in 1:(numDrugs))
    {
      avgDrugCoVariDataSet[i,y] = cov((allWbarsAllIterations[i,y, ]), (allWbarsAllIterations[i,enviroEvolveID, ])) ## all analysis occurs in log transformed data, so take variance of the log trans data
    }
  }
  
  allCovars[,,yay] = avgDrugCoVariDataSet
  
  # CORREALTIONS 
  avgDrugCoriDataSet = array(data = 0, dim = c(numGens, numDrugs))
  for( i in 1:numGens)
  {
    for(y in 1:(numDrugs))
    {
      avgDrugCoriDataSet[i,y] = cor((allWbarsAllIterations[i,y, ]), (allWbarsAllIterations[i,enviroEvolveID, ])) ## all analysis occurs in log transformed data, so take variance of the log trans data
    }
  }
  
  
  
  
  # average Wbar across all iterations for this evolve enviro in this data set 
  averageAllWbar = array(data = 0, dim = c(numGens, (numDrugs)))
  for( i in 1:numGens)
  {
    for(y in 1:(numDrugs))
    {
      averageAllWbar[i,y] = mean(allWbarsAllIterations[i,y, ])
    }
  }
  
  
  
  
  
  
  ##### PLot results for this EnviroEvolveID ######
  WbarData = averageAllWbar
  varData = avgDrugVariDataSet_noLog
  covarData = avgDrugCoVariDataSet
  corData = avgDrugCoriDataSet
  evo = enviroEvolveID
  
  
  
  
  genVec = c()
  nameVec = c()
  fitnessVec = c()
  sdVec = c()
  
  for( resp in 1:numDrugs)
  {
    
    # make drug column have drugs actual names
    
    if(resp == 1){
      name = "No Drug"
    }else if(resp == 2){
      name = "CHL"
    }else if(resp == 3){
      name = "CPR"
    }else if(resp ==4){
      name = "MEC"
    }else if(resp == 5){
      name = "NIT"
    }else if(resp == 6){
      name = "TET"
    }else if(resp == 7){
      name = "TMP"
    }
    
    
    
    
    genVec = c(genVec, 1:numGens)
    nameVec = c(nameVec,rep(name, times = numGens) )
    fitnessVec = c(fitnessVec, (WbarData[,resp]))
    sdVec = c(sdVec, sqrt(varData[,resp]))
    
    
    
    
  }
  
  minVec = fitnessVec - sdVec
  maxVec= fitnessVec + sdVec
  
  
  
  
  # create data frame from the matrix
  WbarDf = data.frame(genVec, nameVec, fitnessVec, minVec, maxVec)
  colnames(WbarDf) = c("time", "drug", "fitness", 'min', 'max')
  
  
  
  # only use a fewe points to make the graphs so dont have obscene number error bars
  WbarDf = WbarDf[seq(from = 1, to = numDrugs*numGens, by = 99), ]
  
  
  
  # then plot using ggplot specifications have created below 
  
  if(evo == 1){
    title = "No Drug"
  }else if(evo == 2){
    title = "CHL"
  }else if(evo == 3){
    title = "CPR"
  }else if(evo ==4){
    title = "MEC"
  }else if(evo == 5){
    title = "NIT"
  }else if(evo == 6){
    title = "TET"
  }else{
    title = "TMP"
  }
  
  
  
  # for error bar plot
  print( ggplot(WbarDf, aes(x=as.numeric(time),y= as.numeric(fitness), colour = drug )) + geom_line()  +
           geom_errorbar(aes(x= time, ymin=min, ymax=max), width=0.25) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                               panel.background = element_blank(), axis.line = element_line(colour = "black"))+
           ggtitle(as.character(title)) + xlab("time")+ ylab("fitness") )
  
  
  ## save average slope
  for(slope in 1:numDrugs)
  {                            ###NEED TO TAKE  OF FIT VALUES ##
    slopeArray[evo, slope] = ((averageAllWbar[numGens/2, slope]) - (averageAllWbar[numGens, slope]))/(numGens/2 - numGens)
  }
  
  
  
  ##save variance
  for(iEnv in 1:numDrugs)
  {
    varArray[evo, iEnv] = avgDrugVariDataSet[numGens, iEnv]
  }
  
  for(iEnv in 1:numDrugs)
  {
    varArray_noLog[evo, iEnv] = avgDrugVariDataSet_noLog[numGens, iEnv]
  }
  
  ##savecovariance
  for(iEnv in 1:numDrugs)
  {
    covarArray[evo, iEnv] = avgDrugCoVariDataSet[numGens, iEnv]
  }
  
  
  # evo drug                                # resp drug  #time     #  fitness          min  and max at various tp
  errorBarsthisDrug = cbind(rep(title, times = dim(WbarDf)[1]), WbarDf$drug, WbarDf$time , WbarDf$fitness , WbarDf$min, WbarDf$max)
  errorBarsAllDrugs = rbind(errorBarsAllDrugs, errorBarsthisDrug)
  
  
  
  
}


# correct for drugs which focus on 
varArray = varArray[c(3,4,5,6),c(3,4,5,6)]
covarArray = covarArray[c(3,4,5,6),c(3,4,5,6)]
slopeArray = slopeArray[c(3,4,5,6),c(3,4,5,6)]


r2Easy= r2Easy[c(3,4,5,6),c(3,4,5,6)]
D22Easy= D22Easy[c(3,4,5,6),c(3,4,5,6)]
D12Easy= D12Easy[c(3,4,5,6),c(3,4,5,6)]


scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}

N = popsize


UbVec = c()
for(i in c(3,4,5,6)) #home
{
  for(j in c(3,4,5,6)) #resp
  {
    UbVec = c(UbVec,U*( 1 - pnorm(0, WTGrowthRateMeans[i], WTGrowthRateSDs[i])))
  }
}

adjR2Vals = c(r2Easy)*2*N*UbVec
adjD22Vals = c(D22Easy)*2*N*UbVec
adjD12Vals = c(D12Easy)*2*N*UbVec
adjD22endTimeVals= c(D22Easy)*2*N*UbVec*numGens
adjD12endTimeVals =c(D12Easy)*2*N*UbVec*numGens




data_thisNu = data.frame(c(adjR2Vals), c(adjD22endTimeVals), c(adjD12endTimeVals), c(slopeArray), c(varArray), c(covarArray))
colnames(data_thisNu) = c('r2', 'D22', 'D12', 'meanRespSlopes', 'meanVarResp', 'meanCoVarResp')


color = 'hotpink4'
## MEAN SLOPE

assign(paste('r2vSlope_ABRSim_Nu', N*U,sep = ""), ggplot(data_thisNu, aes(x=r2,y=meanRespSlopes)) + geom_point(col = color) +
         geom_smooth(method='lm', formula= y~x, col = color, lwd = 0.5) + geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey')+
         theme_bw()+theme(panel.border = element_rect(colour = 'black', fill = NA, size = 1), 
                          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
                          axis.text.y =element_text(size = 35,family = 'Helvetica'),axis.text.x =element_text(size = 35,family = 'Helvetica', color = c('grey30', 'transparent','grey30', 'transparent','grey30', 'transparent')), axis.title = element_blank()) +
         scale_x_continuous(label = scientific_10) +  scale_y_continuous(label = scientific_10))

ggsave(r2vSlope_ABRSim_Nu1, file = 'r2vSlope_ABRSim_Nu1.pdf')

## SLOPE OF VARIANCE TRAJECTORY

assign(paste('D22vVarSlope_ABRSim_Nu', N*U,sep = ""), ggplot(data_thisNu, aes(x=D22,y=meanVarResp)) + geom_point(col = color) +
         geom_smooth(method='lm', formula= y~x, col = color, lwd = 0.5) + geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey')+
         theme_bw()+theme(panel.border = element_rect(colour = 'black', fill = NA, size = 1), 
                          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
                          axis.text.y =element_text(size = 35,family = 'Helvetica'),axis.text.x =element_text(size = 35,family = 'Helvetica', color = c('grey30', 'transparent','grey30', 'transparent','grey30', 'transparent')), axis.title = element_blank()) +
         scale_x_continuous(label = scientific_10) +  scale_y_continuous(label = scientific_10))


ggsave(D22vVarSlope_ABRSim_Nu1, file = 'D22vVarSlope_ABRSim_Nu1.pdf')





## SLOPE OF COVARIANCE TRAJECTORY
assign(paste('D12vCoVarSlope_ABRSim_Nu', N*U,sep = ""),  ggplot(data_thisNu, aes(x=D12,y=meanCoVarResp)) + geom_point(col = color) +
         geom_smooth(method='lm', formula= y~x, col = color, lwd = 0.5) + geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey')+
         theme_bw()+theme(panel.border = element_rect(colour = 'black', fill = NA, size = 1), 
                          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
                          axis.text.y =element_text(size = 35,family = 'Helvetica'),axis.text.x =element_text(size = 35,family = 'Helvetica', color = c('grey30', 'transparent','grey30', 'transparent','grey30', 'transparent')), axis.title = element_blank()) +
         scale_x_continuous(label = scientific_10) +  scale_y_continuous(label = scientific_10))



ggsave(D12vCoVarSlope_ABRSim_Nu1, file = 'D12vCoVarSlope_ABRSim_Nu1.pdf')


