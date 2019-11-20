
install.packages('MASS')
install.packages("corpcor")
data = read.table(file = 'clipboard', sep = '\t', header = FALSE) #copy with NO headers or row labels 


library("corpcor", lib.loc="~/R/win-library/3.5")
library("MASS", lib.loc="C:/Program Files/R/R-3.5.3/library")




## Simulate 100 full data sets based on error distributions for WT ##
# choose value for all knockouts from distributioon with center = reported growth rate 
# and SD = SD of wt replicate measures
WTGrowthRateMeans = c(1.002672, 0.7145567, 0.7089853, 0.8037542 ,0.7257983, 0.7347584, 0.6161176, 0.7, 0.7)
WTGrowthRateSDs   = c(0.03858625, 0.03253058, 0.03956598, 0.05238393, 0.01608364, 0.0213947, 0.03916092, 0.03, 0.03)
numKnockouts = dim(data)[1]
numDrugs = dim(data)[2]
data = as.matrix(data)

allSimmedData = array(data = 0, dim = c(numKnockouts, numDrugs, 100))
allSimmedData[,,1] = data
for(i in 2:100)
{
  print(i)
  for(iDrug in 1:numDrugs)
  {
    for(iKnockout in 1:numKnockouts)
    {
      allSimmedData[iKnockout, iDrug, i] = rnorm(1, mean = data[iKnockout, iDrug], sd = WTGrowthRateSDs[iDrug])
    }
  }
}



# Generate JDFE 
# Pr(s) is a mutivariante (n+1 dimmensional, n = #drugs) gaussian 
# characterized by mean fitnesses in drug environments and variance-covariance between different environments(sigma)

N = 8 #number of drug environments

# make all muVecs; vector with average fitness in each environment
allmuVecs = array(data = 0, dim = c(100, (numDrugs)))
for( q in 1:100)
{
  muVec = numeric(numDrugs)
  
  for( i in 1:(numDrugs))
  {
    muVec[i] = mean(allSimmedData[ , i, q])
  }
  
  allmuVecs[q, ] = muVec
}


# make all sigmas, (n+1)x(n+1) variance-covariance matrix
allSigmas = array(data = 0, dim = c((N+1), (N+1), 100))
for(q in 1:100)
{
  sigma = array(data = 0, dim = c((N+1), (N+1)))
  for( x in 1:(N+1))
  {
    for(y in 1:(N+1))
    {
      if(x == y)
      {
        sigma[x,y] = var(data[ , x])
        
      }else{
        
        sigma[x,y] = cov(data[,x], data[,y])
      }
    }
  }
  allSigmas[,,q] = sigma
}

# mutivariate normal distribution - characterized by muVec and sigma
# sample with n =1 from mvrnorm(n = 1, muVec, sigma) to get 1 potential s value vector (s in all environments with 1st in drug free)
sample =  mvrnorm(n = 1, muVec, sigma) - c(1, 0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7)
# subtract by 1 and 0.7s to turn relative fitnesses into selection coefficients (negative if worse than WT and positive if better)
# assumes WT fitn = 1 in no drug and 0.7 in all drugs 




#### Simulate Evolution in Switching Environments ####

# function to evolve population in given drug 1 time step 
evolveInDrug <- function(startIDarray, mutFitnesses, popsize, enviroEvolveID, U, muVec, sigma)
{
  numMuts = dim(startIDarray)[1]
  
  ## Selection and Drift ##
  
  Wbar = sum(startIDarray[ , 3] * (startIDarray[ ,2] / popsize))
  
  # generate probability vector
  probvec = numeric(numMuts)
  for( i in 1:numMuts)
  {
    probvec[i] = (startIDarray[i,3] * (startIDarray[i,2] / popsize)) / Wbar
  }
  
  # generate IDarray (with only sel/drift) for next gen
  nextGenIDarray = array(data = 0, dim = c(dim(startIDarray)[1], dim(startIDarray)[2]))
  nextGenOffspring = rmultinom(1, size = popsize, prob = probvec)
  
  ##???just use nextGenIDarray = startIDarray
  # nextGenIDarray[, 2] = nextGenOffspring
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
    # find s vecotr for this mutant
    svec = mvrnorm(n = 1, muVec, sigma) - c(1, 0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7)
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
      
      #make new row for this new mutant
      newRow = c(nextGenIDarray[dim(nextGenIDarray)[1], 1]+ 1, 1, startIDarray[i,3]*(1+ svec[enviroEvolveID]))
      
      #bind newrow with nextgenIDarray
      nextGenIDarray = rbind(nextGenIDarray, newRow)
      
      #add row with fitness of new mutant in all environments to mutFitnesses
      newFit = mutFitnesses[i , ]*(1+ svec)
      mutFitnesses = rbind(mutFitnesses, newFit)
    }
  }
  
  output = list(nextGenIDarray, mutFitnesses)
  
  return(output)
}  



popsize = 10^6
U = 10^-4
N = 8
numGens = 1000
numIterations = 100



allSlopeArrays = array(data = 0, dim = c((N+1), (N+1), 100))
for(numSimData in 2:100)
{
  print("Simmed Data Set #")
  print(numSimData)
  
  allAvgWbarAllEnviroEvolves = array(data = 0, dim = c(numGens, (N+1), (N+1)))
  for( k in 1:(N+1))
  {
    
    enviroEvolveID = k
    
    print("enviro #")
    print(k)
    
    allWbarsAllIterations = array(data = 0, dim = c(numGens, (N+1), numIterations))
    
    for(w in 1:numIterations)
    { 
      #print(w)
      
      startIDarray = array(data = c(1, popsize, 0.7), dim = c(1, 3))
      mutFitnesses = array(data = c(1, 0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7), dim = c(1, 9))
      allWbars = array(data = 0, dim = c(numGens, (N+1)))
      
      for( x in 1:numGens)
      {
        #print(x)
        output =  evolveInDrug(startIDarray , mutFitnesses, popsize, enviroEvolveID, U, muVec = allmuVecs[numSimData, ], sigma = allSigmas[,,numSimData])
        startIDarray = output[[1]]
        mutFitnesses = output[[2]]
        
        
        # find and save Wbars in all enviros 
        
        for( n in 1:(N+1))
        {
          allWbars[x, n] = sum(mutFitnesses[ , n] * (startIDarray[ ,2] / popsize))
        }
        
        # remove all muttypes with frequency of 0 from idarray and mutfitarray
        toRemove = c()
        for(q in 1:(dim(startIDarray)[1]))
        {
          if(startIDarray[q,2] == 0)
          {
            toRemove = c(toRemove, q)
            
          }
        }
        
        
      
        if(length(toRemove) > 0)
        {
          startIDarray = startIDarray[-toRemove, ]
          mutFitnesses = mutFitnesses[ -toRemove, ]
          
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
    
    averageAllWbar = array(data = 0, dim = c(numGens, (N+1)))
    for( i in 1:numGens)
    {
      for(y in 1:(N+1))
      {
        averageAllWbar[i,y] = mean(allWbarsAllIterations[i,y, ])
      }
    }
    
    allAvgWbarAllEnviroEvolves[ , , k] = averageAllWbar
  }
  
  
  slopeArray = array(data = 0, dim = c((N+1), (N+1)))
  for( i in 1:(N+1)) #through each evolve enviro
  { 
    for(q in 1:(N+1)) #through each enviro in each evolve enviro
    {
      #find slope of line q  evolved in enviro i
      slope = (log(allAvgWbarAllEnviroEvolves[ 1000, q, i]) - log(allAvgWbarAllEnviroEvolves[500, q , i])) / (1000 - 500)
      
      slopeArray[i,q] = slope
    }
  }
  
  allSlopeArrays[,,numSimData] = slopeArray
  
}




for( e in 1:(N+1))
{
  lwdVec = numeric(N+1) + 1
  lwdVec[e] = 3
  
  plot(log(allAvgWbarAllEnviroEvolves[ , 1, e]) , type = 'l',lwd = lwdVec[1],  ylim = c(-1,5),  xlab = 'Time', ylab = ' log(Population Mean Fitness)')
  lines(log(allAvgWbarAllEnviroEvolves[ , 2, e]),  col = 'green', lwd = lwdVec[2])
  lines(log(allAvgWbarAllEnviroEvolves[ , 3, e]), col = 'blue', lwd = lwdVec[3])
  lines(log(allAvgWbarAllEnviroEvolves[ , 4, e]), col = 'purple', lwd = lwdVec[4])
  lines(log(allAvgWbarAllEnviroEvolves[ , 5, e]), col = 'orange', lwd = lwdVec[5])
  lines(log(allAvgWbarAllEnviroEvolves[ , 6, e]), col = 'cyan', lwd = lwdVec[6])
  lines(log(allAvgWbarAllEnviroEvolves[ , 7, e]), col = 'red', lwd = lwdVec[7])
  lines(log(allAvgWbarAllEnviroEvolves[ , 8, e]), col = 'grey', lwd = lwdVec[8])
  lines(log(allAvgWbarAllEnviroEvolves[ , 9, e]), col = 'brown', lwd = lwdVec[9])
}




### Find slopes of all lines in the average trajectories in all drug evolve environments

slopeArray = array(data = 0, dim = c((N+1), (N+1)))
for( i in 1:(N+1)) #through each evolve enviro
{ 
  for(q in 1:(N+1)) #through each enviro in each evolve enviro
  {
    #find slope of line q  evolved in enviro i
    slope = (log(allAvgWbarAllEnviroEvolves[ 1000, q, i]) - log(allAvgWbarAllEnviroEvolves[500, q , i])) / (1000 - 500)
    
    slopeArray[i,q] = slope
  }
}

# calculate correlation coefficients of drug pairs and put in array
noDrug = data[ ,1]
CHL = data[ , 2]
CPR = data[ , 3]
MEC = data[ , 4]
NIT = data[ , 5]
TET = data[ , 6]
TMP = data[ , 7]
FOX = data[ , 8]
AMP = data[ , 9]

corrMatrix = array(data = 0, dim = c((N+1), (N+1)))
for(i in 1:(N+1))
{
  for(q in 1:(N+1))
  {
    cor = cor(data[ ,i], data[,q])
    corrMatrix[i,q] = cor
  }
}


# calculatee center of distributions 
# center of distribution of points (x,y) equals (xbar, ybar)

centerarrayX = array(data = 0 , dim = c((N+1), (N+1)))
centerarrayY = array(data = 0 , dim = c((N+1), (N+1)))

for(i in 1:(N+1))
{
  for(q in 1:(N+1))
  {
    midpoint = c(mean(data[ , i]), mean(data[ , q]))
    centerarrayX[i,q] = midpoint[1]
    centerarrayY[i,q] = midpoint[2]
    
  }
}





SlopeCorrData = read.table(file = 'clipboard', sep = '\t', header = FALSE) #copy with NO headers or row labels 

correlation = cor(SlopeCorrData[ ,1], SlopeCorrData[ ,2])
plot(SlopeCorrData[ ,1], SlopeCorrData[ ,2])

lm = array(data = 0, dim = c(length(seq(0.1, 1, by = 0.1)), 2))
for(i in seq(0.1,1, by = 0.1))
{
  lm[i, 1] = i
  lm[i, 2] = i*0.0062 + 0.0022
}
lines(lm[ ,1], lm[ ,2])









