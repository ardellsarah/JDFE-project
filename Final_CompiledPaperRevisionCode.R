########## COMPILED CODE FOR PROJECT REVISION 1 ##########
# NEW TO THIS FILE: 
# 1. new scaling of r2 by UbAdjust (to adress non-monotonicity)
# 2. New epistasis section
#     - no more of the old theory, just use 4 JDFEs from rank order and have 
#       shrink with fitness in home 


#### all necessary packages, download if do not have already 
library('MASS')# matrix stats
library(viridis)# color pallete 
library(ggplot2) # main plotting function
library(dplyr) # data reformatting
library(tidyr) # data reformatting
library(stringr) # string manipulation
library(gridExtra)# for more complex figures 
library(multipanelfigure)#for more complex figures 
library(magrittr)# special operator for vector/array reformatting
library(ggpubr)#for more complex figures 
library("corpcor")# matrix stats
library('matrixStats')# matrix stats
library('readr')#simple reading of rectangular data
library('vctrs') # help in function development
library('colorRamps')# color pallete 
library('fitdistrplus')# advanced distribution fitting
library(pracma) # numeric mathematical functions
library("scales") # help with plotting, determining breaks
library('ggallin') # additions to ggplot 
library(utils) # write csvs
library(cowplot)

#### All core functions for the code ####

scientific_10 <- function(x) { # for plotting with scientific notation scales
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}




###################################### F I G U R E     1  ################################################################
                                 ##     CHEVERAU   JDFES   ##
setwd("~/Dropbox/Antibiotic Resistance JDFE Project copy/Paper Revision 1 outputs/code/testCompleteCode/fig_cheverauJDFES")
# (optional) set working directory to where you'd like to save the outputs of this section, example below 
# setwd("~/Documents/JDFE_Project/Figures/Figure_ABR_JDFES")

### NOTE: WHEREVER YOU SET THE WORKING DIRECTORY MUST CONTAIN THE FILES "KOGrDATAPrtoStop.csv' AND 'WTReplicate.csv' FOR THIS SECTION TO RUN 


####    Import data     ####

# Import KO growth rate data  so can have full data (including KO names)
KOGrowthRates = read.csv('KOGrDATAPrtoStop.csv', header = 1, sep = ",")
# remove KOID No drug and and FOX and AMP (last 2 drugs), as dont have WT replicate data for them
KOPrtoStop = KOGrowthRates[,11]
KOGrowthRates = as.data.frame(KOGrowthRates[,3:8])
colnames(KOGrowthRates) = c('CHL', 'CPR', 'MEC', 'NIT', 'TET', 'TMP')


WTReplicateData = read.csv('WTReplicate.csv', header = 1, sep = ":")
# remove excess rows; keep only  1st 6 drug enviros 
WTReplicateData = as.data.frame(WTReplicateData[,3:8])
colnames(WTReplicateData) = c( 'CHL', 'CPR', 'MEC', 'NIT', 'TET', 'TMP')

# Important Variables #
WTGrowthRateMeans = colMeans(as.matrix(WTReplicateData))
WTGrowthRateSDs   = colSds(as.matrix(WTReplicateData) )
numKnockouts = dim(KOGrowthRates)[1]
numDrugs = dim(KOGrowthRates)[2]



sValAllMuts = KOGrowthRates
count = 0
for(i in colnames(KOGrowthRates))
{
  count = count + 1
  sValAllMuts[[i]] = KOGrowthRates[[i]] - WTGrowthRateMeans[count]
}


####   False Discovery Rate (FDR): Benjamini Hochberg procedure   ####

## STEP 1: Find p-val of all muts in each drug 
  # chance that each of them is actually from null WT error dist
  # assume WT error dist gaussian

## STEP 2: rank order muts by ascending p-val

## STEP 3: Clac BH critical value, (i/m)*Q
   # i = pVal RANK
   # m = number samples tested total
   # Q = desired FDR 

## STEP 4: find largest pVal < corresponding BH value 
   # this mut s-val is the cutoff for calling beneficial mutations with FDR = desired FDR

desiredFDRVec = c(0.25 ) #0.1, , 0.5)

benCutoffVec = c()
deltCutoffVec = c()
countBenMuts = c()
countDeltMuts = c()
sampsvec = c()
FDRID = c()
nameVec = c()


for(desiredFDR in desiredFDRVec)
{
  count = 0
  for(i in colnames(KOGrowthRates)) # go through all drugs
  {
    count = count + 1
    muts = sValAllMuts[[i]]
    
    
    ### BENEFICIAL MUTS
    #step 1
    pVals = 1- pnorm(muts, mean = 0, sd = WTGrowthRateSDs[count] ) # pVal = chance of being equal or HIGHER than each mut val by chance, center of dist at 0 because sVals have been adjusted with the WT GR mean 
    assign(paste0('pVals_ben_home',i), pVals)
    
    # step 2
    pRanks = rank(pVals, ties.method="min")
    df = data.frame(muts, pVals, pRanks)
    df = df[order(pRanks),]
    
    # step 3
    #numSamps = length(muts)
    numSamps = length(muts)
    df$BH_Val = (df$pRanks / numSamps) * desiredFDR
    
    # step 4
    benCutoff = df$muts[which.max(rank(df$pVals < df$BH_Val, ties.method = "first") )] # cuttoff is at highest pval which is leff than its BH, so use ascending order ranking for ties 
    benCutoffVec = c(benCutoffVec, benCutoff)
    countBenMuts = c(countBenMuts ,  sum(muts >= benCutoff))
    
    
    ### DELETERIOUS MUTS
    #step 1
    pVals =  pnorm(muts, mean = 0, sd = WTGrowthRateSDs[count] ) # pVal = chance of being equal or LOWER than each mut val by chance, center of dist at 0 because sVals have been adjusted with the WT GR mean 
    
    # step 2
    pRanks = rank(pVals, ties.method="min")
    df = data.frame(muts, pVals, pRanks)
    df = df[order(pRanks),]
    
    # step 3
    #numSamps = length(muts)
    numSamps = length(muts)
    df$BH_Val = (df$pRanks / numSamps) * desiredFDR
    
    # step 4
    deltCutoff = df$muts[which.max(rank(df$pVals < df$BH_Val, ties.method = "first") )] # cuttoff is at highest pval which is leff than its BH, so use ascending order ranking for ties 
    deltCutoffVec = c(deltCutoffVec, deltCutoff)
    countDeltMuts = c( countDeltMuts , sum(muts <= deltCutoff))
    
    
    
    FDRID = c(FDRID, desiredFDR)
    nameVec = c(nameVec, i)
    sampsvec = c(sampsvec, numSamps)
  }
  
}


allFDR_Results = data.frame(FDRID, nameVec,  benCutoffVec , countBenMuts,  deltCutoffVec , countDeltMuts)
colnames(allFDR_Results) = c('FDR', 'Drug', 'benCutoff_sVal', 'numBenMuts', 'deltCutoff_sVal', 'numDeltMuts')
write.csv(allFDR_Results, file = 'allFDR_Results.csv')




##### ID KOS as 'beneficial', 'deleterious', or 'neutral' in all Envs and all FDRs  ####
getGeneNames =  read.csv('KOGrDATAPrtoStop.csv', header = 1, sep = ",") 
mutTypeDf = cbind(getGeneNames$Common.gene.name, sValAllMuts)
colnames(mutTypeDf) = c('geneName', colnames(sValAllMuts))


for( desiredFDR in desiredFDRVec)
{
  count = 0 
  for(home in colnames(sValAllMuts))
  {
    count = count + 1
    
    benCutoff = subset(allFDR_Results, FDR == desiredFDR & Drug == home )$benCutoff_sVal
    deltCutoff = subset(allFDR_Results, FDR == desiredFDR & Drug == home )$deltCutoff_sVal
    
    mutFitEffects = mutTypeDf[[home]] # select all mut fit effects in this home
    
    mutTypes_thisFDR_thisHome = as.character(numeric(length(mutFitEffects)))
    mutTypes_thisFDR_thisHome[mutFitEffects >= benCutoff] = "Beneficial"
    mutTypes_thisFDR_thisHome[mutFitEffects <= deltCutoff] = "Deleterious"
    mutTypes_thisFDR_thisHome[(mutFitEffects < benCutoff) & (mutFitEffects > deltCutoff)] = "Neutral"
    
    nameForDF  = paste0("mutType_", home ,"_FDR_", desiredFDR)
    mutTypeDf[[nameForDF]] = mutTypes_thisFDR_thisHome
    

  }
  

}


## save Data 
write.csv(mutTypeDf, "mutTypeDf.csv")


 
#### Collateral Resistance and Sensivitity  ####
desiredFDRVec_colEffect = c( 0.1 )#, 0.05, 0.25, 0.5) # what FDR want for calling collateral effect
drugNames = colnames(sValAllMuts)[-c(1,6)] # remove CHL and TMP because they find so few beneficial mutants it is not possible to do col res or sens test effectively 
WTGrowthRateSDs_focalDrugs  =  WTGrowthRateSDs[-c(1,6)]

colResCutoffVec = c()
countCRMuts = c()
colSensCutoffVec = c()
countCSMuts = c()
countCNMuts = c()
FDR_collEffect = c()
FDR_home = c()
homeDrug= c()
nonHomeDrug = c()



for(desiredFDR_collateralEffect in desiredFDRVec_colEffect)
{
  
  for(desiredFDR_home in desiredFDRVec) # fdr for calling beneficial and delterious in home
  {
    for(i in 1:length(drugNames)) 
    {
      for(j in 1:length(drugNames)) 
      {
        
        if(i!=j)
        {
          
          title = drugNames[i]
          title2 = drugNames[j]
          
          benCutoff_i = subset(allFDR_Results, FDR == desiredFDR_home & Drug == title )$benCutoff_sVal
          deltCutoff_i = subset(allFDR_Results, FDR == desiredFDR_home & Drug == title )$deltCutoff_sVal
          
          
          mutFitEffects_i = sValAllMuts[[title]] # select all mut fit effects in this home
          mutFitEffects_j_benin_i = sValAllMuts[[title2]][mutFitEffects_i >= benCutoff_i] # select all mut fit effects in this Non-home that are called beneficial in the nonHome
          
          
          ### COLLATERAL RESISTANCE
          #step 1
          pVals = 1- pnorm(mutFitEffects_j_benin_i , mean = 0, sd = WTGrowthRateSDs_focalDrugs[j] )
          # pVal = chance of being equal or greater than each mut val by chance in env2 for muts called ben in env1
          
          # savepval calcd for all genes
          assign(paste0('pVals_CR_home',title, '_nonhome_',title2 ),( 1- pnorm(sValAllMuts[[title2]] , mean = 0, sd = WTGrowthRateSDs_focalDrugs[j]) ))
          
          
          # step 2
          pRanks = rank(pVals, ties.method="min") # min ties gives conservative estimate because higher rank = better change > p for same pval 
          df = data.frame( mutFitEffects_j_benin_i, pVals, pRanks)
          df = df[order(pRanks),]
          
          # step 3
          numSamps = length(  mutFitEffects_j_benin_i)
          df$BH_Val = (df$pRanks / numSamps) * desiredFDR_collateralEffect
          
          # step 4
          colResCutoff = df$mutFitEffects_j_benin_i[which.max(rank(df$pVals < df$BH_Val, ties.method = "first") )] # cuttoff is at highest pval which is leff than its BH, so use ascending order ranking for ties 
          # if col res cutoff < 0, set so col res cutoff = 0 
          if(colResCutoff < 0){colResCutoff = 0}
          colResCutoffVec = c(colResCutoffVec, colResCutoff)
          countCRMuts = c( countCRMuts , sum( mutFitEffects_j_benin_i >= colResCutoff))
          
          
          ### COLLATERAL SENSITIVITY
          #step 1
          pVals = pnorm(mutFitEffects_j_benin_i , mean = 0, sd = WTGrowthRateSDs_focalDrugs[j] )
          # pVal = chance of being equal or LESSER than each mut val by chance in env2 for muts called ben in env1
          
          # savepval calcd for all genes
          assign(paste0('pVals_CS_home',title, '_nonhome_',title2 ),(  pnorm(sValAllMuts[[title2]] , mean = 0, sd = WTGrowthRateSDs_focalDrugs[j]) ))
          
          # step 2
          pRanks = rank(pVals, ties.method="min") # min ties gives conservative estimate because higher rank = better change > p for same pval 
          df = data.frame( mutFitEffects_j_benin_i, pVals, pRanks)
          df = df[order(pRanks),]
          
          # step 3
          numSamps = length(mutFitEffects_j_benin_i)
          df$BH_Val = (df$pRanks / numSamps) * desiredFDR_collateralEffect
          
          # step 4
          colSensCutoff = df$mutFitEffects_j_benin_i[which.max(rank(df$pVals < df$BH_Val, ties.method = "first") )] # cuttoff is at highest pval which is leff than its BH, so use ascending order ranking for ties 
          
          
          # if col sens cutoff > 0, set so col sens cutoff = 0 
          if(colSensCutoff >0){colSensCutoff = 0}
          colSensCutoffVec = c(colSensCutoffVec, colSensCutoff)
          countCSMuts = c( countCSMuts , sum( mutFitEffects_j_benin_i <= colSensCutoff))
          
  
          countCNMuts = c(countCNMuts,  length(mutFitEffects_j_benin_i) - sum( mutFitEffects_j_benin_i <= colSensCutoff) - sum( mutFitEffects_j_benin_i >= colResCutoff) )
          FDR_collEffect = c(FDR_collEffect, desiredFDR_collateralEffect)
          FDR_home = c(FDR_home, desiredFDR_home)
          homeDrug= c(homeDrug, title)
          nonHomeDrug = c(nonHomeDrug, title2)
          
          
          
        }
      }
    }
  }
}


FDR_results_colRes_colSens_df = data.frame(homeDrug, nonHomeDrug, FDR_home, FDR_collEffect, colResCutoffVec, countCRMuts, colSensCutoffVec, countCSMuts, countCNMuts)
FDR_results_colRes_colSens_df$expNumFalseCR = FDR_results_colRes_colSens_df$FDR_collEffect * FDR_results_colRes_colSens_df$countCRMuts
FDR_results_colRes_colSens_df$expNumFalseCS = FDR_results_colRes_colSens_df$FDR_collEffect * FDR_results_colRes_colSens_df$countCSMuts


# save Data
write.csv(FDR_results_colRes_colSens_df, file = 'FDR_results_colRes_colSens_df.csv' )


### Make CSV of the ID of the coll effect of each mutant in each durg pair ####

focalDrugs = c("CPR", "MEC", "NIT", "TET")


CRCS_ID_df = data.frame(getGeneNames$Common.gene.name)
colnames(CRCS_ID_df) = 'GeneName'

for(home in focalDrugs) 
{
  for(nonhome in focalDrugs) 
  {
    
    if(home!=nonhome)
    {
    
      mutFitEffects_home = mutTypeDf[[home]]
      mutFitEffects_nonhome = mutTypeDf[[nonhome]]
      
     # save cutoffs of relevance 
      benCutoff_home = subset(allFDR_Results, FDR == 0.25 & Drug == home)$benCutoff_sVal
      colSens_nonHome = subset(FDR_results_colRes_colSens_df, homeDrug == home & nonHomeDrug == nonhome)$colSensCutoffVec
      colRes_nonHome = subset(FDR_results_colRes_colSens_df, homeDrug == home & nonHomeDrug == nonhome)$colResCutoffVec
      
      # if mut beneficial in home with FDR 0.25, ID coll effect with FDR 0.1
      mutIDS_thispair = as.character(numeric(length(mutFitEffects_home)))
      
      beninHome = mutFitEffects_home>=benCutoff_home
      
      # ID Col Sens
      mutIDS_thispair[(mutFitEffects_home >= benCutoff_home &mutFitEffects_nonhome <= colSens_nonHome )] = 'Collaterally Sensitive'
      
      # ID Col Res
      mutIDS_thispair[(mutFitEffects_home >= benCutoff_home & mutFitEffects_nonhome >= colRes_nonHome )] = 'Collaterally Resistant'
      
      # ID CN
      mutIDS_thispair[(mutFitEffects_home >= benCutoff_home & (mutFitEffects_nonhome < colRes_nonHome) & (mutFitEffects_nonhome > colSens_nonHome))] = 'Collaterally Neutral'
      
      
      CRCS_ID_df[[paste0('Home_', home, '_NonHome_', nonhome)]] = mutIDS_thispair
      
    }
  }
}

write.csv(CRCS_ID_df, file = 'CRCS_ID_df_cheverauJDFES_full.csv')



## Now, go through and only save data on genes with a call in at least one pair 

# which gene indexes/names have significance in at least one drug
indexes_haveCall = c()
geneNames_haveCall = c()
for(i in 1:length(getGeneNames$Common.gene.name) )
{
  calls_thisGene = as.character(CRCS_ID_df[i, c(-1)])
  if(sum(calls_thisGene == '0')==12)
  {
    
  }else{
    indexes_haveCall = c(indexes_haveCall,i )
    geneNames_haveCall = c( geneNames_haveCall,getGeneNames$Common.gene.name[i])

  }
  
}


# go through these genes and make vectors for all outputs

callDF = data.frame(matrix(ncol = 41, nrow = 1))
colnames(callDF) = c('geneName', 
                     'call_CPR_MEC', 'call_CPR_NIT',  'call_CPR_TET',
                     'call_MEC_CPR', 'call_MEC_NIT',  'call_MEC_TET',
                     'call_NIT_CPR', 'call_NIT_MEC',  'call_NIT_TET',
                     'call_TET_CPR', 'call_TET_MEC',  'call_TET_NIT',
                     'pResistantCPR', 'pResistantMEC',  'pResistantNIT',  'pResistantTET', 
                     'pCR_CPR_MEC', 'pCR_CPR_NIT',  'pCR_CPR_TET',
                     'pCR_MEC_CPR', 'pCR_MEC_NIT',  'pCR_MEC_TET',
                     'pCR_NIT_CPR', 'pCR_NIT_MEC',  'pCR_NIT_TET',
                     'pCR_TET_CPR', 'pCR_TET_MEC',  'pCR_TET_NIT',
                     'pCS_CPR_MEC', 'pCS_CPR_NIT',  'pCS_CPR_TET',
                     'pCS_MEC_CPR', 'pCS_MEC_NIT',  'pCS_MEC_TET',
                     'pCS_NIT_CPR', 'pCS_NIT_MEC',  'pCS_NIT_TET',
                     'pCS_TET_CPR', 'pCS_TET_MEC',  'pCS_TET_NIT')

for(gene in indexes_haveCall)
{

      geneName = as.character(getGeneNames$Common.gene.name[gene])
      call_CPR_MEC = CRCS_ID_df[[paste0('Home_', 'CPR', '_NonHome_', 'MEC')]][gene]
      call_CPR_NIT  = CRCS_ID_df[[paste0('Home_', 'CPR', '_NonHome_', 'NIT')]][gene]
      call_CPR_TET  = CRCS_ID_df[[paste0('Home_', 'CPR', '_NonHome_', 'TET')]][gene]
      call_MEC_CPR  = CRCS_ID_df[[paste0('Home_', 'MEC', '_NonHome_', 'CPR')]][gene]
      call_MEC_NIT  = CRCS_ID_df[[paste0('Home_', 'MEC', '_NonHome_', 'NIT')]][gene]
      call_MEC_TET =   CRCS_ID_df[[paste0('Home_', 'MEC', '_NonHome_', 'TET')]][gene]
      call_NIT_CPR  = CRCS_ID_df[[paste0('Home_', 'NIT', '_NonHome_', 'CPR')]][gene]
      call_NIT_MEC   = CRCS_ID_df[[paste0('Home_', 'NIT', '_NonHome_', 'MEC')]][gene]
      call_NIT_TET = CRCS_ID_df[[paste0('Home_', 'NIT', '_NonHome_', 'TET')]][gene]
      call_TET_CPR  = CRCS_ID_df[[paste0('Home_', 'TET', '_NonHome_', 'CPR')]][gene]
      call_TET_MEC   = CRCS_ID_df[[paste0('Home_', 'TET', '_NonHome_', 'MEC')]][gene]
      call_TET_NIT = CRCS_ID_df[[paste0('Home_', 'TET', '_NonHome_', 'NIT')]][gene]
      pResistantCPR = pVals_ben_homeCPR[gene]
      pResistantMEC = pVals_ben_homeMEC[gene]
      pResistantNIT  =  pVals_ben_homeNIT[gene]
      pResistantTET = pVals_ben_homeTET[gene]
      pCR_CPR_MEC = pVals_CR_homeCPR_nonhome_MEC[gene]
      pCR_CPR_NIT = pVals_CR_homeCPR_nonhome_MEC[gene]
      pCR_CPR_TET= pVals_CR_homeCPR_nonhome_TET[gene]
      pCR_MEC_CPR= pVals_CR_homeMEC_nonhome_CPR[gene]
      pCR_MEC_NIT = pVals_CR_homeMEC_nonhome_NIT[gene]
      pCR_MEC_TET= pVals_CR_homeMEC_nonhome_TET[gene]
      pCR_NIT_CPR = pVals_CR_homeNIT_nonhome_CPR[gene]
      pCR_NIT_MEC = pVals_CR_homeNIT_nonhome_MEC[gene]
      pCR_NIT_TET= pVals_CR_homeNIT_nonhome_TET[gene]
      pCR_TET_CPR = pVals_CR_homeTET_nonhome_CPR[gene]
      pCR_TET_MEC = pVals_CR_homeTET_nonhome_MEC[gene]
      pCR_TET_NIT= pVals_CR_homeTET_nonhome_NIT[gene]
      pCS_CPR_MEC = pVals_CS_homeCPR_nonhome_MEC[gene]
      pCS_CPR_NIT = pVals_CS_homeCPR_nonhome_MEC[gene]
      pCS_CPR_TET= pVals_CS_homeCPR_nonhome_MEC[gene]
      pCS_MEC_CPR= pVals_CS_homeMEC_nonhome_CPR[gene]
      pCS_MEC_NIT = pVals_CS_homeMEC_nonhome_NIT[gene]
      pCS_MEC_TET= pVals_CS_homeMEC_nonhome_TET[gene]
      pCS_NIT_CPR = pVals_CS_homeNIT_nonhome_CPR[gene]
      pCS_NIT_MEC = pVals_CS_homeNIT_nonhome_MEC[gene]
      pCS_NIT_TET= pVals_CS_homeNIT_nonhome_TET[gene]
      pCS_TET_CPR = pVals_CS_homeTET_nonhome_CPR[gene]
      pCS_TET_MEC = pVals_CS_homeTET_nonhome_MEC[gene]
      pCS_TET_NIT= pVals_CS_homeTET_nonhome_NIT[gene]
      
      
      NEWROW = c(geneName, 
                 call_CPR_MEC, call_CPR_NIT,  call_CPR_TET,
                 call_MEC_CPR, call_MEC_NIT,  call_MEC_TET,
                 call_NIT_CPR, call_NIT_MEC,  call_NIT_TET,
                 call_TET_CPR, call_TET_MEC,  call_TET_NIT,
                 pResistantCPR, pResistantMEC,  pResistantNIT,  pResistantTET, 
                 pCR_CPR_MEC, pCR_CPR_NIT,  pCR_CPR_TET,
                 pCR_MEC_CPR, pCR_MEC_NIT,  pCR_MEC_TET,
                 pCR_NIT_CPR, pCR_NIT_MEC,  pCR_NIT_TET,
                 pCR_TET_CPR, pCR_TET_MEC,  pCR_TET_NIT,
                 pCS_CPR_MEC, pCS_CPR_NIT,  pCS_CPR_TET,
                 pCS_MEC_CPR, pCS_MEC_NIT,  pCS_MEC_TET,
                 pCS_NIT_CPR, pCS_NIT_MEC,  pCS_NIT_TET,
                 pCS_TET_CPR, pCS_TET_MEC,  pCS_TET_NIT)
      
      callDF  = rbind(callDF, NEWROW)

}

callDF_complete = callDF[c(-1), ]
write.csv(callDF_complete, file ='callDF_complete_cheverau.csv' )



 
#### Graph the JDFEs####

# Choose the FDRs will use for Homes and Coll effects from tested ones
desiredFDR_home = 0.25
desiredFDR_collateralEffect = 0.1
focalDrugs = c("CPR", "MEC", "NIT", "TET") # only drugs with actual beneficial mutants for this FDR
WTGrowthRateSDs_focalDrugs  =  WTGrowthRateSDs[-c(1,6)]


## Save scatter plot color coded with collateral effects for all pairs
 # for diagonal (eg home = nonhome), save histogram with DFE + error dist + indicator of ben and delt cutoffs 

for(i in 1:length(focalDrugs) )
{
  for(j in 1:length(focalDrugs)) 
  {
    
    if(i!=j)
    {
      
      title = focalDrugs[i]
      title2 = focalDrugs[j]
      
      # home cutoffs
      benCutoff_home = subset(allFDR_Results, FDR == desiredFDR_home & Drug == title )$benCutoff_sVal
      deltCutoff_home = subset(allFDR_Results, FDR == desiredFDR_home & Drug == title )$deltCutoff_sVal
      
      # collateral cutoffs, CR = collateral resistance, CS = collateral sensitivity , CN = collateral neutrality 
      CR_cutoff = subset(FDR_results_colRes_colSens_df, FDR_home == desiredFDR_home & FDR_collEffect == desiredFDR_collateralEffect & homeDrug == title & nonHomeDrug ==title2)$colResCutoffVec
      CS_cutoff = subset(FDR_results_colRes_colSens_df, FDR_home == desiredFDR_home & FDR_collEffect == desiredFDR_collateralEffect & homeDrug == title & nonHomeDrug ==title2)$colSensCutoffVec
      
      
      # make vector ID-ing every mut by collateral effect
      muts_home = sValAllMuts[[title]]
      muts_nonHome = sValAllMuts[[title2]]
      
      
      colID_muts = as.character(numeric(length(muts_home)))
      colID_muts[muts_home >= benCutoff_home & muts_nonHome >= CR_cutoff] = 'CR'  # muts greater than cutoffs in both home are collaterally resistant 
      colID_muts[muts_home >= benCutoff_home & muts_nonHome <= CS_cutoff] = 'CS'  # muts greater than ben cutoff in  home and less than col sens in nonhome are collaterally sensitive
      colID_muts[muts_home >= benCutoff_home & muts_nonHome > CS_cutoff & muts_nonHome < CR_cutoff] = 'CN'  # muts greater than ben cutoff in  home and in between nonhome cutoffs are collateral neutral 
      colID_muts[muts_home < benCutoff_home] = 'Other'
      
      df = data.frame(muts_home, muts_nonHome, colID_muts)
      colnames(df) = c('HomeS', 'NonHomeS', 'mutType')
      
      
      
      JDFE = ggplot(df, aes(x=HomeS, y=NonHomeS, color = mutType)) +  xlim(-1,0.3) +ylim(-1, 0.3)+
        geom_point( size = 3.5, alpha = 0.5)+ scale_color_manual(values = c( '#13944a','orange', 'blue','grey'), drop = 'FALSE') +
        theme_bw()+  xlab(title) + ylab(title2)+
        theme(panel.border = element_rect(color = 'black', fill = NA, size = 3), panel.grid.major = element_blank(),
              legend.position = 'none',axis.text.x = element_text(size = 25, color = 'black'), axis.text.y = element_text(size = 25, color = 'black'),panel.grid.minor = element_blank(), 
              axis.line = element_line(colour = "black"), axis.title.x = element_blank(), axis.title.y = element_blank()) +
        geom_hline(yintercept =  0, linetype = 'dashed') + 
        geom_vline(xintercept =  0, linetype = 'dashed')
         # green = CN. Orange = CR, Blue = CS, Grey = Other 
    
      ggsave(paste0('JDFE_home_', title, '_nonHome_', title2, '_FDRhome_', desiredFDR_home,'_FDRnonHome_', desiredFDR_collateralEffect, '.pdf' ),JDFE , width = 4.71, height = 4.65)
      
      
      
    }else
    {
      
      title = focalDrugs[i]
      
      # cutoffs
      benCutoff_home = subset(allFDR_Results, FDR == desiredFDR_home & Drug == title )$benCutoff_sVal
      deltCutoff_home = subset(allFDR_Results, FDR == desiredFDR_home & Drug == title )$deltCutoff_sVal
    
      muts_home = sValAllMuts[[title]]
      simError_dist = rnorm(length(muts_home), mean = 0, sd = WTGrowthRateSDs_focalDrugs[i])
      
      df = data.frame(muts_home, simError_dist)
      
      DFE = ggplot() + geom_histogram(data = df, aes(x = muts_home , y =..count../sum(..count..)),fill = 'grey40',bins = 50) +
        geom_histogram( aes(x = simError_dist, y =..count../sum(..count..)),fill = '#f80101',alpha = 0.5,bins = 50) +
        theme_bw() +  scale_fill_manual(values = c('grey40', '#f80101'))+  
        xlab(title)+ ylab("Density") +
        theme(panel.border = element_rect(color = 'black', fill = NA, size = 3), panel.grid.major = element_blank(),axis.text.x = element_text(size = 25, color = 'black'), axis.text.y = element_text(size = 25, color = 'black'),
              panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = 'none') +
        geom_vline(xintercept =  0, linetype = 'dashed')+ 
        coord_cartesian(ylim=c(0,0.3), xlim = c(-1,0.3)) +
        scale_y_continuous(expand = c(0,0)) 
      
      ggsave(paste0('DFE_', title,  '_FDRhome_', desiredFDR_home, '.pdf' ),DFE , width = 4.71, height = 4.65)
      
      
      
    }# end if else
    
  } # j loop
  
} # i loop


#### SUPP TO FIG 1: Stifler et al JDFEs ####

#setwd("~/FolderForOutputs")

## Must have the following file in working directory


stifflerData = read.csv("Single_Gene_JDFE_forR.csv")


### Find 90% CI for effect of every mutation for which there is replication   ####
# the AMP concentrations for which there is replication are  39, 156, 625 and 2500 
# CEF has no reported replication

z_90CI = 1.645 # z-score to use in calculation

ampConcVec = c('Amp_39', 'Amp_156', 'Amp_625', 'Amp_2500')


concID = c()
countBenMuts = c()
countDeltMuts = c()
countNeutralMuts = c()
countisNA = c()

countCRmuts = c()
countCSmuts = c()
countCNmuts = c()


for(conc in ampConcVec)
{
  
  t1_mutEffects = stifflerData[[paste0(conc, "_t1")]]
  t2_mutEffects = stifflerData[[paste0(conc, "_t2")]]
  array_mutEffects = cbind(t1_mutEffects, t2_mutEffects)
  
  meanEffects = rowMeans(array_mutEffects)
  sds = rowSds(array_mutEffects)
  
  # calc CI
  CI_lowerBound = meanEffects - z_90CI*(sds / sqrt(2)) # n = 2
  CI_upperBound = meanEffects + z_90CI*(sds / sqrt(2)) # n = 2
  
  
  # assign p-vals to vec
  assign(paste0('pVals_CR', conc), 1- pnorm(meanEffects, mean = 0, sd = sds) )
  assign(paste0('pVals_CS', conc),  pnorm(meanEffects, mean = 0, sd = sds) )
  
  
  
  # add to data frame
  stifflerData[[paste0(conc, '_CI_lowerBound' )]] = CI_lowerBound
  stifflerData[[paste0(conc, '_CI_upperBound' )]] = CI_upperBound
  
  # make mut effect ID, where identify all muts as beneficial, deleterious or neutral 
  mutIDVec = as.character(numeric(length(meanEffects)))
  mutIDVec[CI_lowerBound > 0] = 'Beneficial' # if the lower bound above 0, we're confident is actually beneficial (where confidence ~90%)
  mutIDVec[CI_upperBound < 0] = 'Deleterious' # if the upper bound is below 0, confident is deleterious
  mutIDVec[CI_upperBound > 0 & CI_lowerBound < 0] = 'Neutral' # if the CI crosses 0, will assume neutral 
  mutIDVec[is.na(meanEffects)] = NA # anything where arent measured values will have NA mut catagory 
  
  # add to datafram
  stifflerData[[paste0(conc, '_mutID')]] = mutIDVec
  
  
  # save summary data for these envs 
  
  # first, need to remove NAs from the mutIDvec
  mutIDVec_nonNA = mutIDVec[!is.na(mutIDVec)]
  
  concID = c(concID, conc)
  countBenMuts = c(countBenMuts, sum(mutIDVec_nonNA == 'Beneficial'))
  countDeltMuts = c(countDeltMuts, sum(mutIDVec_nonNA == 'Deleterious'))
  countNeutralMuts = c(countNeutralMuts, sum(mutIDVec_nonNA == "Neutral"))
  countisNA = c(countisNA, sum(is.na(mutIDVec)))
  
  # Identify collateral effects 
  muts_CEF = stifflerData$cef_0.2ngmL # no error data for CEF, so assume all mutations correctly identified by cutoff at 0 
  
  mutIDVec_col = as.character(numeric(length(meanEffects)))
  mutIDVec_col[ muts_CEF > 0 &  CI_lowerBound > 0] = 'Collaterally Resistant' # beneficial in both , collateral resistance
  mutIDVec_col[muts_CEF > 0 &  CI_upperBound < 0] = 'Collaterally Sensitive' # deleterious in non home ben in home, collateral sensitivity 
  mutIDVec_col[muts_CEF > 0 & CI_upperBound > 0 & CI_lowerBound < 0] = 'Collaterally Neutral' # ben in home, neutral in non home, collateral neutrality 
  mutIDVec_col[muts_CEF < 0 ] = 'Other' # is not ben in home, is other 
  mutIDVec_col[is.na(meanEffects) | is.na(muts_CEF)] = NA # anything where arent measured values will have NA mut catagory 
  
  
  # add to datafram
  stifflerData[[paste0(conc, '_collateralMutID')]] = mutIDVec_col
  assign(paste0(conc, '_collateralMutID'),mutIDVec_col )
  # save summary data
  mutIDVec_col_nonNA = mutIDVec_col[!is.na(mutIDVec_col)]
  
  countCRmuts = c(countCRmuts, sum(mutIDVec_col_nonNA == 'Collaterally Resistant'))
  countCSmuts = c(countCSmuts, sum(mutIDVec_col_nonNA == 'Collaterally Sensitive'))
  countCNmuts = c(countCNmuts, sum(mutIDVec_col_nonNA == 'Collaterally Neutral'))
  
  
  
  
  
  #### Plot and Save the JDFE 
  
  df = data.frame(muts_CEF, meanEffects, mutIDVec_col)
  colnames(df) = c('HomeS', 'NonHomeS', 'mutType')
  
  
  df$mutType = factor(df$mutType, levels = c('Other','Collaterally Neutral', 'Collaterally Sensitive','Collaterally Resistant' ))
  df = df[order(df$mutType), ]
  
  
  JDFE = ggplot(df, aes(x=HomeS, y=NonHomeS, color = mutType)) + 
    geom_point( size = 3.5, alpha = 0.5)+ scale_color_manual(values = c('grey', '#13944a','blue','orange'), drop = 'FALSE') +
    theme_bw()+  xlab("CEF") + ylab(conc)+
    theme(panel.border = element_rect(color = 'black', fill = NA, size = 3), panel.grid.major = element_blank(),
          legend.position = 'none',axis.text.x = element_text(size = 25, color = 'black'), axis.text.y = element_text(size = 25, color = 'black'),panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"), axis.title.x = element_blank(), axis.title.y = element_blank()) +
    geom_hline(yintercept =  0, linetype = 'dashed') + 
    geom_vline(xintercept =  0, linetype = 'dashed') +
    scale_y_continuous(limits = c(-4,1), breaks = c(0, -2, -4)) +
    scale_x_continuous(limits = c(-0.5, 1.5), breaks = c(0,1))
  # green = CN. Orange = CR, Blue = CS, Grey = Other 
  
  ggsave(paste0('JDFE_home_CEF', '_nonHome_', conc, '_90CINonHome.pdf' ),JDFE , width = 3.81, height = 3.48)
  
  
}

# make summary DF
summary_df = data.frame(concID, countBenMuts, countDeltMuts, countNeutralMuts, countisNA, countCRmuts, countCSmuts, countCNmuts)

write.csv(summary_df, file = 'summary_mutEfffects_StifflerData.csv')
write.csv(stifflerData, file = 'fullData_mutEffects_and_types.csv')


## Go through all muts, identify ones with calls
IDs_benInHome = stifflerData$cef_0.2ngmL > 0

geneNames = stifflerData$Mutation_amp1[IDs_benInHome]
call_AMP156 = Amp_156_collateralMutID[IDs_benInHome]
call_AMP2500 = Amp_2500_collateralMutID[IDs_benInHome]
call_AMP39 = Amp_39_collateralMutID[IDs_benInHome]
call_AMP625 = Amp_625_collateralMutID[IDs_benInHome]
p_CRAmp_156 = pVals_CRAmp_156[IDs_benInHome]
p_CRAmp_2500= pVals_CRAmp_2500[IDs_benInHome]
p_CRAmp_39= pVals_CRAmp_39[IDs_benInHome]
p_CRAmp_625= pVals_CRAmp_625[IDs_benInHome]
p_CSAmp_156 = pVals_CSAmp_156[IDs_benInHome]
p_CSAmp_2500= pVals_CSAmp_2500[IDs_benInHome]
p_CSAmp_39= pVals_CSAmp_39[IDs_benInHome]
p_CSAmp_625= pVals_CSAmp_625[IDs_benInHome]

calls_andPvals_stiffler = data.frame(geneNames , call_AMP39, call_AMP156,call_AMP625, call_AMP2500, p_CRAmp_39, p_CRAmp_156, p_CRAmp_625, p_CRAmp_2500, p_CSAmp_39, p_CSAmp_156, p_CSAmp_625, p_CSAmp_2500)

write.csv(calls_andPvals_stiffler, 'calls_andPvals_stiffler.csv')


###############################################################################################################
######################################### F  I  G  U  R  E        2  ##########################################
###############################################################################################################

# (optional) set working directory to where you'd like to save the outputs of this section, example below 
# setwd("~/")

## function for running the evolution sim
# inputs: uVec, Sigma, popsize, U, numiterations, numGens, wantPlot, yplotmin, yplotmax
   # uVec = c(home DFE mean , nonHome DFE mean)
   # Sigma = variance covariance matrix of JDFE
   # popsize = population size
   # U = per genome per generation mutation rate
   # numiterations = number of replicate evolution simulations to run
   # numGens = number of gens per simulation; !!NOTE: must be AT LEAST 1000 to accomodate variation accumulation
   # wantPlot = TRUE or FALSE if want plot of the fitness trajectory in home and non home
   # yplotmin and yplotmax = boundaries of the fitness trajectory plot, require input even if wantPLot is FALSE

# outputs: save graph of fit traj in home/non-home (if wantPLot == TRUE), data frame of Home fitness trajectory slope, NonHome fitness trajectory slope, Home fitness variance slope, NonHome fitness variance slope, covarSlope,r1,r2,D11,D12,D22
## NOTE: only written to accomodate simple 2D gaussian JDFE
simEvo_WrightFisher_savePlots <- function(uVec, Sigma, popsize, U, numiterations, numGens, wantPlot, yplotmin, yplotmax) 
{
  
  uHome = uVec[1]  # mean of home DFE
  uNonHome = uVec[2] # mean non home DFE
  sigHome = Sigma[1,1] # SD of Home DFE
  sigNonHome = Sigma[2,2] # SD of Nonhome DFE
  covar = Sigma[1,2] # covar of JDFE
  
  ################## M A K E     T H E      J D F E
  sValallTheoMuts = mvrnorm(10000, mu =uVec, Sigma = Sigma)
  DFE1 = sValallTheoMuts[,1]
  DFE2 = sValallTheoMuts[,2]
  

  allMeanFitTrajHome = array(0, dim = c(numGens, numiterations))
  allMeanFitTrajNonHome = array(0, dim = c(numGens, numiterations))
  
  for( it in 1:numiterations)
  {
    print(it)
    ID_df  = as.data.frame( array( data = c(1, 0,0,popsize), dim = c(1,4)))
    colnames(ID_df) = c('Strain', 'HomeGR', 'NonHomeGR', 'Count')
    
    
    
    allIDArrayList = list()
    
    
    meanFitnessTrajHome = c()
    meanFitnessTrajNonHome = c()
    
    for(gen in 1:numGens)
    {
      #####    S  E  L  E  C   T  I  O  N    
      #---------------------------------------------------------------------------------------_#
      # select next gen offspring
      meanFitHome = (1/popsize)*sum(ID_df$Count * ID_df$HomeGR)
      meanFitNonHome = (1/popsize)*sum(ID_df$Count * ID_df$NonHomeGR)
      
      probVec = ID_df$Count / popsize + (ID_df$Count / popsize)*(ID_df$HomeGR -meanFitHome)
      # all neg probs = 0 prob
      probVec[probVec<0] = 0
      offspring = rmultinom(1, size = popsize, prob = probVec)
      #   
      
  
      # save important parameters
      allIDArrayList[[gen]] = ID_df
      meanFitnessTrajHome = c(meanFitnessTrajHome, meanFitHome)
      meanFitnessTrajNonHome = c(meanFitnessTrajNonHome, meanFitNonHome)
      
      
      
      # update ID df with new count
      ID_df$Count = offspring 
      
      # remove any extinct lineages
      ID_df = ID_df[!ID_df$Count == 0, ]
      
      # re-do counting so dont have issues calling indexes of parents
      ID_df$Strain = 1:length(ID_df$Strain)
      
      
      
      #### M  U  T  A  T  I  O  N  S
      #---------------------------------------------#
      
      # Sample number of muts that will happen and then choose parent IDs
      M = rpois(1, popsize*U)
      
      if(M > 0 )
      {
        parentIDs = sample(ID_df$Strain, M, prob = (ID_df$Count/popsize), replace = TRUE)
        
        
        # sample M mutation effects from JDFE, usimg Sigma (covariance matrix ) from above fit
        mutEffects = as.array(mvrnorm(M, uVec, Sigma))
        
        if(M ==1) # to correct for a single sample frommvrnorm being read as a vector
        {
          mutEffects = array(mutEffects, dim = c(1,2))
        }
        
        # make M new rows for ID_df, one for each mut (4 cols, ID#, homeFit, respFit (!!!!!which not doing anything with RN!!!!))
        newMutArray = as.data.frame(array(data = 0, dim = c(M, 4)))
        colnames(newMutArray) = c('Strain', 'HomeGR', 'NonHomeGR', 'Count')
        
        # make new mutant ID's ,, just count up from highest current ID
        newMutArray$Strain = (max(ID_df$Strain)+1):(max(ID_df$Strain) + M )
        
        # update new Mut Fitness values by adding sampled fit Effects
        newMutArray$HomeGR = ID_df$HomeGR[parentIDs] + mutEffects[,1]
        newMutArray$NonHomeGR = ID_df$NonHomeGR[parentIDs] + mutEffects[,2]
        
        # add count of ONE to each new mut, subtract offspring from parent
        newMutArray$Count = rep(1, times = M)
        
        
        ####??? faster way to do this
        for(rem in parentIDs)
        {
          ID_df$Count[rem] = ID_df$Count[rem] -1
        }
        
        # collate old Id_df and new muts
        ID_df = rbind(ID_df, newMutArray)
        
      }
      # start again with next gen :) 
    }
    
    allMeanFitTrajHome[,it] = meanFitnessTrajHome
    allMeanFitTrajNonHome[,it] = meanFitnessTrajNonHome
    
  }
  

  time = 1:numGens
  slopeHome = summary(lm(rowMeans(allMeanFitTrajHome)[800:numGens]~time[800:numGens]))$coefficients[2,1]
  slopeNonHome = summary(lm(rowMeans(allMeanFitTrajNonHome)[800:numGens]~time[800:numGens]))$coefficients[2,1]
  
  varHome_eachGen = rowVars(allMeanFitTrajHome)
  varSlopeHome = summary(lm(rowVars(allMeanFitTrajHome)[800:numGens]~time[800:numGens]))$coefficients[2,1]
  
  varNonHome_eachGen = rowVars(allMeanFitTrajNonHome)
  varSlopeNonHome = summary(lm(rowVars(allMeanFitTrajNonHome)[800:numGens]~time[800:numGens]))$coefficients[2,1]
  
  covar_eachGen = c()
  for(covT in 1:numGens)
  {
    covar_eachGen = c(covar_eachGen, cov(allMeanFitTrajHome[covT, ], allMeanFitTrajNonHome[covT, ]))
  }
  
  covarSlope = summary(lm(covar_eachGen[800:numGens]~time[800:numGens]))$coefficients[2,1]
  
  
  
  
  UbAdjust = 1 - pnorm(0, uHome, sigHome)
  homeSvals = sValallTheoMuts[,1]*(sValallTheoMuts[,1]>=0)
  r1 = 2*mean((homeSvals)^2)/ UbAdjust
  r2 = 2*mean((homeSvals)*(sValallTheoMuts[,2]))/ UbAdjust
  D11 = 2*mean((homeSvals)^3)/ UbAdjust
  D12 = 2*mean((homeSvals)^2*(sValallTheoMuts[,2]))/ UbAdjust
  D22 = 2*mean((homeSvals)*(sValallTheoMuts[,2])^2)/ UbAdjust
  

  if(wantPlot == TRUE)
  {
    meanFitTrajHome = rowMeans(allMeanFitTrajHome)
    
    meanFitTrajNonHome = rowMeans(allMeanFitTrajNonHome)
    
    fitness = c(meanFitTrajHome, meanFitTrajNonHome)
    min = c(meanFitTrajHome - sqrt(varHome_eachGen), meanFitTrajNonHome - sqrt(varNonHome_eachGen)) # min of ribbon is -1 SD
    max = c(meanFitTrajHome + sqrt(varHome_eachGen), meanFitTrajNonHome + sqrt(varNonHome_eachGen)) # max of ribbon is +1 SD
    time = rep(1:numGens, times = 2)
    Env = as.factor(c(rep(1, times = numGens), rep(2, times = numGens)))
    WbarDf = data.frame(time, fitness, min, max, Env)
    
    
   plot =( ggplot(WbarDf, aes(x=as.numeric(time),y= fitness, colour = Env )) + geom_line(size = 2)  +
             geom_ribbon(aes(x= time, ymin=min, ymax=max,fill = Env,group = Env), alpha = 0.3, linetype = 0) + 
             scale_color_manual(values=c('chocolate1','dodgerblue2'))+  geom_hline(yintercept =  0, linetype = 'dashed')+
             theme_bw()+theme(panel.border = element_rect(colour = 'black', fill = NA, size = 2), 
                              panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
                              axis.text =element_text(size = 20,family = 'Helvetica'), axis.title = element_blank(), legend.position = 'none')+ 
             xlab("Time")+ ylab("Fitness")+ylim(yplotmin,yplotmax))
   ggsave(paste0('SimTraj_JDFEnum_', JDFEnum, '.pdf'), plot, width = 3.82, height = 3.7)
   
  }
  
 
  output = data.frame(slopeHome, slopeNonHome,varSlopeHome, varSlopeNonHome,  covarSlope,r1,r2,D11,D12,D22)
  return(output)
  
  
  
}


########### Simple 2D Gaussian (main text) ###############


# vars that are changed for text figures
respMeanVec = c(0.08, 0.145, 0,-0.145,-0.08)
corVec = c(-0.8,-0.5, 0, 0.5, 0.8)

# vars that can be changed, but are constant for text figures 
homeMean = -0.05 
homeSd = 0.1
respSd = 0.1
xplotmin = -0.35
xplotmax = 0.15
yplotmin = -0.38
yplotmax = 0.38
numiterations = 100
numGens = 1000
popsize = 10^4
U = 10^-4
Ubadjust = 1 - pnorm(0, homeMean, homeSd)
UB = U*Ubadjust

for(JDFEnum in 1:length(respMeanVec))
{
  print(JDFEnum)
  
  respMean = respMeanVec[JDFEnum]
  corr = corVec[JDFEnum]
  xplotmin = -0.35
  xplotmax = 0.25
  yplotmin = -0.42
  yplotmax = 0.42
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
 
  JDFE_fig  = ggplot(df, aes(x=vals1, y=vals2, z=prob)) + 
    geom_contour_filled(breaks = c(0.45,3,10,100))+
    scale_fill_manual(values = c('#d0cece', '#5b5958', '#333231')) +
    xlim(xplotmin, xplotmax)+ ylim(yplotmin, yplotmax) +
    theme_bw()+theme(panel.border = element_rect(colour = 'black', fill = NA, size = 2), 
                     panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
                     axis.text =element_text(size = 20,family = 'Helvetica'), axis.title = element_blank(), legend.position = 'none')+ 
    geom_hline(yintercept =  0, linetype = 'dashed') + 
    geom_vline(xintercept =  0, linetype = 'dashed')+
    geom_point(aes(weighted.mean(vals1, prob), weighted.mean(vals2, prob)), size = 5, pch = 4, color = 'white')
  
  # Save JDFE 
  ggsave(paste0('JDFE_fig_',JDFEnum,'.pdf'), plot = JDFE_fig, width = 3.82, height = 3.70)
 
  
   # RUN EVOLUTION FUNCTION  #
  # if wantPLot TRUE, will SAVE to working directory graph of mean fitness trajectory in home and non-home; ribbons show +/- 1sd
  # also saves results in data fram (in order):
  # slopeHome, slopeNonHome,covarSlope,r1,r2,D11,D12,D22
  # prints number each iteration so can keep track of progress, is slow for NU > ~100

  uVec = meanVec # from above
  Sigma = sigma # from above
 # output = simEvo_WrightFisher_savePlots(uVec, Sigma, popsize , U , numiterations , numGens , wantPlot =  TRUE, yplotmin = -2, yplotmax = 2)
  
  
  
}




########### More Complex Double Gaussian (Supplement) ###############
setwd("~/")
## function for running the evolution sim
# inputs: uVec, Sigma, popsize, U, numiterations, numGens, wantPlot, yplotmin, yplotmax
# uVec1 = c(home DFE mean , nonHome1 DFE mean)
# uVec2 = c(home DFE mean , nonHome2 DFE mean)
# Sigma1 = variance covariance matrix of 1st gaussian
# Sigma2 = variance covariance matrix of 2nd gaussian
# popsize = population size
# U = per genome per generation mutation rate
# numiterations = number of replicate evolution simulations to run
# numGens = number of gens per simulation; !!NOTE: must be AT LEAST 1000 to accomodate variation accumulation
# wantPlot = TRUE or FALSE if want plot of the fitness trajectory in home and non home
# yplotmin and yplotmax = boundaries of the fitness trajectory plot, require input even if wantPLot is FALSE

# outputs: SAVES graph of fit traj in home/non-home (if wantPLot == TRUE), data frame of Home fitness trajectory slope, NonHome fitness trajectory slope, Home fitness variance slope, NonHome fitness variance slope, covarSlope,r1,r2,D11,D12,D22
## NOTE: only written to accomodate double 2D gaussian JDFE as used in this section
simEvo_doubleGauss <- function(uVec1, Sigma1, uVec2, Sigma2, popsize, U, numiterations, numGens, wantPlot, yplotmin, yplotmax) 
{
  
  # 1st Gaussian Params
  uHome1 = uVec1[1]  # mean of home 1 DFE
  uNonHome1 = uVec1[2] # mean non home 1 DFE
  sigHome1 = Sigma1[1,1] # SD of Home 1 DFE
  sigNonHome1 = Sigma1[2,2] # SD of Nonhome 1 DFE
  covar1 = Sigma1[1,2] # covar of JDFE 1 
  
  # 2nd Gaussian Params
  uHome2 = uVec2[1]  # mean of home 2 DFE
  uNonHome2 = uVec2[2] # mean non home 2 DFE
  sigHome2 = Sigma2[1,1] # SD of Home 2 DFE
  sigNonHome2 = Sigma2[2,2] # SD of Nonhome 2 DFE
  covar2 = Sigma2[1,2] # covar of JDFE 12
  
  
  
  ################## M A K E     T H E      J D F E
  # do this to help calc r and d parameters 
  ## 1st env ##
  meanVec_1 <- c(uHome1, uNonHome1)
  sigma_1 <- matrix(c(sigHome1,covar1,covar1,sigNonHome1), nrow=2)
  sValallTheoMuts_1 = mvrnorm(10000, mu =meanVec_1, Sigma = sigma_1)
  DFE1_1 = sValallTheoMuts_1[,1]
  DFE2_1 = sValallTheoMuts_1[,2]
  
  
  ## 2nd env ##
  meanVec_2 <- c(uHome2, uNonHome2)
  sigma_2 <- matrix(c(sigHome2,covar2,covar2,sigNonHome2), nrow=2)
  sValallTheoMuts_2 = mvrnorm(10000, mu =meanVec_2, Sigma = sigma_2)
  DFE1_2 = sValallTheoMuts_2[,1]
  DFE2_2 = sValallTheoMuts_2[,2]
  
  
  ## together ##
  sValallTheoMuts = array(0, dim= c(20000, 2))
  sValallTheoMuts[,1] = c(DFE1_1, DFE1_2)
  sValallTheoMuts[,2] = c(DFE2_1, DFE2_2)
  
  
  allMeanFitTrajHome = array(0, dim = c(numGens, numiterations))
  allMeanFitTrajNonHome = array(0, dim = c(numGens, numiterations))
  
  for( it in 1:numiterations)
  {
    print(it)
    ID_df  = as.data.frame( array( data = c(1, uHome+0.001, uNonHome+0.001,popsize), dim = c(1,4)))
    colnames(ID_df) = c('Strain', 'HomeGR', 'NonHomeGR', 'Count')
    
    
    
    allIDArrayList = list()
    
    
    meanFitnessTrajHome = c()
    meanFitnessTrajNonHome = c()
    
    for(gen in 1:numGens)
    {
      #####    S  E  L  E  C   T  I  O  N    
      #---------------------------------------------------------------------------------------_#
      # select next gen offspring
      meanFitHome = (1/popsize)*sum(ID_df$Count * ID_df$HomeGR)
      meanFitNonHome = (1/popsize)*sum(ID_df$Count * ID_df$NonHomeGR)
      
      probVec = ID_df$Count / popsize + (ID_df$Count / popsize)*(ID_df$HomeGR -meanFitHome)
      # all neg probs = 0 prob
      probVec[probVec<0] = 0
      offspring = rmultinom(1, size = popsize, prob = probVec)
      #   
      
      
      # save important parameters
      allIDArrayList[[gen]] = ID_df
      meanFitnessTrajHome = c(meanFitnessTrajHome, meanFitHome)
      meanFitnessTrajNonHome = c(meanFitnessTrajNonHome, meanFitNonHome)
      
      
      
      # update ID df with new count
      ID_df$Count = offspring 
      
      # remove any extinct lineages
      ID_df = ID_df[!ID_df$Count == 0, ]
      
      # re-do counting so dont have issues calling indexes of parents
      ID_df$Strain = 1:length(ID_df$Strain)
      
      
      
      #### M  U  T  A  T  I  O  N  S
      #---------------------------------------------#
      
      # Sample number of muts that will happen and then choose parent IDs
      M = rpois(1, popsize*U)
      
      if(M > 0 )
      {
        parentIDs = sample(ID_df$Strain, M, prob = (ID_df$Count/popsize), replace = TRUE)
        
        
        ### sample M mutation effects from JDFE, 
        # first determine which gaussian each mut effect sampled from
        whichGuass = sample(c(1,2), M, replace = TRUE, prob = c(1,1))
        howMany1 = sum(whichGuass == 1)
        howMany2 = sum(whichGuass == 2)
        
        mutEffects = array(data = NA, dim = c(M, 2))
        if(howMany1 != 0)
        {
          mutEffects1 = as.array(mvrnorm(howMany1, uVec1, Sigma1))
          mutEffects[(whichGuass == 1), ] = mutEffects1
        }
        
        if(howMany2 != 0)
        {
          mutEffects2 = as.array(mvrnorm(howMany2, uVec2, Sigma2))
          mutEffects[(whichGuass == 2), ] = mutEffects2
        }
        
        
        if(M ==1) # to correct for a single sample frommvrnorm being read as a vector
        {
          mutEffects = array(mutEffects, dim = c(1,2))
        }
        
        # make M new rows for ID_df, one for each mut (4 cols, ID#, homeFit, respFit (!!!!!which not doing anything with RN!!!!))
        newMutArray = as.data.frame(array(data = 0, dim = c(M, 4)))
        colnames(newMutArray) = c('Strain', 'HomeGR', 'NonHomeGR', 'Count')
        
        # make new mutant ID's ,, just count up from highest current ID
        newMutArray$Strain = (max(ID_df$Strain)+1):(max(ID_df$Strain) + M )
        
        # update new Mut Fitness values by adding sampled fit Effects
        newMutArray$HomeGR = ID_df$HomeGR[parentIDs] + mutEffects[,1]
        newMutArray$NonHomeGR = ID_df$NonHomeGR[parentIDs] + mutEffects[,2]
        
        # add count of ONE to each new mut, subtract offspring from parent
        newMutArray$Count = rep(1, times = M)
        
        
        ####??? faster way to do this
        for(rem in parentIDs)
        {
          ID_df$Count[rem] = ID_df$Count[rem] -1
        }
        
        # collate old Id_df and new muts
        ID_df = rbind(ID_df, newMutArray)
        
      }
      # start again with next gen :) 
    }
    
    allMeanFitTrajHome[,it] = meanFitnessTrajHome
    allMeanFitTrajNonHome[,it] = meanFitnessTrajNonHome
    
  }
  
  
  time = 1:numGens
  slopeHome = summary(lm(rowMeans(allMeanFitTrajHome)[800:numGens]~time[800:numGens]))$coefficients[2,1]
  slopeNonHome = summary(lm(rowMeans(allMeanFitTrajNonHome)[800:numGens]~time[800:numGens]))$coefficients[2,1]
  
  varHome_eachGen = rowVars(allMeanFitTrajHome)
  varSlopeHome = summary(lm(rowVars(allMeanFitTrajHome)[800:numGens]~time[800:numGens]))$coefficients[2,1]
  
  varNonHome_eachGen = rowVars(allMeanFitTrajNonHome)
  varSlopeNonHome = summary(lm(rowVars(allMeanFitTrajNonHome)[800:numGens]~time[800:numGens]))$coefficients[2,1]
  
  covar_eachGen = c()
  for(covT in 1:numGens)
  {
    covar_eachGen = c(covar_eachGen, cov(allMeanFitTrajHome[covT, ], allMeanFitTrajNonHome[covT, ]))
  }
  
  covarSlope = summary(lm(covar_eachGen[800:numGens]~time[800:numGens]))$coefficients[2,1]
  
  
  
  
  UbAdjust = mean(1 - pnorm(0, uHome1, sigHome1), 1 - pnorm(0, uHome2, sigHome2))
  homeSvals = sValallTheoMuts[,1]*(sValallTheoMuts[,1]>=0)
  r1 = 2*mean((homeSvals)^2)/UbAdjust 
  r2 = 2*mean((homeSvals)*(sValallTheoMuts[,2]))/UbAdjust 
  D11 = 2*mean((homeSvals)^3)/UbAdjust 
  D12 = 2*mean((homeSvals)^2*(sValallTheoMuts[,2]))/UbAdjust 
  D22 = 2*mean((homeSvals)*(sValallTheoMuts[,2])^2)/UbAdjust 
  
  
  if(wantPlot == TRUE)
  {
    meanFitTrajHome = rowMeans(allMeanFitTrajHome)
    
    meanFitTrajNonHome = rowMeans(allMeanFitTrajNonHome)
    
    fitness = c(meanFitTrajHome, meanFitTrajNonHome)
    min = c(meanFitTrajHome - sqrt(varHome_eachGen), meanFitTrajNonHome - sqrt(varNonHome_eachGen)) # min of ribbon is -1 SD
    max = c(meanFitTrajHome + sqrt(varHome_eachGen), meanFitTrajNonHome + sqrt(varNonHome_eachGen)) # max of ribbon is +1 SD
    time = rep(1:numGens, times = 2)
    Env = as.factor(c(rep(1, times = numGens), rep(2, times = numGens)))
    WbarDf = data.frame(time, fitness, min, max, Env)
    
    
    plot =( ggplot(WbarDf, aes(x=as.numeric(time),y= fitness, colour = Env )) + geom_line(size = 2)  +
              geom_ribbon(aes(x= time, ymin=min, ymax=max,fill = Env,group = Env), alpha = 0.3, linetype = 0) + 
              scale_color_manual(values=c('chocolate1','dodgerblue2'))+  geom_hline(yintercept =  0, linetype = 'dashed')+
              theme_bw()+theme(panel.border = element_rect(colour = 'black', fill = NA, size = 2), 
                               panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
                               axis.text =element_text(size = 20,family = 'Helvetica'), axis.title = element_blank(), legend.position = 'none')+ 
              xlab("Time")+ ylab("Fitness")+ylim(yplotmin,yplotmax))
    ggsave(paste0('SimTraj_DoubleGauss_JDFEnum_', JDFEnum, '.pdf'), plot, width = 3.82, height = 3.7)
  }
  
  
  output = data.frame(slopeHome, slopeNonHome,varSlopeHome, varSlopeNonHome,  covarSlope,r1,r2,D11,D12,D22)
  return(output)
  
  
  
}

#### vals that are changed
## !!! NOTE:: For plotting technicality, much have uNonHome1 > uNonHome2
uNonHome1Vec = c(0.1,0.5, 0.17, 0.5) # means of DFE of Nonhome for 1st gaussian
uNonHome2Vec = c(-0.1,-0.5, -0.5, -0.17) # means of DFE of Nonhome for 2nd gaussian
corr1 = 0 # corr of 1st gaussian for all
corr2 = 0 # corr of 2nd gaussian for all
x1plotmin = -1  
x1plotmax = 1 
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
uHome1 = 0.4 # mean of DFE of home for 1st gaussian
homeSd_1 = 0.1 # SD of DFE of home for 1st gaussian
uHome2 = 0.4 # mean of DFE of home for 2nd gaussian
homeSd_2 = 0.1 # SD of DFE of home for 2nd gaussian

nonHomeSd_1 = 0.1 # SD of DFE of Nonhome for 1st gaussian
nonHomeSd_2 = 0.1 # SD of DFE of Nonhome for 2nd gaussian

numiterations = 100
numGens = 1000
popsize = 10^5
U = 10^-4


for(JDFEnum in 1:length(uNonHome1Vec))
{
  uNonHome1 = uNonHome1Vec[JDFEnum]
  uNonHome2 = uNonHome2Vec[JDFEnum]
  
  
  # Make the JDFE
  #-------------#
  homeVar_1 = homeSd_1^2
  homeVar_2 = homeSd_2^2
  nonHomeVar_1 = nonHomeSd_1^2
  nonHomeVar_2 = nonHomeSd_2^2
  covar1 = corr1*homeSd_1*nonHomeSd_1
  covar2 = corr2*homeSd_2*nonHomeSd_2
  
  ## NonHome 1
  uVec1 <- c(uHome1, uNonHome1)
  sigma1 <- matrix(c(homeVar_1,covar1,covar1,nonHomeVar_1), nrow=2)
  data.grid <- expand.grid(vals1 = seq(x1plotmin,x1plotmax, length.out=200), vals2 = seq(y1plotmin, y1plotmax, length.out=200))
  df1 <- cbind(data.grid, prob = mvtnorm::dmvnorm(data.grid, mean = uVec1, sigma = sigma1))
  
  ## NonHome 2
  uVec2 <-  c(uHome2, uNonHome2)
  sigma2 <- matrix(c(homeVar_2,covar2,covar2,nonHomeVar_2), nrow=2)
  data.grid <- expand.grid(vals1 = seq(x2plotmin,x2plotmax, length.out=200), vals2 = seq(y2plotmin, y2plotmax, length.out=200))
  df2 <- cbind(data.grid, prob = mvtnorm::dmvnorm(data.grid, mean = uVec2, sigma = sigma2))
  
  ## together
  df = data.frame(c(df1$vals1, df2$vals1), c(df1$vals2, df2$vals2),c(df1$prob, df2$prob))
  colnames(df) = c('vals1', 'vals2', 'prob')
  
  ## plot the JDFE
  JDFE_fig = (ggplot(df, aes(x=vals1, y=vals2, z=prob)) + 
                geom_contour_filled(breaks = c(4,8, 12,100))+
                scale_fill_manual(values = c('#d0cece', '#5b5958', '#333231')) +
                xlim(xTOTplotmin, xTOTplotmax)+ ylim(yTOTplotmin, yTOTplotmax) +
                theme_bw()+theme(panel.border = element_rect(colour = 'black', fill = NA, size = 2), 
                                 panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
                                 axis.text =element_text(size = 20,family = 'Helvetica'), axis.title = element_blank(), legend.position = 'none')+ 
                geom_hline(yintercept =  0, linetype = 'dashed') + 
                geom_vline(xintercept =  0, linetype = 'dashed')+
                geom_point(aes(weighted.mean(vals1, prob), weighted.mean(vals2, prob)), pch = 4, size = 3))
  
  
  # Save JDFE 
  ggsave(paste0('JDFE_DoubleGauss_fig_',JDFEnum,'.pdf'), plot = JDFE_fig, width = 3.82, height = 3.70)
  
  
  
  #### EXAMPLE OF EVOLUTION FUNCTION CALL ###
  # outputs graph of mean fitness trajectory in home and non-home; ribbons show +/- 1sd
  # also saves results in DF slopeHome, slopeNonHome,varSlopeHome, varSlopeNonHome,  covarSlope,r1,r2,D11,D12,D22
  
  #from plotting section above
  uVec1 = uVec1
  uVec2 = uVec2
  Sigma1 = sigma1
  Sigma2 = sigma2 
  
  output = simEvo_doubleGauss(uVec1, Sigma1, uVec2, Sigma2, popsize = popsize, U = U, numiterations= numiterations, numGens= numGens, wantPlot= TRUE, yplotmin = -50, yplotmax = 50) 
  
}



#### END FIGURE 2 ####





######################################### F  I  G  U  R  E        3  ##########################################
                                       # non epistatic JDFE  vs sim #


# (optional) set working directory to where you'd like to save the outputs of this section, example below 
# setwd("~/Documents/JDFE Project/Figures/Figure gaussian JDFES")


###################   SSWM SECTION ###################
## function for running the evolution sim
# inputs: uVec, Sigma, popsize, U, numiterations, numGens, wantPlot, yplotmin, yplotmax
# uVec = c(home DFE mean , nonHome DFE mean)
# Sigma = variance covariance matrix of JDFE
# popsize = population size
# U = per genome per generation mutation rate
# numiterations = number of replicate evolution simulations to run
# numGens = number of gens per simulation; !!NOTE: must be AT LEAST 1000 to accomodate variation accumulation
# wantPlot = TRUE or FALSE if want plot of the fitness trajectory in home and non home

# outputs: graph of fit traj in home/non-home (if wantPLot == TRUE), data frame of Home fitness trajectory slope, NonHome fitness trajectory slope, Home fitness variance slope, NonHome fitness variance slope, covarSlope,r1,r2,D11,D12,D22
## NOTE: only written to accomodate simple 2D gaussian JDFE
simEvo_SSWM <- function(uVec, Sigma, popsize, U, numiterations, numGens, wantPlot ) # dont have call to ustomize yplot axes
{
  
  uHome = uVec[1]  # mean of home DFE
  uNonHome = uVec[2] # mean non home DFE
  sigHome = sqrt(Sigma[1,1]) # SD of Home DFE
  sigNonHome = sqrt(Sigma[2,2]) # SD of Nonhome DFE
  covar = Sigma[1,2] # covar of JDFE
  
  ################## M A K E     T H E      J D F E
  sValallTheoMuts = mvrnorm(100000, mu =uVec, Sigma = Sigma)
  DFE1 = sValallTheoMuts[,1]
  DFE2 = sValallTheoMuts[,2]
  UbAdjust = 1 - pnorm(0, uHome, sigHome)
  
  
  
  
  ### Find R and D-values ####
  
  homeSvals = sValallTheoMuts[,1]*(sValallTheoMuts[,1]>=0)
  r1 = 2*mean((homeSvals)^2)/UbAdjust 
  r2 = 2*mean((homeSvals)*(sValallTheoMuts[,2]))/UbAdjust 
  D11 = 2*mean((homeSvals)^3)/UbAdjust 
  D12 = 2*mean((homeSvals)^2*(sValallTheoMuts[,2]))/UbAdjust 
  D22 = 2*mean((homeSvals)*(sValallTheoMuts[,2])^2)/UbAdjust 
  
  
  
  
  numDrugs = length(uVec) # how many envs
  FitnessTrajectoryHome = array(0, dim= c(numGens, numiterations))
  FitnessTrajectoryNonHome= array(0, dim = c(numGens, numiterations))
  
  for(it in 1:numiterations)
  {
    wBarHome = 0
    wBarResp = 0
   
    
    ### make vec of times in which get mutations
    end = 0
    current = 0
    mutTimesVec = c()
    Ub = U*(1 - pnorm(0, uHome, sigHome)) 
    while(end < numGens)
    {
      current = end
      nextMutTime = ceiling(rexp(1, rate = popsize*Ub))
      end = current + nextMutTime
      mutTimesVec = c(mutTimesVec, end)
      
    }
    mutTimesVec = mutTimesVec[-length(mutTimesVec)] # remoe the last entry, because is greater than the numGens have
    
    for(gen in 1:numGens)
    {
      if(sum(gen==mutTimesVec)!=0) # if at a generation with a mutation, sample a mutant
      {
        # get home S-val
        homeS = -1 
        while(homeS < 0 )
        {
          samp =  mvrnorm(1, uVec, Sigma)
          homeS = samp[1]
          respS = samp[2]
        }
          
        
        
        
        prob = 2*homeS # probability of fixation ~2s
        
        # adjust for impossible probabilities 
        if(prob > 1)
        {
          prob = 1
        }else if(prob < 0)
        {
          prob = 0
        }
        
        # 'flip coin' to see if mutation fixes 
        flip = sample(c(0,1), 1, prob = c(1-prob,prob))
        
        if(flip == 1) ## if mut does fix
        {
          respS =samp[2] # nonHome s value
          wBarHome = (homeS)+wBarHome # additive fitness effect to home Mean fitness
          wBarResp = (respS)+wBarResp # additive fitness effect to Nonhome Mean fitness
          
          FitnessTrajectoryHome[gen, it] = wBarHome
          FitnessTrajectoryNonHome[gen, it] = wBarResp
          
        }else{ # if mut does not fix 
          # just propogate last Wbar in both envs to this timepoint
          FitnessTrajectoryHome[gen,it] = wBarHome
          FitnessTrajectoryNonHome[gen,it] = wBarResp
          
        }
        
        
        
        
      }else{ # at a gen where dont get mut
        
        # just propogate last Wbar in both envs to this timepoint
        FitnessTrajectoryHome[gen,it] = wBarHome
        FitnessTrajectoryNonHome[gen,it] = wBarResp
        
      }
    }# end gen loop
   
  }#end iteration loop 
  
  time = 1:numGens
  meanFitTrajHome = rowMeans(FitnessTrajectoryHome)
  meanFitTrajNonHome = rowMeans(FitnessTrajectoryNonHome)
  slopeHome = summary(lm(meanFitTrajHome[500:numGens]~time[500:numGens]))$coefficients[2,1] # SSWM needs very little initial time to build variation but start LM at gen 500 to be safe
  slopeNonHome = summary(lm(meanFitTrajNonHome[500:numGens]~time[500:numGens]))$coefficients[2,1]
  
  
  varHome_eachGen = rowVars(FitnessTrajectoryHome)
  varSlopeHome = summary(lm(varHome_eachGen[500:numGens]~time[500:numGens]))$coefficients[2,1]
  
  varNonHome_eachGen = rowVars(FitnessTrajectoryNonHome)
  varSlopeNonHome = summary(lm(varNonHome_eachGen[500:numGens]~time[500:numGens]))$coefficients[2,1]
  
  covar_eachGen = c()
  for(covT in 1:numGens)
  {
    homeFitsThisGen_allits = FitnessTrajectoryHome[covT, ]
    nonHomeFitsThisGen_allits = FitnessTrajectoryNonHome[covT, ]
    covar_eachGen = c(covar_eachGen, cov(homeFitsThisGen_allits, nonHomeFitsThisGen_allits))
  }
  
  covarSlope = summary(lm(covar_eachGen[500:numGens]~time[500:numGens]))$coefficients[2,1]
  
  
  
  if(wantPlot == TRUE)
  {

    fitness = c(meanFitTrajHome, meanFitTrajNonHome)
    min = c(meanFitTrajHome - sqrt(varHome_eachGen), meanFitTrajNonHome - sqrt(varNonHome_eachGen)) # min of ribbon is -1 SD
    max = c(meanFitTrajHome + sqrt(varHome_eachGen), meanFitTrajNonHome + sqrt(varNonHome_eachGen)) # max of ribbon is +1 SD
    time = rep(1:numGens, times = 2)
    Env = as.factor(c(rep(1, times = numGens), rep(2, times = numGens)))
    WbarDf = data.frame(time, fitness, min, max, Env)
    
    
    print( ggplot(WbarDf, aes(x=as.numeric(time),y= fitness, colour = Env )) + geom_line()  +
             geom_ribbon(aes(x= time, ymin=min, ymax=max,fill = Env,group = Env), alpha = 0.3, linetype = 0) + 
             scale_color_manual(values=c('chocolate1','dodgerblue2'))+  geom_hline(yintercept =  0, linetype = 'dashed')+
             theme_bw()+theme(panel.border = element_rect(colour = 'black', fill = NA, size = 2), 
                              panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
                              axis.text =element_text(size = 20,family = 'Helvetica'), axis.title = element_blank(), legend.position = 'none')+ 
             xlab("Time")+ ylab("Fitness")) 
  }
  
  
  output = list(meanFitTrajHome, meanFitTrajNonHome, varHome_eachGen, varNonHome_eachGen, covar_eachGen, data.frame(slopeHome, slopeNonHome,varSlopeHome, varSlopeNonHome,  covarSlope,r1,r2,D11,D12,D22))
  return(output)
}


## have GAUSSIAN JDFE with home parameters
sigHome = 0.01 #home DFE SD 
uHome = - 0.001 # mean home DFE

# Non Home Parameters
sdVec = seq(0.0001, 0.01, length.out =5)
m2Vec = seq(0.0001, 0.01, length.out = 5)
corVec = seq(-0.9, 0.9, length.out = 5)
numiterations = 300
numDrugs = 2
numGens = 1000
wantPlot = FALSE # dont want sim evo to make plot for all 125 JDFEs

# Empty vectors for saving results
meanHomeSlopes = c()
meanNonHomeSlopes = c()
meanVarSlopeHome =  c()
meanVarSlopeNonHome = c()
meanCoVarSlope =  c()

r1Vals = c()
r2Vals = c()
D11Vals = c()
D22Vals= c()
D12Vals = c()

allM2 = c()
allsd = c()
allCor = c()

popsize = 10^3
N = popsize
U = 10^-4
Ub = U*(1 - pnorm(0, uHome, sigHome)) # prob of beneficial mut (as definined in theory)

count = 0
for(sd in sdVec) # non Home DFE SD
{
  for(mu2 in m2Vec) # non Home DFE mean
  {
    for(c0 in corVec) # JDFE correlation
    {
      count = count+1
      print(count)
      

      covar = c0*sqrt(sigHome^2 * sd^2)
      uVec = c(uHome, mu2)
      Sigma =  matrix(c(sigHome^2,covar,covar,sd^2), nrow=2)
     
      
      
      # save vars for this it
      allM2 = c(allM2, mu2)
      allsd = c(allsd, sd)
      allCor = c(allCor, c0)
      UbAdjust = 1 - pnorm(0, uHome, sigHome)
      
      
      output = simEvo_SSWM(uVec, Sigma, popsize, U, numiterations ,numGens , wantPlot = FALSE )
          # list(meanFitTrajHome, meanFitTrajNonHome, varHome_eachGen, varNonHome_eachGen, covar_eachGen, data.frame(slopeHome, slopeNonHome,varSlopeHome, varSlopeNonHome,  covarSlope,r1,r2,D11,D12,D22))

      outputDF  = output[[6]] 
      
      
      
      # save results
      meanHomeSlopes = c(meanHomeSlopes, outputDF$slopeHome)
      meanNonHomeSlopes = c(meanNonHomeSlopes,outputDF$slopeNonHome)
      meanVarSlopeHome = c(meanVarSlopeHome, outputDF$varSlopeNonHome)
      meanVarSlopeNonHome = c(meanVarSlopeNonHome, outputDF$varSlopeNonHome)
      meanCoVarSlope = c(meanCoVarSlope, outputDF$covarSlope)
      
      r1Vals = c(r1Vals, outputDF$r1)
      r2Vals = c(r2Vals,outputDF$r2 )
      D11Vals = c(D11Vals, outputDF$D11)
      D22Vals= c(D22Vals, outputDF$D22)
      D12Vals = c(D12Vals, outputDF$D12)
      
    }
  }
}


## PLOT RESULTS 
color = 'dodgerblue4' # for sswm plots use this color 


# adjust all sim results by 2NUb because theory is scaled by this factor
scaled_meanNonHomeSlopes= meanNonHomeSlopes/(N*Ub)
scaled_meanHomeSlopes = meanHomeSlopes/(N*Ub)
scaled_meanVarSlopeHome = meanVarSlopeHome/(N*Ub)
scaled_meanVarSlopeNonHome = meanVarSlopeNonHome/(N*Ub)
scaled_meanCoVarSlope = meanCoVarSlope/(N*Ub)



data_thisNu_SSWM = data.frame(r2Vals,D22Vals,D12Vals,scaled_meanHomeSlopes, scaled_meanNonHomeSlopes, scaled_meanVarSlopeHome, scaled_meanVarSlopeNonHome, scaled_meanCoVarSlope, allM2, allsd, allCor,rep( uHome, times = length(allsd)), rep(sigHome, times = length(allsd)), r1Vals)
colnames(data_thisNu_SSWM) = c('r2', 'D22', 'D12', 'scaled_slopeHome', 'scaled_slopeNonHome', 'scaled_varSlopeHome', 'scaled_varSlopeNonHome', 'scaled_covarSlope', 'uDFE_nonHome', 'SD_DFE_nonHome', 'corJDFE', 'uDFE_Home', 'SD_DFE_Home', 'r1')
assign(paste("data_thisNu_",U*popsize,"_SSWM", sep = ""), data_thisNu_SSWM) # have saved in global environmnet with unique name

## save CSV of data in working directory
  write.csv(data_thisNu_SSWM, paste("data_thisNu_",U*popsize,"_SSWM_numGens_", numGens, "_numIts_", numiterations, ".csv", sep = ""))




## PLOT R2 vs MEAN SLOPE

# Note, breaks are chosen manually for the data range here
r2_SSWM_v_sim_Linear =  ggplot(data_thisNu_SSWM, aes(x=r2,y=scaled_slopeNonHome)) + 
  geom_hline(yintercept = 0 , color = 'grey', size = 1.5) +geom_vline(xintercept = 0, color = 'grey', size= 1.5) +
  geom_point(col = color, size = 7, alpha = 0.6) +
  geom_smooth(method='lm', formula= y~x, col = color, lwd = 0.8) + geom_abline(size  = 1.5, intercept = 0, slope = 1, linetype = 'dashed', col = 'grey')+
  theme_bw()+theme(panel.border = element_rect(colour = 'black', fill = NA, size = 1.5), 
                   panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
                   axis.text.x =element_text(size = 35,family = 'Helvetica', color = c('black')),axis.text.y =element_text(size = 35,family = 'Helvetica', color = c('black')), axis.title = element_blank()) +
  scale_y_continuous(label=scientific_10, breaks = c(-5*10^-5,0*10^-5, 5*10^-5, 10*10^-5)) + scale_x_continuous(label=scientific_10, breaks = c(0*10^-5, 5*10^-5, 10*10^-5, -5*10^-5))


#linear  model for above data
LM_SSWM_r2_vs_SimNonHomeSlope = summary(lm(scaled_meanNonHomeSlopes~ r2Vals))
rSq_SSWM_LM_r2_v_SimSlopeNonHome = LM_SSWM_r2_vs_SimNonHomeSlope$r.squared
Coeffs_SSWM_LM_r2_v_SimSlopeNonHome  = LM_SSWM_r2_vs_SimNonHomeSlope$coefficients


# save PDF of the graph
ggsave(paste0("r2_SSWM_v_sim_Linear_NU_",U*popsize,"_SSWM_numGens_", numGens, "_numIts_", numiterations, ".pdf"), r2_SSWM_v_sim_Linear, width = 5.26, height = 3.82)


## D22 vs SLOPE OF VARIANCE TRAJECTORY, on log scale because goes across several orders of magnitude 
D22_SSWM_v_sim_LogScale =  ggplot(data_thisNu_SSWM[-c(1:5), ], aes(x=D22,y=scaled_varSlopeNonHome)) +   geom_point(col = color, size = 7, alpha = 0.6) +
  geom_smooth(method='lm', formula= y~x, col = color, lwd = 0.8) + geom_abline(size  = 1.5, intercept = 0, slope = 1, linetype = 'dashed', col = 'grey')+
  theme_bw()+theme(panel.border = element_rect(colour = 'black', fill = NA, size = 1.5), 
                   panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
                   axis.text.x =element_text(size = 35,family = 'Helvetica', color = c('black')),axis.text.y =element_text(size = 35,family = 'Helvetica', color = c('black')), axis.title = element_blank()) +
  scale_y_continuous(trans = 'log10' ,breaks =  c(10^-8, 10^-7, 10^-6),
                     labels = scientific_10) + scale_x_continuous(trans = 'log10' ,breaks = c(10^-8, 10^-7, 10^-6),
                                                                  labels = scientific_10)  

LM_SSWM_D22_vs_SimNonHomeVarSlope = summary(lm(scaled_meanVarSlopeNonHome~ D22Vals))
rSq_SSWM_LM_D22_v_SimNonHomeVarSlope = LM_SSWM_D22_vs_SimNonHomeVarSlope$r.squared
Coeffs_SSWM_LM_D22_v_SimNonHomeVarSlope = LM_SSWM_D22_vs_SimNonHomeVarSlope$coefficients


# save PDF of the graph
ggsave(paste0("D22_SSWM_v_sim_LogScale_NU_",U*popsize,"_SSWM_numGens_", numGens, "_numIts_", numiterations, ".pdf"), D22_SSWM_v_sim_LogScale, width = 5.26, height = 3.82)






### Save linear model Rsq, coefficients and P-vals
rSqVals_SSWM_LMS = data.frame(rSq_SSWM_LM_r2_v_SimSlopeNonHome, rSq_SSWM_LM_D22_v_SimNonHomeVarSlope)
colnames(rSqVals_SSWM_LMS) = c('R2_v_SimNonHomeSlope', 'D22_v_SimNonHomeVarSlope')

coefficients_SSWM_LMS = rbind(Coeffs_SSWM_LM_r2_v_SimSlopeNonHome, Coeffs_SSWM_LM_D22_v_SimNonHomeVarSlope)

write.csv(rSqVals_SSWM_LMS,  paste0('rSqVals_SSWM_LMS_numGens',numGens, '_numIts_', numiterations ,'.csv'))
write.csv(coefficients_SSWM_LMS, paste0('coefficients_SSWM_LMS_numGens',numGens, '_numIts_', numiterations ,'.csv'))






################  WRIGHT  FISHER   SECTION   ###################

#### Functions ####

# for calculating Clonal Interference R and D values, based on Good et al PNAS

w <- function(x) #  background function from Good et al PNAS ,, depends on xc and v saved in environment
{
  result_w = numeric(length(x))
  
  result_w[x>0 & x<=xc] = xc*exp((x[x>0 & x<=xc]^2 - xc^2)/(2*v)) 
  result_w[x>xc] = x[x>xc] 
  return(result_w)
}

f <- function(x) # from Good et al PNAS ,, depends on  v saved in environment
{
  result_f = 2*(1/sqrt(2*pi*v))*exp(-x^2/(2*v))
  return(result_f)
}

fixP_num <- function(sVec)  #numerical implementation of Good et al PNAS fixation probability,, depends on xc, v, popsize and U being in global env 
{
  allfixP = c()
  for(sVal in sVec)
  {
    if(sVal < xc)
    {
      integrand = function(x)
      {
        return(w(x)*f(x-sVal))
      }
      
      x = seq(0,0.5, length.out = 1000)
      result = sum(integrand(x)*(x[2]-x[1])) 
      
      allfixP = c(allfixP, result)
    }else{
      result = sVal # if s > xc then fix prob is just = s
      allfixP = c(allfixP, result)
    }
    
  }
  
  
  return(allfixP)
}



## function for running the evolution sim
# inputs: uVec, Sigma, popsize, U, numiterations, numGens, wantPlot, yplotmin, yplotmax
# uVec = c(home DFE mean , nonHome DFE mean)
# Sigma = variance covariance matrix of JDFE
# popsize = population size
# U = per genome per generation mutation rate
# numiterations = number of replicate evolution simulations to run
# numGens = number of gens per simulation; !!NOTE: must be AT LEAST 1000 to accomodate variation accumulation
# wantPlot = TRUE or FALSE if want plot of the fitness trajectory in home and non home
# yplotmin and yplotmax = boundaries of the fitness trajectory plot, require input even if wantPLot is FALSE

# outputs: graph of fit traj in home/non-home (if wantPLot == TRUE),  
   #output = list(allMeanFitTrajNonHome, allMeanFitTrajHome, meanFitTrajHome, meanFitTrajNonHome, varHome_eachGen, varNonHome_eachGen, covar_eachGen, 
             # data.frame(slopeHome, slopeNonHome,varSlopeHome, varSlopeNonHome,  covarSlope,r1,r2,D11,D12,D22))

## NOTE: only written to accomodate simple 2D gaussian JDFE
simEvo_WrightFisher <- function(uVec, Sigma, popsize, U, numiterations, numGens, wantPlot) 
{
  
  uHome = uVec[1]  # mean of home DFE
  uNonHome = uVec[2] # mean non home DFE
  sigHome = sqrt(Sigma[1,1]) # SD of Home DFE
  sigNonHome = sqrt(Sigma[2,2]) # SD of Nonhome DFE
  covar = Sigma[1,2] # covar of JDFE
  UbAdjust = 1 - pnorm(0, uHome, sigHome)
  ################## M A K E     T H E      J D F E
  sValallTheoMuts = mvrnorm(100000, mu =uVec, Sigma = Sigma)
  DFE1 = sValallTheoMuts[,1]
  DFE2 = sValallTheoMuts[,2]
  
  
  allMeanFitTrajHome = array(0, dim = c(numGens, numiterations))
  allMeanFitTrajNonHome = array(0, dim = c(numGens, numiterations))
  
  for( it in 1:numiterations)
  {
    #print(it)
    ID_df  = as.data.frame( array( data = c(1, uHome, uHome,popsize), dim = c(1,4)))
    colnames(ID_df) = c('Strain', 'HomeGR', 'NonHomeGR', 'Count')
    
    
    
    allIDArrayList = list()
    
    
    meanFitnessTrajHome = c()
    meanFitnessTrajNonHome = c()
    
    for(gen in 1:numGens)
    {
      #####    S  E  L  E  C   T  I  O  N    
      #---------------------------------------------------------------------------------------_#
      # select next gen offspring
      meanFitHome = (1/popsize)*sum(ID_df$Count * ID_df$HomeGR)
      meanFitNonHome = (1/popsize)*sum(ID_df$Count * ID_df$NonHomeGR)
      
      probVec = ID_df$Count / popsize + (ID_df$Count / popsize)*(ID_df$HomeGR -meanFitHome)
      # all neg probs = 0 prob
      probVec[probVec<0] = 0
      offspring = rmultinom(1, size = popsize, prob = probVec)
      #   
      
      
      # save important parameters
      allIDArrayList[[gen]] = ID_df
      meanFitnessTrajHome = c(meanFitnessTrajHome, meanFitHome)
      meanFitnessTrajNonHome = c(meanFitnessTrajNonHome, meanFitNonHome)
      
      
      
      # update ID df with new count
      ID_df$Count = offspring 
      
      # remove any extinct lineages
      ID_df = ID_df[!ID_df$Count == 0, ]
      
      # re-do counting so dont have issues calling indexes of parents
      ID_df$Strain = 1:length(ID_df$Strain)
      
      
      
      #### M  U  T  A  T  I  O  N  S
      #---------------------------------------------#
      
      # Sample number of muts that will happen and then choose parent IDs
      M = rpois(1, popsize*U)
      
      if(M > 0 )
      {
        parentIDs = sample(ID_df$Strain, M, prob = (ID_df$Count/popsize), replace = TRUE)
        
        
        # sample M mutation effects from JDFE, usimg Sigma (covariance matrix ) from above fit
        mutEffects = as.array(mvrnorm(M, uVec, Sigma))
        
        if(M ==1) # to correct for a single sample frommvrnorm being read as a vector
        {
          mutEffects = array(mutEffects, dim = c(1,2))
        }
        
        # make M new rows for ID_df, one for each mut (4 cols, ID#, homeFit, respFit (!!!!!which not doing anything with RN!!!!))
        newMutArray = as.data.frame(array(data = 0, dim = c(M, 4)))
        colnames(newMutArray) = c('Strain', 'HomeGR', 'NonHomeGR', 'Count')
        
        # make new mutant ID's ,, just count up from highest current ID
        newMutArray$Strain = (max(ID_df$Strain)+1):(max(ID_df$Strain) + M )
        
        # update new Mut Fitness values by adding sampled fit Effects
        newMutArray$HomeGR = ID_df$HomeGR[parentIDs] + mutEffects[,1]
        newMutArray$NonHomeGR = ID_df$NonHomeGR[parentIDs] + mutEffects[,2]
        
        # add count of ONE to each new mut, subtract offspring from parent
        newMutArray$Count = rep(1, times = M)
        
        
        ####??? faster way to do this
        for(rem in parentIDs)
        {
          ID_df$Count[rem] = ID_df$Count[rem] -1
        }
        
        # collate old Id_df and new muts
        ID_df = rbind(ID_df, newMutArray)
        
      }
      # start again with next gen :) 
    }
    
    allMeanFitTrajHome[,it] = meanFitnessTrajHome
    allMeanFitTrajNonHome[,it] = meanFitnessTrajNonHome
    
  }
  
  
  time = 1:numGens
  meanFitTrajHome = rowMeans(allMeanFitTrajHome)
  meanFitTrajNonHome = rowMeans(allMeanFitTrajNonHome)
  
  slopeHome = summary(lm(rowMeans(allMeanFitTrajHome)[800:numGens]~time[800:numGens]))$coefficients[2,1]
  slopeNonHome = summary(lm(rowMeans(allMeanFitTrajNonHome)[800:numGens]~time[800:numGens]))$coefficients[2,1]
  
  varHome_eachGen = rowVars(allMeanFitTrajHome)
  varSlopeHome = summary(lm(rowVars(allMeanFitTrajHome)[800:numGens]~time[800:numGens]))$coefficients[2,1]
  
  varNonHome_eachGen = rowVars(allMeanFitTrajNonHome)
  varSlopeNonHome = summary(lm(rowVars(allMeanFitTrajNonHome)[800:numGens]~time[800:numGens]))$coefficients[2,1]
  
  covar_eachGen = c()
  for(covT in 1:numGens)
  {
    covar_eachGen = c(covar_eachGen, cov(allMeanFitTrajHome[covT, ], allMeanFitTrajNonHome[covT, ]))
  }
  
  covarSlope = summary(lm(covar_eachGen[800:numGens]~time[800:numGens]))$coefficients[2,1]
  
  
  
  
  UbAdjust = 1 - pnorm(0, uHome, sigHome)
  homeSvals = sValallTheoMuts[,1]*(sValallTheoMuts[,1]>=0)
  r1 = 2*mean((homeSvals)^2) / UbAdjust 
  r2 = 2*mean((homeSvals)*(sValallTheoMuts[,2])) / UbAdjust 
  D11 = 2*mean((homeSvals)^3)/UbAdjust 
  D12 = 2*mean((homeSvals)^2*(sValallTheoMuts[,2]))/UbAdjust 
  D22 = 2*mean((homeSvals)*(sValallTheoMuts[,2])^2)/UbAdjust 
  
  
  if(wantPlot == TRUE)
  {
   
    fitness = c(meanFitTrajHome, meanFitTrajNonHome)
    min = c(meanFitTrajHome - sqrt(varHome_eachGen), meanFitTrajNonHome - sqrt(varNonHome_eachGen)) # min of ribbon is -1 SD
    max = c(meanFitTrajHome + sqrt(varHome_eachGen), meanFitTrajNonHome + sqrt(varNonHome_eachGen)) # max of ribbon is +1 SD
    time = rep(1:numGens, times = 2)
    Env = as.factor(c(rep(1, times = numGens), rep(2, times = numGens)))
    WbarDf = data.frame(time, fitness, min, max, Env)
    
    
    print( ggplot(WbarDf, aes(x=as.numeric(time),y= fitness, colour = Env )) + geom_line()  +
             geom_ribbon(aes(x= time, ymin=min, ymax=max,fill = Env,group = Env), alpha = 0.3, linetype = 0) + 
             scale_color_manual(values=c('chocolate1','dodgerblue2'))+  geom_hline(yintercept =  0, linetype = 'dashed')+
             theme_bw()+theme(panel.border = element_rect(colour = 'black', fill = NA, size = 2), 
                              panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
                              axis.text =element_text(size = 20,family = 'Helvetica'), axis.title = element_blank(), legend.position = 'none')+ 
             xlab("Time")+ ylab("Fitness"))
  }
  
  
  output = list(allMeanFitTrajNonHome, allMeanFitTrajHome, meanFitTrajHome, meanFitTrajNonHome, varHome_eachGen, varNonHome_eachGen, covar_eachGen, data.frame(slopeHome, slopeNonHome,varSlopeHome, varSlopeNonHome,  covarSlope,r1,r2,D11,D12,D22))
  
  
  return(output)
  
  
  
}

# easily changed parameters
numGens = 1000 # mut be at least 1000 to get accurate estimation of slope, can be more
numiterations = 300 # do 300 for paper results, is quite computationally intensive, recommend running on a computing cluster
wantPlot = FALSE # printing a plot for all 125 JDFEs is very computationally intensive


# foundational parameters, can change only if recalculate xc and v values as well
popsizeVec = c(10^4, 10^5, 10^6) 
U = 10^-4
sigHome = 0.01 #home DFE SD 
uHome = - 0.001 # mean home DFE

## Values from numerically solving eq 18 and 19 from Good et al PNAS with values of home DFE and beta = 2 (because gaussian)
xc1 = 0.015296	#xc for Nu1 (popsize = 10^4, U = 10^-4) and out uHome and sdHome
v1 = 1.96423E-05
xc10 = 0.0286986
v10 = 4.57933E-05
xc100 = 0.0398821	
v100 = 6.81092E-05


##### Run Evolution Simulation ####

meanHomeSlopes = c()
meanNonHomeSlopes = c()
meanVarSlopeHome =  c()
meanVarSlopeNonHome = c()
meanCoVarSlope =  c()

r1Vals = c()
r2Vals = c()
D11Vals = c()
D22Vals= c()
D12Vals = c()

fullMeanFitTrajs_Home = list()
fullMeanFitTrajs_nonHome = list()
percent_colSens_lastGen = c()

sigNonHomeList = c()
uNonHomeList = c()
corrList = c()
allNu = c()
allpopsize = c()
allUb = c()




for(popsize in popsizeVec)
{
  print('popsize')
  print(popsize)
  JDFEnum = 0
  for(sigNonHome in seq(0.0001, 0.01, length.out =5))
  {
    for(uNonHome in seq(0.0001, 0.01, length.out = 5))
    {
      for(corr in seq(-0.9, 0.9, length.out = 5))
      {
        JDFEnum = JDFEnum + 1
        print(JDFEnum)
        
        
        covar = corr*sqrt(sigHome^2 * sigNonHome^2)
        uVec = c(uHome, uNonHome)
        Sigma =  matrix(c(sigHome^2,covar,covar,sigNonHome^2), nrow=2)
        UbAdjust = 1 - pnorm(0, uHome, sigHome)
        
        # save parameters of JDFE 
        sigNonHomeList = c(sigNonHomeList, sigNonHome)
        uNonHomeList = c(uNonHomeList, uNonHome)
        corrList = c(corrList, corr)
        allNu = c(allNu, popsize*U)
        allpopsize= c(allpopsize, popsize)
        allUb = c(allUb, U*UbAdjust)
        
        # run evo sim
        output = simEvo_WrightFisher(uVec, Sigma, popsize, U, numiterations, numGens, wantPlot) 
        outputDF = output[[8]]
        
         # save results
        meanHomeSlopes = c(meanHomeSlopes, outputDF$slopeHome)
        meanNonHomeSlopes = c(meanNonHomeSlopes,outputDF$slopeNonHome)
        meanVarSlopeHome = c(meanVarSlopeHome, outputDF$varSlopeNonHome)
        meanVarSlopeNonHome = c(meanVarSlopeNonHome, outputDF$varSlopeNonHome)
        meanCoVarSlope = c(meanCoVarSlope, outputDF$covarSlope)
        
        r1Vals = c(r1Vals, outputDF$r1)
        r2Vals = c(r2Vals,outputDF$r2 )
        D11Vals = c(D11Vals, outputDF$D11)
        D22Vals= c(D22Vals, outputDF$D22)
        D12Vals = c(D12Vals, outputDF$D12)
       
        fullMeanFitTrajs_Home[[JDFEnum]] = output[[2]]
        fullMeanFitTrajs_nonHome[[JDFEnum]] = output[[1]]
        
        lastgen_nonHomeFit = output[[1]][200, ] # changed to200 because in all NU pattern is established but not all the way to all 0s and 1s yet
        percent_colSens_lastGen = c(percent_colSens_lastGen, sum(lastgen_nonHomeFit<0) / length(lastgen_nonHomeFit)) 
   
        
      }
    }
    
  }
  
}



#### Calculate Clonal Interference R and D Values  ####
    # contains variables which rely on running above code, at least the JDFEs parameters, not the evolution results
    
UbAdjust = 1 - pnorm(0, uHome, sigHome)
Ub = U*UbAdjust

#numJDFEs = length(r1Vals) / length(popsizeVec)
allr1_new = c()
allr2_new = c()
allD11_new = c()
allD12_new = c()
allD22_new = c()

end =125
for(popsize in popsizeVec)
{
  if(popsize*U == 1)
  {
    xc = xc1
    v = v1
    
  }else if(popsize*U == 10){
    xc = xc10
    v = v10
    
  }else if(popsize*U == 100)
  {
    xc = xc100
    v = v100
  }
    
  
  for(i in 1:end) 
  {
    print(i)
    uNonHome = uNonHomeList[i]
    sigNonHome = sigNonHomeList[i]
    corr =corrList[i]
    
    UbAdjust = 1 - pnorm(0, uHome, sigHome)
    cov = corr*sqrt(sigHome^2 * sigNonHome^2)
    
    ### FIT MVRNORM TO DATA 
    Sigma = array(data = c(sigHome^2, cov, cov, sigNonHome^2), dim = c(2,2))
    sValallTheoMuts= mvrnorm(n = 10000, c(uHome,uNonHome), Sigma)
    
    homeSvals = sValallTheoMuts[,1]
    respSvals = sValallTheoMuts[,2]
 
    r1 = mean((homeSvals)*fixP_num(homeSvals))
    r2 = mean((respSvals)*fixP_num(homeSvals))
    D11 = mean(((homeSvals)^2) * fixP_num(homeSvals))
    D12 = mean((homeSvals)*(respSvals)*fixP_num(homeSvals))
    D22 = mean(fixP_num(homeSvals)*(respSvals)^2)
    
    
    allr1_new = c(allr1_new, r1)
    allr2_new = c(allr2_new, r2)
    allD11_new = c(allD11_new, D11)
    allD12_new = c(allD12_new, D12)
    allD22_new = c(allD22_new, D22)
    
  }
  
}

## Save All Data
data_WrightFisherSims_allNu = data.frame(allNu, allpopsize, allUb,rep(uHome, times = length(allNu)),rep(sigHome, times = length(allNu)),  uNonHomeList, sigNonHomeList, corrList, allr1_new, allr2_new, allD11_new, allD22_new, allD12_new, meanHomeSlopes, meanNonHomeSlopes, meanVarSlopeHome, meanVarSlopeNonHome, meanCoVarSlope, r1Vals, r2Vals, D11Vals, D22Vals, D12Vals,  percent_colSens_lastGen  )
colnames(data_WrightFisherSims_allNu) = c('NU', "PopSize", "Ub", "uHome", 'sdHome', 'uNonHome', 'sdNonHome', 'corr', 'r1_CI', 'r2_CI', 'D11_CI', "D22_CI", "D12_CI", 'homeSlope', 'nonHomeSlope', 'homeVarSlope', 'nonHomeVarSlope', 'covarSlope', 'r1_SSWM', 'r2_SSWM', 'D11_SSWM', 'D22_SSWM', 'D12_SSWM', 'percentColSens_lastgen')

write.csv(data_WrightFisherSims_allNu, paste0('data_WrightFisherSims_allNu_numGens_',numGens, '_numIts_', numiterations ,'.csv'))



data_thisNu_1   = subset(data_WrightFisherSims_allNu , NU == 1)
data_thisNu_10  = subset(data_WrightFisherSims_allNu, NU == 10)
data_thisNu_100 = subset(data_WrightFisherSims_allNu, NU == 100)



####   Nu1 Plots and LMs
color = 'palegreen4'

r2BenGood_v_sim_LINEAR_NU1 = ggplot(data_thisNu_1, aes(x=(r2_CI) ,y= (nonHomeSlope/(10^4*Ub)))) + 
  geom_hline(yintercept = 0 , color = 'grey', size = 1.5) +geom_vline(xintercept = 0, color = 'grey', size= 1.5) +
  geom_point(col = color, size = 7, alpha = 0.6) +
  geom_smooth(method='lm', formula= y~x, col = color, lwd = 0.8) + geom_abline(size  = 1.5, intercept = 0, slope = 1, linetype = 'dashed', col = 'grey')+
  theme_bw()+theme(panel.border = element_rect(colour = 'black', fill = NA, size = 1.5), 
                   panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
                   axis.text.x =element_text(size = 35,family = 'Helvetica', color = c('black')),axis.text.y =element_text(size = 35,family = 'Helvetica', color = c('black')), axis.title = element_blank()) +
  scale_y_continuous(label=scientific_10) + scale_x_continuous(label=scientific_10, breaks = c(0*10^-5, 2*10^-5, -2*10^-5, 4*10^-5))  


# LMs r2 v slope, Nu1
nonHomeSlope_nu1 = data_thisNu_1$nonHomeSlope
r2_CI_nu1 = data_thisNu_1$r2_CI

LM_Nu1_r2_vs_SimNonHomeSlope = summary(lm(nonHomeSlope_nu1 ~ r2_CI_nu1 ))
rSq_Nu1_LM_r2_v_SimNonHomeSlope = LM_Nu1_r2_vs_SimNonHomeSlope$r.squared
Coeffs_Nu1_LM_r2_v_SimNonHomeSlope = LM_Nu1_r2_vs_SimNonHomeSlope$coefficients

                            # remove 1:5 for ease of viewing, as they have much lower D22 than other JDFEs
D22BenGood_v_sim_Log_NU1 = ggplot(data_thisNu_1[-c(1:5), ], aes(x=D22_CI ,y= nonHomeVarSlope/(10^4*Ub))) + 
  geom_point(col = color, size = 7, alpha = 0.6) +
  geom_smooth(method='lm', formula= y~x, col = color, lwd = 0.8) + geom_abline(size  = 1.5, intercept = 0, slope = 1, linetype = 'dashed', col = 'grey')+
  theme_bw()+theme(panel.border = element_rect(colour = 'black', fill = NA, size = 1.5), 
                   panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
                   axis.text.x =element_text(size = 35,family = 'Helvetica', color = c('black')),axis.text.y =element_text(size = 35,family = 'Helvetica', color = c('black')), axis.title = element_blank()) +
  scale_y_continuous(trans = 'log10' ,breaks =  c(3*10^-9, 3*10^-8, 3*10^-7),
                     labels = scientific_10) + scale_x_continuous(trans = 'log10' ,breaks = c(3*10^-9, 3*10^-8, 3*10^-7),
                                                                  labels = scientific_10)  

# LMs D22 v Varslope, Nu1
varSlope_log_nu1 = log(data_thisNu_1$nonHomeVarSlope)
D22_CI_log_nu1 = log(data_thisNu_1$D22_CI)

LM_Nu1_D22_vs_SimNonHomeVarSlope = summary(lm(varSlope_log_nu1~D22_CI_log_nu1))
rSq_Nu1_LM_D22_v_SimNonHomeVarSlope = LM_Nu1_D22_vs_SimNonHomeVarSlope$r.squared
Coeffs_Nu1_LM_D22_v_SimNonHomeVarSlope = LM_Nu1_D22_vs_SimNonHomeVarSlope $coefficients




# save plots
ggsave('r2BenGood_v_sim_LINEAR_NU1.pdf', r2BenGood_v_sim_LINEAR_NU1, width = 5.26, height = 3.82)
ggsave('D22BenGood_v_sim_LINEAR_NU1.pdf', D22BenGood_v_sim_Log_NU1, width = 5.26, height = 3.82)

# save ALL LMs
rSqVals_Nu1_LMS = data.frame(rSq_Nu1_LM_r2_v_SimNonHomeSlope, rSq_Nu1_LM_D22_v_SimNonHomeVarSlope)
colnames(rSqVals_Nu1_LMS) = c('R2_v_SimNonHomeSlope', 'D22_v_SimNonHomeVarSlope')
coefficients_Nu1_LMS = rbind(Coeffs_Nu1_LM_r2_v_SimNonHomeSlope, Coeffs_Nu1_LM_D22_v_SimNonHomeVarSlope)










####     Nu10 Plots and LMS
color = 'darkgoldenrod'
r2BenGood_v_sim_LINEAR_NU10 = ggplot(data_thisNu_10, aes(x=(r2_CI) ,y= (nonHomeSlope/(10^5*Ub)))) + 
  geom_hline(yintercept = 0 , color = 'grey', size = 1.5) +geom_vline(xintercept = 0, color = 'grey', size= 1.5) +
  geom_point(col = color, size = 7, alpha = 0.6) +
  geom_smooth(method='lm', formula= y~x, col = color, lwd = 0.8) + geom_abline(size  = 1.5, intercept = 0, slope = 1, linetype = 'dashed', col = 'grey')+
  theme_bw()+theme(panel.border = element_rect(colour = 'black', fill = NA, size = 1.5), 
                   panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
                   axis.text.x =element_text(size = 35,family = 'Helvetica', color = c('black')),axis.text.y =element_text(size = 35,family = 'Helvetica', color = c('black')), axis.title = element_blank()) +
  scale_y_continuous(label=scientific_10) + scale_x_continuous(label=scientific_10, breaks = c(0*10^-5, -1*10^-5, 1*10^-5))  

# LMs r2 v slope, Nu10
nonHomeSlope_nu10 = data_thisNu_10$nonHomeSlope
r2_CI_nu10 = data_thisNu_10$r2_CI

LM_Nu10_r2_vs_SimNonHomeSlope = summary(lm(nonHomeSlope_nu10  ~ r2_CI_nu10 ))
rSq_Nu10_LM_r2_v_SimNonHomeSlope = LM_Nu10_r2_vs_SimNonHomeSlope$r.squared
Coeffs_Nu10_LM_r2_v_SimNonHomeSlope = LM_Nu10_r2_vs_SimNonHomeSlope$coefficients


D22BenGood_v_sim_Log_NU10 = ggplot(data_thisNu_10[-c(1:5), ], aes(x=D22_CI ,y= nonHomeVarSlope/(10^5*Ub))) + 
  geom_point(col = color, size = 7, alpha = 0.6) +
  geom_smooth(method='lm', formula= y~x, col = color, lwd = 0.8) + geom_abline(size  = 1.5, intercept = 0, slope = 1, linetype = 'dashed', col = 'grey')+
  theme_bw()+theme(panel.border = element_rect(colour = 'black', fill = NA, size = 1.5), 
                   panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
                   axis.text.x =element_text(size = 35,family = 'Helvetica', color = c('black')),axis.text.y =element_text(size = 35,family = 'Helvetica', color = c('black')), axis.title = element_blank()) +
  scale_y_continuous(trans = 'log10' ,breaks =  c(3*10^-9, 3*10^-8, 3*10^-7),
                     labels = scientific_10) + scale_x_continuous(trans = 'log10' ,breaks = c(3*10^-10, 3*10^-9, 3*10^-8, 3*10^-7),
                                                                  labels = scientific_10)  


# LMs D22 v Varslope, Nu10
varSlope_log_nu10 = log(data_thisNu_10$nonHomeVarSlope)
D22_CI_log_nu10 = log(data_thisNu_10$D22_CI)

LM_Nu10_D22_vs_SimNonHomeVarSlope = summary(lm(varSlope_log_nu10~ D22_CI_log_nu10))
rSq_Nu10_LM_D22_v_SimNonHomeVarSlope = LM_Nu10_D22_vs_SimNonHomeVarSlope$r.squared
Coeffs_Nu10_LM_D22_v_SimNonHomeVarSlope = LM_Nu10_D22_vs_SimNonHomeVarSlope $coefficients


# save plots
ggsave('r2BenGood_v_sim_LINEAR_NU10.pdf', r2BenGood_v_sim_LINEAR_NU10, width = 5.26, height = 3.82)
ggsave('D22BenGood_v_sim_LINEAR_NU10.pdf', D22BenGood_v_sim_Log_NU10, width = 5.26, height = 3.82)

# save LMs
rSqVals_Nu10_LMS = data.frame(rSq_Nu10_LM_r2_v_SimNonHomeSlope, rSq_Nu10_LM_D22_v_SimNonHomeVarSlope)
colnames(rSqVals_Nu10_LMS) = c('R2_v_SimNonHomeSlope', 'D22_v_SimNonHomeVarSlope')
coefficients_Nu10_LMS = rbind(Coeffs_Nu10_LM_r2_v_SimNonHomeSlope, Coeffs_Nu10_LM_D22_v_SimNonHomeVarSlope)




####    Nu100 Plots and LMS
color = 'hotpink4'
r2BenGood_v_sim_LINEAR_NU100 = ggplot(data_thisNu_100,aes(x=(r2_CI) ,y= (nonHomeSlope/(10^6*Ub)))) + 
  geom_hline(yintercept = 0 , color = 'grey', size = 1.5) +geom_vline(xintercept = 0, color = 'grey', size= 1.5) +
  geom_point(col = color, size = 7, alpha = 0.6) +
  geom_smooth(method='lm', formula= y~x, col = color, lwd = 0.8) + geom_abline(size  = 1.5, intercept = 0, slope = 1, linetype = 'dashed', col = 'grey')+
  theme_bw()+theme(panel.border = element_rect(colour = 'black', fill = NA, size = 1.5), 
                   panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
                   axis.text.x =element_text(size = 35,family = 'Helvetica', color = c('black')),axis.text.y =element_text(size = 35,family = 'Helvetica', color = c('black')), axis.title = element_blank()) +
  scale_y_continuous(label=scientific_10) + scale_x_continuous(label=scientific_10, breaks = c(0*10^-5, 2*10^-6,-2*10^-6))  

# LMs r2 v slope, Nu1
nonHomeSlope_nu100 = data_thisNu_100$nonHomeSlope
r2_CI_nu100 = data_thisNu_100$r2_CI

LM_Nu100_r2_vs_SimNonHomeSlope = summary(lm(nonHomeSlope_nu100 ~ r2_CI_nu100))
rSq_Nu100_LM_r2_v_SimNonHomeSlope = LM_Nu100_r2_vs_SimNonHomeSlope$r.squared
Coeffs_Nu100_LM_r2_v_SimNonHomeSlope = LM_Nu100_r2_vs_SimNonHomeSlope$coefficients


D22BenGood_v_sim_Log_NU100 = ggplot(data_thisNu_100[-c(1:5), ], aes(x=D22_CI ,y= nonHomeVarSlope/(10^6*Ub))) + 
  geom_point(col = color, size = 7, alpha = 0.6) +
  geom_smooth(method='lm', formula= y~x, col = color, lwd = 0.8) + geom_abline(size  = 1.5, intercept = 0, slope = 1, linetype = 'dashed', col = 'grey')+
  theme_bw()+theme(panel.border = element_rect(colour = 'black', fill = NA, size = 1.5), 
                   panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
                   axis.text.x =element_text(size = 35,family = 'Helvetica', color = c('black')),axis.text.y =element_text(size = 35,family = 'Helvetica', color = c('black')), axis.title = element_blank()) +
  scale_y_continuous(trans = 'log10' ,breaks =  c(3*10^-11, 3*10^-10, 3*10^-9),
                     labels = scientific_10) + scale_x_continuous(trans = 'log10' ,breaks = c(3*10^-10, 3*10^-9, 3*10^-8),
                                                                  labels = scientific_10)  



# LMs D22 v Varslope, Nu10
varSlope_log_nu100 = log(data_thisNu_100$nonHomeVarSlope)
D22_CI_log_nu100 = log(data_thisNu_100$D22_CI)

LM_Nu100_D22_vs_SimNonHomeVarSlope = summary(lm(varSlope_log_nu100 ~ D22_CI_log_nu100))
rSq_Nu100_LM_D22_v_SimNonHomeVarSlope = LM_Nu100_D22_vs_SimNonHomeVarSlope$r.squared
Coeffs_Nu100_LM_D22_v_SimNonHomeVarSlope = LM_Nu100_D22_vs_SimNonHomeVarSlope $coefficients


# save plots
ggsave('r2BenGood_v_sim_LINEAR_NU100.pdf', r2BenGood_v_sim_LINEAR_NU100, width = 5.26, height = 3.82)
ggsave('D22BenGood_v_sim_LINEAR_NU100.pdf', D22BenGood_v_sim_Log_NU100, width = 5.26, height = 3.82)

# save LMs
rSqVals_Nu100_LMS = data.frame(rSq_Nu100_LM_r2_v_SimNonHomeSlope, rSq_Nu100_LM_D22_v_SimNonHomeVarSlope)
colnames(rSqVals_Nu100_LMS) = c('R2_v_SimNonHomeSlope', 'D22_v_SimNonHomeVarSlope')
coefficients_Nu100_LMS = rbind(Coeffs_Nu100_LM_r2_v_SimNonHomeSlope, Coeffs_Nu100_LM_D22_v_SimNonHomeVarSlope)





## Compile all LMS

coefficients_allNu_LMS_WrightFisher = rbind(coefficients_Nu1_LMS, coefficients_Nu10_LMS, coefficients_Nu100_LMS)

rSqVals_WrightFisher_LMS = data.frame(c(rSq_Nu1_LM_r2_v_SimNonHomeSlope, rSq_Nu10_LM_r2_v_SimNonHomeSlope, rSq_Nu100_LM_r2_v_SimNonHomeSlope ), c(rSq_Nu1_LM_D22_v_SimNonHomeVarSlope, rSq_Nu10_LM_D22_v_SimNonHomeVarSlope, rSq_Nu100_LM_D22_v_SimNonHomeVarSlope))
colnames(rSqVals_WrightFisher_LMS) = c('R2_v_SimNonHomeSlope', 'D22_v_SimNonHomeVarSlope')

write.csv(rSqVals_WrightFisher_LMS, paste0('rSqVals_WrightFisher_LMS_numGens',numGens, '_numIts_', numiterations ,'.csv'))
write.csv(coefficients_allNu_LMS_WrightFisher,paste0('coefficients_allNu_LMS_WrightFisher_numGens',numGens, '_numIts_', numiterations ,'.csv'))




####  Plot CI slope results vs SSWM predictions  ####
####   Nu1 Plots and LMs
color = 'palegreen4'


## R2 
r2SSWMM_sim_LINEAR_NU1 = ggplot(data_thisNu_1, aes(x=(r2_SSWM/10^-5) ,y= (nonHomeSlope/(10^4*Ub)/10^-5))) + 
  geom_hline(yintercept = 0 , color = 'grey', size = 1.5) +geom_vline(xintercept = 0, color = 'grey', size= 1.5) +
  geom_point(col = color, size = 7, alpha = 0.6) +
  labs(title = 'x10^-5', x = 'x10^-5')+
  geom_smooth(method='lm', formula= y~x, col = color, lwd = 0.8) + geom_abline(size  = 1.5, intercept = 0, slope = 1, linetype = 'dashed', col = 'grey')+
  theme_bw()+theme(panel.border = element_rect(colour = 'black', fill = NA, size = 1.5), 
                   panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
                   axis.text.x =element_text(size = 16,family = 'Helvetica', color = c('black')),axis.text.y =element_text(size = 16,family = 'Helvetica', color = c('black')),
                   axis.title.y = element_blank(),
                   axis.title.x = element_text(size = 16,family = 'Helvetica', color = c('black'), hjust = 1)) +
  scale_y_continuous(  breaks = c(-5,0, 5)) + scale_x_continuous(breaks = c(0, -10,10))  


# LMs r2 v slope, Nu1
nonHomeSlope_nu1 = data_thisNu_1$nonHomeSlope
r2_SSWM_nu1 = data_thisNu_1$r2_SSWM

LM_Nu1_r2SSWM_vs_SimNonHomeSlope = summary(lm(nonHomeSlope_nu1 ~ r2_SSWM_nu1 ))
rSq_Nu1_LM_r2SSWM_v_SimNonHomeSlope = LM_Nu1_r2SSWM_vs_SimNonHomeSlope$r.squared
Coeffs_Nu1_LM_r2SSWM_v_SimNonHomeSlope = LM_Nu1_r2SSWM_vs_SimNonHomeSlope$coefficients

## D22
D22SSWM_v_sim_Log_NU1 = ggplot(data_thisNu_1[-c(1:5), ], aes(x=D22_SSWM/ 10^-8 ,y= nonHomeVarSlope/(10^4*Ub)/10^-8)) + 
  geom_point(col = color, size = 7, alpha = 0.6) +
  labs(title = 'x10^-8', x = 'x10^-8')+
  geom_smooth(method='lm', formula= y~x, col = color, lwd = 0.8) + geom_abline(size  = 1.5, intercept = 0, slope = 1, linetype = 'dashed', col = 'grey')+
  theme_bw()+theme(panel.border = element_rect(colour = 'black', fill = NA, size = 1.5), 
                   panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
                   axis.text.x =element_text(size = 16,family = 'Helvetica', color = c('black')),axis.text.y =element_text(size = 16,family = 'Helvetica', color = c('black')), 
                   axis.title.y = element_blank(),
                   axis.title.x = element_text(size = 16,family = 'Helvetica', color = c('black'), hjust = 1)) +
  scale_y_continuous(trans = 'log10' ,breaks =  c(3, 30)) + scale_x_continuous(trans = 'log10' ,breaks = c(3, 30, 300)) 

# LMs D22 v Varslope, Nu1
varSlope_log_nu1 = log(data_thisNu_1$nonHomeVarSlope)
D22_SSWM_log_nu1 = log(data_thisNu_1$D22_SSWM)

LM_Nu1_D22SSWM_vs_SimNonHomeVarSlope = summary(lm(varSlope_log_nu1~D22_SSWM_log_nu1))
rSq_Nu1_LM_D22SSWM_v_SimNonHomeVarSlope = LM_Nu1_D22SSWM_vs_SimNonHomeVarSlope$r.squared
Coeffs_Nu1_LM_D22SSWM_v_SimNonHomeVarSlope = LM_Nu1_D22SSWM_vs_SimNonHomeVarSlope$coefficients




# save plots
ggsave('WFSim_r2SSWM_v_sim_LINEAR_NU1.pdf', r2SSWMM_sim_LINEAR_NU1, width = 5.26, height = 3.82)
ggsave('WFSim_D22SSWM_v_sim_LINEAR_NU1.pdf', D22SSWM_v_sim_Log_NU1, width = 5.26, height = 3.82)

# save LMs
rSqVals_Nu1_LMS_SSWM = data.frame(rSq_Nu1_LM_r2SSWM_v_SimNonHomeSlope, rSq_Nu1_LM_D22SSWM_v_SimNonHomeVarSlope)
colnames(rSqVals_Nu1_LMS_SSWM) = c('R2_v_SimNonHomeSlope', 'D22_v_SimNonHomeVarSlope')
coefficients_Nu1_LMS_SSWM = rbind(Coeffs_Nu1_LM_r2SSWM_v_SimNonHomeSlope, Coeffs_Nu1_LM_D22SSWM_v_SimNonHomeVarSlope)









####     Nu10 Plots and LMS
color = 'darkgoldenrod'    # scale  by 10^-5 on y and x
r2SSWM_v_sim_LINEAR_NU10 = ggplot(data_thisNu_10, aes(x=(r2_SSWM/ 10^-5) ,y= (nonHomeSlope/(10^5*Ub)/10^-5))) + 
  geom_hline(yintercept = 0 , color = 'grey', size = 1.5) +geom_vline(xintercept = 0, color = 'grey', size= 1.5) +
  geom_point(col = color, size = 7, alpha = 0.6) +
  geom_smooth(method='lm', formula= y~x, col = color, lwd = 0.8) + geom_abline(size  = 1.5, intercept = 0, slope = 1, linetype = 'dashed', col = 'grey')+
  labs(title = 'x10^-5', x = 'x10^-5')+
  theme_bw()+theme(panel.border = element_rect(colour = 'black', fill = NA, size = 1.5), 
                   panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
                   axis.text.x =element_text(size = 16,family = 'Helvetica', color = c('black')),axis.text.y =element_text(size = 16,family = 'Helvetica', color = c('black')), 
                   axis.title.y = element_blank(),
                   axis.title.x = element_text(size = 16,family = 'Helvetica', color = c('black'), hjust = 1)) +
  scale_y_continuous(breaks = c(-1, 0, 2)) + scale_x_continuous( breaks = c(0, -10,20))  

# LMs r2 v slope, Nu10
nonHomeSlope_nu10 = data_thisNu_10$nonHomeSlope
r2_SSWM_nu10 = data_thisNu_10$r2_SSWM

LM_Nu10_r2SSWM_vs_SimNonHomeSlope = summary(lm(nonHomeSlope_nu10  ~ r2_SSWM_nu10 ))
rSq_Nu10_LM_r2SSWM_v_SimNonHomeSlope = LM_Nu10_r2SSWM_vs_SimNonHomeSlope$r.squared
Coeffs_Nu10_LM_r2SSWM_v_SimNonHomeSlope = LM_Nu10_r2SSWM_vs_SimNonHomeSlope$coefficients


# scale by 10^-10 on y and 10^-8 on x for ease of reading
D22SSWM_v_sim_Log_NU10 = ggplot(data_thisNu_10[-c(1:5), ], aes(x=D22_SSWM /10^-10 ,y= nonHomeVarSlope/(10^5*Ub)/10^-10)) + 
  geom_point(col = color, size = 7, alpha = 0.6) +
  labs(title = 'x10^-10', x = 'x10^-10')+
  geom_smooth(method='lm', formula= y~x, col = color, lwd = 0.8) + geom_abline(size  = 1.5, intercept = 0, slope = 1, linetype = 'dashed', col = 'grey')+
  theme_bw()+theme(panel.border = element_rect(colour = 'black', fill = NA, size = 1.5), 
                   panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
                   axis.text.x =element_text(size = 16,family = 'Helvetica', color = c('black')),axis.text.y =element_text(size = 16,family = 'Helvetica', color = c('black')), 
                   axis.title.y = element_blank(),
                   axis.title.x = element_text(size = 16,family = 'Helvetica', color = c('black'), hjust = 1)) +
  scale_y_continuous(trans = 'log10' ,breaks =  c(3, 30, 300)) + scale_x_continuous(trans = 'log10' ,breaks = c(300, 3000, 30000))  


# LMs D22 v Varslope, Nu10
varSlope_log_nu10 = log(data_thisNu_10$nonHomeVarSlope)
D22_SSWM_log_nu10 = log(data_thisNu_10$D22_SSWM)

LM_Nu10_D22SSWM_vs_SimNonHomeVarSlope = summary(lm(varSlope_log_nu10~ D22_SSWM_log_nu10))
rSq_Nu10_LM_D22SSWM_v_SimNonHomeVarSlope = LM_Nu10_D22SSWM_vs_SimNonHomeVarSlope$r.squared
Coeffs_Nu10_LM_D22SSWM_v_SimNonHomeVarSlope = LM_Nu10_D22SSWM_vs_SimNonHomeVarSlope$coefficients


# save plots
ggsave('WFSim_r2SSWM_v_sim_LINEAR_NU10.pdf', r2SSWM_v_sim_LINEAR_NU10, width = 5.26, height = 3.82)
ggsave('WFSim_D22SSWM_v_sim_LINEAR_NU10.pdf', D22SSWM_v_sim_Log_NU10, width = 5.26, height = 3.82)

# save LMs
rSqVals_Nu10_LMS_SSWM = data.frame(rSq_Nu10_LM_r2SSWM_v_SimNonHomeSlope, rSq_Nu10_LM_D22SSWM_v_SimNonHomeVarSlope)
colnames(rSqVals_Nu10_LMS_SSWM) = c('R2_v_SimNonHomeSlope', 'D22_v_SimNonHomeVarSlope')
coefficients_Nu10_LMS_SSWM = rbind(Coeffs_Nu10_LM_r2SSWM_v_SimNonHomeSlope, Coeffs_Nu10_LM_D22SSWM_v_SimNonHomeVarSlope)




####    Nu100 Plots and LMS
color = 'hotpink4'
r2SSWM_v_sim_LINEAR_NU100 = ggplot(data_thisNu_100,aes(x=(r2_SSWM/10^-6) ,y= (nonHomeSlope/(10^6*Ub)/10^-6))) + 
  geom_hline(yintercept = 0 , color = 'grey', size = 1.5) +geom_vline(xintercept = 0, color = 'grey', size= 1.5) +
  geom_point(col = color, size = 7, alpha = 0.6) +
  labs(title = 'x10^-6', x = 'x10^-6')+
  geom_smooth(method='lm', formula= y~x, col = color, lwd = 0.8) + geom_abline(size  = 1.5, intercept = 0, slope = 1, linetype = 'dashed', col = 'grey')+
  theme_bw()+theme(panel.border = element_rect(colour = 'black', fill = NA, size = 1.5), 
                   panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
                   axis.text.x =element_text(size = 16,family = 'Helvetica', color = c('black')),axis.text.y =element_text(size = 16,family = 'Helvetica', color = c('black')), 
                   axis.title.y = element_blank(),
                   axis.title.x = element_text(size = 16,family = 'Helvetica', color = c('black'), hjust = 1)) +
  scale_y_continuous( breaks = c(-1,0,2) )+ scale_x_continuous( breaks = c(-100,0,100))  

# LMs r2 v slope, Nu1
nonHomeSlope_nu100 = data_thisNu_100$nonHomeSlope
r2_SSWM_nu100 = data_thisNu_100$r2_SSWM

LM_Nu100_r2SSWM_vs_SimNonHomeSlope = summary(lm(nonHomeSlope_nu100 ~ r2_SSWM_nu100))
rSq_Nu100_LM_r2SSWM_v_SimNonHomeSlope = LM_Nu100_r2SSWM_vs_SimNonHomeSlope$r.squared
Coeffs_Nu100_LM_r2SSWM_v_SimNonHomeSlope = LM_Nu100_r2SSWM_vs_SimNonHomeSlope$coefficients


D22SSWM_v_sim_Log_NU100 = ggplot(data_thisNu_100[-c(1:5), ], aes(x=D22_SSWM /10^-10 ,y= nonHomeVarSlope/(10^6*Ub)/10^-10)) + 
  geom_point(col = color, size = 7, alpha = 0.6) +
  labs(title = 'x10^-10', x = 'x10^-10')+
  geom_smooth(method='lm', formula= y~x, col = color, lwd = 0.8) + geom_abline(size  = 1.5, intercept = 0, slope = 1, linetype = 'dashed', col = 'grey')+
  theme_bw()+theme(panel.border = element_rect(colour = 'black', fill = NA, size = 1.5), 
                   panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
                   axis.text.x =element_text(size = 16,family = 'Helvetica', color = c('black')),axis.text.y =element_text(size = 16,family = 'Helvetica', color = c('black')), 
                   axis.title.y = element_blank(),
                   axis.title.x = element_text(size = 16,family = 'Helvetica', color = c('black'), hjust = 1)) +
  scale_y_continuous(trans = 'log10' ,breaks =  c(3, 30,300)) + scale_x_continuous(trans = 'log10' ,breaks = c(3000,30000))  



# LMs D22 v Varslope, Nu100
varSlope_log_nu100 = log(data_thisNu_100$nonHomeVarSlope)
D22_SSWM_log_nu100 = log(data_thisNu_100$D22_SSWM)

LM_Nu100_D22SSWM_vs_SimNonHomeVarSlope = summary(lm(varSlope_log_nu100 ~ D22_SSWM_log_nu100))
rSq_Nu100_LM_D22SSWM_v_SimNonHomeVarSlope = LM_Nu100_D22SSWM_vs_SimNonHomeVarSlope$r.squared
Coeffs_Nu100_LM_D22SSWM_v_SimNonHomeVarSlope = LM_Nu100_D22SSWM_vs_SimNonHomeVarSlope$coefficients


# save plots
ggsave('r2SSWM_v_sim_LINEAR_NU100.pdf', r2SSWM_v_sim_LINEAR_NU100, width = 5.26, height = 3.82)
ggsave('D22SSWM_v_sim_LINEAR_NU100.pdf', D22SSWM_v_sim_Log_NU100, width = 5.26, height = 3.82)

# save LMs
rSqVals_Nu100_LMS_SSWM = data.frame(rSq_Nu100_LM_r2SSWM_v_SimNonHomeSlope, rSq_Nu100_LM_D22SSWM_v_SimNonHomeVarSlope)
colnames(rSqVals_Nu100_LMS_SSWM) = c('R2_v_SimNonHomeSlope', 'D22_v_SimNonHomeVarSlope')
coefficients_Nu100_LMS_SSWM = rbind(Coeffs_Nu100_LM_r2SSWM_v_SimNonHomeSlope, Coeffs_Nu100_LM_D22SSWM_v_SimNonHomeVarSlope)




## Compile Plots
compiled_SSWMPredictorsResultsPLot = plot_grid(r2SSWMM_sim_LINEAR_NU1, r2SSWM_v_sim_LINEAR_NU10, r2SSWM_v_sim_LINEAR_NU100, D22SSWM_v_sim_Log_NU1,D22SSWM_v_sim_Log_NU10,D22SSWM_v_sim_Log_NU100,labels = c('A', 'C', 'E','B','D','F'), label_size = 20, nrow = 2, ncol = 3, align = 'v')
ggsave('compiled_SSWMPredictorsResultsPLot.pdf',compiled_SSWMPredictorsResultsPLot, width  = 8.67, height = 6.17)




## Compile all LMS

coefficients_allNu_LMS_WrightFisher_SSWMPredictors = rbind(coefficients_Nu1_LMS_SSWM, coefficients_Nu10_LMS_SSWM, coefficients_Nu100_LMS_SSWM)

rSqVals_WrightFisher_SSWMPredictors_LMS = data.frame(c(rSq_Nu1_LM_r2SSWM_v_SimNonHomeSlope, rSq_Nu10_LM_r2SSWM_v_SimNonHomeSlope, rSq_Nu100_LM_r2SSWM_v_SimNonHomeSlope ), c(rSq_Nu1_LM_D22SSWM_v_SimNonHomeVarSlope, rSq_Nu10_LM_D22SSWM_v_SimNonHomeVarSlope, rSq_Nu100_LM_D22SSWM_v_SimNonHomeVarSlope))
colnames(rSqVals_WrightFisher_LMS) = c('R2_v_SimNonHomeSlope', 'D22_v_SimNonHomeVarSlope')

write.csv(rSqVals_WrightFisher_SSWMPredictors_LMS, paste0('rSqVals_WrightFisher_SSWMPredictors_LMS_numGens',numGens, '_numIts_', numiterations ,'.csv'))
write.csv(coefficients_allNu_LMS_WrightFisher_SSWMPredictors,paste0('coefficients_allNu_LMS_WrightFisher_SSWMPredictors_numGens',numGens, '_numIts_', numiterations ,'.csv'))




#### end figure 3   #####






#####################################################################################################################################
###################      F I G U R E   4:      R A N K   O R D E R  ##################################################################
#####################################################################################################################################

# choose directory to save outputs
# setwd("~/")

## NOTE:: THIS SECTION REQUIRES RUNNING THE FIGURE 3 SECTION SINCE WILL NEED TO COMPARE OUTPUTS FROM THAT
   # SECTION WITH THE RANK PREDICTIONS




####  M A K E   J D F E ####
# use same home for all of them ,, same one use throughout the paper

sigHome = 0.01 #home DFE SD 
uHome = - 0.001 # mean home DFE 



## CHOOSE from set of JDFEs from non epistatic section

chosenJDFEs = c(116,86,81,51) 
numDrugs = length(chosenJDFEs)
sigNonHomeVec = data_WrightFisherSims_allNu$sdNonHome[chosenJDFEs] # make sure data_WrightFisherSims_allNu and others is properly saved from the SSWM sims fig 3
uNonHomeVec =data_WrightFisherSims_allNu$uNonHome[chosenJDFEs]
corVec = data_WrightFisherSims_allNu$corr[chosenJDFEs]

r1Vec = data_WrightFisherSims_allNu$r1_SSWM[chosenJDFEs]
r2Vec = data_WrightFisherSims_allNu$r2_SSWM[chosenJDFEs]
D11Vec = data_WrightFisherSims_allNu$D11_SSWM[chosenJDFEs]
D12Vec = data_WrightFisherSims_allNu$D12_SSWM[chosenJDFEs]
D22Vec = data_WrightFisherSims_allNu$D22_SSWM[chosenJDFEs]

homeID = rep(1, times = numDrugs)
nonHomeID = 1:numDrugs
homeMeanVec = rep(uHome, times = numDrugs)
homeVarVec = rep(sigHome, times = numDrugs)

cVec = r2Vec/sqrt(D22Vec)
hypotheticalDrugJDFEData = data.frame( homeID, nonHomeID, homeMeanVec,uNonHomeVec, homeVarVec,sigNonHomeVec, corVec, r1Vec, r2Vec, D11Vec, D22Vec, D12Vec, cVec)
colnames(hypotheticalDrugJDFEData) = c('Home', 'NonHome','meanHome', 'meanNonHome', 'varHome', 'varNonHome', 'cor', 'r1','r2','D11', 'D22', 'D12', 'C_val')



##### A: plot and SAVE the focal JDFEs #####
colors =  c('#fe0059', '#fea500', '#0059fe','#00e594' ) # for line color
colorPalletes = list()
colorPalletes[[1]] = c('white', '#ffc3d8', '#ff88b2', '#ff3a7f')
colorPalletes[[2]] = c('white','#ffe3b0', '#ffcf75','#ffba3a')
colorPalletes[[3]] = c('white', '#b0cbff', '#75a5ff', '#2672ff')
colorPalletes[[4]] = c('white', '#beffe8', '#83ffd3', '#34ffb7')


hypDrugJDFEplts = list()
for(d in 1:numDrugs)
{
  
  corr = corVec[d]
  respMean = uNonHomeVec[d]
  color = colors[d]
  pallete = colorPalletes[[d]]
  
  
  homeMean = uHome
  homeSd = sigHome
  respSd = sigNonHomeVec[d]
  xplotmin = -0.025
  xplotmax = 0.02
  yplotmin = -0.015
  yplotmax = 0.028
  
  
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
  hypDrugJDFEplts[[d]] = ggplot(df, aes(x=vals1, y=vals2, z=prob)) + 
    geom_contour_filled(bins = 4)+
    scale_fill_manual(values = pallete) + 
    theme_bw()+
    theme(panel.border = element_rect(colour = 'black', fill = NA, size = 2), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
          axis.text =element_text(size = 20,family = 'Helvetica', color ='black'), axis.title = element_blank(), legend.position = 'none')+ 
    geom_point(aes(weighted.mean(vals1, prob), weighted.mean(vals2, prob)), pch = 4, col = 'black', size = 3) + 
    scale_x_continuous( limits = c(xplotmin,xplotmax), breaks = c(-0.02,0,0.02))+
    scale_y_continuous(limits =c(yplotmin, yplotmax), breaks = c(-0.01,0,0.01,0.02,0.03) ) +
    geom_hline(yintercept =  0, linetype = 'dashed') +   geom_vline(xintercept =  0, linetype = 'dashed')
  
  
  
  
  
  
  
  ggsave(paste0('JDFEnum_', chosenJDFEs[d], '.pdf'), hypDrugJDFEplts[[d]], width =3.89, height = 3.67)
}





#### B. Predict and Test Rank Order of the chosen JDFES ####

numGens = 300 # factor of 2 diff from in sim for NUb 1 since this is scale of 2Nub gend
chanceColRes = 1-pnorm(0, (hypotheticalDrugJDFEData$r2)*numGens, sqrt(hypotheticalDrugJDFEData$D22*numGens)) # this is exactly percent chance of collateral resistance, lowest value of this is best (rank 1), highest is worst
hypotheticalDrugJDFEData$rankValue = rank(hypotheticalDrugJDFEData$C_val) # because biggest value of pnorm should be rank 1
hypotheticalDrugJDFEData$chanceColresGen300 = chanceColRes 

write.csv(hypotheticalDrugJDFEData , file = 'HypotheticalDrugJDFEData.csv')

### RUN SIM 
numGens = 1000
numiterations = 1000
popsize = 10^4
U = 10^-4

##### R U N     S I M 
simEvo_WrightFisher_rank <- function(uVec, Sigma, popsize, U, numiterations, numGens, wantPlot) 
{
  
  uHome = uVec[1]  # mean of home DFE
  uNonHome = uVec[2] # mean non home DFE
  sigHome = sqrt(Sigma[1,1]) # SD of Home DFE
  sigNonHome = sqrt(Sigma[2,2]) # SD of Nonhome DFE
  covar = Sigma[1,2] # covar of JDFE
  
  ################## M A K E     T H E      J D F E
  sValallTheoMuts = mvrnorm(100000, mu =uVec, Sigma = Sigma)
  DFE1 = sValallTheoMuts[,1]
  DFE2 = sValallTheoMuts[,2]
  
  
  allMeanFitTrajHome = array(0, dim = c(numGens, numiterations))
  allMeanFitTrajNonHome = array(0, dim = c(numGens, numiterations))
  
  for( it in 1:numiterations)
  {
    print(it)
    ID_df  = as.data.frame( array( data = c(1, 0, 0,popsize), dim = c(1,4)))
    colnames(ID_df) = c('Strain', 'HomeGR', 'NonHomeGR', 'Count')
    
    
    
    allIDArrayList = list()
    
    
    meanFitnessTrajHome = c()
    meanFitnessTrajNonHome = c()
    
    for(gen in 1:numGens)
    {
      #####    S  E  L  E  C   T  I  O  N    
      #---------------------------------------------------------------------------------------_#
      # select next gen offspring
      meanFitHome = (1/popsize)*sum(ID_df$Count * ID_df$HomeGR)
      meanFitNonHome = (1/popsize)*sum(ID_df$Count * ID_df$NonHomeGR)
      
      probVec = ID_df$Count / popsize + (ID_df$Count / popsize)*(ID_df$HomeGR -meanFitHome)
      # all neg probs = 0 prob
      probVec[probVec<0] = 0
      offspring = rmultinom(1, size = popsize, prob = probVec)
      #   
      
      
      # save important parameters
      allIDArrayList[[gen]] = ID_df
      meanFitnessTrajHome = c(meanFitnessTrajHome, meanFitHome)
      meanFitnessTrajNonHome = c(meanFitnessTrajNonHome, meanFitNonHome)
      
      
      
      # update ID df with new count
      ID_df$Count = offspring 
      
      # remove any extinct lineages
      ID_df = ID_df[!ID_df$Count == 0, ]
      
      # re-do counting so dont have issues calling indexes of parents
      ID_df$Strain = 1:length(ID_df$Strain)
      
      
      
      #### M  U  T  A  T  I  O  N  S
      #---------------------------------------------#
      
      # Sample number of muts that will happen and then choose parent IDs
      M = rpois(1, popsize*U)
      
      if(M > 0 )
      {
        parentIDs = sample(ID_df$Strain, M, prob = (ID_df$Count/popsize), replace = TRUE)
        
        
        # sample M mutation effects from JDFE, usimg Sigma (covariance matrix ) from above fit
        mutEffects = as.array(mvrnorm(M, uVec, Sigma))
        
        if(M ==1) # to correct for a single sample frommvrnorm being read as a vector
        {
          mutEffects = array(mutEffects, dim = c(1,2))
        }
        
        # make M new rows for ID_df, one for each mut (4 cols, ID#, homeFit, respFit (!!!!!which not doing anything with RN!!!!))
        newMutArray = as.data.frame(array(data = 0, dim = c(M, 4)))
        colnames(newMutArray) = c('Strain', 'HomeGR', 'NonHomeGR', 'Count')
        
        # make new mutant ID's ,, just count up from highest current ID
        newMutArray$Strain = (max(ID_df$Strain)+1):(max(ID_df$Strain) + M )
        
        # update new Mut Fitness values by adding sampled fit Effects
        newMutArray$HomeGR = ID_df$HomeGR[parentIDs] + mutEffects[,1]
        newMutArray$NonHomeGR = ID_df$NonHomeGR[parentIDs] + mutEffects[,2]
        
        # add count of ONE to each new mut, subtract offspring from parent
        newMutArray$Count = rep(1, times = M)
        
        
        ####??? faster way to do this
        for(rem in parentIDs)
        {
          ID_df$Count[rem] = ID_df$Count[rem] -1
        }
        
        # collate old Id_df and new muts
        ID_df = rbind(ID_df, newMutArray)
        
      }
      # start again with next gen :) 
    }
    
    allMeanFitTrajHome[,it] = meanFitnessTrajHome
    allMeanFitTrajNonHome[,it] = meanFitnessTrajNonHome
    
  }
  
  
  time = 1:numGens
  meanFitTrajHome = rowMeans(allMeanFitTrajHome)
  meanFitTrajNonHome = rowMeans(allMeanFitTrajNonHome)
  
  slopeHome = summary(lm(rowMeans(allMeanFitTrajHome)[800:numGens]~time[800:numGens]))$coefficients[2,1]
  slopeNonHome = summary(lm(rowMeans(allMeanFitTrajNonHome)[800:numGens]~time[800:numGens]))$coefficients[2,1]
  
  varHome_eachGen = rowVars(allMeanFitTrajHome)
  varSlopeHome = summary(lm(rowVars(allMeanFitTrajHome)[800:numGens]~time[800:numGens]))$coefficients[2,1]
  
  varNonHome_eachGen = rowVars(allMeanFitTrajNonHome)
  varSlopeNonHome = summary(lm(rowVars(allMeanFitTrajNonHome)[800:numGens]~time[800:numGens]))$coefficients[2,1]
  
  covar_eachGen = c()
  for(covT in 1:numGens)
  {
    covar_eachGen = c(covar_eachGen, cov(allMeanFitTrajHome[covT, ], allMeanFitTrajNonHome[covT, ]))
  }
  
  covarSlope = summary(lm(covar_eachGen[800:numGens]~time[800:numGens]))$coefficients[2,1]
  
  
  
  
  
  homeSvals = sValallTheoMuts[,1]*(sValallTheoMuts[,1]>=0)
  r1 = mean((homeSvals)^2)/UbAdjust
  r2 = mean((homeSvals)*(sValallTheoMuts[,2]))/UbAdjust
  D11 = mean((homeSvals)^3)/UbAdjust
  D12 = mean((homeSvals)^2*(sValallTheoMuts[,2]))/UbAdjust
  D22 = mean((homeSvals)*(sValallTheoMuts[,2])^2)/UbAdjust
  
  
  if(wantPlot == TRUE)
  {
    
    fitness = c(meanFitTrajHome, meanFitTrajNonHome)
    min = c(meanFitTrajHome - sqrt(varHome_eachGen), meanFitTrajNonHome - sqrt(varNonHome_eachGen)) # min of ribbon is -1 SD
    max = c(meanFitTrajHome + sqrt(varHome_eachGen), meanFitTrajNonHome + sqrt(varNonHome_eachGen)) # max of ribbon is +1 SD
    time = rep(1:numGens, times = 2)
    Env = as.factor(c(rep(1, times = numGens), rep(2, times = numGens)))
    WbarDf = data.frame(time, fitness, min, max, Env)
    
    
    print( ggplot(WbarDf, aes(x=as.numeric(time),y= fitness, colour = Env )) + geom_line()  +
             geom_ribbon(aes(x= time, ymin=min, ymax=max,fill = Env,group = Env), alpha = 0.3, linetype = 0) + 
             scale_color_manual(values=c('chocolate1','dodgerblue2'))+  geom_hline(yintercept =  0, linetype = 'dashed')+
             theme_bw()+theme(panel.border = element_rect(colour = 'black', fill = NA, size = 2), 
                              panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
                              axis.text =element_text(size = 20,family = 'Helvetica'), axis.title = element_blank(), legend.position = 'none')+ 
             xlab("Time")+ ylab("Fitness"))
  }
  
  
  colResID_allTimePoints =  (allMeanFitTrajNonHome > 0) # will give array of same dim as allmeanfit traj with 1 for CR and 0 for CS
  
  CR_overtime = c()
  for(check in 1:numGens)
  {
    CR_overtime = c(CR_overtime, sum(colResID_allTimePoints[check,] == 1)/ length(colResID_allTimePoints[check,]))
  }
  
  output = list(allMeanFitTrajNonHome, allMeanFitTrajHome, meanFitTrajHome, meanFitTrajNonHome, varHome_eachGen, varNonHome_eachGen, covar_eachGen, data.frame(slopeHome, slopeNonHome,varSlopeHome, varSlopeNonHome,  covarSlope,r1,r2,D11,D12,D22), CR_overtime)
  
  
  return(output)
  
  
  
}


## have GAUSSIAN JDFE with home parameters
sigHome = sigHome #home DFE SD 
uHome =  uHome# mean home DFE
UbAdjust =  1 - pnorm(0, uHome, sigHome)


ColRes_proportion_overTime = c()
drugID_fordf = c()
r2_all = c()
D22_all = c()
homeSlope_vec = c()
nonHomeSlope_vec = c()
varslopeHome_vec = c()
varSlopeNonHome_vec = c()
timeVec = c()

for(drug in 1:numDrugs)
{
  # JDFE Parameters
  sigNonHome = sigNonHomeVec[drug]
  uNonHome =  uNonHomeVec[drug]
  corr = corVec[drug]
  
  covar = corr*sqrt(sigHome^2 * sigNonHome^2)
  uVec = c(uHome, uNonHome)
  Sigma =  matrix(c(sigHome^2,covar,covar,sigNonHome^2), nrow=2)
  UbAdjust = 1 - pnorm(0, uHome, sigHome)
  
  # Run evo Sim 
  output = simEvo_WrightFisher_rank(uVec, Sigma, popsize, U, numiterations, numGens, wantPlot = TRUE) 
  # list(allMeanFitTrajNonHome, allMeanFitTrajHome, meanFitTrajHome, meanFitTrajNonHome, varHome_eachGen, varNonHome_eachGen, covar_eachGen, data.frame(slopeHome, slopeNonHome,varSlopeHome, varSlopeNonHome,  covarSlope,r1,r2,D11,D12,D22), CR_overtime)
  data_frame_out = output[[8]]
  
  
  ColRes_proportion_overTime = c(ColRes_proportion_overTime, output[[9]])
  drugID_fordf = c(drugID_fordf, rep(drug, times = length( output[[9]])))
  r2_all = c(r2_all, rep(data_frame_out$r2, times = length( output[[9]])))
  D22_all = c(D22_all, rep(data_frame_out$D22, times = length( output[[9]])))
  homeSlope_vec = c(homeSlope_vec, rep(data_frame_out$slopeHome, times = length( output[[9]])))
  nonHomeSlope_vec = c(nonHomeSlope_vec, rep(data_frame_out$slopeNonHome, times = length( output[[9]])))
  varslopeHome_vec = c(varslopeHome_vec, rep(data_frame_out$varSlopeHome, times = length( output[[9]])))
  varSlopeNonHome_vec = c(varSlopeNonHome_vec, rep(data_frame_out$varSlopeNonHome, times = length( output[[9]])))
  timeVec = c(timeVec, c(1:numGens))
}

c_all = r2_all/sqrt(D22_all)


rankOrder_chosenJDFEs_out = data.frame(timeVec, as.factor(drugID_fordf), r2_all, D22_all, homeSlope_vec, nonHomeSlope_vec, varslopeHome_vec, varSlopeNonHome_vec,  ColRes_proportion_overTime)
colnames(rankOrder_chosenJDFEs_out) = c('Time', "Drug", "r2", "D22", "homeSlope", 'nonHomeSlope', "varSlopeHome", 'varSlopeNonHome', 'CR_proportion')
write.csv(rankOrder_chosenJDFEs_out, paste0('rankOrder_chosenJDFEs_df_NU', popsize*U, '_numGens_', numGens, '_numiterations_', numiterations, '.csv'))


ColResOverTime_plot = ggplot(rankOrder_chosenJDFEs_out, aes(x = Time, y = CR_proportion, color = Drug, group = Drug )) +
  geom_line(size = 2)+ xlim(c(1,numGens))+
  scale_color_manual(values = c("#fe0059" ,"#fea500" ,"#0059fe" ,"#00e594"))+
  theme_bw()+theme(panel.border = element_rect(colour = 'black', fill = NA, size = 1), legend.position = 'none',
                   panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
                   axis.text =element_text(size = 15,family = 'Helvetica', color = 'black'), axis.title = element_text(size = 20, family = 'Helvetica', colour = 'black'))



ggsave(paste0('ColResOverTime_plot_selectedJDFEs_NU', popsize*U, '_numGens_', numGens, '_numiterations_', numiterations, '.pdf'), ColResOverTime_plot, width =3.89, height = 3.67)




###### Plot Rank Order from all JDFEs #######

time = 300 # time in sim to test 

df1 = subset(data_WrightFisherSims_allNu, NU == 1) 
df10 =subset(data_WrightFisherSims_allNu, NU == 10) 
df100 = subset(data_WrightFisherSims_allNu, NU == 100) 
allR2 = subset(data_WrightFisherSims_allNu, NU == 1)$r2_SSWM # NU doesnt matter
allD22 = subset(data_WrightFisherSims_allNu, NU == 1)$D22_SSWM # NU doesnt matter

# calculate the observed percent col sens at the given time
chanceCollSens_SimOutput_allJDFEs_nu1 = pnorm(0, mean= (df1$nonHomeSlope )*time , sd  =  sqrt((df1$nonHomeVarSlope)*time ))
chanceCollSens_SimOutput_allJDFEs_nu10 =  pnorm(0, mean= (df10$nonHomeSlope )*time , sd  =  sqrt((df10$nonHomeVarSlope)*time ))
chanceCollSens_SimOutput_allJDFEs_nu100 =  pnorm(0, mean= (df100$nonHomeSlope )*time , sd  =  sqrt((df100$nonHomeVarSlope)*time ))

# calculate the predicted percent col sens and c value
chanceCollSens_Pred_allJDFEs = pnorm(0, mean= (allR2)*time , sd  =  sqrt((allD22)*time ))
c_allJDFEs = df1$r2_SSWM / sqrt(df1$D22_SSWM)


# calculate the rank of JDFEs in sim 
rank_sim_1 = rank(1/chanceCollSens_SimOutput_allJDFEs_nu1)
rank_sim_10 = rank(1/chanceCollSens_SimOutput_allJDFEs_nu10)
rank_sim_100 = rank(1/chanceCollSens_SimOutput_allJDFEs_nu100)
allrank_sim = c( rank_sim_1, rank_sim_10, rank_sim_100)

# calculate the predicted rank 
rank_pred = rep(rank(c_allJDFEs), times = 3)

NuID = as.character(c(rep(1, times = length(df1$NU)) ,rep(10, times = length(df1$NU)) ,rep(100, times = length(df1$NU))    ) )


##plot for rank
plotDFrank = data.frame(allrank_sim, rank_pred, NuID)

rankPlot_all125JDFEs = ggplot(plotDFrank, aes(x = rank_pred, y = allrank_sim, color = NuID)) +
  geom_point(size = 7, alpha = 0.6) + geom_abline(size  = 1.5, intercept = 0, slope = 1, linetype = 'dashed', col = 'grey')+
  theme_bw()+theme(panel.border = element_rect(colour = 'black', fill = NA, size = 1.5), 
                   panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
                   axis.text.x =element_text(size = 30,family = 'Helvetica', color = c('black')),axis.text.y =element_text(size = 30,family = 'Helvetica', color = c('black')), axis.title = element_blank()) +
  scale_color_manual(values =c( 'palegreen4', 'darkgoldenrod', 'hotpink4'))

ggsave(paste0("rankPlot_all125JDFEs_gen", time, ".pdf") ,rankPlot_all125JDFEs,width= 5.91, height=4.42)


allpercDif = c()
for(i in 1:(length(test)-1))
{
  c1 = abs(test[i])
  c2 = abs(test[i+1])
  percDiff = abs(c2-c1)/ ((c1+c2)/2) * 100
  allpercDif = c(allpercDif,percDiff)
}

###########################################################################################
######################  FIGURE 6: HOW    TO   MEASURE    JDFEs    ######################
########################################################################################
## this section does use results from fgure 3 because needs to call actual sim slope 
# choose working directory 
  # setwd("~/")


#####    Luria Delbruk   ####
# Random sampling above a given fitness threshold 


####  All JDFEs ####
  
# home parameters
sigHome = 0.01 #home DFE SD 
uHome = - 0.001 # mean home DFE
UbAdjust = 1 - pnorm(0, uHome, sigHome)
focalNU = 100 # only affects which slope compare to 
sigThresh = 0.5
sampThreshold = sigThresh*sigHome # in main text show 0.5 sigma, in supp have  1, 2 and 3

data_thisNu = subset(data_WrightFisherSims_allNu, NU == focalNU)


numSamples = c(1,10,100,1000)
numReps = 100
allTestr2 = c()
allTestD22 = c()
allTestD12 = c()
sampleVec = c()
allr2 = c()
allD22= c()
allD12 = c()
jdfeID = c()
allsimSlope = c()
allsimVar = c()

for(samp in numSamples)
{
  count = 0
  for(sigNonHome in seq(0.0001, 0.01, length.out =5))
  {
    for(uNonHome in seq(0.0001, 0.01, length.out = 5))
    {
      for(corr in seq(-0.9, 0.9, length.out = 5))
      {
        count = count + 1
        print(count)
        #for( rep in 1:numReps)
       # {
          
          
          jdfeID = c(jdfeID, count)
          sampleVec = c(sampleVec, samp)
          ################## M A K E     T H E      J D F E
          cov = corr*sqrt(sigHome^2 * sigNonHome^2)
          Sigma = array(data = c(sigHome^2, cov, cov, sigNonHome^2), dim = c(2,2))
          
          
          
          # Only sample "SAMP" number of muts
          
          #start with initial try 
          
          #how many samples to take to try to get right number in 1 try
          percentile_toreach =  pnorm(sigThresh*sigHome + uHome , mean = uHome, sd = sigHome)
          howManySampsNeeded_for1correct = 1 / (1-percentile_toreach )
          howManySamps_toTake = ceiling(howManySampsNeeded_for1correct)*samp*2 # do x2 to be extra sure
          
          trial = as.array(mvrnorm(howManySamps_toTake, c(uHome,uNonHome), Sigma))
          
          homeS_all = trial[,1]
          nonHomeS_all = trial[,2]
          
          homeS_toKeep = homeS_all[homeS_all > sampThreshold]
          nonHomeS_toKeep = nonHomeS_all[homeS_all > sampThreshold]
          
          homeSvals =  homeS_toKeep
          nonHomeSvals = nonHomeS_toKeep
          
          
          # If not enough, Run while loop until get right num muts
          num = 0 
          while(length(homeSvals) < samp)
          {
            trial = as.array(mvrnorm(1, c(uHome,uNonHome), Sigma))
            
            if(trial[1] > sampThreshold)
            {
              homeSvals = c(homeSvals, trial[1])
              nonHomeSvals = c(nonHomeSvals, trial[2])
            }
          }
          
          
          #make sure are of exact right length
          homeSvals = homeSvals[1:samp]
          nonHomeSvals = nonHomeSvals[1:samp]
          
          
          homeSvals =  homeSvals*( homeSvals>0)
          r1_test = mean((homeSvals)^2)
          r2_test = mean((homeSvals)*(nonHomeSvals))
          D11_test = mean((homeSvals)^3)
          D12_test = mean((homeSvals)^2*(nonHomeSvals))
          D22_test = mean((homeSvals)*(nonHomeSvals)^2)
          
          allTestr2 = c(allTestr2, r2_test)
          allTestD22 = c(allTestD22, D22_test)
          allTestD12 = c(allTestD12, D12_test)
          
          
          #### full approx, 100,000 muts sampled, gives 'true' r2
          test = as.array(mvrnorm(100000, c(uHome,uNonHome), Sigma))
          homeSvals = test[,1]*(test[,1]>0)
          nonHomeSvals = test[,2]
          r1_test = mean((homeSvals)^2)
          r2_test = mean((homeSvals)*(nonHomeSvals))
          D11_test = mean((homeSvals)^3)
          D12_test = mean((homeSvals)^2*(nonHomeSvals))
          D22_test = mean((homeSvals)*(nonHomeSvals)^2)
          
          allr2 = c(allr2, r2_test)
          allD22 = c(allD22, D22_test)
          allD12 = c(allD12, D12_test)
          
          allsimSlope = c(allsimSlope, data_thisNu$nonHomeSlope[count])
          allsimVar = c(allsimVar, data_thisNu$nonHomeVarSlope[count])
        #}
        
      }
    }
  }
}

# show plot OF MEASURED VS PTRUE R2 FOR ALL JDFES
    # keep in mind, this plot is only showing 1 itaeration of measuring for 
    # each JDFE in each sample number, so will change somewhat each time is run
    
df = data.frame(as.factor(sampleVec), allTestr2, allTestD12, allTestD22, allr2, allD22, allD12, jdfeID, allsimSlope, allsimVar)
colnames(df) = c('Sample', 'r2_approx', 'D12_approx', 'D22_approx', 'r2', 'D22', 'D12', 'JDFE', 'slopeSim', 'varSlopeSim')
df$c_approx = df$r2_approx/ sqrt(df$D22_approx)
df$c_real = df$r2 / sqrt(df$D22)

df$rank_approx = df$c_approx
df$rank_real = df$c_real

lmRsq_out  = c()
numSampID = c()
for(i in numSamples)
{
  df[df$Sample == i, ]$rank_approx = rank(df[df$Sample == i, ]$c_approx)
  df[df$Sample == i, ]$rank_real = rank(df[df$Sample == i, ]$c_real)

  real_rank = df[df$Sample == i, ]$rank_real
  approx_rank =  df[df$Sample == i, ]$rank_approx
  lm = summary(lm(real_rank ~ approx_rank))
  
  lmRsq_out = c(lmRsq_out, lm$adj.r.squared)
}
lm_df = data.frame(numSamples, lmRsq_out)
colnames(lm_df) = c('NumSamples_LD', 'Rsq_rankRealvEstimated')


write.csv(df, file = paste0("allJDFE_LuriaDelbrukOutputs_sig_", sigThresh, ".csv"))
write.csv(lm_df, file =paste0("LuriaDelbruk_Rsq_lms_real_v_estimated_rank_sig_", sigThresh, ".csv"))



all125JDFEs_r2simvPred_randomMuts = ggplot(df, aes(x = r2_approx , y = r2, color = Sample)) + 
  geom_hline(yintercept = 0 , color = 'grey', size = 1.5) +geom_vline(xintercept = 0, color = 'grey', size= 1.5) +
  geom_point( size = 7, alpha = 0.6) +
  geom_abline(size  = 1.5, intercept = 0, slope = 1, linetype = 'dashed', col = 'black')+
  theme_bw()+theme(panel.border = element_rect(colour = 'black', fill = NA, size = 1.5), 
                   panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.position = 'none',
                   axis.text.x =element_text(size = 20,family = 'Helvetica', color = c('black','transparent')),axis.text.y =element_text(size = 20,family = 'Helvetica', color = c('black')), axis.title = element_blank()) +
  scale_y_continuous(label=scientific_10) + scale_x_continuous(label=scientific_10)

ggsave( paste0("all125JDFEs_r2simvPred_LuriaDelbruck_sig_", sigThresh, ".pdf"), all125JDFEs_r2simvPred_randomMuts,width = 4.43, height = 3.38 )


all125JDFEs_D22simvPred_randomMuts = ggplot(df, aes(x = D22_approx , y = D22, color = Sample)) + 
  geom_hline(yintercept = 0 , color = 'grey', size = 1.5) +geom_vline(xintercept = 0, color = 'grey', size= 1.5) +
  geom_point( size = 7, alpha = 0.6) +
  geom_abline(size  = 1.5, intercept = 0, slope = 1, linetype = 'dashed', col = 'black')+
  theme_bw()+theme(panel.border = element_rect(colour = 'black', fill = NA, size = 1.5), 
                   panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.position = 'none',
                   axis.text.x =element_text(size = 20,family = 'Helvetica', color = c('black','transparent')),axis.text.y =element_text(size = 20,family = 'Helvetica', color = c('black')), axis.title = element_blank()) +
  scale_y_continuous(label=scientific_10) + scale_x_continuous(label=scientific_10)

ggsave( paste0("all125JDFEs_D22simvPred_LuriaDelbruck_sig_", sigThresh, ".pdf"), all125JDFEs_D22simvPred_randomMuts, width = 4.43, height = 3.38)



# rank for all JDFEs from NU 100 sim, uses the slope + var in sim at gen 300 to rank
rank_sim_1 = rank(1/chanceCollSens_SimOutput_allJDFEs_nu1)
df$rank_sim = rep(rank_sim_100, times = 4)

all125JDFEs_rankApprox_v_rankReal_randomMuts = ggplot(df, aes(x = rank_approx , y =rank_sim, color = as.factor(Sample))) + 
  geom_point( size = 7, alpha = 0.6) +
  geom_abline(size  = 2.5, intercept = 0, slope = 1, linetype = 'dashed', col = 'black')+
  theme_bw()+theme(panel.border = element_rect(colour = 'black', fill = NA, size = 1.5), 
                   panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                   axis.line = element_line(colour = "black"), legend.position = 'none',
                   axis.text.x =element_text(size = 20,family = 'Helvetica', color = c('black')),
                   axis.text.y =element_text(size = 20,family = 'Helvetica', color = c('black')), 
                   axis.title = element_blank()) 

ggsave( paste0("all125JDFEs_rankApprox_v_rankReal_LuriaDelbruck_sig_", sigThresh, ".pdf"), all125JDFEs_rankApprox_v_rankReal_randomMuts, width = 3.65, height = 3.45)





####  Lineage Tracking ####



### All JDFES  ####

# this section relies on findings from figure 3 and that data frame being saved
focalNU = 100 # choose which NU (from those in fig 3) to use for real slope and sims
numMutsToSampleVec = c(1, 10, 100,1000 )# vector with each number of mutations to sample 

###
sigHome = 0.01 #home DFE SD 
uHome = - 0.001 # mean home DFE
U = 10^-4
popsize = focalNU/U
UbAdjust = 1 - pnorm(0, uHome, sigHome)
numGens = 250
#mutThreshold = 22*(1/popsize) # threshold necessary for mutation to reach to be sampled

data_thisNu = subset(data_WrightFisherSims_allNu, NU ==focalNU)


r2_evoGuess_total = c()
D22_evoGuess_total = c()
D12_evoGuess_total = c()
r2_evoGuess_mostrecentMut = c()
D22_evoGuess_mostrecentMut = c()
D12_evoGuess_mostrecentMut = c()
realr2Vec =c()
realD22Vec = c()
realSimSlopeVec = c()
numMutsTosampleVec = c()
jdfeID = c()

averageNonHomeEffect_forsampledHomeBenmuts_total = c()
averageNonHomeEffect_forsampledHomeBenmuts_mostrecent = c()
var_sampledbenMuts_total = c()
var_sampledbenMuts_mostrecent = c()

count = 0
for(sigNonHome in seq(0.0001, 0.01, length.out =5))
{
  for(uNonHome in seq(0.0001, 0.01, length.out = 5))
  {
    for(corr in seq(-0.9, 0.9, length.out = 5))
    {
      count = count + 1
      
      print('count')
      print(count)
      
      ### M A K E     T H E      J D F E
      sigNonHome = data_thisNu$sdNonHome[count]
      uNonHome =  data_thisNu$uNonHome[count]
      corr = data_thisNu$corr[count]
      
      simSlope_real = data_thisNu$nonHomeSlope[count]
      truer2 = data_thisNu$r2_SSWM[count]
      trueD22 = data_thisNu$D22_SSWM[count]
      
      UbAdjust = 1 - pnorm(0, uHome, sigHome)
      cov = corr*sqrt(sigHome^2 * sigNonHome^2)
      
      ### FIT MVRNORM TO DATA 
      Sigma = array(data = c(sigHome^2, cov, cov, sigNonHome^2), dim = c(2,2))
      sValallTheoMuts= mvrnorm(n = 100000, c(uHome,uNonHome), Sigma)
      
      for(numMuts in numMutsToSampleVec)
      {
        
        print(numMuts)
        
        sampledMutHomeGR_total = c()
        sampledMutHomeGR_mostrecentMut = c()
        
        sampledMutNonHomeGR_total = c()
        sampledMutNonHomeGR_mostrecentMut = c()
        
        
        
        
        
        allIDArrayList = list()
        
        
        sampMutCount = 0
        gen = 0
        
        
        numiterations = 1
        numToSamp_perIt = numMuts
        if(numMuts > 50 ) # were going to do independent trials of the sims and sample max 50 each time
        {
          numiterations = numMuts/100
          numToSamp_perIt = 100
        }
           # when samp w replace
       # for(iteration in 1:numiterations)
        while(length(sampledMutHomeGR_total) < numMuts)
        {
          
          ID_df  = as.data.frame( array( data = c(1, uHome, uNonHome,popsize, 0, 0, 0, uHome, uNonHome), dim = c(1,9)))
          colnames(ID_df) = c('Strain', 'HomeGR', 'NonHomeGR', 'Count', 'isNewMut', 'numMutsThisLine', 'sampled', 'parentFitHome', 'parentFitNonHome')
          
          
          
          for(gen in 1:numGens)
          {
            # gen = gen + 1
            #####    S  E  L  E  C   T  I  O  N    
            #---------------------------------------------------------------------------------------_#
            # select next gen offspring
            meanFitHome = (1/popsize)*sum(ID_df$Count * ID_df$HomeGR)
            meanFitNonHome = (1/popsize)*sum(ID_df$Count * ID_df$NonHomeGR)
            
            
            
            ####    r muultinom
            probVec = ID_df$Count / popsize + (ID_df$Count / popsize)*(ID_df$HomeGR -meanFitHome)
            #  probVec = (1 + ID_df$HomeGR/meanFit)*ID_df$Count # from ben good
            # all neg probs = 0 prob
            probVec[probVec<0] = 0
            offspring = rmultinom(1, size = popsize, prob = probVec)
            #   
            
            
            
            # update ID df with new count
            ID_df$Count = offspring 
            
            # remove any extinct lineages
            ID_df = ID_df[!ID_df$Count == 0, ]
            
            # re-do counting so dont have issues calling indexes of parents
            ID_df$Strain = 1:length(ID_df$Strain)
            
            
            
            #### M  U  T  A  T  I  O  N  S
            #---------------------------------------------#
            
            # Sample number of muts that will happen and then choose parent IDs
            M = rpois(1, popsize*U)
            
            if(M > 0 )
            {
              parentIDs = sample(ID_df$Strain, M, prob = (ID_df$Count/popsize), replace = TRUE)
              
              
              # sample M mutation effects from JDFE, usimg Sigma (covariance matrix ) from above fit
              mutEffects = as.array(mvrnorm(M, c(uHome,uNonHome), Sigma))
              
              if(M ==1) # to correct for a single sample frommvrnorm being read as a vector
              {
                mutEffects = array(mutEffects, dim = c(1,2))
              }
              
              # make M new rows for ID_df, one for each mut (4 cols, ID#, homeFit, respFit (!!!!!which not doing anything with RN!!!!))
              newMutArray = as.data.frame(array(data = 0, dim = c(M, 9)))
              colnames(newMutArray) = c('Strain', 'HomeGR', 'NonHomeGR', 'Count', 'isNewMut', 'numMutsThisLine', 'sampled', 'parentFitHome', 'parentFitNonHome')
              
              # make new mutant ID's ,, just count up from highest current ID
              newMutArray$Strain = (max(ID_df$Strain)+1):(max(ID_df$Strain) + M )
              
              # update new Mut Fitness values by adding sampled fit Effects
              newMutArray$HomeGR = ID_df$HomeGR[parentIDs] + mutEffects[,1]
              newMutArray$NonHomeGR = ID_df$NonHomeGR[parentIDs] + mutEffects[,2]
              newMutArray$parentFitHome = ID_df$HomeGR[parentIDs]
              newMutArray$parentFitNonHome = ID_df$NonHomeGR[parentIDs]
              
              
              # add count of ONE to each new mut, subtract offspring from parent
              newMutArray$Count = rep(1, times = M)
              newMutArray$isNewMut = 1
              newMutArray$numMutsThisLine = ID_df$numMutsThisLine[parentIDs] + 1
              newMutArray$sampled = rep(0, times = M)
              
              ####??? faster way to do this
              for(rem in parentIDs)
              {
                ID_df$Count[rem] = ID_df$Count[rem] -1
              }
              
              # collate old Id_df and new muts
              ID_df = rbind(ID_df, newMutArray)
              
            }# M>0 if statement
            
            
             
            hasBeenSampled = numeric(length(ID_df$Strain))
         
            
            ID_df$sampled = hasBeenSampled
            
            
             
            
          } # gen loop 
          
          #### sample mutants from last gen 
          
          #  find ben Mut IDs 
          benMutIDs = (ID_df$HomeGR - ID_df$parentFitHome)  > 0
          
          # sample as many muts as supposed to , up to 50 
             
          # without replacement
          if(sum(benMutIDs)> numToSamp_perIt)
          {
            sampleIDS = sample(ID_df$Strain[benMutIDs], numToSamp_perIt, replace = FALSE, prob = ID_df$Count[benMutIDs] )
            
          }else{
            sampleIDS = sample(ID_df$Strain[benMutIDs], sum(benMutIDs), replace = FALSE, prob = ID_df$Count[benMutIDs] )
            
          }
          
          sampledMutHomeGR_total = c(sampledMutHomeGR_total, ID_df$HomeGR[sampleIDS])
          sampledMutHomeGR_mostrecentMut = c( sampledMutHomeGR_mostrecentMut , ( ID_df$HomeGR[sampleIDS] - ID_df$parentFitHome[sampleIDS]))
          
          sampledMutNonHomeGR_total = c(sampledMutNonHomeGR_total, ID_df$NonHomeGR[sampleIDS])
          sampledMutNonHomeGR_mostrecentMut = c(sampledMutNonHomeGR_mostrecentMut, ( ID_df$NonHomeGR[sampleIDS] - ID_df$parentFitNonHome[sampleIDS]))
          
          
          
           
          
          
        } # iterationloop
        
       
        
        # Save for ouput
        realr2Vec =c(realr2Vec, truer2)
        realD22Vec = c(realD22Vec, trueD22)
        realSimSlopeVec = c(realSimSlopeVec, simSlope_real)
        
        
        ## calc the pleiotropy parameters based on these saved muts USING TOTALS
        homeSvals = sampledMutHomeGR_total*(sampledMutHomeGR_total>=0) 
        nonHomeSvals = sampledMutNonHomeGR_total
        r1 = mean((homeSvals)^2)
        r2 = mean((homeSvals)*(nonHomeSvals))
        D11 = mean((homeSvals)^3)
        D12 = mean((homeSvals)^2*(nonHomeSvals))
        D22 = mean((homeSvals)*(nonHomeSvals)^2)
        
        
        r2_evoGuess_total = c(r2_evoGuess_total, r2)
        D22_evoGuess_total = c(D22_evoGuess_total, D22)
        D12_evoGuess_total = c(D12_evoGuess_total, D12)
        
        
        ## calc the pleiotropy parameters based on these saved muts USING MOST RECENT MUTS
        homeSvals = sampledMutHomeGR_mostrecentMut*(sampledMutHomeGR_mostrecentMut>=0) 
        nonHomeSvals = sampledMutNonHomeGR_mostrecentMut
        r1 = mean((homeSvals)^2)
        r2 = mean((homeSvals)*(nonHomeSvals))
        D11 = mean((homeSvals)^3)
        D12 = mean((homeSvals)^2*(nonHomeSvals))
        D22 = mean((homeSvals)*(nonHomeSvals)^2)
        
        r2_evoGuess_mostrecentMut = c(r2_evoGuess_mostrecentMut, r2)
        D22_evoGuess_mostrecentMut = c(D22_evoGuess_mostrecentMut, D22)
        D12_evoGuess_mostrecentMut = c(D12_evoGuess_mostrecentMut, D12)
        
        
        
        jdfeID = c(jdfeID, count)
        numMutsTosampleVec = c( numMutsTosampleVec , numMuts)
        
        
        ### 0ther parameters to try 

        
         averageNonHomeEffect_forsampledHomeBenmuts_total = c(averageNonHomeEffect_forsampledHomeBenmuts_total, mean(sampledMutNonHomeGR_total[sampledMutHomeGR_total>0]))
         averageNonHomeEffect_forsampledHomeBenmuts_mostrecent = c(averageNonHomeEffect_forsampledHomeBenmuts_mostrecent, mean(sampledMutNonHomeGR_mostrecentMut[sampledMutHomeGR_mostrecentMut>0]))
         
         var_sampledbenMuts_total = c(var_sampledbenMuts_total , var(sampledMutNonHomeGR_total[sampledMutHomeGR_total>0]))
         var_sampledbenMuts_mostrecent = c(var_sampledbenMuts_mostrecent, var(sampledMutNonHomeGR_mostrecentMut[sampledMutHomeGR_mostrecentMut>0]))
        
      } # numMuts Vec loop 
    } # corr loop
  } # uNonHome loop
} # sigNonHome loop


# Save data

df_lineageTracking_ALLJDFEs =  data.frame(jdfeID, numMutsTosampleVec, r2_evoGuess_total, D22_evoGuess_total, r2_evoGuess_mostrecentMut, D22_evoGuess_mostrecentMut, realr2Vec, realD22Vec ,realSimSlopeVec,  averageNonHomeEffect_forsampledHomeBenmuts_total,  averageNonHomeEffect_forsampledHomeBenmuts_mostrecent, var_sampledbenMuts_total, var_sampledbenMuts_mostrecent )
colnames(df_lineageTracking_ALLJDFEs) = c(  'JDFE', 'numMuts', 'r2_guess_total', 'D22_guess_total','r2_guess_mostrecentMut', 'D22_guess_mostrecentMut', 'r2_SSWM_real', 'D22_SWWM_real', 'realSimSlope', 'averageNonHomeEffect_forsampledHomeBenmuts_total', 'averageNonHomeEffect_forsampledHomeBenmuts_mostrecent', 'var_sampledbenMuts_total', 'var_sampledbenMuts_mostrecent')
df_lineageTracking_ALLJDFEs$numMuts = as.factor(df_lineageTracking_ALLJDFEs$numMuts)

df_lineageTracking_ALLJDFEs$c_approx_total = df_lineageTracking_ALLJDFEs$r2_guess_total / sqrt(df_lineageTracking_ALLJDFEs$D22_guess_total )
df_lineageTracking_ALLJDFEs$c_approx_mostrecent = df_lineageTracking_ALLJDFEs$r2_guess_mostrecentMut / sqrt(df_lineageTracking_ALLJDFEs$D22_guess_mostrecentMut )
df_lineageTracking_ALLJDFEs$c_real = df_lineageTracking_ALLJDFEs$r2_SSWM_real / sqrt(df_lineageTracking_ALLJDFEs$D22_SWWM_real) 
df_lineageTracking_ALLJDFEs$c_apriori = df_lineageTracking_ALLJDFEs$averageNonHomeEffect_forsampledHomeBenmuts_mostrecent / sqrt(df_lineageTracking_ALLJDFEs$var_sampledbenMuts_mostrecent) 


df_lineageTracking_ALLJDFEs$rank_approx_total = rank(df_lineageTracking_ALLJDFEs$c_approx_total)
df_lineageTracking_ALLJDFEs$rank_approx_mostrecent = rank(df_lineageTracking_ALLJDFEs$c_approx_mostrecent)
df_lineageTracking_ALLJDFEs$rank_real = rank(df_lineageTracking_ALLJDFEs$c_real)
df_lineageTracking_ALLJDFEs$rank_apriori = rank(df_lineageTracking_ALLJDFEs$c_apriori)


lm_rsq_mostrecent  = c()
lm_rsq_total = c()
lm_rsq_apriori = c()
for(i in numMutsToSampleVec)
{
  df_lineageTracking_ALLJDFEs[df_lineageTracking_ALLJDFEs$numMuts == i, ]$rank_approx_total = rank(df_lineageTracking_ALLJDFEs[df_lineageTracking_ALLJDFEs$numMuts == i, ]$c_approx_total)
  df_lineageTracking_ALLJDFEs[df_lineageTracking_ALLJDFEs$numMuts == i, ]$rank_approx_mostrecent = rank(df_lineageTracking_ALLJDFEs[df_lineageTracking_ALLJDFEs$numMuts == i, ]$c_approx_mostrecent)
  df_lineageTracking_ALLJDFEs[df_lineageTracking_ALLJDFEs$numMuts == i, ]$rank_real = rank(df_lineageTracking_ALLJDFEs[df_lineageTracking_ALLJDFEs$numMuts == i, ]$c_real)
  df_lineageTracking_ALLJDFEs[df_lineageTracking_ALLJDFEs$numMuts == i, ]$rank_apriori = rank(df_lineageTracking_ALLJDFEs[df_lineageTracking_ALLJDFEs$numMuts == i, ]$c_apriori)
  
  # total 
  real_rank = df_lineageTracking_ALLJDFEs[df_lineageTracking_ALLJDFEs$numMuts == i, ]$rank_real
  approx_rank_total =  df_lineageTracking_ALLJDFEs[df_lineageTracking_ALLJDFEs$numMuts == i, ]$rank_approx_total
  lm = summary(lm(real_rank ~ approx_rank_total))
  
  lm_rsq_total= c(lm_rsq_total, lm$adj.r.squared)
  
  
  # most recent
  real_rank = df_lineageTracking_ALLJDFEs[df_lineageTracking_ALLJDFEs$numMuts == i, ]$rank_real
  approx_rank_mostrecent =  df_lineageTracking_ALLJDFEs[df_lineageTracking_ALLJDFEs$numMuts == i, ]$rank_approx_mostrecent
  lm = summary(lm(real_rank ~ approx_rank_mostrecent))
  
  lm_rsq_mostrecent= c(lm_rsq_mostrecent, lm$adj.r.squared)
  
  # estimate - avg for r2 and var for d22
  real_rank = df_lineageTracking_ALLJDFEs[df_lineageTracking_ALLJDFEs$numMuts == i, ]$rank_real
  approx_rank_apriori =  df_lineageTracking_ALLJDFEs[df_lineageTracking_ALLJDFEs$numMuts == i, ]$rank_apriori
  lm = summary(lm(real_rank ~ approx_rank_apriori))
  
  lm_rsq_apriori= c(lm_rsq_apriori, lm$adj.r.squared)
  
  
  
}
lm_df_LT = data.frame(numMutsToSampleVec, lm_rsq_total, lm_rsq_mostrecent, lm_rsq_apriori)
colnames(lm_df_LT) = c('NumSamples_LD', 'total_Rsq_rankRealvEstimated', 'mostrecent_Rsq_rankRealvEstimated','apriori_Rsq_rankRealvEstimated' )


write.csv(lm_df_LT, file = 'LineageTracking_Rsq_lms_real_v_estimated_rank.csv')
write.csv( df_lineageTracking_ALLJDFEs ,paste0('df_lineageTracking_ALLJDFEs_NU', focalNU,'_allSampNums',  '.csv'))



#### Plots


# Real vs approx rank 
#dft = subset(df_lineageTracking_ALLJDFEs, numMuts = 10000)
all125JDFEs_rankApprox_v_rankReal_lineageTrack = ggplot(df_lineageTracking_ALLJDFEs, aes(x = rank_approx_mostrecent , y =rank_real, color = as.factor(numMuts))) + 
  geom_point( size = 7, alpha = 0.6) +
  geom_abline(size  = 2.5, intercept = 0, slope = 1, linetype = 'dashed', col = 'black')+
  theme_bw()+theme(panel.border = element_rect(colour = 'black', fill = NA, size = 1.5), 
                   panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                   axis.line = element_line(colour = "black"), legend.position = 'none',
                   axis.text.x =element_text(size = 20,family = 'Helvetica', color = c('black')),
                   axis.text.y =element_text(size = 20,family = 'Helvetica', color = c('black')), 
                   axis.title = element_blank()) #+
#scale_y_continuous(label=scientific_10) + scale_x_continuous(label=scientific_10)



ggsave( "all125JDFEs_rankApprox_v_rankReal_lineageTrack.pdf",all125JDFEs_rankApprox_v_rankReal_lineageTrack, width = 3.65, height = 3.45 )




###############################################################################################################
######################################### APPENDIX : EPISTASIS  ##########################################
###############################################################################################################


# (optional) set working directory to where you'd like to save the outputs of this section, example below 
# setwd("~/Documents/JDFE Project/Figures/Figure Epistasis")

#### First, Set rules on how JDFEs change with fitness  #####
# We will have JDFE mean stay same, but variance decrease with fitness in both enviornments 


# Function to determine how jdfe SD changes , input curren mean fitness
  # plus initial DFEmean, gamma (rate of decrease) and xmax (max fitness )
newSDVec<- function(meanFitHome, meanFitNonHome,  sigHomeInitial,sigNonHomeInitial, gamma1, gamma2)
{
  
  newSDhome = sigHomeInitial - gamma1*meanFitHome
  newSDNonhome = sigNonHomeInitial  - gamma2*meanFitNonHome
  return(c( newSDhome, newSDNonhome))
}

## Set parameters
gamma1 = 0.5
gamma2 = 0.5
popsize = 10^4
U = 10^-4
numiterations = 1500
numGens = 1000
wantPlot = TRUE

# test how the function behaves
sdOverFitness = c()
fitVec = c()
for(i in seq(0,0.01, length.out = 30))
{
  meanFitHome = i
  meanFitNonHome = i
  
  sdHome = newSDVec(meanFitHome, meanFitNonHome,  sigHomeInitial,sigNonHomeInitial, gamma1, gamma2)[1]
  sdOverFitness = c(sdOverFitness, sdHome)
  fitVec = c(fitVec, i)
}
plot(fitVec, sdOverFitness)



#### Evolution Simulation  ####
simEvo_WrightFisher_plusEpistasis <- function(uVec, Sigma, popsize, U, numiterations, numGens, wantPlot, gamma1, gamma2) 
{
  
  uHome = uVec[1]  # mean of home DFE
  uNonHome = uVec[2] # mean non home DFE
  sigHomeInitial = sqrt(Sigma[1,1]) # SD of Home DFE
  sigNonHomeInitial = sqrt(Sigma[2,2]) # SD of Nonhome DFE
  covar = Sigma[1,2] # covar of JDFE
  UbAdjust = 1 - pnorm(0, uHome, sigHome)
  ################## M A K E     T H E      J D F E
  sValallTheoMuts = mvrnorm(100000, mu =uVec, Sigma = Sigma)
  DFE1 = sValallTheoMuts[,1]
  DFE2 = sValallTheoMuts[,2]
  
  
  allMeanFitTrajHome = array(0, dim = c(numGens, numiterations))
  allMeanFitTrajNonHome = array(0, dim = c(numGens, numiterations))
  
  for( it in 1:numiterations)
  {
    print(it)
    ID_df  = as.data.frame( array( data = c(1,0, 0,popsize), dim = c(1,4)))
    colnames(ID_df) = c('Strain', 'HomeGR', 'NonHomeGR', 'Count')
    
    
    
    allIDArrayList = list()
    
    
    meanFitnessTrajHome = c()
    meanFitnessTrajNonHome = c()
    
    for(gen in 1:numGens)
    {
      #####    S  E  L  E  C   T  I  O  N    
      #---------------------------------------------------------------------------------------_#
      # select next gen offspring
      meanFitHome = (1/popsize)*sum(ID_df$Count * ID_df$HomeGR)
      meanFitNonHome = (1/popsize)*sum(ID_df$Count * ID_df$NonHomeGR)
      
      probVec = ID_df$Count / popsize + (ID_df$Count / popsize)*(ID_df$HomeGR -meanFitHome)
      # all neg probs = 0 prob
      probVec[probVec<0] = 0
      offspring = rmultinom(1, size = popsize, prob = probVec)
      #   
      
      
      # save important parameters
      allIDArrayList[[gen]] = ID_df
      meanFitnessTrajHome = c(meanFitnessTrajHome, meanFitHome)
      meanFitnessTrajNonHome = c(meanFitnessTrajNonHome, meanFitNonHome)
      
      
      
      # update ID df with new count
      ID_df$Count = offspring 
      
      # remove any extinct lineages
      ID_df = ID_df[!ID_df$Count == 0, ]
      
      # re-do counting so dont have issues calling indexes of parents
      ID_df$Strain = 1:length(ID_df$Strain)
      
      
      
      #### M  U  T  A  T  I  O  N  S
      #---------------------------------------------#
      
      # Sample number of muts that will happen and then choose parent IDs
      M = rpois(1, popsize*U)
      
      if(M > 0 )
      {
        parentIDs = sample(ID_df$Strain, M, prob = (ID_df$Count/popsize), replace = TRUE)
        
        mutEffects = array(data = NA, dim = c(M, 2))
        count = 0 
        for(parent in parentIDs)
       {
          count = count + 1
          # Determine parent's JDFE shape
          sdVec_thisparent = newSDVec(meanFitHome = ID_df$HomeGR[parent], meanFitNonHome =  ID_df$NonHomeGR[parent],sigHomeInitial,sigNonHomeInitial, gamma1, gamma2 )
          covar_thisparent = corr*sqrt(sdVec_thisparent[1]^2 * sdVec_thisparent[2]^2)
          UbAdjust = 1 - pnorm(0, uHome, sigHome)
          Sigma_thisparent =  matrix(c(sdVec_thisparent[1]^2,covar_thisparent ,covar_thisparent ,sdVec_thisparent[2]^2), nrow=2)
          
          # sample mut effect 
          mutEffects[count, ] = as.array(mvrnorm(1, uVec, Sigma_thisparent))
          
          # correct for 0 vars, these ones will make an incorrect Sigma, 
          if(sdVec_thisparent[1] < 0) # sd has reached max, make the muthome effect = uVec[1]
          {
            mutEffects[count, 1] = uVec[1]
          }
          if(sdVec_thisparent[2] < 0)
          {
            mutEffects[count, 2] = uVec[2]
          }
          
          
          
          
        }
        
      
        if(M ==1) # to correct for a single sample frommvrnorm being read as a vector
        {
          mutEffects = array(mutEffects, dim = c(1,2))
        }
        
        # make M new rows for ID_df, one for each mut (4 cols, ID#, homeFit, respFit (!!!!!which not doing anything with RN!!!!))
        newMutArray = as.data.frame(array(data = 0, dim = c(M, 4)))
        colnames(newMutArray) = c('Strain', 'HomeGR', 'NonHomeGR', 'Count')
        
        # make new mutant ID's ,, just count up from highest current ID
        newMutArray$Strain = (max(ID_df$Strain)+1):(max(ID_df$Strain) + M )
        
        # update new Mut Fitness values by adding sampled fit Effects
        newMutArray$HomeGR = ID_df$HomeGR[parentIDs] + mutEffects[,1]
        newMutArray$NonHomeGR = ID_df$NonHomeGR[parentIDs] + mutEffects[,2]
        
        # add count of ONE to each new mut, subtract offspring from parent
        newMutArray$Count = rep(1, times = M)
        
        
        ####??? faster way to do this
        for(rem in parentIDs)
        {
          ID_df$Count[rem] = ID_df$Count[rem] -1
        }
        
        # collate old Id_df and new muts
        ID_df = rbind(ID_df, newMutArray)
        
      }
      # start again with next gen :) 
    }
    
    allMeanFitTrajHome[,it] = meanFitnessTrajHome
    allMeanFitTrajNonHome[,it] = meanFitnessTrajNonHome
    
  }
  
  
  time = 1:numGens
  meanFitTrajHome = rowMeans(allMeanFitTrajHome)
  meanFitTrajNonHome = rowMeans(allMeanFitTrajNonHome)
  
  slopeHome = summary(lm(rowMeans(allMeanFitTrajHome)[800:numGens]~time[800:numGens]))$coefficients[2,1]
  slopeNonHome = summary(lm(rowMeans(allMeanFitTrajNonHome)[800:numGens]~time[800:numGens]))$coefficients[2,1]
  
  varHome_eachGen = rowVars(allMeanFitTrajHome)
  varSlopeHome = summary(lm(rowVars(allMeanFitTrajHome)[800:numGens]~time[800:numGens]))$coefficients[2,1]
  
  varNonHome_eachGen = rowVars(allMeanFitTrajNonHome)
  varSlopeNonHome = summary(lm(rowVars(allMeanFitTrajNonHome)[800:numGens]~time[800:numGens]))$coefficients[2,1]
  
  covar_eachGen = c()
  for(covT in 1:numGens)
  {
    covar_eachGen = c(covar_eachGen, cov(allMeanFitTrajHome[covT, ], allMeanFitTrajNonHome[covT, ]))
  }
  
  covarSlope = summary(lm(covar_eachGen[800:numGens]~time[800:numGens]))$coefficients[2,1]
  
  
  
  
  UbAdjust = 1 - pnorm(0, uHome, sigHome)
  homeSvals = sValallTheoMuts[,1]*(sValallTheoMuts[,1]>=0)
  r1 = 2*mean((homeSvals)^2) / UbAdjust 
  r2 = 2*mean((homeSvals)*(sValallTheoMuts[,2])) / UbAdjust 
  D11 = 2*mean((homeSvals)^3)/UbAdjust 
  D12 = 2*mean((homeSvals)^2*(sValallTheoMuts[,2]))/UbAdjust 
  D22 = 2*mean((homeSvals)*(sValallTheoMuts[,2])^2)/UbAdjust 
  
  
  if(wantPlot == TRUE)
  {
    
    fitness = c(meanFitTrajHome, meanFitTrajNonHome)
    min = c(meanFitTrajHome - sqrt(varHome_eachGen), meanFitTrajNonHome - sqrt(varNonHome_eachGen)) # min of ribbon is -1 SD
    max = c(meanFitTrajHome + sqrt(varHome_eachGen), meanFitTrajNonHome + sqrt(varNonHome_eachGen)) # max of ribbon is +1 SD
    time = rep(1:numGens, times = 2)
    Env = as.factor(c(rep(1, times = numGens), rep(2, times = numGens)))
    WbarDf = data.frame(time, fitness, min, max, Env)
    
    # for adding points at certain gens
   # point_Xvals = c(1,500,750,1000)
  #  point_Yvals_1 = c(meanFitTrajHome[point_Xvals])
   # point_Yvals_2 = c( meanFitTrajNonHome[point_Xvals])
   # df_points = data.frame(point_Xvals, point_Yvals_1, point_Yvals_2)
    
  print( ggplot() + geom_line(data = WbarDf, aes(x=as.numeric(time),y= fitness, colour = Env ))  +
             geom_ribbon(data = WbarDf, aes(x= time, ymin=min, ymax=max,fill = Env,group = Env), alpha = 0.3, linetype = 0) + 
             scale_color_manual(values=c('chocolate1','dodgerblue2'))+  geom_hline(yintercept =  0, linetype = 'dashed')+
            # geom_point(data= df_points, aes(x = point_Xvals, y = point_Yvals_1), size = 3, shape = 24, fill  = 'black')+
             #geom_point(data= df_points, aes(x = point_Xvals, y = point_Yvals_2), size = 3, shape = 25, fill  = 'black')+
             theme_bw()+theme(panel.border = element_rect(colour = 'black', fill = NA, size = 2), 
                              panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
                              axis.text =element_text(size = 20,family = 'Helvetica', color = 'black'), axis.title = element_blank(), legend.position = 'none')+ 
             xlab("Time")+ ylab("Fitness"))
  }
  
  
  percent_colRes_overTime = numeric(numGens)
  for(gen in 1:numGens)
  {
    percent_colRes_overTime[gen] = sum(allMeanFitTrajNonHome[gen, ] > 0) / numiterations
  }
  
  output = list(allMeanFitTrajNonHome, allMeanFitTrajHome, meanFitTrajHome, meanFitTrajNonHome, varHome_eachGen, varNonHome_eachGen, covar_eachGen, percent_colRes_overTime)
  
  
  return(output)
  
  
  
}


# output = list(allMeanFitTrajNonHome, allMeanFitTrajHome, meanFitTrajHome, meanFitTrajNonHome, varHome_eachGen, varNonHome_eachGen, covar_eachGen, percent_colRes_overTime)

sigHome = 0.01 #home DFE SD 
uHome = - 0.001 # mean home DFE 

chosenJDFEs = c(116,86,81,51) # same as rank order section
numDrugs = length(chosenJDFEs)
sigNonHomeVec = data_WrightFisherSims_allNu$sdNonHome[chosenJDFEs] # make sure data_WrightFisherSims_allNu and others is properly saved from the SSWM sims fig 3
uNonHomeVec =data_WrightFisherSims_allNu$uNonHome[chosenJDFEs]
corVec = data_WrightFisherSims_allNu$corr[chosenJDFEs]

r1Vec = data_WrightFisherSims_allNu$r1_SSWM[chosenJDFEs]
r2Vec = data_WrightFisherSims_allNu$r2_SSWM[chosenJDFEs]
D11Vec = data_WrightFisherSims_allNu$D11_SSWM[chosenJDFEs]
D12Vec = data_WrightFisherSims_allNu$D12_SSWM[chosenJDFEs]
D22Vec = data_WrightFisherSims_allNu$D22_SSWM[chosenJDFEs]

r1_BGVec = data_WrightFisherSims_allNu$r1_CI
r2_BGVec = data_WrightFisherSims_allNu$r2_CI


homeID = rep(1, times = numDrugs)
nonHomeID = 1:numDrugs
homeMeanVec = rep(uHome, times = numDrugs)
homeVarVec = rep(sigHome, times = numDrugs)

cVec = r2Vec/sqrt(D22Vec)
hypotheticalDrugJDFEData = data.frame( homeID, nonHomeID, homeMeanVec,uNonHomeVec, homeVarVec,sigNonHomeVec, corVec, r1Vec, r2Vec, D11Vec, D22Vec, D12Vec, cVec)
colnames(hypotheticalDrugJDFEData) = c('Home', 'NonHome','meanHome', 'meanNonHome', 'varHome', 'varNonHome', 'cor', 'r1','r2','D11', 'D22', 'D12', 'C_val')

# run evolution simulation for each JDFE, and plot perc col res over time
rankVec = rank(cVec)

jdfeID = c()
percColRes_overTime_allJDFEs = c()
timeVec = c()
rankID = c()
for(JDFEnum in 1:numDrugs)
{
  print('JDFE')
  print(JDFEnum)
  

  sigNonHome = sigNonHomeVec[JDFEnum]
  uNonHome = uNonHomeVec[JDFEnum]
  corr = corVec[JDFEnum]
  covar = corr*sigHome*sigNonHome
  
  uVec <- c(uHome, uNonHome)
  Sigma <- matrix(c(sigHome^2,covar,covar,sigNonHome^2), nrow=2)
  
  output = simEvo_WrightFisher_plusEpistasis(uVec, Sigma, popsize, U, numiterations, numGens, wantPlot, gamma1 , gamma2 ) 
    
  crOverTime = output[[8]]
  
  jdfeID = c(jdfeID, rep(JDFEnum, times = numGens))
  rankID = c(rankID,rep(rankVec[JDFEnum], times = numGens) )
  percColRes_overTime_allJDFEs = c(percColRes_overTime_allJDFEs , crOverTime  )
  timeVec = c(timeVec, 1:numGens)
}

df_epistasis_cr = data.frame( as.factor(rankID) , percColRes_overTime_allJDFEs, timeVec)
colnames(df_epistasis_cr) = c('rank', 'CR', 'Time')

df_epistasis_cr_resort =df_epistasis_cr
df_epistasis_cr_resort$rank = factor(rankID, levels = c(4,3,2,1))

ColResOverTime_wEpistasis_plot = ggplot(df_epistasis_cr_resort, aes(x = Time, y = CR, color =rank, group =rank)) +
  geom_line(size = 2)+ xlim(c(1,numGens))+
  scale_color_manual(values = c("#fe0059","#fea500", "#0059fe", "#00e594"))+
  theme_bw()+theme(panel.border = element_rect(colour = 'black', fill = NA, size = 2), legend.position = 'none',
                   panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
                   axis.text =element_text(size = 20,family = 'Helvetica', color = 'black'), axis.title = element_blank())

ggsave('ColResOverTime_wEpistasis_plot_1500its.pdf', ColResOverTime_wEpistasis_plot, width = 4.79, height = 4.34 )


#### Plot focal JDFEs change with Fitness ####
# Use example JDFEs from rank order section

colors =  c('#fe0059', '#fea500', '#0059fe','#00e594' ) # for line color
colorPalletes = list()
colorPalletes[[1]] = c('white', '#ffc3d8', '#ff88b2', '#ff3a7f')
colorPalletes[[2]] = c('white','#ffe3b0', '#ffcf75','#ffba3a')
colorPalletes[[3]] = c('white', '#b0cbff', '#75a5ff', '#2672ff')
colorPalletes[[4]] = c('white', '#beffe8', '#83ffd3', '#34ffb7')


timeVec = c(1,500,750,1000) # takes about 300 gens to build diversity 
meanFitHome_fromSims_forabovetimes = c(0,0.008416299, 0.01462998, 0.01731425)
meanFitNonHome_fromSims_forabovetimes = c(0, -0.004725939, -0.007722825, -0.009679302)
hypDrugJDFEplts = list()
tot = 0 
for(d in c(4,3,2,1))
{
  sigNonHomeInitial =  sigNonHomeVec[d]
  sigHomeInitial = sigHome
  timeCount = 0
  for(t in timeVec)
  {
    tot = tot + 1
    UbAdjust = 1 - pnorm(0, uHome, sigHome)
    timeCount = timeCount +1

    corr = corVec[d]
    respMean = uNonHomeVec[d]
    color = colors[d]
    pallete = colorPalletes[[d]]
    homeMean = uHome
    
    sdThisFitness = newSDVec(meanFitHome = meanFitHome_fromSims_forabovetimes[timeCount], meanFitNonHome_fromSims_forabovetimes[timeCount],sigHomeInitial,sigNonHomeInitial, gamma1, gamma2 )
    homeSd = sdThisFitness[1]
    respSd = sdThisFitness[2]
    xplotmin = -0.025
    xplotmax =  0.02
    yplotmin = -0.018
    yplotmax = 0.032
    
    
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
    hypDrugJDFEplts[[tot]] = ggplot(df, aes(x=vals1, y=vals2, z=prob)) + 
      geom_contour_filled(bins = 4)+
      scale_fill_manual(values = pallete) + 
      theme_bw()+
      theme(panel.border = element_rect(colour = 'black', fill = NA, size = 2), 
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
            axis.text =element_text() , axis.title = element_blank(), legend.position = 'none', 
            axis.ticks = element_line(size = 1.3), axis.ticks.length = unit(0.25,'cm') )+ 
      geom_point(aes(weighted.mean(vals1, prob), weighted.mean(vals2, prob)), pch = 4, col = 'black', size = 3) + 
      scale_x_continuous( limits = c(xplotmin,xplotmax), breaks = c(-0.02,0,0.02))+
      scale_y_continuous(limits =c(yplotmin, yplotmax), breaks = c(-0.01,0,0.01,0.02,0.03) ) +
      geom_hline(yintercept =  0, linetype = 'dashed') +   geom_vline(xintercept =  0, linetype = 'dashed')
    
    
    
    
    
    
    
    ggsave(paste0('JDFEnum_', chosenJDFEs[d], '_timeNum_', t, '.pdf'),  hypDrugJDFEplts[[tot]], width =3.89, height = 3.67)
  }
}

allJDFEoverFitPlot = plot_grid(plotlist = hypDrugJDFEplts, ncol = 4, nrow = 4, align = 'v') 
ggsave('allJDFEoverFitPlot.pdf', allJDFEoverFitPlot, width = 7.88, height = 6.69)


