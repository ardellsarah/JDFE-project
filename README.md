# JDFE-project
A model of evolution of collateral resistance and sensitivity. The code is compiled in a way that, once all packages indicated a top of file are downloaded and loaded for use, the entirety of the code can be sourced and it will produce the majority of the figures 2-5 for the main text and supplement and save them to the specified working directory. Figure 1 requires some manual saving. As sourcing the whole file is computationally intensive and as there are several parameters which the user can change to explore results, the code is broken up into sections which can be sourced separately to produce different figures. I detail each section below.


#####################################################################

FIGURE 1:
All code needed to produce figure 1 and the supplement to figure 1. This section can be run independently of all other sections. 

Subsection: Simple 2D Gaussian (main text)
- for generating the each JDFE graph (FIG 1 A-E) and the array containing the DFEs in both environments, simply change the correlation ('corr') and non-home DFE mean ('respMean') to the value specified in caption to figure 1 and run the code through the 'example of evolution function call' line. 
- Then, for generating the corresponding fitness trajectory, run the 3 lines of code beneath the 'example of evolution function call', feeling free to change any variable other than 'sValallTheoMuts', as this is the array of the JDFE data generated above. To produce the same results as Figure 1, use variables specified in the figure caption 
- Must be done separately for each JDFE in figure 1

Subsection: More Complex 'Double Gaussian' for Supplement
- for generating the JDFE graphs (SUPP FIG 1 A-D) and the array containing the DFEs in both environments, change the variables 'resp1mean' and 'resp2mean' to the values specified in caption to supplemental figure 1 and run the code through the 'example of evolution function call' line. 
- As before, for generating the corresponding fitness trajectory, run the 3 lines of code beneath the 'example of evolution function call', feeling free to change any variable other than 'sValallTheoMuts', as this is the array of the JDFE data generated above. To produce the same results as Supplemental Figure 1, use variables specified in the figure caption 
- Must be done separately for each JDFE in Supplemental figure 1
#######################################################################

FIGURE 2:
All code needed to produce figure 2. This section can be run independently of all other sections. 

Subsection: 'Non-SSWM'
This subsection is written to automatically produce all of the graphs (other than the first column) of figure 2. It will automatically run through all JDFEs produced by looping through all parameters for mean in environment 2 (m2Vec), sd of environment 2 (sdVec) and correlation between home and non-home (corVec) and all mutation rate values used for the figure (specified in U loop, line 811). For each JDFE, it will run the wright-fisher evolution simulation desicribed in the Methods of the text for the number of iterations specified by the 'numIterations' variable (line 801), use 100 for the same results as the main text figure. 

The graphs for figure 2 will all be saved as variables with names r2VSlope_Nu_<INSERTNU>, D22vVarSlope_Nu_<INSERTNU>, and D12vCoVarSlope_Nu_<INSERTNU>, depending on the contents of the plot and the current N (population size) * U (mutation rate). They are saved to the working directory and printed at the end of the full section. The results of these simulations are also saved in a .csv file 'allData_nonSSWM_TheoJDFEs', which is exported to the working directory as well.
  
  
Subsection: 'SSWM Theo JDFE'
This subsection is written to automatically produce all of the graphs for the first column of figure 2. It will automatically run through all JDFEs produced by looping through all parameters for mean in environment 2 (m2Vec), sd of environment 2 (sdVec) and correlation between home and non-home (corVec) for the number of iterations specified by 'numIterations', use 300 for the same results as the main text figure. 
  
  As above, the graphs for figure 2 will all be saved to the environment as variables with names r2VSlope_Nu_0.1, D22vVarSlope_Nu_0.1, and D12vCoVarSlope_Nu_0.1, depending on the contents of the plot and the current N (population size) * U (mutation rate). They are saved to the working directory and printed at the end of the full section, as indicated in code comments. 
  
#######################################################################
  
FIGURE 3:
All code needed to produce figure 3 and the supplement to figure 3 (Supp Fig 2). This section can be run independently of all other sections.

Subsection: SSWM With Epistasis Sim
Will produce the 1st graphs in panels B and C in main text Figure 3. To generate each fitness trajectory panel and the corresponding r2 vs deltaY panel, use the U and gamma1 parameters as specified in figure 3 caption. Graphs with be automatically saved to working directory with names speficing the Nu of the simulation they are from and the correlation of the corresponding JDFE. Fitness trajectories for all correlations will be saved to the working directory along with the panel for r2 vs deltaY.

Subsection:  EPISTASIS SIM, NOT SSWM 
Will produce the 2nd - 4th graphs in panels B and C in main text Figure 3. Running this whole section will cycle through all U and corresponding gamma1 values used for the main text figure. Fitness trajectories for all correlations will be saved to the working directory along with the panel for r2 vs deltaY.
  

Subsection: PANEL A: how JDFES change over time
Running this section as is will produce panel A for figure 3. The code automatically generates the JDFEs expected at times 1, 100 and 1000 for a population following the theoretical trajectories along with a graph of the pairwise fitnesses (fit in home, fit in non-home) for all 3 JDFEs. 

#######################################################################

FIGURE 4:
All code needed to produce figure 4 and the supplementary tables 1a and 1b. This section can be run independently of all other sections. NOTE: in order for the section to work, the working directory must contain files 'KOGrDATAPrtoStop.csv' (from Cheverau et al 2015) and 'WTReplicate.csv' (sent to us by authors of Cheverau et al 2015), both of which are also available in this repository. 

Running the entirety of this section will automatically produce the graphs comprising figure 4 and save them to the working directory. It will also save 4 .csv files: 

'mutTypeDf.csv' - contains the identity (2 - beneficial ; 1 - deleterious; 0 - neutral) we assign to all KO mutants in the Cheverau et al 2015 data set. See methods section on identifying benefical, deleterious and neutral mutants for more

'mutTypeCounts.csv' - contains the information pertinent to the false discovery rate of beneficial mutations in each environment 

'allABRJDFEStats.csv' - contains the important statistics of all the antibiotic resistance JDFEs of interest. Including, the mean, variance and covariance of the home and non-home and the pleiotropy parameters r1, r2,d11,d22,d12 and c calculated using all values reported in the Cheverau et al 2015 data set

'allABRJDFEStats_sigVals.csv' - contains the important statistics of all the antibiotic resistance JDFEs of interest. Including, the mean, variance and covariance of the home and non-home and the pleiotropy parameters r1, r2,d11,d22,d12 and c calculated using only the 'significant' (non-neutral) values we called from the Cheverau et al 2015 data set

#######################################################################

FIGURE 5:
All code needed to produce figure 5 and the supplementary figure to figure 5 (supp fig 4). This section CANNOT be run independently of all other sections. It relies on results from the figure 4 section, so must be run after figure 4 section is run. 

Subsection: ALL VALUES
When the entirety of this section is run, it produces all panels of the Supplemental figure 4 (supplement to figure 5) and saves them to the working directory. It also produces several extra graphs which may be of interest to the user (explained in comments in code). These extra graphs are not automatically saved to the working directory but can be by un-commenting their corresponding ggsave lines. 

Subsection: Main text - only sig values 
When the entirety of this section is run, it produces all panels of the main text figure 5  and saves them to the working directory. It also produces several extra graphs which may be of interest to the user (explained in comments in code). These extra graphs are not automatically saved to the working directory but can be by un-commenting their corresponding ggsave lines. 
#######################################################################

SUPP - ANTIBIOTIC RESISTANCE DATA SIMULATION
All code needed to produce the supplmentary figure to figure 4/5 (supp fig 5) as well as supplementary table 1 and supplementary file of all Weibull bin fit parameters. This section can be run independently of all other sections. NOTE: in order for the section to work, the working directory must contain files 'KOGrDATAPrtoStop.csv' (from Cheverau et al 2015) and 'WTReplicate.csv' (sent to us by authors of Cheverau et al 2015), both of which are also available in this repository. 


This section runs the full wright-fisher simulation for the antibiotic resistance JDFES. Parameters that can be changed for the simulation are 'popsize', 'U', 'NumGens', and 'numIterations' (lines 4215-4219). When the entirety of this section is run, all panels for supplemental figure 5 are saved to the working directory. The files DFE_DistFit_Data.csv  and 'allABRWeibullBinParameters.csv' are also saved to the working directory. This former file contains the log liklihood goodness of fit estimates for the home DFEs for the normal, eponential, and weibull distributions, and the later contains the shape and scale parameters for the weibull distributions fit to all home environmnet 'bins'. NA values in this file represent bins for certain home/non-home pairs which had no mutation representation in the Cheverau et al 2015 data set (ie bins for which a distribution could not be fit/bins in which it is not possible for a mutation to occur in the simulation). 







