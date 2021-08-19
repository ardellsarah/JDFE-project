# JDFE-project
A model of evolution of collateral resistance and sensitivity. The code is compiled in a way that, once all packages indicated a top of file are downloaded and loaded for use, the entirety of the code can be sourced and it will produce the majority of the figure components and files for the main text and supplement and save them to the specified working directory.  As sourcing the whole file is computationally intensive and as there are several parameters which the user can change to explore results, the code is broken up into sections, some of which can be sourced separately to produce different figures. I detail each section below.


#####################################################################
FIGURE 1:
All code needed to produce figure 1. This section can be run independently of all other sections. NOTE: in order for the section to work, the working directory must contain files 'KOGrDATAPrtoStop.csv' (from Cheverau et al 2015) and 'WTReplicate.csv' (sent to us by authors of Cheverau et al 2015), both of which are also available in this repository. 

Running the entirety of this section will automatically produce the graphs comprising the figure and save them to the working directory. It will also save 5 .csv files: 

'mutTypeDf.csv' - contains the identity (2 - beneficial ; 1 - deleterious; 0 - neutral) we assign to all KO mutants in the Cheverau et al 2015 data set. See methods section on identifying benefical, deleterious and neutral mutants for more

'mutTypeCounts.csv' - contains the information pertinent to the false discovery rate of beneficial mutations in each environment 

'allABRJDFEStats.csv' - contains the important statistics of all the antibiotic resistance JDFEs of interest. Including, the mean, variance and covariance of the home and non-home and the pleiotropy parameters r1, r2,d11,d22,d12 and c calculated using all values reported in the Cheverau et al 2015 data set

'allABRJDFEStats_sigVals.csv' - contains the important statistics of all the antibiotic resistance JDFEs of interest. Including, the mean, variance and covariance of the home and non-home and the pleiotropy parameters r1, r2,d11,d22,d12 and c calculated using only the 'significant' (non-neutral) values we called from the Cheverau et al 2015 data set

'probCollRes_CollSens_df.csv' - contains the expected number of mutations which would meet the criteria for collateral resistance and sensitivity in each antibiotic pair based off the noise distribution alone. 

#######################################################################

SUPP FIGURE 1 :
Running this entire section will produce all components of the Supplementary figure 1. Note, the working directory must contain file "Single_Gene_JDFE_forR.csv", which is available in this repository 


#######################################################################

FIGURE 2 : 5 example JDFEs + Trajectories
All code needed to produce figure 2 and the supplement to figure 2. This section can be run independently of all other sections. 

Subsection: Simple 2D Gaussian (main text)
- running this whole section will create and save all components of figure 2

Subsection: More Complex 'Double Gaussian' for Supplement
- running this whole section will create and save all components of supplementary figure 2


#######################################################################

FIGURE 3: Gaussian JDFEs
All code needed to produce figure 3. This section can be run independently of all other sections. 


Subsection: 'SSWM Section'
This subsection is written to automatically produce all of the graphs for the first column of figure 3. It will automatically run through all JDFEs produced by looping through all parameters for mean in environment 2 (m2Vec), sd of environment 2 (sdVec) and correlation between home and non-home (corVec) for the number of iterations specified by 'numIterations', use 1000 for the same results as the main text figure. 
  
The graphs for figure 2 will all be saved to the to the working directory 
  
  
Subsection: 'Wright-Fisher Section'
This subsection is written to automatically produce all of the graphs (other than the first column) of figure 2. It will automatically run through all JDFEs produced by looping through all parameters for mean in environment 2 (m2Vec), sd of environment 2 (sdVec) and correlation between home and non-home (corVec) and all population size values used for the figure . For each JDFE, it will run the wright-fisher evolution simulation desicribed in the Methods of the text for the number of iterations specified by the 'numIterations' variable , use 300 for the same results as the main text figure. 

The graphs for figure 2 will all be saved a to the working directory and printed at the end of the full section. The results of these simulations are also saved in a .csv file 'data_WrightFisherSims_allNu', which is exported to the working directory as well.
  
  
 


#######################################################################

FIGURE 4: Rank Order
All code needed to produce figure 4 and the supplementary figure to figure 4  This section CANNOT be run independently of all other sections. It relies on results from the figure 3 section, so must be run after figure 3 section is run.  All Figure components will be saved to the working directory.



  
#######################################################################

FIGURE 5: Measuring JDFEs
All code needed to produce figure 5 and the supplementary figure to figure 5.  This section CANNOT be run independently of all other sections. It relies on results from the figure 3 section, so must be run after figure 3 section is run. All Figure components will be saved to the working directory.

In order to generate the supplementary figures for the LD simulations, simply change the parameter sigThresh. This parameter takes the value 0.5 in the main text and 1,2,and 3 in the supplement
  
  
#######################################################################

APPENDIX: JDFEs with Epistasis
All code needed to produce the appendix figure. This section CANNOT be run independently of all other sections. It relies on results from the figure 3 section, so must be run after figure 3 section is run.  All Figure components will be saved to the working directory.


  



