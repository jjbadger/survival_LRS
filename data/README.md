# `data`
This directory houses all data used in the analysis portion of this repository. 

## jags_dataPP.RData
This is a list of 16 data objects needed to run the model describing variation in survival due to individual variation in lifetime provisioning performance (PP). 


 1) A list of data for JAGS survival model.

         nind = number of individuals, 1x1.
         
         nyears = time series length,  1x1.
         
         ageb = age, discretized into age groups 1-3, for each individual and time step, 273 x 32.
         
         afr = age at first reproduction for each individual, 273 x 1
         
         ch = capture histories for each individual over time series. Coded as 1 for sighted and 0 for not sighted, 273 x 32. 
         
          
 2) List of data for embedded regression model on pup weaning mass. 
 
        Variables with a "v" after them are in the form of vectors for the regression 
        on pup weaning mass because we did not have data for every individual every year 
        (due to skip breeding or failed to collect mass) and loop through the vectors instead of each individual and year. So,  
        
              n.obs = number of observations, 1x1.
              
              pwm_v = pup weaning masses, 2368 x 1, i.e. n.obs x 1.
              
              age_v = standardized age, 2368 x 1.
              
              age_v2 = standardized age^2, 2368 x 1.
              
              par_v = parity, discretized into 1, 2, 3, and 4+ parities, 2368 x 1.
              
              sex_v =  pup sex, coded 0 for male and 1 for female, 2368 x 1.
              
              afr_v = age at first reproduction, 2368 x 1.
              
              amo_v = Atlantic Multidecadal Oscillation, averaged over preceeding 3 years, 2368 x 1.
              
              ind_v = individual identifier, 2368 x 1.
              
              year_v = year, 2368 x 1.


You can also find this description commented-out in the analysis script. 

## pwm.txt 
This is a large dataframe of 
