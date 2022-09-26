# `data`
This directory houses all data used in the analysis portion of this repository. 

## jags_dataPP.RData
This is a list of 16 data objects needed to run the model describing variation in survival due to individual variation in lifetime provisioning performance (PP). 


 A list of data for JAGS survival model.

         nind = number of individuals, 1x1.
         
         nyears = time series length,  1x1.
         
         ageb = age, discretized into age groups 1-3, for each individual and time step, 273 x 32.
         
         afr = age at first reproduction for each individual, 273 x 1
         
         ch = capture histories for each individual over time series. Coded as 1 for sighted and 0 for not sighted, 273 x 32. 
         
          
 And, a list of data for embedded regression model on pup weaning mass. 
 
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

## Methods in data collection 

We used a 32-year mark-resighting dataset (1989-2020) of known-aged grey seal females to determine how lifetime reproductive performance influences survival. Individuals were marked as pups shortly after weaning with unique alpha-numeric hot-iron brands in 1985, 1986, 1987, and 1989. These permanent brands allowed reliable identification of individuals over the course of their lives. Females can recruit to the breeding population as early as 4 years old, but this is uncommon, and the average age of first reproduction is 5.6 +- 0.12 SE years for these cohorts (den Heyer et. al 2013). Each breeding season 1989-2020, teams of researchers conducted 5-7 roughly weekly censuses of branded females returning to the island to give birth and mate. Once sighted, branded individuals with dependent pups were visited daily from a roughly 30m distance. Prior to weaning, pups were sexed and marked with semipermanent, uniquely numbered tags in the hind flipper to ensure accurate identification after the marked female ended lactation and returned to sea, leaving her pup in the colony. Females attend their pups continuously throughout lactation. Therefore, once a pup was sighted alone, it was considered weaned and weighed to the nearest 0.5-kg.

The probability of observing a marked female during any given year includes both the probability the female is present, and the probability that she is detected given presence at the breeding colony. Females not seen in any of the weekly brand-resighting censuses in a given year could arise via one of three events: 1) the female died, 2) the female was not present at the breeding colony due to temporary emigration ("reproductive skipping") or they quickly abandoned their pup, or 3) was present in the breeding colony, gave birth and lactated normally, but was not detected during weekly censuses either because of her position on the island or unreadable brand markings. For individuals branded from 1998-2002, 3.7% or 170 of 4569 sightings (from breeding seasons 2002 to 2012) were not readable either because of a poor-quality brand or poor conditions (Bowen et al. 2015). Low sighting probability associated with poor brands would obscure female quality, so we removed individuals whose brands were sighted as unreadable for more than 25\% of their total sightings. Of the remaining brands, a recent CMR analysis of this population indicated that, if a female rears a pup on the island, there is less than a 5% chance researchers will fail to detect her in at least one resighting census (Badger et al. 2020). Grey seals are highly site philopatric, and once recruited to a breeding colony, will very rarely pup elsewhere (Bowen et al. 2015). Thus, we are able to reliably follow the reproductive history of individuals, and do not expect permanent emigration to other colonies to be a significant source of sighting error. We removed sightings of females without pups from this analysis. Individuals that are not rearing pups can be skittish and may flee to the water, resulting in a lower sighting probability than females nursing and defending young. Pregnant females generally arrive at the breeding colony five days before parturition and can travel large distances before settling to give birth (median: 2.4 km, Weitzman et al. 2017). The distance moved by nursing females from one day to the next averages about 5 m (Boness1979, Weitzman et al. 2017), so movement by females during lactation is not a source of sighting error. Nursing typically lasts 2-3 weeks, so females are usually available to be resighted during two or more of the weekly re-sighting censuses in a given year.   

Individual sighting histories were collected from age at first reproduction (first sighting in breeding colony) until death (or 2020 for animals still living). Sighting histories of individuals were scored as a 0 (not sighted) or 1 (sighted) for each year from 1989 to 2020. Females sighted in only one breeding season were omitted from this analysis (n = 27) to ensure that they had in fact recruited to the Sable Island breeding population and we have adequate data to estimate trends in reproductive performance. 


Environmental Data

The Atlantic Multidecadal Oscillation (AMO) index describes multidecadal atmosphere and sea variability in the Atlantic, with the warm phase associated with positive SST anomalies over most of the North Atlantic. Here, we use the AMO to account for large-scale environmental variability that can influence  both survival and individual reproductive investment and success. We used the annual mean of AMO unsmoothed from the Kaplan SST V2 index calculated at NOAA/ESRL/PSD1 (http://www.esrl.noaa.gov/psd/data/timeseries/AMO/downloaded from NOAA on 22 August 2019). The time series was detrended with 10- year low-pass filtered annual mean area-averaged SST anomalies over the North Atlantic basin. For this analysis, we used the mean annual index in the preceding 3 years to investigate environmental effects on survival rates (as in den Heyer et al. 2021). 




