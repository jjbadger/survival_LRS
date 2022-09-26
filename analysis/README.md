## Analysis

This script fits a Cormack Jolly Seber open population model to grey seal capture histories. 






# Methods 




## Defining lifetime and age class-specific provisioning performance 

Provisioning performance, *PP*, was estimated as the individual's effect (intercept) on weaning mass during the course of our study, relative to other females, after accounting for confounding effects of age at first reproduction, age, parity, AMO, and offspring sex. We modeled the weaning mass of pup $j$ born to of female $i$ in year $t$ ($mass_{j,t}$) as a linear mixed-effects model with random individual and year intercepts: 

$mass_{j,t} = \pi_1 \cdot age_{i,t} + \pi_2 \cdot age^2_{i,t} + \pi_{3, m} + \pi_4 \cdot \textbf{I}(sex_{i,t} = female) +  \pi_5 \cdot AFR_i + \pi_6 \cdot AMO + PP_i + \eta_t + \upsilon_{j,t}$

Where parameters ${\pi}$ = $\{\pi_1, \pi_2, \pi_3, \pi_4, \pi_5, \pi_6\}$ reflect the quadratic age effect, effect of parity, pup sex, age at first reproduction (primiparity), and AMO, respectively, where **I** signifies an indicator variable, and $m$ denotes the parity group (1, 2, 3, or 4+) of female $i$ in year $t$ so $m \in \{1, 2, 3, 4\}$, and parameters $\pi_{3, m}$ sum to zero. $PP_i$ is the random effect of individual such that $PP_i \sim N(0,\sigma^2_{PP})$, $\eta_t$ reflects the random year effect, where $\eta_{t} \sim N(0,\sigma^2_{\eta})$, and $\upsilon_{j,t}$ is the error term where $\upsilon_{j,t} \sim N(0,\sigma^2_{\upsilon})$. Fitted estimates of the random effect of each female (i.e. individual intercepts $PP_i$) provide an estimate of that individual $i$'s performance relative to other females in the study after accounting for fixed effects. For example, a positive fitted intercept $PP_i$ indicates that individual $i$  on average produces larger pups at weaning than the population mean after accounting for all effects of age, sex, parity, etc. We fit this model to the weaning masses of pups produced over the course of the female's lives and also separately to the pups produced in the three discretized periods of the female's lives, yielding performance metrics $PP_{life}$, $PP_{early}$, $PP_{peak}$, and $PP_{late}$, respectively.  We standardized these metrics for easier interpretation of results, where the standardized metric for each of the individuals $i = {1,...,n}$ is defined as $PP_i = \frac{PP_i - \mu}{\sigma}$ where $\mu = \frac{1}{n} \sum _1^n PP_i$ and $\sigma = \sqrt{\frac{\Sigma(PP_i - \mu)^2}{n}}$.



## Defining lifetime and age-class specific reproductive frequency *RF*

We estimated a female's reproductive frequency, *RF*, by modeling her reproductive history as a Markov chain in a multi-state open robust design capture recapture modeling framework (MSORD) to account for imperfect detection in reproductive rate. Between her first and last sightings on the island during our study, a female transitions among three reproductive states: initially a first time breeder $F$, then switching between a breeder state $B$, or non-breeder state $NB$. Reproductive frequency is then defined as the probability of transition $\psi^{kB}$ into a breeding state $B$ from any reproductive state $k$. An individual's state transition from year $t$ to $t+1$ is modeled as a categorical trial with probabilities of transition  $\psi^{ks}$ from any state $k$ to any state $s$, $k,s \in \{F, B, NB\}$. We used mixed-effects logistic regression embedded in the MSORD to account for standardized female age, parity, previous breeding state, age at first reproduction, AMO, and random individual and year effects in probability of breeding $\psi^{kB}$ to estimate individual intercepts:


$\psi^{kB} = \lambda_1 \cdot_{i,t} + \lambda_2 \cdot age^2_{i,t} + \lambda_{3, m} +  \lambda_{4, k}  + \lambda_5 \cdot AFR_i + \lambda_6 \cdot AMO_t + RF_i + \theta_{t} + \omega_{i,t}$

Where parameters $\lambda = \{\lambda_1, \lambda_2, \lambda_3, \lambda_4, \lambda_5, \lambda_6\}$ represent the quadratic age effect, effect of parity, the effect of the previous state $k$, age at first reproduction (primiparity), and AMO, respectively, where $m$ denotes the parity group (1, 2, 3, or 4+) of female $i$ in year $t$ so $m \in \{1, 2, 3, 4\}$, parameters $\lambda_{4, k}$ sum to zero. $RF_i$ is the random effect of individual such that $RF_i \sim N(0,\sigma^2_{RF})$, $\theta_{t}$ reflects the random year effect, where $\theta_{t} \sim N(0,\sigma^2_{\theta})$, and $\omega_{i,t}$ is the error term where $\omega_{i,t} \sim N(0,\sigma^2_{\omega})$. As for provisioning performance, the fitted estimate $RF_i$ provides an estimate of individual $i$'s reproductive frequency relative to other females in the study. So, if a female $i$'s fitted intercept $RF_i$ is negative, she reproduces less frequently than the population average after accounting for all fixed effects. We fit this model with capture histories over the course of the females' lives and for the three age-classes separately, yielding performance metrics $RF_{life}$, $RF_{early}$, $RF_{peak}$, and $RF_{late}$, respectively.  We standardized these metrics for easier interpretation of results similarly to above. 



## Survival Model 


 We determined whether $PP$ and $RF$ influenced an individual's survival probability using a state-space formulation of the Cormack-Jolly-Seber open-population model fit to female sighting histories. Our capture history matrix $ch$ has dimension $I x T$, where $I$ is the total number of marked and recruited individuals and $T$ is the number of years in this study. We condition on first capture (defined for individuals in vector $f$), so entrance into the breeding population and age at first reproduction is not estimated here. We define the latent state variable $z$ which takes the value 1 if individual $i$ is alive at time $t$ (*t* = 1,...,$T$) and value 0 if the individual is dead, thus $z$ defines the true state of individual $i$ at time $t$. The individual is necessarily alive at its first sighting ($f_i$), but states of subsequent occasions are modeled as Bernoulli trials. Conditional on being alive at occasion $t$, individual $i$ may survive to occasion $t+1$ with probability $\phi_{i,t}$ ($t$ = 1,..., $T$-1). The Bernoulli success parameter is the product of the survival probability $\phi_{i,t}$ and the state of the individual $i$ at time *t*, ensuring dead individuals ($z_{i,t}$ = 0) remain dead. 
