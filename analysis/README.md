## Analysis

This script fits a Cormack Jolly Seber open population model to grey seal capture histories. 






Methods 




Defining lifetime and age class-specific provisioning performance 

During lactation, grey seal pups consume only milk provided by the female, and as capital breeders, females fast for the entire lactation period. Therefore, in our study, the body mass of a pup at weaning is a reasonable estimate of the energy (i.e. nutrients) transferred to young. While provisioning performance or maternal investment may be better measured by quantifying the change in offspring mass and/or the change in maternal mass from birth to weaning, such an endeavor is logistically limited for a longitudinal study of this scale. By using pup weaning mass as a proxy for maternal investment, we are assuming that the bulk of the variation in pup weaning masses stems from nursing and maternal behavior, and that the variation from other sources is comparatively negligible. However, as maternal performance is determined largely by the resulting offspring size, we determined this assumption reasonable for the question at hand. \textit{PP} was estimated as the individual's effect (intercept) on weaning mass during the course of our study, relative to other females, after accounting for confounding effects of age at first reproduction, age, parity, AMO, and offspring sex. We modeled the weaning mass of pup $j$ born to of female $i$ in year $t$ ($mass_{j,t}$) as a linear mixed-effects model with random individual and year intercepts: 

$mass_{j,t} = \pi_1 \cdot age_{i,t} + \pi_2 \cdot age^2_{i,t} + \pi_{3, m} + \pi_4 \cdot \textbf{I}(sex_{i,t} = female) +  \pi_5 \cdot AFR_i + \pi_6 \cdot AMO + PP_i + \eta_t + \upsilon_{j,t}$

Where parameters **${\pi}$** = $\{\pi_1, \pi_2, \pi_3, \pi_4, \pi_5, \pi_6\}$ reflect the quadratic age effect, effect of parity, pup sex, age at first reproduction (primiparity), and AMO, respectively, where \textbf{I} signifies an indicator variable, and $m$ denotes the parity group (1, 2, 3, or 4+) of female $i$ in year $t$ so $m \in \{1, 2, 3, 4\}$, and parameters $\pi_{3, m}$ sum to zero. $PP_i$ is the random effect of individual such that $PP_i \sim N(0,\sigma^2_{PP})$, $\eta_t$ reflects the random year effect, where $\eta_{t} \sim N(0,\sigma^2_{\eta})$, and $\upsilon_{j,t}$ is the error term where $\upsilon_{j,t} \sim N(0,\sigma^2_{\upsilon})$. Fitted estimates of the random effect of each female (i.e. individual intercepts $PP_i$) provide an estimate of that individual $i$'s performance relative to other females in the study after accounting for fixed effects. For example, a positive fitted intercept $PP_i$ indicates that individual $i$  on average produces larger pups at weaning than the population mean after accounting for all effects of age, sex, parity, etc. We fit this model to the weaning masses of pups produced over the course of the female's lives and also separately to the pups produced in the three discretized periods of the female's lives, yielding performance metrics $PP_{life}$, $PP_{early}$, $PP_{peak}$, and $PP_{late}$, respectively.  We standardized these metrics for easier interpretation of results, where the standardized metric for each of the individuals $i = {1,...,n}$ is defined as $PP_i = \frac{PP_i - \mu}{\sigma}$ where $\mu = \frac{1}{n} \sum _1^n PP_i$ and $\sigma = \sqrt{\frac{\Sigma(PP_i - \mu)^2}{n}}$.
