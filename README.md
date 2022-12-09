# Bayesian Length Catch Curve 
This model has been replaced by fishblicc which has better performance and is more flexible.  

This is the Bayesian Length Catch Curve project to estimate mortality on single length frequency samples. This is explicitly a "data-limited" method designed to be used as a risk assessment approach rather than provide definitive information on stock status. It is suitable for end-of-project evaluation of length frequency data that may have been collected over a short time period (e.g. a year) or to monitor status of many species that may take up a small proportion of the catches.  
In common with such methods, there are significant assumptions that will not be met in practice, therefore the focus of this method is robustness and looking for ways to reduce the impact of model structural errors or at least spot when such errors invalidate the results. Although the aim is to make the estimation method as accurate as possible, because the method is risk based, the decision has been made to make the results more precautionary rather than less biased.  
Catch curves and related interpretations (e.g. spawner potential ratio) look at the ratio of large (older, mature) fish in a sample compared to small (immature young) fish. If the ratio is too small compared to what might be expected in a sample of a population that is unfished, the popualtion will be at high risk of overexploitation.  
The estimation method works by assuming von Bertalanffy growth model describes the mean length at age in a gamma probability density function (constant CV) and then intergrating over age and length to estimate the expected proportion of the population in each length bin. The model is flexible enough to include all data, so for example observed lengths above L_inf are not rejected.
The model is fitted to a single length frequency sample with 6 parameters:  
  - S50% and Ssp for the selectivity at 50% and selectivity steepness for a logistic selectivity curve.  
  - L_inf the asymptotic mean length  
  - Z_K the total mortality in units of K, the growth rate.  
  - Gs the growth model error parameter.  
  - phi the over-dispersion parameter for the negative binomial.  

The model is fitted using Stan (mc-stan.org). 
 
