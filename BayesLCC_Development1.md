---
title: "Developing Bayes LCC"
author: "Paul A H Medley"
date: "2 July 2019"
output: 
  html_document:
    toc: true
    code_folding: hide
    toc_depth: 2.0
    keep_md: true
---



# Introduction

The probability that a randomly sampled fish is a particular length depends upon its age, and other assorted parameters for growth, selectivity and so on. The problem is that for most length frequencies, the ages of any of the fish will be unknown. Age can be effectively removed by integrating over age to obtain the probability a randomly sampled fish will be in length bin $i$:

$$
p_i = \int_{L_i}^{L_{i+1}} \int_0^\infty \  P(L | a, \theta) {\ \ \ \ } P(a | \theta) \ da \ dL
$$

For the above equation, we can neglect the prior for age. Although it may be possible to construct a meaningful prior based on an exponential pdf, it will still be difficult to do this without applying some sort of empirical Bayes procedure. Therefore, it will only be considered if found necessary to obtain sensible results. 

If we assume $P(a | \theta)$ is uniform, the prior will be uninformative. If no limit is placed on age, the prior is improper. However, any arbitrarily large age which is above any possible maximum age (e.g. 10000 years) produces a proper prior with the same effect. Therefore, we note that this makes $p_i$ proportional rather than equal to the integral. This leaves us with the length-given-age function to define.  

# Lognormal Growth Error Model  

Most information on length is derived from age through a growth model. For this we can assume the mean length at age follows a von Bertalanffy (vB) model. The error model for the growth will need to limit sizes to be greater than zero (sizes cannot be zero or negative) and allow variance to increase with size. One option for this which meets these criteria is a log-normal error model:  

$$   
\int_{L_i}^{L_{i+1}} P(L\ |\ a,\, \theta) \ dL 
= \int_{L_i}^{L_{i+1}} Lognormal(L,\, log(L_{\infty}(1-{\large e}^{-K(a-a_0)})),\, \sigma_g) \ dL 
$$

If length is lognormal, length-based selectivity is logistic and the mortality based on age is linear, the likelihood for an observed length taken at random from the catch would be:  

$$ 
P(L | a, \theta) \ \  P(a | \theta) \propto  {\large exp}
  \left ( - \left (
  {log(L)-log(L_{\infty}(1-{\large e}^{-Ka})) \over \sqrt 2 \sigma_g} 
  \right )^2-Za \right )
\  {1 \over {\sigma_g (1+e^{- l_s (L-l_{50})} ) }} 
$$

This needs to be integrated over the length bin and over all ages. Integrating over age removes the age variable and therefore accounts for age-length conversion. However, there is no analytical solution for the age integral, so some numerical approximation is required.

Note that because we only have lengths, we do not need $a_0$ because the ages are relatve ages and their start age is arbitrary. The integral is over age, so the logistic selectivity is also removed for the moment, but needs to be included in the integral over age.  

## Gaussian Quadrature

Numerical integration techniques are all based on the assumption that there is an underlying polynomial that provides a good approximation to the function over some range. Gaussian quadrature provides a very rapid and accurate approximation to integrals for specific functions. It may be particularly appropriate for Stan MCMC because the integral is a weight sum of function evaluations, so the function is fast and smooth (i.e. suitable for auto-differentiation). A common approach used for the normal distribution integral is Gauss-Hermite quadrature.  

The aim is to carry out Gauss–Hermite quadrature for $e^{-x^2}$. To do this variables are changed so that the normal distribution is simplified and defined in terms of $x$ integrating between $+\infty$ and $-\infty$.  

$$
x   =  {log(L)-log(L_{\infty}(1-{\large e}^{-Ka})) \over \sqrt {2} \sigma_g} \\
a   = -log(1-(L / L_{\infty}) \ {\large e}^{-\sqrt {2} \sigma_g x})/K   \\
da  =
{{\sqrt 2 \sigma_g}  \over K } 
{{(L / L_{\infty}) \ {\large e}^{-\sqrt {2} \sigma_g x}} \over 
  {1- (L / L_{\infty}) \ {\large e}^{-\sqrt {2} \sigma_g x }}}      \ dx
$$

The variable $a$ can then be substituted:  

$$
\int_{0}^{\infty}  P(a | L, \theta) \ da \  \propto   
\int_{-\infty}^{\infty} e^{-x^2 \  +(Z/K) \ log(1-(L / L_{\infty}) \ {\large e}^{-\sqrt {2} \sigma_g \ x})}
{{\sqrt 2 }  \over K } 
{{(L / L_{\infty}) \ {\large e}^{-\sqrt {2} \sigma_g \ x}} \over 
{1 - (L / L_{\infty}) \ {\large e}^{-\sqrt {2} \sigma_g \ x}}}      \ dx
$$

and further simplified:

$$
\int_{0}^{\infty}  P(a | L, \theta) \ da \  \propto   
\int_{-\infty}^{\infty} 
e^{-x^2} 
 \ (1-(L / L_{\infty}) \ {\large e}^{-\sqrt {2} \sigma_g x})^{(Z/K - 1)}  \ 
{{ \sqrt 2 }  \over {K } } 
{{(L / L_{\infty}) \ {\large e}^{-\sqrt {2} \sigma_g x}} }  \ dx
$$

In this form, it should be possible to obtain estimates of the integral using Gaussian quadrature. This can be tested using the R function "integrate()" which uses an adaptive method that should be accurate.  


```r
gr <- function(a, pL, K, sg, Z) {
  #Calculates the log-normal probability for L/Linf
  mL <- log((1-exp(-K*a)))
  mL[is.na(mL)] <- 0
  x <- (log(pL) - mL)/(sqrt(2)*sg)
  r <- exp(-(x^2) - a*Z)/sg
  return(r)
}

grint <- function(x, pL, K, sg, Z) {
  #Function for lognormal pdf with change of variable to x
  sx <- sqrt(2)*x*sg
  const <- pL*sqrt(2)/K
  r <- exp(-x^2-sx) * (1-pL*exp(-sx))^(Z/K-1)
  r[is.na(r)] <- 0
  r <- r * const 
  return(r)
}

L <- 5
Linf <- 10
K <- 0.2
sg <- 0.1
Z <- 0.3

n <- 15
le <- double(n)
for (i in 1:n) {
  le[i] <- integrate(gr, 0, Inf, i/Linf, K, sg, Z)$value
}
plot(y=le, x=1:n, type="l")
```

![](BayesLCC_Development1_files/figure-html/CheckIntegration-1.png)<!-- -->

```r
x <- 1
a <- -log(1-(L/Linf)*exp(-sqrt(2)*x*sg))/K

"Should be the same:"
```

```
## [1] "Should be the same:"
```

```r
integrate(gr, 0, Inf, L/Linf, K, sg, Z)
```

```
## 4.413506 with absolute error < 0.00037
```

```r
integrate(grint, -Inf, Inf, L/Linf, K, sg, Z)
```

```
## 4.413506 with absolute error < 1.2e-05
```

Note that simulated size frequency is already dome-shaped without a selectivity function because of the growth meaning fish spend less time in smaller size categories.  

A random draw of parameters allows a test of the integration procedures. In this case, we use the R "integrate()" function and a Gauss-Hermite integral. The latter uses a function with the $exp(-x^2)$ removed as this is approximated by the quadrature polynomial. This is done because the numerical behaviour of the function may make the Gaussian quadrature inaccurate over some ranges of parameters.  


```r
grintGQ <- function(x, pL, K, sg, Z) {
  #Function for lognormal pdf grint with exp(-x^2) removed for Gauss quadrature
  sx <- sqrt(2)*x*sg
  const <- pL*sqrt(2)/K
  r <- exp(-sx) * (1-pL*exp(-sx))^(Z/K-1)
  r[is.na(r)] <- 0
  r <- r * const 
  return(r)
}

nv <- 20
dg1 <- gauss.quad(nv, kind="hermite")
df <- tibble(L=double(), Linf=double(), pL=double(), K=double(), sg=double(), Z=double(), 
                 int_a = double(), int_b = double(), int_d=double())

for (i in 1:2000) {
  repeat {
    Linf <- 5.0+95*runif(1)
    L <- 1.0+Linf*runif(1)
    K <- 0.05+1.0*runif(1)
    sg <- 0.01+0.3*runif(1)
    Z <- 0.05 + 1.0*runif(1)
    int <- grintGQ(dg1$`nodes`, L/Linf, K, sg, Z)
    int_a <- integrate(gr, 0, Inf, pL=L/Linf, K=K, sg=sg, Z=Z)$value
    int_b <- sum(int*dg1$weights)
    if (int_a>0.000001) break
  }
  df <- rbind(df, tibble(L, Linf, pL=L/Linf, K, sg, Z, 
                         int_a, int_b, int_d=(int_a - int_b)))
}

summary(df)
```

```
##        L               Linf              pL                K          
##  Min.   : 1.025   Min.   : 5.055   Min.   :0.01604   Min.   :0.05131  
##  1st Qu.: 9.409   1st Qu.:29.648   1st Qu.:0.27570   1st Qu.:0.32532  
##  Median :21.393   Median :51.729   Median :0.53296   Median :0.57142  
##  Mean   :27.171   Mean   :52.352   Mean   :0.53369   Mean   :0.56294  
##  3rd Qu.:40.460   3rd Qu.:76.166   3rd Qu.:0.78647   3rd Qu.:0.80938  
##  Max.   :98.060   Max.   :99.866   Max.   :1.14401   Max.   :1.04980  
##        sg                Z              int_a              int_b         
##  Min.   :0.01022   Min.   :0.0502   Min.   :  0.0000   Min.   :  0.0000  
##  1st Qu.:0.09206   1st Qu.:0.2879   1st Qu.:  0.8105   1st Qu.:  0.8095  
##  Median :0.16029   Median :0.5493   Median :  1.5537   Median :  1.5468  
##  Mean   :0.16124   Mean   :0.5447   Mean   :  3.4675   Mean   :  2.9468  
##  3rd Qu.:0.23188   3rd Qu.:0.7901   3rd Qu.:  2.9612   3rd Qu.:  2.9032  
##  Max.   :0.31000   Max.   :1.0491   Max.   :193.7660   Max.   :198.3108  
##      int_d          
##  Min.   :-63.74442  
##  1st Qu.:  0.00000  
##  Median :  0.00000  
##  Mean   :  0.52076  
##  3rd Qu.:  0.00007  
##  Max.   : 86.30420
```

```r
sum(df$int_d^2)
```

```
## [1] 41357.65
```

int_a and int_b should be equal, so int_d should be zero. They are generally close, but there are clearly significant differences in some cases.  


```r
plot(x=df$int_a, y=df$int_b)
```

![](BayesLCC_Development1_files/figure-html/Plots1-1.png)<!-- -->

Significant differences occur when the integrals are large. While many calculations give reasonable results, others are unacceptible. Plotting the difference (int_d) against various parameters estimates indicates the problem.   


```r
plot(x=df$int_d, y=df$pL)
```

![](BayesLCC_Development1_files/figure-html/Plots2-1.png)<!-- -->

```r
plot(x=df$int_d, y=df$Z/df$K)
```

![](BayesLCC_Development1_files/figure-html/Plots2-2.png)<!-- -->

```r
plot(x=df$int_d, y=df$K)
```

![](BayesLCC_Development1_files/figure-html/Plots2-3.png)<!-- -->

```r
plot(x=df$int_d, y=df$Z)
```

![](BayesLCC_Development1_files/figure-html/Plots2-4.png)<!-- -->

```r
plot(x=df$int_d, y=df$sg)
```

![](BayesLCC_Development1_files/figure-html/Plots2-5.png)<!-- -->

```r
plot(x=df$int_d, y=(df$Z*df$sg*df$K))
```

![](BayesLCC_Development1_files/figure-html/Plots2-6.png)<!-- -->

The outliers occur when size is close to the asymptote and Z/K is very small.  

## Summary of Findings

Lognormal is a reasonable error model for growth, but does not lead to a simple closed form for the integral over age. 

Gaussian quadrature is very efficient and seems to work well, but only robustly where Z>K. Errors for the more general function are potentially very high even though the majority of function evaluations are reasonably accurate even in these cases.

Problems occur when $Z << K$ and $L / L_{\infty} > 0.60$. Within this region the term $(1-L / L_{\infty} e^{-\sqrt 2 \sigma_g x})^{-1}$ becomes more significant and causes the model to depart from the assumed polynomial used to estimate the integral. Unfortunately estimation of $Z$ is the central objective of this exercise, so it is necessary that any function is robust to the full range of these values. While some improvement may be possible using generalized Gauss-Hermite quadrature, it is not clear that a particularly efficient procedure will be found, so an alternative error model was considered.

# Gamma Growth Error Model

Instead of the lognormal, size at age variation may be modelled using the gamma distribution. This function has similar attributes to the log-normal, with increasing size variance with mean size and probability mass only on the positive range. The functional form is also simpler, so there is an opportunity to derive a simpler form for the integration. There is no reason to prefer the lognormal over the gamma distribution.  

The gamma distribution is given as:

$$
 f(L; \alpha, \beta) = {{\beta^\alpha L^{\alpha-1} {\large e}^{-\beta L}} \over \Gamma(\alpha)}
$$

The mean for the gamma distribution is $\alpha / \beta$, so we can incorporate growth by defining $\beta$ based on mean length defined by age :

$$
  L_{\infty}(1-{\large e}^{-K(a-a_0)} ) = L_{\infty}(1 - b_0 {\large e}^{-Ka} ) \\
  \beta = {\alpha \over L_{\infty}(1 - b_0 e^{-Ka})}
$$

and then substituting $\beta$ in the above equation in the integral across possible ages:

$$
\int_{0}^{\infty}  P(a | L, \theta) \ da \  \propto   
\int_{0}^{\infty} 
{{\left (  {\alpha \over L_{\infty}(1- b_0 e^{-Ka})} \right )^\alpha \ L^{\alpha-1} \ 
{\large e}^{-{\alpha L \over L_{\infty}(1- b_0 e^{-Ka})}-Za}} \over \Gamma(\alpha)}
  \ da
$$

Before proceeding, we provide an R function that can be integrated.  


```r
gam1 <- function(a, alpha, L, Linf, K, Z, b0) {
  AmL <- alpha/(Linf*(1 - b0*exp(-K*a)))
  return((AmL^alpha)*(L^(alpha-1))*exp(-AmL*L-Z*a)/gamma(alpha))
}
alpha <- 10.0
L <- 35
Linf <- 70
K <- 0.2
Z <- 0.1
b0 <- exp(-K*1.0)

integrate(gam1, 0, Inf, alpha, L, Linf, K, Z, b0)
```

```
## 0.1470326 with absolute error < 5.9e-06
```

For this function, the appropriate Gauss quadrature polynomial is the generalized Laguerre-Gauss, which approximates functions containing $e^{-x} x^{-\alpha}; \ \ 0 \leq  x < \infty$ 

To develop such a function, an appropriate substitution is to define a variable $p$ such that:  

$$
\begin{align}
  1+p & = {1 \over 1 - b_0 {\large e}^{-Ka}} \\
  a & = -log \left ( {1-{1 \over (1+p)}} \right ) / K  + log(b_0)/K \\
  da & = {(1+p)^{-2} \over \left ( {1-{1 \over (1+p)}} \right ) K } dp = {1 \over p(1+p)K} dp
\end{align}
$$

Substituting this back into the equation yields:  

$$
\int_{0}^{\infty}  P(a | L, \theta) \ da \  \propto   
\int_{0}^{\infty} 
{{({\alpha (1+p) / L_{\infty} })^\alpha \ L^{\alpha-1} \ 
{\large e}^{-\alpha (1+p) L / L_{\infty} }} \over \Gamma(\alpha) p (1+p) K}
\left ( {{p \over 1+p}} \right )^{Z / K} b_0^{-Z / K}
  \ dp
$$

A little reorganising of this function returns it close to the gamma form:  

$$
\int_{0}^{\infty}  P(a | L, \theta) \ da \  \propto  
{ \alpha L^{\alpha-1} \ b_0^{-Z / K} \over K \ L_{\infty} \ \Gamma(\alpha)}
\int_{0}^{\infty} 
{({\alpha (1+p) / L_{\infty} })^{\alpha-1} \  
{\large e}^{- L \ \alpha (1+p) / L_{\infty} }} {1 \over  \  p }
\left ( {{p \over (1+p)}} \right )^{Z / K}
  \ dp
$$

This function is the most likely to provide a good basis for estimation. A R function to check this produces the same result:    


```r
gam3 <- function(p, alpha, L, Linf, K, Z, b0) {
  AmL <- alpha*(1+p)/Linf
  const <- alpha*L^(alpha-1) * b0^(-Z/K) / (K*Linf*gamma(alpha))
  return(const*(AmL^(alpha-1))*exp(-AmL*L)*(1/p)*((p/(1+p))^(Z/K)))
}
integrate(gam3, 0, Inf, alpha, L, Linf, K, Z, b0)
```

```
## 0.1470338 with absolute error < 1.1e-05
```

We further re-arrange the equation to extract all values that are not required in the integral, and further simplify the equation.  

$$
\int_{0}^{\infty}  P(a | L, \theta) \ da \  \propto  
{ \ b_0^{-Z / K} \over L K \ \Gamma(\alpha)}
\ \left ( {\alpha L \over L_{\infty}} \right ) ^{\alpha}
{\large e}^{- \alpha L / L_{\infty}}
\ \int_{0}^{\infty} 
\ {(1+p)^{\alpha-Z/K-1} \  
\ p^{Z/K-1}
\ {\large e}^{-p \ \alpha L / L_{\infty} }} 
\ dp
$$

Again, we check the function remains consistent:  


```r
gam4 <- function(p, alpha, L, Linf, K, Z, b0) {
  AmL <- alpha*L/Linf
  const <- (AmL^alpha) * exp(-AmL) * b0^(-Z/K) / (K*L*gamma(alpha))
  return(const*((1+p)^(alpha-Z/K-1)) * (p^(Z/K-1)) * exp(-p*AmL))
}
integrate(gam4, 0, Inf, alpha, L, Linf, K, Z, b0)
```

```
## 0.1470338 with absolute error < 1.1e-05
```

A final simplification reduces the integral to a dimensionless values of length as a proportion of $L_\infty$ and $Z/K$.   

$$
\begin{align}
c & = \alpha L / L_\infty \\
m & = Z/K \\
x & = cp     \\
p & = x/c   \\
dx & = 1/c \ dp \\
\int_{0}^{\infty}  P(a | L, \theta) \ da \  & \propto  
{ \alpha \ b_0^{-m} \ {\large e}^{-c} \over L_\infty \ K \ \Gamma(\alpha)}

\ \int_{0}^{\infty} 
\ {(c+x)^{\alpha-m-1} \  
\ x^{m-1}
\ {\large e}^{-x }} 
\ dx
\end{align}
$$

A check for the new function that it produces the same integral value:  


```r
gam5 <- function(x, alpha, L, Linf, K, Z, b0) {
  c <- alpha*L/Linf
  m <- Z/K
  const <-  c * b0^(-m) * exp(-c) / (L*K*gamma(alpha))
  return(const*((c+x)^(alpha-m-1)) * (x^(m-1)) * exp(-x))
}
integrate(gam5, 0, Inf, alpha, L, Linf, K, Z, b0)
```

```
## 0.1470341 with absolute error < 3.1e-05
```

We now run a simulation test with widely varying parameters to see whether the standard integration routine and Laguerre-Gauss integration yield the same results. We firstly check how many nodes are required to minimise the error for the generalised Laguerre-Gauss function. The function itself is now adjusted so only the non-Laguerre-Gauss term is retained.  


```r
V_nv <- seq(10, 100, by=5)
err <- double()
dev <- double()

gam5 <- function(x, alpha, L, Linf, K, Z, b0) {
  c <- alpha*L/Linf
  m <- Z/K
  const <-  c * b0^(-m) * exp(-c) / (L*K*gamma(alpha))
  return(const*((c+x)^(alpha-m-1))) #(x^(m-1)) * exp(-x) has been removed as this is in the LG polynomial
}

for (nv in V_nv) {
  df <- tibble(L=double(), Linf=double(), pL=double(), K=double(), Z=double(), 
               int_a = double(), int_b = double(), int_p=double(), dev=double())

  for (i in 1:20000) {
    repeat {
      Linf <- 5.0+95*runif(1)
      L <- 1.0+Linf*runif(1)
      K <- 0.05+1.0*runif(1)
      Z <- 0.05 + 1.0*runif(1)
      b0 <- exp(-K*runif(1))
      int_aI <- integrate(gam1, 0, Inf, alpha, L, Linf, K=K, Z=Z, b0=b0)
      dg1 <- gauss.quad(nv, kind="laguerre", alpha=Z/K-1)
      int <- gam5(dg1$`nodes`, alpha, L, Linf, K, Z, b0)
      int_b <- sum(int*dg1$weights)
      if (int_aI$value>0.000001) break
    }
    df <- rbind(df, tibble(L, Linf, pL=L/Linf, K, Z, 
                           int_a=int_aI$value, int_b, int_p=(int_aI$value - int_b), 
                           dev = max(0, int_aI$value - int_b - int_aI$abs.error)))
  }
    #summary(df)
  err <- c(err, sum(df$dev^2))
  dev <- c(dev, sum(df$int_p^2))
}
plot(x=V_nv, y=dev)
```

![](BayesLCC_Development1_files/figure-html/GaussQuad5-1.png)<!-- -->

```r
plot(x=V_nv, y=err)
```

![](BayesLCC_Development1_files/figure-html/GaussQuad5-2.png)<!-- -->

The errors decrease with increasing numbers of nodes, with diminishing returns. The differences are minimised around 20 nodes and there appears little improvement above this level. Bearing in mind the adaptive Gauss-Kronrod algorithm used in "integrate" also has an associated error, this appears to provide a reasonable basis for integrating over age. However, in considering the full parameter range errors increase at extremes.    


```r
nv <- 60L
df <- tibble(L=double(), Linf=double(), pL=double(), K=double(), Z=double(), 
             int_a = double(), int_b = double(), int_p=double(), dev=double())

for (i in 1:20000) {
  repeat {
    Linf <- 5.0+95*runif(1)
    L <- 1.0+Linf*runif(1)
    K <- 0.05+1.0*runif(1)
    Z <- 0.05 + 1.0*runif(1)
    int_x <- integrate(gam1, 0, Inf, alpha, L, Linf, K=K, Z=Z, b0=b0)
    dg1 <- gauss.quad(nv, kind="laguerre", alpha=Z/K-1)
    int <- gam5(dg1$`nodes`, alpha, L, Linf, K, Z, b0)
    int_b <- sum(int*dg1$weights)
    if (int_x$value > 0.000001) break
  }
  df <- rbind(df, tibble(L, Linf, pL=L/Linf, K, Z,
                         int_a=int_x$value, int_b, int_p=(int_x$value - int_b), 
                         dev = max(0, int_x$value - int_b - int_x$abs.error)))
}

plot(x=df$int_a, y=df$int_b)
```

![](BayesLCC_Development1_files/figure-html/Plots3-1.png)<!-- -->

```r
plot(x=df$int_p, y=df$pL)
```

![](BayesLCC_Development1_files/figure-html/Plots3-2.png)<!-- -->

```r
plot(x=df$int_p, y=df$K/df$Z)
```

![](BayesLCC_Development1_files/figure-html/Plots3-3.png)<!-- -->

```r
plot(x=df$int_p, y=df$Z/df$K)
```

![](BayesLCC_Development1_files/figure-html/Plots3-4.png)<!-- -->

```r
plot(x=df$int_p, y=df$K)
```

![](BayesLCC_Development1_files/figure-html/Plots3-5.png)<!-- -->

```r
plot(x=df$int_p, y=df$Z)
```

![](BayesLCC_Development1_files/figure-html/Plots3-6.png)<!-- -->

Generally the largest errors occur at higher $Z/K$ values. This must be born in mind when setting the prior for $Z/K$ and in the analysis of very slow growing species. In addition, errors were only significant for sizes less than 30% of $L_\infty$. Generally for these, further improvements in accuracy may be necessary, for example including more nodes at smaller sizes. However, any inaccuracies can be checked after analysis.  

Further attempted changes in variables did not yield improved results. Errors occur at the function margins.   

Unlike the log-normal above, the Gaussian quadrature nodes and weights will need to be re-estimated for each change in $Z/K$. This will mean recalculating the nodes and weights for each parameter draw once, then these can be used for all lengths to be evaluated.  

Stan is introducing a function for 1D integration in version 2.2 (next version). This would be an adaptive integration method. This may prove to be more suitable, but a bespoke integration may provide a faster calculation, since an adaptive method may require a large number of calculations close to function boundaries and it may be difficult to take advantage of using the same node weights over the whole length composition. 

## Fixed Gaussian Quadrature Weight Version

The main problem found with the method above is that the calculation of the Gaussian quadrature is slow. It makes more sense to try to do this once to provide as good an approximation as possible before the MCMC run and accept that there will be some increased uncertainity. Any uncertainty in the numerical integration calculations may be dwarfed by the structural uncertainty and model assumptions.  

To manage this, we introduce a new parameter $d$, which is the $\alpha$ parameter in the original Generalized Gauss–Laguerre quadrature formula. This is set at a fixed value based on the problem, and remains fixed through the simulation.  

$$
\begin{align}
c & = \alpha L_i / L_\infty \\
m & = Z/K \\
x & = cp     \\
p & = x/c   \\
dx & = 1/c \ dp \\
\int_{0}^{\infty}  P(a | L_i, \theta) \ da \  & \propto  
{ \alpha \ b_0^{-m} \ {\large e}^{-c} \over L_\infty K \ \Gamma(\alpha)}

\ \int_{0}^{\infty} 
\ {(c+x)^{\alpha-m-1} \  
\ x^{m-1-d}  \ x^{d}
\ {\large e}^{-x }} 
\ dx
\end{align}
$$

In this form, with a fixed $d$ parameter, the integral does not perform well and can be highly inaccurate across the range of realistic parameters. The $m$ parameter, which the $d$ replaces in the integral, is the important parameter of interest and this is significantly affected by the change in variables. Therefore, this was not pursued further.  

# Length Bin Integration

For integration over lengths in each length bin, "Simpson's Rule" should provide adequate accuracy. This only requires evaluation of the density of fish at each length on either end and in the centre of each bin, so if there are N bins, 2N+1 evaluations are required. Again, because this only requires function evaluations at fixed intervals, the resulting function should be smooth and quick to calculate, making it suitable for Stan MCMC.  

## Simpson's Rule

The rule tabulates the function across the bin from the lower (0) to upper (2) bound. 

$$S_N(f) = \frac{\Delta x}{3} \sum_{i=1}^{N/2} \left( f(x_{2i-2}) + 4 f(x_{2i-1}) + f(x_{2i}) \right)$$
where N is even number of divisions between end points $a$ and $b$, $\Delta x = (b-a)/N $. Splitting each bin into 2 divisions, each 1cm bins would be 0.5cm. The expected density of fish is then estimated at each interval and mid-interval for each bin and the calculation applied to obtain the expected proportion of fish in each bin. The final values will need to be normalised by dividing by the sum of the values across all bins, including bins at either end which might be extended to ensure non-negligible probabilities are included.   

This last issue, ensuring that zeroes are represented so that the normalisation does not eliminate significant likelihoods for zero observations. Given that the growth parameters will need highly informative priors, it makes sense to set zero observations at either end of the range based on these, so they would be defined in the length frequencies in each case based on $L_\infty$.

In this form, the integral calculated for each length bin needs to be normalised to obtain the expected proportion of the sample. Therefore the constants which are fixed over length and age could be neglected.   

$$
S_i = { \alpha \over L_\infty \ K \ \Gamma(\alpha)}
\int_{L_i}^{L_{i+1}} 
{{\large e}^{-(\alpha L / L_\infty)} \over 1 + {\large e}^{-S_s (L - S_{50}) }}
\ \int_{0}^{\infty} \ {((\alpha L / L_\infty)+x)^{(\alpha-Z/K-1)} \  
\ x^{(Z/K-1)} \ {\large e}^{-x }} \  dx
$$

$$
p_i = {S_i \over \sum_{i=1}^{n} S_i }
$$

However, they are useful in providing the scale and therefore avoiding overflow/underflow, so are retained in the calculations, except the parameter $K$ and $b_0$ are removed. The final form of the model is:

$$
p_i \propto { \alpha \over L_\infty \ \Gamma(\alpha)}
\int_{L_i}^{L_{i+1}} 
{{\large e}^{-(\alpha L / L_\infty)} \over 1 + {\large e}^{-S_s (L - S_{50}) }}
\ \int_{0}^{\infty} \ {((\alpha L / L_\infty)+x)^{(\alpha-Z/K-1)} \  
\ x^{(Z/K-1)} \ {\large e}^{-x }} \  dx \ dL
$$

Calculations are carried out on a log-scale where possible also to avoid overflow/undeflow.  

# Likelihood

The observations are counts within each length bin. An appropriate likelihood is the poisson distribution.

$$
P(x) = {\large e^{-\mu} \mu^x \over x!}
$$

# Overdispersion

There are two options for dealing with overdispersion in the observations. Overdispersion in observation error can be accounted for either using an alternative likelihood to account for overdispersion in observation error or to model process errors directly. Observation error might well be higher than expected in the poisson because samples are not independent. Fishing by its very nature will tend to catch fish in groups which may share similar sizes. Process error may occur where the abundance of different ages is not due to a constant mortality. This is mostly likely the result of recruitment variation and, as demonstrated in the simulations above, this may have a profound effect on the estimates. It will not be possible to estimate the processes (e.g. recruitment variation) without a time series of data. Instead the process error might be accounted for as a random effect, which might improve the estimates of model uncertainty rather than dismiss errors that can be neglected.

$$
P(x) = {\beta^\alpha \over \Gamma(\alpha) }

x^{\alpha-1} \large e^{-\beta x}

$$


$$
P(x) = \int_0^\infty {\beta^\alpha \over \Gamma(\alpha) }

\mu^{\alpha-1} \large e^{-\beta \mu} \ \ 

\small { \large e^{-\mu} \ \small \mu^x \over \small x!} d\mu \ = \ 

\int_0^\infty {\beta^\alpha \over \Gamma(\alpha) \ x! }

\mu^{\alpha+x-1} \ \large e^{-\beta \mu -\mu} \ \ 

 d\mu \ = \ { \Gamma(\alpha+x) \ \beta^\alpha \over \Gamma(\alpha) \ \Gamma(x) \ (\beta+1)^{\alpha+x}}

$$

This is a negative binomial distribution. Stan provides a negative binomial probability mass function with an additional term to represent over-dispersion relative to the binomial distribution. This was used to fit the model to see whether the performance improved using this likelihood.   

# Spawning Potential Ratio

One option is to calculate the SPR based on the proportions with alternative mortality rates:   

$$
W_i(Z) \propto 
\int_{L_i}^{L_{i+1}} 
{L^b \ {\large e}^{-(\alpha L / L_\infty)} \over 1 + {\large e}^{-M_s (L - M_{50}) }}
\ \int_{0}^{\infty} \ {((\alpha L / L_\infty)+x)^{(\alpha-Z/K-1)} \  
\ x^{(Z/K-1)} \ {\large e}^{-x }} \  dx
$$

where $M_s$ and$M_{50}$ are the maturity ogive parameters and $b$ is the length-weight exponential parmeter (usually $b$ ~ 3).  The SPR is then approximated by:  

$$
  SPR = {\sum_{i=1}^{n-1} W_i(Z) \over  \sum_{i=1}^{n-1} W_i(1.5)}
$$

However, the problem here is because we are integrating over age, the same mortality is being applied across all ages, regardless of selectivity. Selectivity is an important part of the model, as it offers a way to attempt to increase SPR without necessarily reducing the mortality (i.e. fishing effort) applied. Therefore, it makes more sense to return to explicit age to account for explicit age-length effects, specifically selectivity which affects which ages/lengths fishing mortality applies. Therefore, it makes sense to limit this simple calculation from the fully selected to maximum length and avoid using lengths from the period selectivity is less than 1.0.  
