 //><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>
 // ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>
 //  
 //###########    Bayes LCC Stock Assessment    #############
 //###########    PAUL MEDLEY                   #############
 //###########    paulahmedley@gmail.com        #############
 //###########    December 2019                 #############
 //><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>
 // ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>

//###########################################
//####  FUNCTIONS   #########################
//###########################################


functions {

  vector[] gauss_laguerre_quad(int nt, real alpha) {
    //  Compute a generalized Gauss-Laguerre rule for approximating
    //    Integral ( a <= x < +oo ) |x-a|^ALPHA exp(-B*x(x-a)) f(x) dx
    //  of order "nt".
    //  Order nt is the number of points in the rule:
    //  ALPHA is the exponent of |X|:
    real r8_epsilon = 2.220446049250313E-016;
    int itn  = 30;
    vector[nt] aj;
    vector[nt] bj;
    vector[nt] wt = rep_vector(0, nt);
    real       zemu = tgamma(alpha+1.0);
    vector[nt] knot_wt[2];
    int  m;
    int  mml;
    int  k;
    real b;
    real c;
    real f;
    real g;
    real sign_g;
    real p;
    real r;
    real s;

    {
      vector[nt] i = cumulative_sum(rep_vector(1, nt));
      aj = 2.0*i - 1.0 + alpha;
      bj = sqrt(i .* (i+alpha));
    }
    wt[1] = sqrt(zemu);
    //  Diagonalize the Jacobi matrix.
    if (nt == 1) {
      knot_wt[1] = aj;
      knot_wt[2] = square(wt);
      return knot_wt;
    }
    
    bj[nt] = 0.0;
    
    for (ll in 1:nt) {
      int j = 0;
      while (1) {
        m = ll;
        while (m<nt) { //set m
          //if (m == nt) break;
          if (fabs(bj[m]) <= r8_epsilon * (fabs(aj[m]) + fabs(aj[m+1])) ) break;
          m += 1;
        }
        if (m == ll) break;  //ends while(1) loop

        p = aj[ll];
        if (itn <= j) {
          reject("gauss_laguerre_quad - Fatal error! Iteration limit exceeded")
        }
        j += 1;
        g = (aj[ll+1]-p) / (2.0*bj[ll]);
        if (g>=0) sign_g = 1.0; else sign_g = -1.0; 
        r =  sqrt(square(g) + 1.0);
        s = g + fabs(r)*sign_g;
        g = aj[m] - p + (bj[ll] / s);
        s = 1.0;
        c = 1.0;
        p = 0.0;
        mml = m - ll;
        
        for (ii in 1:mml) {
          int i = m - ii;
          f = s * bj[i];
          b = c * bj[i];
          if (fabs(g) <= fabs(f)) {
            c = g / f;
            r =  sqrt(square(c) + 1.0);
            bj[i+1] = f * r;
            s = 1.0 / r;
            c *= s;
          } else {
            s = f / g;
            r =  sqrt(square(s) + 1.0);
            bj[i+1] = g * r;
            c = 1.0 / r;
            s *= c;
          }
          g = aj[i+1] - p;
          r = (aj[i]-g)*s + 2.0 * c * b;
          p = s * r;
          aj[i+1] = g + p;
          g = c * r - b;
          f = wt[i+1];
          b = wt[i];
          wt[i+1] = s * b + c * f;
          wt[i] = c * b - s * f;
        }
        aj[ll] -= p;
        bj[ll] = g;
        bj[m] = 0.0;
      }  //while
    }
    //  Sorting
    for (ii in 2:m)  {
      int i = ii - 1;
      k = i;
      p = aj[i];
      
      for (j in ii:nt) {
        if (aj[j] < p) {
          k = j;
          p = aj[j];
        }
      }
      
      if (k != i) {
        g = aj[i];
        aj[k] = g;
        aj[i] = p;
        p = wt[i];
        wt[i] = wt[k];
        wt[k] = p;
      }
    }
    knot_wt[1] = aj;
    knot_wt[2] = square(wt);
    return knot_wt;
  }


  vector prob_length(vector L, real Gamma_alpha, real Linf, real Z_K, 
                  real Sels, real Sel50, 
                  vector LG_x, vector LG_wt) {
    int nl = rows(L);
    int ng = rows(LG_x); 
    vector[ng] ss;
    vector[nl] v;
    real       ga_l = Gamma_alpha/Linf;
    vector[nl] lSel = log_inv_logit(Sels * (L - Sel50));
    vector[nl] c = ga_l*L;
    vector[nl] cn =  lSel - c + log(ga_l) - lgamma(Gamma_alpha);
    real       expn1 = Gamma_alpha-Z_K-1.0;
    vector[ng] logwt = log(LG_wt);
    
    //Gauss-Laguerre Quadrature
    for (i in 1:nl) {
      ss = log(c[i] + LG_x)*expn1 + logwt;
      v[i] = log_sum_exp(ss);
    }
    return exp(cn + v); //(x^(m-1)) * exp(-x) has been removed as this is in the LG polynomial
  }  //prob_length


  vector prob_mass(vector x, vector y) {
    //Uses mid point bin values to obtain quadratic approximation to probability mass in bin
    // x are bin bounds (N Bins + 1)
    // y bn heights including mid-points (2N + 1)
    int BN = rows(x)-1;
    vector[BN] sv;
    real norm;
    for (i in 1:BN) {
      sv[i] = (y[2*i-1] + 4*y[2*i] + y[2*i+1])*(x[i+1]-x[i])/3.0;  //Simpson's rule
    }
    norm = sum(sv);
    return sv / norm;
  }

}  //FUNCTIONS


//###########################################
//####  DATA   #########################
//###########################################


data {
  //Likelihood are strictly applied to the bins supplied. 
  //It is important to include all zero observations with non-negligible probability
  int<lower=2>   NK;
  //Assume non-overlapping length bins
  int<lower=1>   oBN;          //Number of length bins

  vector[oBN+1]  Bnd;        //Bin lower bounds plus final upper bound
  int            fq[oBN];      //Frequency in each bin

  //Constant Hyperparameters for priors
  real poLinfm;
  real poLinfs;
  real polZ_Km;
  real polZ_Ks;
  real polGam;
  real polGas;
  real poSel50m;
  real poSel50s;
  real poSelsm;
  real poSelss;
  real polNB_phim; 
  real polNB_phis;
}

transformed data {
  int BN_1 = oBN+1;       //Includes an upper bound
  int NV = sum(fq);
  int BN2_1 = 2*oBN+1;
  vector[BN2_1] Len;      //All length bin points that need to be evaluated (includes mid point for each bin)

  for (i in 1:oBN) {
    Len[2*i-1] = Bnd[i];
    Len[2*i] = (Bnd[i]+Bnd[i+1])*0.5;
  }
  Len[BN2_1] = Bnd[BN_1];
}


//###########################################
//####  PARAMETERS  #########################
//###########################################


parameters {  //modelled param
  real<lower = -polZ_Ks>           nlZ_K;
  real<lower = -poLinfm/poLinfs>   nLinf;  
  real                             nGalpha;
  real<lower = -poSel50m/poSel50s> nSel50;
  real<lower = 0>                  nSels;
  real                             nNB_phi;
}

transformed parameters {
  real Z_K  = exp(polZ_Km + nlZ_K*polZ_Ks);
  real Linf = poLinfm + nLinf*poLinfs;  
  real Galpha = exp(polGam + nGalpha*polGas);
  real Sel50  = poSel50m + nSel50*poSel50s;
  real Sels   = nSels*poSelsm;
  real NB_phi = exp(polNB_phim + nNB_phi*polNB_phis);
}


//###########################################
//####    MODEL     #########################
//###########################################

model {
  vector[oBN] efq; 
  //<><  ><>  <><  ><>  <><  ><>  <><  ><>  <><  ><>  
  //###   PRIORS   ###  <><  ><>  <><  ><>  <><  ><>
  //<><  ><>  <><  ><>  <><  ><>  <><  ><>  <><  ><>
  
  //see https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations
  target += normal_lpdf(nlZ_K | 0, 1);
  target += normal_lpdf(nLinf | 0, 1);           
  target += normal_lpdf(nSel50 | 0, 1);
  target += normal_lpdf(nGalpha | 0, 1);
  target += exponential_lpdf(nSels | 1);
  target += normal_lpdf(nNB_phi | 0, 1);

  //<><  ><>  <><  ><>  <><  ><>  <><  ><>  <><  ><>  
  //### Expecteds ###  <><  ><>  <><  ><>  <><  ><>
  //<><  ><>  <><  ><>  <><  ><>  <><  ><>  <><  ><>  
  {
    //Prepare Gauss-Laguerre Quadrature
    vector[NK]     kn_wt[2] = gauss_laguerre_quad(NK, Z_K-1.0);
    //calculate the expected density at each length point integrating over age
    vector[BN2_1]  PL = prob_length(Len, Galpha, Linf, Z_K, Sels, Sel50, 
                                    kn_wt[1], kn_wt[2]);
    //sum over each bin to get the normalised probability mass and raise to an expected number
    efq = prob_mass(Bnd, PL)*NV;
  }
  
  //<><  ><>  <><  ><>  <><  ><>  <><  ><>  <><  ><>  
  //### LIKELIHOOD ###  <><  ><>  <><  ><>  <><  ><>
  //<><  ><>  <><  ><>  <><  ><>  <><  ><>  <><  ><>  

  {
    target += neg_binomial_2_lpmf(fq | efq, NB_phi);
  }
} //model


generated quantities {
}


