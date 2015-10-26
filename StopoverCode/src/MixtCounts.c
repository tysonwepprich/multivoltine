// Built on command line using:
// R CMD SHLIB -o MixtCounts.so MixtCounts.c

#include <R.h>
#include <Rmath.h>


void MixtCountsLL(
    int    *counts,//counts 
    int    *SS, //number of sites
    int    *TT, //number of samples
    double *N,     // Parameter : superpopulation size   
    double *p,     // Parameter : detection probabilities 
    double *Beta,  // Parameter : entry proportions
    double *phi,   // Parameter : retention probabilities 
    double *LL,    //Poisson loglikelihood
    int    *resultcode)
  
{
  int  i, j, b, m;    
  double lambda, temp; //expected counts
  
  int Tcounts [*SS][*TT]; //counts matrix, element counts[i,j] is the count obtained for site i at time j
  for (j=0; j<(*TT); j++)
    for (i=0; i<(*SS); i++) 
      Tcounts[i][j] = counts[(*SS) * j + i];
  
  double Tp [*SS][*TT]; //p matrix, element p[i,j] is the probability of detecting an individual at site i, occasion j 
  for (j=0; j<(*TT); j++)
    for (i=0; i<(*SS); i++) 
      Tp[i][j] = p[(*SS) * j + i];   
  
  double TBeta [*SS][*TT]; //beta matrix, element beta[i,b] is the proportion of individuals at site i that first became available for detection at time b
  for (b=0; b<(*TT); b++)
    for (i=0; i<(*SS); i++) 
      TBeta[i][b] = Beta[(*SS) * b + i];      
  
  double Tphi [*TT-1][*TT-1][*SS]; //phi array, element phi[b,j,i] is the probability for individuals of cohort b surviving from time j to time j+1 on site i
  for(i=0; i<(*SS); i++) 
    for (j=0; j<(*TT-1); j++)
      for (b=0; b<(*TT-1); b++) 
        Tphi[b][j][i] = phi[i* (*TT-1) * (*TT-1) + (*TT-1) * j + b];
  
  
  *resultcode = 9;   // no result
  
  *LL = 0;
  
  for (i=0; i< *SS; i++)
  {
    for (j=0; j< *TT; j++)
    {
      lambda = 0;
      
      for(b=0; b<(j+1); b++)
      {
        temp = N[i]*Tp[i][j];
        
        temp *= TBeta[i][b];
        
        for (m=b; m<j; m++)
          temp *= Tphi[b][m][i];
        
        lambda += temp;
        
      }  
      
      if(Tcounts[i][j]> -1) *LL -= -lambda + Tcounts[i][j]*log(lambda) - lgamma(Tcounts[i][j]+1); 
      
    } 
  }
  
  
  
  
  *resultcode = 0;   // successful completion
}
