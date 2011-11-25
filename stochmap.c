/*stochmap.c: stochastic mapping on phylogenetic trees

 Copyright 2009 Matthew Spencer, Simon Whelan

 This file is part of stochmap.

    stochmap is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    stochmap is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with stochmap.  If not, see <http://www.gnu.org/licenses/>.

 **********************************

 This program uses substantial amounts of code from MrBayes 3.1.2:
 mbmath_sub.c (just added a few things from bayes.c, command.c, model.c, mb.c to mbmath.c, to avoid having to compile the whole of MrBayes)
 mb.h
 mbmath_sub.h
 globals.h

 **********************************
 This program uses my_getopt, for which the license is given below:

my_getopt - a command-line argument parser
Copyright 1997-2001, Benjamin Sittler

Permission is hereby granted, free of charge, to any person
obtaining a copy of this software and associated documentation
files (the "Software"), to deal in the Software without
restriction, including without limitation the rights to use, copy,
modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT.  IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.

*******************************************/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <stdarg.h>
#include "mb.h"
#include "mbmath_sub.h"
#include "globals.h"
#include "getopt.h" /*uses my_getopt because GNU getopt isn't available on all systems*/
#include "stochmap.h"

/* prototypes */
MrBFlt ***AllocatecondE(int nbranch,int nproc,int nsite);
MrBFlt *****AllocatePartials(int nbranch,int nproc,int nsite,int nstate);
MrBFlt ***AllocateQset(int nproc, int nstate);
int ****Allocatescalefact(int nbranch,int nproc, int nsite);
void computeEHD(MrBFlt **Fp, MrBFlt **Fc, MrBFlt **ENLtD, MrBFlt *L, MrBFlt *E, int m, int n, int **scalefact);
void computeENLt(int n, MrBFlt ** QL, MrBFlt t, MrBFlt **ENLt, MrBFlt **Eigvecs, MrBFlt **inverseEigvecs, MrBFlt *EigenValues);
void computeENLtD(int n, MrBFlt **ENLt, MrBFlt **Pt, MrBFlt **ENLtD);
MrBFlt computeKi(MrBFlt *EigenValues, MrBFlt t, int i);
MrBFlt computeIijt(int i, int j, MrBFlt t, MrBFlt *EigenValues);
MrBFlt computemu2(MrBFlt *pi_i,MrBFlt **QL,MrBFlt t,int m,MrBFlt **Eigvecs,MrBFlt **inverseEigvecs,MrBFlt *EigenValues);
void computepi(int m, MrBFlt **inverseEigvecs, MrBFlt *EigenValues, MrBFlt *pi_i);
void computepimix(int m, MrBFlt **inverseEigvecs, MrBFlt *EigenValues, MrBFlt *pi_i, int k, MrBFlt *mixprobs);
MrBFlt computepriorE(MrBFlt *pi_i, MrBFlt **QL, MrBFlt t, int m);
void computePt(int n, MrBFlt t, MrBFlt **Eigvecs, MrBFlt **inverseEigvecs, MrBFlt *EigenValues, MrBFlt **Pt);
MrBFlt computepriorV(MrBFlt priorE, MrBFlt *pi_i, MrBFlt **QL, MrBFlt t, int m, MrBFlt **Eigvecs, MrBFlt **inverseEigvecs, MrBFlt *EigenValues);
MrBFlt computepriorVE(MrBFlt **ENLtD, MrBFlt priorE,MrBFlt *pi_i,MrBFlt **Pt,int n);
void computeSi(int n, int ii, MrBFlt **Eigvecs, MrBFlt **inverseEigvecs, MrBFlt **Si);
MrBFlt computeZscore(MrBFlt xi,MrBFlt mu, MrBFlt V);
int CountBranches(char *partialfilename);
int CountCols( char *partialfilename);
int CountProc(char *partialfilename);
int CountSites(char *partialfilename);
int CountStates(char *partialfilename);
void DebugPartials(MrBFlt **Fp,MrBFlt **Fc, MrBFlt **Pt, int nstate, int nsite);
void EigenDecomp(int n, MrBFlt **Q, MrBFlt *eigenValues, MrBFlt **eigvecs, MrBFlt **inverseEigvecs);
void FreecondE(MrBFlt ***condE, int nbranch);
void FreePartials(MrBFlt *****partials,int nbranch,int nproc);
void FreeQset(MrBFlt ***Qset, int nproc);
void Freescalefact(int ****scalefact, int nbranch, int nproc, int nsite);
int GetQMatrixSize(char *qfilename);
void makeQL (MrBFlt **Q, int **L, MrBFlt **QL, int dim);
void printHelp(void);
void ReadLMatrix(char *lfilename, int **L, int n);
void ReadMap(int *sitemap,int ncols, char *line);
void ReadOneSite(FILE *partial_fv, int branch, int nstate, int site, MrBFlt *****partials, int ****scalefact);
void ReadParams(int argc, char **argv, settings *sets, char *infilename, char *outfilename, char *lfilename);
void ReadPartialLine(char *line, MrBFlt *****partials, int branch, int proc, int which, int site, int nstate, int ****scalefact);
void ReadPartials(char *partialfilename, MrBFlt ***Qset, MrBFlt *****partials, MrBFlt **sitelikes, MrBFlt *tbranch, MrBFlt *mixprobs, int nstate, int nproc, int nsite, int *multiplicities, int ****scalefact, int *sitemap, int ncols, MrBFlt **pi_i);
void ReadPi(FILE *fv, int nstate, MrBFlt **pi_i, int proc);
void ReadQMatrix(FILE *fv, MrBFlt **Q, int nstate);
int ReadSiteLogLikes(char *line, MrBFlt **sitelikes, int nproc, int *multiplicities);
void ReadSites(FILE *partial_fv, MrBFlt **sitelikes, int branch, int nstate, int nsite, MrBFlt *****partials, int nproc, int *multiplicities, int ****scalefact);
void testEHD(MrBFlt **Fp,MrBFlt **Fc, MrBFlt *L);
void WeightByPi(MrBFlt **F,MrBFlt *pi_i, int nsite, int nstate);
void WriteResults(FILE *outfile,MrBFlt *****partials,int nbranch, int nproc, int nsite, MrBFlt ***condE, MrBFlt **priorE, MrBFlt **priorV, int *multiplicities, int *sitemap, int ncols, MrBFlt *tbranch);

/*allocate space for conditional expectations*/
MrBFlt ***AllocatecondE(int nbranch,int nproc,int nsite)
{
  int i;
  MrBFlt ***condE;

  condE=(MrBFlt ***)SafeMalloc((size_t)((nbranch+1)*sizeof(MrBFlt**)));
  if (!condE)
    {
      MrBayesPrint ("%s   Error: Problem allocating conditional expectations.\n", spacer);
      exit(1);
    }
  for(i=0;i<nbranch;i++){
    condE[i]=AllocateDoubleMatrix(nproc+1,nsite);/*extra row is for sum*/
  }
  return condE;
}

/*allocate space for partial likelihoods*/
MrBFlt *****AllocatePartials(int nbranch,int nproc,int nsite,int nstate)
{
  int i,j,k;
  MrBFlt *****partials;
  
  partials=(MrBFlt *****)SafeMalloc((size_t)((nbranch)*sizeof(MrBFlt****)));
  if (!partials)
    {
      MrBayesPrint ("%s   Error: Problem allocating partials.\n", spacer);
      exit(1);
    }
  for(i=0;i<nbranch;i++){/*each branch*/
    partials[i]=(MrBFlt ****)SafeMalloc((size_t)((nproc)*sizeof(MrBFlt***)));
    if (!partials[i])
      {
	MrBayesPrint ("%s   Error: Problem allocating partials.\n", spacer);
	exit(1);
      }    
    for(j=0;j<nproc;j++){/*each process*/
      partials[i][j]=(MrBFlt ***)SafeMalloc((size_t)((2)*sizeof(MrBFlt**)));
      if (!partials[i][j])
	{
	  MrBayesPrint ("%s   Error: Problem allocating partials.\n", spacer);
	  exit(1);
	}
      for(k=0;k<2;k++){/*left and right*/
	partials[i][j][k]=AllocateDoubleMatrix(nsite,nstate);/*rows are sites, cols are states*/
      }
    }
  }
  return partials;
}

/*allocate space for a Q matrix for each process*/
MrBFlt ***AllocateQset(int nproc, int nstate)
{
  int i;
  MrBFlt ***Qset;

  Qset=(MrBFlt ***)SafeMalloc((size_t)((nproc)*sizeof(MrBFlt**)));
  if (!Qset)
    {
      MrBayesPrint ("%s   Error: Problem allocating Qset.\n", spacer);
      exit(1);
    }
  for(i=0;i<nproc;i++) Qset[i] = AllocateSquareDoubleMatrix (nstate);
  return Qset;
}

/*allocate space for scale factors*/
int ****Allocatescalefact(int nbranch,int nproc, int nsite)
{
  int i,j;
  int ****scalefact;

  scalefact=(int ****)SafeMalloc((size_t)((nbranch)*sizeof(int***)));
  if (!scalefact)
    {
      MrBayesPrint ("%s   Error: Problem allocating scale factors.\n", spacer);
      exit(1);
    }
  for(i=0;i<nbranch;i++){
    scalefact[i]=(int ***)SafeMalloc((size_t)((nproc)*sizeof(int**)));
    if (!scalefact)
      {
	MrBayesPrint ("%s   Error: Problem allocating scale factors.\n", spacer);
	exit(1);
      }
    for(j=0;j<nproc;j++){
      scalefact[i][j]=AllocateIntegerMatrix(2,nsite);
    }
  }
  return scalefact;
}



/*Equations 3.8 and 3.4 in Minin and Suchard 2008 Phil Trans B 363:3985-3995*/
/*Fp and Fc are partial likelihoods at each end of a branch of interest (rows are SITES, cols are STATES*/
/*ENLtD is expected number of labelled changes given states i,j at each end of branch (m x m matrix)*/
/*L is a vector of site likelihoods (1 x n): scaled by the values in scalefact*/
/*E is a vector in which the expected number of labelled changes on each site on this edge, conditional on the data in the tree, is stored*/
/*m is the number of states (observable x hidden)*/
/*n is the number of sites*/
/*scalefact is scaling factors for likelihoods for left and right partials at each site, such that we multiply the numbers in L by 10^{-scaling factor} */
void computeEHD(MrBFlt **Fp, MrBFlt **Fc, MrBFlt **ENLtD, MrBFlt *L, MrBFlt *E, int m, int n, int **scalefact)
{
  int i,j,s;
  
  for(s=0;s<n;s++){/*each site*/
    E[s]=0.;
    /*compute 3.8 for this site*/
    for(i=0;i<m;i++){
      for(j=0;j<m;j++){
	E[s]+=Fp[s][i]*ENLtD[i][j]*Fc[s][j];
      }
    }
    E[s]=exp(-(scalefact[0][s]+scalefact[1][s])*LOGTEN+log(E[s])-log(L[s]));/*multiply by scale factors, divide by site likelihood*/
    /*E[s]=exp(-(scalefact[0][s]+scalefact[1][s])*LOGTEN+log(E[s]));*//*TEST multiply by scale factors but don't divide by site likelihood: should get log likelihood out if we put partials in*/
    /*printf("Warning: not dividing by L=%f\n",L[s]);*/
  }
}


/*Equation 2.4 in Minin and Suchard 2008 Phil Trans B 363:3985-3995
  The result is put in ENLt. This version matches the way joint.mean.markov.jumps in the R package markovjumps does this*/
void computeENLt(int n, MrBFlt ** QL, MrBFlt t, MrBFlt **ENLt, MrBFlt **Eigvecs, MrBFlt **inverseEigvecs, MrBFlt *EigenValues)
{
  int i,j;
  MrBFlt Iijt;

  MrBFlt **A, **B;

  /*temporary storage*/
  A=AllocateSquareDoubleMatrix(n);
  B=AllocateSquareDoubleMatrix(n);

  /* we want U * [I_{ij} .* (U^{-1} * Q_L * U)] * U^{-1} */
  MultiplyMatrices(n,QL,Eigvecs,A);/*A = Q_L * U */
  MultiplyMatrices(n,inverseEigvecs,A,B);/*B = U^{-1} * A */
  
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      Iijt=computeIijt(i,j,t,EigenValues);
      A[i][j]=Iijt*B[i][j]; /*C = I_{ij} .* B */
    }
  }
  
  MultiplyMatrices(n,A,inverseEigvecs,B);/*D = C * U^{-1} */
  MultiplyMatrices(n,Eigvecs,B,ENLt);/*E = U * D */

  FreeSquareDoubleMatrix(A);
  FreeSquareDoubleMatrix(B);
}

/*conditional expectation of number of labelled transitions, given start and end states and time*/
void computeENLtD(int n, MrBFlt **ENLt, MrBFlt **Pt, MrBFlt **ENLtD)
{
  int i,j;

  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      ENLtD[i][j]=ENLt[i][j]/Pt[i][j];
}


/* Equation 2.5 in Minin and Suchard 2008 Phil Trans B 363:3985-3995*/
MrBFlt computeIijt(int i, int j, MrBFlt t, MrBFlt *EigenValues)
{
  MrBFlt tol=1e-15;

  if(AreDoublesEqual(EigenValues[i],EigenValues[j],tol))
    return(t*exp(EigenValues[i]*t));
  else
    return((exp(EigenValues[i]*t)-exp(EigenValues[j]*t))/(EigenValues[i]-EigenValues[j]));
}

/*used to get second factorial moment of prior number of events*/
MrBFlt computeKi(MrBFlt *EigenValues, MrBFlt t, int i)
{
  MrBFlt K;
  MrBFlt tol=1e-15;

  if(AreDoublesEqual(EigenValues[i],0.,tol))
    K=t*t/2.;
  else
    K=(exp(EigenValues[i]*t)-1)/(EigenValues[i]*EigenValues[i])-t/EigenValues[i];
  return(K);
}

/*second factorial moment of prior number of labelled events, summed over initial states, weighted by stationary probabilities*/
MrBFlt computemu2(MrBFlt *pi_i,MrBFlt **QL,MrBFlt t,int m,MrBFlt **Eigvecs,MrBFlt **inverseEigvecs,MrBFlt *EigenValues)
{
  MrBFlt mu2,K;
  MrBFlt **A, **B;

  int i,j;

  /*temporary storage*/
  A=AllocateSquareDoubleMatrix(m);
  B=AllocateSquareDoubleMatrix(m);

  for(j=0;j<m;j++){                       /*B = U*K */
    K=computeKi(EigenValues,t,j);
    for(i=0;i<m;i++)
      B[i][j]=Eigvecs[i][j]*K;
  }
  MultiplyMatrices(m,B,inverseEigvecs,A); /*A = B*U^{-1} */
  MultiplyMatrices(m,A,QL,B);             /*B = A*QL */
  MultiplyMatrices(m,QL,B,A);             /*A = QL*B */

  mu2=0.;
  for(i=0;i<m;i++)                        /*mu2 = pi' * A * 1 */
    for(j=0;j<m;j++)
      mu2+=pi_i[i]*A[i][j];
  mu2=mu2*2.;
  FreeSquareDoubleMatrix(A);
  FreeSquareDoubleMatrix(B);

  return(mu2);

}

/*compute stationary distribution given eigen decomposition, assuming there is one zero eigenvalue*/
void computepi(int m, MrBFlt **inverseEigvecs, MrBFlt *EigenValues, MrBFlt *pi_i)
{
  int i,j;
  MrBFlt sum;
  MrBFlt tol=1e-6;

  printf("Warning: eigenvalue tol=%e\n",tol);
  i=-1;
  for(j=0;j<m;j++){
    if(AreDoublesEqual(EigenValues[j],0.,tol)){/*approximately zero eigenvalue*/
      if(i==-1)
	i=j;
      else/*we already found a zero eigenvalue, and there should only be one*/
	i=-2;
    }
  }
  if(i==-1){/*we didn't find a zero eigenvalue*/
    MrBayesPrint ("%s   Error: no zero eigenvalue found. Check tolerance, currently %g\n", spacer, tol);
    exit(1);	
  }
  else if(i==-2){/*we found several zero eigenvalues*/
    MrBayesPrint ("%s   Error: more than one zero eigenvalue found. Check tolerance, currently %g\n", spacer, tol);
    exit(1);	
  }
  else{/*we found one zero eigenvalue*/
    /*find the sum of the corresponding eigenvector (note we want the complex conjugate of the corresponding row of the inverse eigenvector matrix)*/
    sum=0.;
    for(j=0;j<m;j++)
      sum+=inverseEigvecs[i][j];/*assumed this eigenvector is real*/
    
    /*rescale the eigenvector*/
    for(j=0;j<m;j++)
      pi_i[j]=inverseEigvecs[i][j]/sum;
  }
}

/*compute stationary distribution given eigen decomposition, assuming there is more than one zero eigenvalue in a block-diagonal rate matrix. mixprobs should give the mixing probabilities in the order in which the zero eigenvalues corresponding to the categories are found (don't know for sure that this will match the block order)*/
void computepimix(int m, MrBFlt **inverseEigvecs, MrBFlt *EigenValues, MrBFlt *pi_i, int k, MrBFlt *mixprobs)
{
  int i,j,h;
  MrBFlt sum;
  MrBFlt tol=1e-9;

  for(j=0;j<m;j++)
    pi_i[j]=0.;
  i=0;
  for(j=0;j<m;j++){
    if(AreDoublesEqual(EigenValues[j],0.,tol)){/*approximately zero eigenvalue*/
      if(i<k){
	sum=0.;
	for(h=0;h<m;h++)/*sum of this eigenvector*/
	  sum+=inverseEigvecs[j][h];/*in general we want the complex conjugate of this row, but we assumed this eigenvector is real*/
	for(h=0;h<m;h++)/*rescale the eigenvector and weight by mix prob*/
	  pi_i[h]+=mixprobs[i]*inverseEigvecs[j][h]/sum;
      }
      i++;
    }
  }
  if(i!=k){/*we didn't find the right number of eigenvalues*/
    MrBayesPrint ("%s   Error: %d zero eigenvalues found, %d expected. Check tolerance (currently %g) and model\n", spacer,i,k,tol);
    exit(1);	
  }
}

/*compute prior expectation of number of labelled changes on edge (for a single site)*/
MrBFlt computepriorE(MrBFlt *pi_i, MrBFlt **QL, MrBFlt t, int m)
{
  int i,j;
  MrBFlt Eprior;

  Eprior=0.;
  for(j=0;j<m;j++)
    for(i=0;i<m;i++)
      Eprior+=pi_i[i]*QL[i][j];
  Eprior*=t;
  return(Eprior);
}

/*compute transition probability, given eigen decomposition*/
/*The result is in Pt*/
void computePt(int n, MrBFlt t, MrBFlt **Eigvecs, MrBFlt **inverseEigvecs, MrBFlt *EigenValues, MrBFlt **Pt)
{
  int i,j;

  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
	Pt[i][j]=inverseEigvecs[i][j]*exp(EigenValues[i]*t);
  MultiplyMatrices(n,Eigvecs,Pt,Pt);
}

/*compute prior variance of number of labelled changes on edge (for a single site)*/
MrBFlt computepriorV(MrBFlt priorE, MrBFlt *pi_i, MrBFlt **QL, MrBFlt t, int m, MrBFlt **Eigvecs, MrBFlt **inverseEigvecs, MrBFlt *EigenValues)
{
  MrBFlt Vprior, mu2;

  mu2=computemu2(pi_i,QL,t,m,Eigvecs,inverseEigvecs,EigenValues);/*second factorial moment summed over initial states, weighted by stationary probabilities*/
  Vprior=mu2+priorE-priorE*priorE;/*variance from factorial moments*/
  return(Vprior);
}

/*compute prior variance (over start and end states) of MEAN number of labelled changes on edge (for a single site)*/
MrBFlt computepriorVE(MrBFlt **ENLtD, MrBFlt priorE,MrBFlt *pi_i,MrBFlt **Pt,int n)
{
  int i,j;
  MrBFlt priorVE;

  priorVE=0.;

  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      priorVE+=pi_i[i]*Pt[i][j]*ENLtD[i][j]*ENLtD[i][j];
  priorVE-=priorE*priorE;
  return(priorVE);
}

/*compute S_i = U*E_i*inv(U), where U is eigenvectors and E_i is a matrix of zeros, except entry e_ii=1*/
void computeSi(int n, int ii, MrBFlt **Eigvecs, MrBFlt **inverseEigvecs, MrBFlt **Si)
{
  int i,j;

  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      Si[i][j]=Eigvecs[i][ii]*inverseEigvecs[ii][j];
}

/*z score: xi is observation, mu is mean, V is variance*/
MrBFlt computeZscore(MrBFlt xi,MrBFlt mu, MrBFlt V)
{
  return((xi-mu)/sqrt(V));
}

/*count the number of branches in the partial likelihood file*/
int CountBranches(char *partialfilename)
{
  FILE *partial_fv;
  char line[LINELEN];
  int nbranch=0;

  /*open the partial file*/
  partial_fv=fopen(partialfilename,"r");
  if(partial_fv==NULL){
    MrBayesPrint ("%s   Error: Problem opening partial likelihood file.\n", spacer);
    exit(1);
  }
  while(fgets(line,sizeof(line),partial_fv)!=NULL){
    if(strstr(line,"Branch[")!=0)/*we found the start of a branch block*/
      nbranch++;
  }
  fclose(partial_fv);
  return(nbranch);  
}

/*count the number of columns in the original alignment*/
int CountCols(char *partialfilename)
{
  int c,f;
  FILE *partial_fv;
  char line[LINELEN];
  char *token=NULL, *delim=" \t\r\n";

  /*open the partial file*/
  partial_fv=fopen(partialfilename,"r");
  if(partial_fv==NULL){
    MrBayesPrint ("%s   Error: Problem opening partial likelihood file.\n", spacer);
    exit(1);
  }
  c=0;
  f=0;
  while(fgets(line,sizeof(line),partial_fv)!=NULL){
    if(strstr(line,"Mapping to real data:")!=0){/*we found the mapping line*/
      f=1;
      token=strtok(line,":");/*get rid of "Mapping to real data:" token, which we know exists*/
      token=strtok(NULL,delim);
      while(token!=NULL){
	c++;
	token=strtok(NULL,delim);
      }
    }
  }
  if(f==0){
    MrBayesPrint ("%s   Error: no site mapping in partial likelihood file.\n", spacer);
    exit(1);
  }
  fclose(partial_fv);
  return(c);
}

/*count the number of processes*/
int CountProc(char *partialfilename)
{
  FILE *partial_fv;
  char line[LINELEN];
  int nproc=0;
  int p,args_assigned;

  /*open the partial file*/
  partial_fv=fopen(partialfilename,"r");
  if(partial_fv==NULL){
    MrBayesPrint ("%s   Error: Problem opening partial likelihood file.\n", spacer);
    exit(1);
  }
  while(fgets(line,sizeof(line),partial_fv)!=NULL){
    if((args_assigned=sscanf(line,"Process[%d] QMat:",&p))==1){/*we found the start of a Q matrix block*/
	nproc++;
      }
    if(strstr(line,"Branch[")!=0){/*we found the start of a branch block, should be no more Q matrices*/
      break;
    }
  }
  fclose(partial_fv);
  if(nproc>1){
    MrBayesPrint ("%s   Error: read %d processes. Multiple processes not yet implemented.\n", spacer, nproc);
    exit(1);
  }
  return(nproc);
}

/*count the number of sites in the partial likelihood file*/
int CountSites(char *partialfilename)
{
  FILE *partial_fv;
  char line[LINELEN];
  int nsite=0;

  /*open the partial file*/
  partial_fv=fopen(partialfilename,"r");
  if(partial_fv==NULL){
    MrBayesPrint ("%s   Error: Problem opening partial likelihood file.\n", spacer);
    exit(1);
  }
  while(fgets(line,sizeof(line),partial_fv)!=NULL){
    if(strstr(line,"Branch[0")!=0){/*we found the start of branch 0*/
      while(fgets(line,sizeof(line),partial_fv)!=NULL){
	if(strstr(line,"//")) break;/*end of branch block*/
	if(strstr(line,"Site[")) nsite++;/*start of site*/
      }
    }
  }
  fclose(partial_fv);
  return(nsite);
}

/*count the number of states in the Q matrix for process 0*/
int CountStates(char *partialfilename)
{
  FILE *partial_fv;
  char line[LINELEN];
  int nstates=0;

  /*open the partial file*/
  partial_fv=fopen(partialfilename,"r");
  if(partial_fv==NULL){
    MrBayesPrint ("%s   Error: Problem opening partial likelihood file.\n", spacer);
    exit(1);
  }
  while(fgets(line,sizeof(line),partial_fv)!=NULL){
    if(strstr(line,"Process[0] QMat:")!=0){/*we found the start of the process 0 matrix*/
      while(fgets(line,sizeof(line),partial_fv)!=NULL){
	if(strstr(line,"Eqm")) break;/*end of Q matrix block*/
	else nstates++;/*assuming there are no blank lines in the block*/
      }
    }
  }
  fclose(partial_fv);
  return(nstates);
}

/*are we doing multiplications in the right way?*/
void DebugPartials(MrBFlt **Fp,MrBFlt **Fc, MrBFlt **Pt, int nstate, int nsite)
{
  int i,j,s;
  MrBFlt **Dp, **Dc;

  Dp=AllocateDoubleMatrix(nsite,nstate);
  Dc=AllocateDoubleMatrix(nsite,nstate);
  
  for(s=0;s<nsite;s++){/*each site*/
    for(i=0;i<nstate;i++){
      Dp[s][i]=0.;
      Dc[s][i]=0.;
    }
    for(i=0;i<nstate;i++){
      for(j=0;j<nstate;j++){
	Dp[s][i]+=Fp[s][j]*Pt[i][j];
	Dc[s][i]+=Fc[s][j]*Pt[i][j];
      }
    }  
  }

/*   for(s=0;s<n;s++){/\*each site*\/ */
/*     E[s]=0.; */
/*     /\*compute 3.8 for this site*\/ */
/*     for(i=0;i<m;i++){ */
/*       for(j=0;j<m;j++){ */
/* 	E[s]+=Fp[s][i]*ENLtD[i][j]*Fc[s][j]; */
/*       } */
/*     } */

  printf("Fp\n");
  for(s=0;s<nsite;s++){
    for(i=0;i<nstate;i++)
      printf("%f ",Fp[s][i]);
    printf("\n");
  }
  printf("Fc\n");
  for(s=0;s<nsite;s++){
    for(i=0;i<nstate;i++)
      printf("%f ",Fc[s][i]);
    printf("\n");
  }

  FreeDoubleMatrix(Dp);
  FreeDoubleMatrix(Dc);
}

/*eigen decomposition of rate matrix*/
void EigenDecomp(int n, MrBFlt **Q, MrBFlt *EigenValues, MrBFlt **Eigvecs, MrBFlt **inverseEigvecs)
{
  MrBFlt **Qw;
  int isComplex,i,j;
  complex **Ceigvecs, **CinverseEigvecs;
  MrBFlt *eigvalsImag;

  Ceigvecs = AllocateSquareComplexMatrix (n);
  CinverseEigvecs = AllocateSquareComplexMatrix (n);
  Qw=AllocateSquareDoubleMatrix(n);
  eigvalsImag=AllocateDoubleVector(n);
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      Qw[i][j]=Q[i][j];/*Q gets destroyed in decomposition*/
  isComplex = GetEigens (n, Qw, EigenValues, eigvalsImag, Eigvecs, inverseEigvecs, Ceigvecs, CinverseEigvecs);
  free(eigvalsImag);
  FreeSquareComplexMatrix(Ceigvecs);	
  FreeSquareComplexMatrix(CinverseEigvecs);	
  FreeSquareDoubleMatrix(Qw);

  if(isComplex){/*we didn't implement anything that will deal with complex eigenvectors*/
    MrBayesPrint ("%s   Error: Q matrix has complex eigenvectors. Is the substitution model reversible?\n", spacer);
    exit(1);
  }
}

/*free the conditional expectations*/
void FreecondE(MrBFlt ***condE, int nbranch)
{
  int i;

  for(i=0;i<nbranch;i++) FreeDoubleMatrix(condE[i]);
  free(condE);
}


/*free space for partial likelihoods*/
void FreePartials(MrBFlt *****partials,int nbranch,int nproc)
{
  int i,j,k;

  for(i=0;i<nbranch;i++){
    for(j=0;j<nproc;j++){
      for(k=0;k<2;k++)
	FreeDoubleMatrix(partials[i][j][k]);
      free(partials[i][j]);
    }
    free(partials[i]);
  }
  free(partials);
}

/*free the Q matrices*/
void FreeQset(MrBFlt ***Qset, int nproc)
{
  int i;

  for(i=0;i<nproc;i++) FreeSquareDoubleMatrix(Qset[i]);
  free(Qset);
}

/*free the scale factors*/
void Freescalefact(int ****scalefact, int nbranch, int nproc, int nsite)
{
  int i,j;

  for(i=0;i<nbranch;i++){
    for(j=0;j<nproc;j++){
      FreeIntegerMatrix(scalefact[i][j]);
    }
    free(scalefact[i]);
  }
  free(scalefact);
}


/*size of Q matrix*/
int GetQMatrixSize(char *qfilename)
{
  int n;
  FILE *fp;
  double x;
  
  n=0;
  fp=fopen(qfilename,"r");
  if(fp==NULL){
    MrBayesPrint ("%s   Error: Problem opening Q matrix file.\n", spacer);
    exit(1);
  }
  else{
    while(!feof(fp)){
      fscanf(fp,"%lf",&x);
      n++;
    }
    n=sqrt(n);/*assuming the matrix is square*/
    fclose(fp);
  }
  return(n);
}

/* construct the QL matrix, where ql_ij = q_ij * l_ij (l_ij is 0 or 1)*/
void makeQL (MrBFlt **Q, int **L, MrBFlt **QL, int dim)
{
  int i,j;
  for(i=0;i<dim;i++)
    for(j=0;j<dim;j++)
      {
	if(L[i][j]==1){
	  QL[i][j]=Q[i][j];
	}
	else{
	  QL[i][j]=0.0;
	}
      }
}


/*display help message*/
void printHelp(){
  printf(PROGRAM_NAME);
  printf(" ");
  printf(PROGRAM_VERSION);

  printf("\nUsage: \n");
  printf(PROGRAM_NAME);
  printf(" [-h] [[infile] [outfile] [lfile]] \n");
  printf("Options: \n");
  printf("-h: this help message\n");
  printf("Input file format: a .partL file generated by Leaphy.\n");
  printf("Output file format: tab-separated text file containing conditional and prior expectations of numbers of labelled events on site and edges. Described in detail in the README.txt file accompanying this code.\n");
  printf("L file format: a whitespace-delimited square matrix of 1s and 0s (same size as the instantaneous rate matrix for the model used to generate the input file) with 1s indicating transitions that we want to count. Should be symmetrical, with 0 on the diagonal.\n");
  printf(COPYRIGHT);

}

/*read the L matrix from a text file.*/
void ReadLMatrix(char *lfilename, int **L, int n)
{
  FILE *fp;
  int i,j;

  fp=fopen(lfilename,"r");
  if(fp==NULL){
    MrBayesPrint ("%s   Error: Problem reading L matrix.\n", spacer);
    exit(1);
  }
  else{
    i=0;
    j=0;
    while(i<n && (fscanf(fp,"%d",&L[i][j]))!= EOF){
      j++;
      if(j==n){
	j=0;
	i++;
      }
    }
    if(i<n){
      MrBayesPrint ("%s   Error: could not read all elements of L matrix. Failed at i=%d, j=%d\n", spacer,i,j);
      exit(1);		
    }
  }
  fclose(fp);
}

/*read the site map*/
void ReadMap(int *sitemap,int ncols, char *line)
{
  int i;
  char *token=NULL, *delim=" \t\r\n";

  i=0;
  token=strtok(line,":");/*get rid of "Mapping to real data:" token, which we know exists*/
  token=strtok(NULL,delim);
  while(token!=NULL){
    sitemap[i] = atoi(token);
    if(i>ncols){
      MrBayesPrint ("%s   Error: too many entries in site map.\n", spacer);
      exit(1);
    }
    i++;
    token=strtok(NULL,delim);
  }
}

/*read partial likelihood information for a single site*/
void ReadOneSite(FILE *partial_fv, int branch, int nstate, int site, MrBFlt *****partials, int ****scalefact)
{
  char line[LINELEN];
  int args_assigned, proc;

  /*read each process*/
  if(fgets(line,sizeof(line),partial_fv)!=NULL){
    if((args_assigned=sscanf(line," Process[%d] Left",&proc))==1){
      if(fgets(line,sizeof(line),partial_fv)!=NULL)/*read the left partial vector*/
	ReadPartialLine(line, partials, branch, proc, 0, site, nstate,scalefact);
      else{
	MrBayesPrint ("%s   Error: couldn't read left partial likelihood line for branch %d, process %d, site %d.\n", spacer, branch, proc, site);
	exit(1);		
      }
    }
    else{
      MrBayesPrint ("%s   Error: expected process number, found %s.\n", spacer, line);
      exit(1);
    }
  }
  else{
    MrBayesPrint ("%s   Error: couldn't read left partial likelihood for branch %d, process %d, site %d.\n", spacer, branch, proc, site);
       exit(1);
  }
  if(fgets(line,sizeof(line),partial_fv)!=NULL){
    if((args_assigned=sscanf(line," Process[%d] Right",&proc))==1){
      if(fgets(line,sizeof(line),partial_fv)!=NULL)/*read the right partial vector*/
	ReadPartialLine(line, partials, branch, proc, 1, site, nstate, scalefact);
      else{
	MrBayesPrint ("%s   Error: couldn't read right partial likelihood line for branch %d, process %d, site %d.\n", spacer, branch, proc, site);
	exit(1);		
      }
    }
    else{
      MrBayesPrint ("%s   Error: expected process number, found %s.\n", spacer, line);
      exit(1);
    }
  }
  else{
    MrBayesPrint ("%s   Error: couldn't read right partial likelihood for branch %d, process %d, site %d.\n", spacer, branch, proc, site);
    exit(1);
  }
}

void ReadParams(int argc, char **argv, settings *sets, char *infilename, char *outfilename, char *lfilename)
{
  int c,offset=0;
  extern char *optarg;
  sets->printHelp=FALSE;

  opterr = 0;
  while ((c = getopt (argc, argv, "h")) != -1) switch (c){
  case 'h':/*help message*/
    sets->printHelp=TRUE;
    return;
    break;
  case '?':
    if (isprint (optopt))
      fprintf (stderr, "Unknown option `-%c'.\n", optopt);
    else
      fprintf (stderr,"Unknown option character `\\x%x'.\n",optopt);
    exit(EXIT_FAILURE);
  default:
    abort ();
}
  if(offset==0){
    if(argc > optind){
      sprintf(infilename,"%s",argv[optind]);
    }
    else{/* ask for the input-file */
      printf("\nEnter input file : ");
      scanf("%s",infilename);
    }
  }
  if(argc>(optind+1+offset)){
    sprintf(outfilename,"%s",argv[optind+1+offset]);
  }
  else{/* ask for output-file.        */
    printf("\nEnter output file : ");
    scanf("%s",outfilename);
  }
  if(argc>(optind+2+offset)){
    sprintf(lfilename,"%s",argv[optind+2+offset]);
  }
  else{/* ask for L matrix file.        */
    printf("\nEnter L matrix file : ");
    scanf("%s",lfilename);
  }
}

/*read a line of partial likelihoods*/
void ReadPartialLine(char *line, MrBFlt *****partials, int branch, int proc, int which, int site, int nstate, int ****scalefact)
{
  int scale, i;
  char *token=NULL, *unconverted, *delim=" \t\r\n";

  if(!sscanf(line," Scale %d: ",&scale)){
    MrBayesPrint ("%s   Error: couldn't read scale from %s.\n", spacer, line);
    exit(1);
  }
  scalefact[branch][proc][which][site]=scale;
  i=0;
  token=strtok(line,":");/*get rid of scale token, which we know exists*/
  token=strtok(NULL,delim);
  while(token!=NULL){
    partials[branch][proc][which][site][i] = strtod(token, &unconverted);/*will rescale this later*/
    if (!isspace(*unconverted) && *unconverted != 0)
      {
	MrBayesPrint ("%s   Error: couldn't read branch %d, proc %d, which %d, site %d, state %d from %s.\n", spacer, branch, proc, which, site, i, line);
	exit(1);
      }
    if(i>nstate){
      MrBayesPrint ("%s   Error: too many states for branch %d, proc %d, which %d, site %d from %s.\n", spacer, branch, proc, which, site, line);
      exit(1);
    }
    i++;
    token=strtok(NULL,delim);
  }
}

/*read the partial likelihoods from a file*/
void ReadPartials(char *partialfilename,MrBFlt ***Qset,MrBFlt *****partials,MrBFlt **sitelikes,MrBFlt *tbranch,MrBFlt*mixprobs, int nstate, int nproc, int nsite, int *multiplicities, int ****scalefact, int *sitemap, int ncols, MrBFlt **pi_i)
{
  FILE *partial_fv;
  char line[LINELEN];
  int args_assigned,branch=0,proc=0,procfound=0,dumi;

  /*open the partial file*/
  partial_fv=fopen(partialfilename,"r");
  if(partial_fv==NULL){
    MrBayesPrint ("%s   Error: Problem opening partial likelihood file.\n", spacer);
    exit(1);
  }
  
  while(fgets(line,sizeof(line),partial_fv)!=NULL){
    if(strstr(line,"Mapping to real data:")!=0)/*we found the mapping line*/
      ReadMap(sitemap,ncols,line);/*read the site map*/

    if((args_assigned=sscanf(line,"Process[%d] QMat:",&proc))==1){/*we found a Q matrix block*/
      procfound++;
      ReadQMatrix(partial_fv, Qset[proc], nstate);/*read the Q matrix*/
      ReadPi(partial_fv,nstate,pi_i,proc);/*look for stationary distribution*/
    }

    /*now look for the branches*/
    if((args_assigned=sscanf(line,"Branch[%d]",&branch))==1){/*found the start of a branch block*/
      while(fgets(line,sizeof(line),partial_fv)!=NULL){/*look for branch length*/     
	if((args_assigned=sscanf(line," Branch %d P(%lf)",&dumi,&tbranch[branch]))==2){/*found the start of a branch block*/
	  break;
	}
      }
      ReadSites(partial_fv, sitelikes, branch, nstate, nsite, partials, nproc, multiplicities,scalefact);/*read the sites for this branch*/
    }
  }
  fclose(partial_fv);
  if(procfound!=nproc){
    MrBayesPrint ("%s   Error: expected %d Q matrices, found %d.\n", spacer, nproc, proc);
    exit(1);
  }
}

/*read the stationary distribution from a text file*/
void ReadPi(FILE *fv, int nstate, MrBFlt **pi_i, int proc)
{
  int i,found=0;
  char line[LINELEN];
  char *token=NULL, *unconverted, *delim=" \t\r\n";

  while(fgets(line,sizeof(line),fv)!=NULL){
    if(strstr(line,"//")!=0){
      break;/*end of block*/
    }
    if(strstr(line,"Eqm:")!=0){/*we found the line*/
      i=0;
      found=1;
      token=strtok(line,":");/*get rid of first token, which we know exists*/
      token=strtok(NULL,delim);
      while(token!=NULL){
	pi_i[proc][i] = strtod(token, &unconverted);
	if (!isspace(*unconverted) && *unconverted != 0){
	  MrBayesPrint ("%s   Error: couldn't read stationary probability %d for process %d from %s.\n", spacer, i, proc, line);
	  exit(1);
	}
	if(i>nstate){
	  MrBayesPrint ("%s   Error: too many states for stationary distribution of process %d from %s.\n", spacer, proc, line);
	  exit(1);
	}
	i++;
	token=strtok(NULL,delim);
      }
      break;
    }
    if(found==0){
      MrBayesPrint ("%s   Error: couldn't find stationary distribution for process %d.\n", spacer, proc);
      exit(1);
    }
  }
}

/*read the Q matrix from a text file*/
void ReadQMatrix(FILE *fv, MrBFlt **Q, int nstate)
{
  int i,j;
  char line[LINELEN];
  char *token=NULL, *unconverted, *delim=" \t\r\n";

  i=0;
  j=0;
  while(fgets(line,sizeof(line),fv)!=NULL){
    token=strtok(line,delim);
    while(token!=NULL){
      Q[i][j] = strtod(token, &unconverted);
       if (!isspace(*unconverted) && *unconverted != 0){
	MrBayesPrint ("%s   Error: couldn't read Q[%d][%d] from %s.\n", spacer,i,j,line);
	exit(1);
      }
      if(j>nstate){
	MrBayesPrint ("%s   Error: too many states for row %d of Q from %s.\n", spacer, i,line);
	exit(1);
      }
      j++;
      token=strtok(NULL,delim);
    }
    i++;
    j=0;
    if(i==nstate)
      break;
  }
}

/*read a line of site log likelihoods: return -1 if we didn't find a site on this line, site number otherwise*/
int ReadSiteLogLikes(char *line, MrBFlt **sitelikes, int nproc, int *multiplicities)
{
  int site,args_assigned,i;
  char *token=NULL, *unconverted, *delim=" []\t\r\n";

  site=-1;
  if((args_assigned=sscanf(line,"Site[%d] ",&site))==1){/*read site number*/
    token=strtok(line,delim);/*get rid of info at start: safe because we know they exist from sscanf*/
    token=strtok(NULL,delim);
    token=strtok(NULL,delim);
    i=0;
    while(token!=NULL){
      if(i<nproc){
	sitelikes[i][site]=exp(strtod(token,&unconverted));/*convert from log likelihood to likelihood*/
	if (!isspace(*unconverted) && *unconverted != 0)
	  {
	    MrBayesPrint ("%s   Error: couldn't read site %d proc %d log likelihood from  %s.\n", spacer, site, i, line);
	    exit(1);
	  }
	i++;
      }
      else
	multiplicities[site]=atoi(token);
      token=strtok(NULL,delim);
    }
    if(i!=nproc){
      MrBayesPrint ("%s   Error: expected %d site log likelihoods, found %d proc for site %s from  %s.\n", spacer, nproc, i, site, line);
      exit(1);
    }
  }
  return site;
}


/*read site partial likelihoods from file*/
void ReadSites(FILE *partial_fv, MrBFlt **sitelikes, int branch, int nstate, int nsite, MrBFlt *****partials, int nproc, int *multiplicities, int ****scalefact)
{
  int site, sitesfound=0;
  char line[LINELEN];
   
  while(fgets(line,sizeof(line),partial_fv)!=NULL){/*look for site*/
    if(strstr(line,"//")!=0) break;/*end of block*/
    if((site=ReadSiteLogLikes(line,sitelikes,nproc,multiplicities))!=-1){
      sitesfound++;
      ReadOneSite(partial_fv,branch,nstate,site,partials,scalefact);
    }
  }
  if(sitesfound!=nsite){/*wrong number of sites read for this branch*/
    MrBayesPrint ("%s   Error: expected %d sites, found %d.\n", spacer, nsite, sitesfound);
    exit(1);
  }
}

/*hard code some test data for m=2 states, n=3 sites*/
/*use with Q=[-0.3 0.3;0.4 -0.4], t=0.9*/
void testEHD(MrBFlt **Fp,MrBFlt **Fc, MrBFlt *L)
{
  Fp[0][0]=0.571428571428571;
  Fp[0][1]=0.571428571428571;
  Fp[0][2]=0.;
  Fp[1][0]=0.;
  Fp[1][1]=0.;
  Fp[1][2]=0.428571428571429;

  Fc[0][0]=0.;
  Fc[0][1]=1.;
  Fc[0][2]=0.;
  Fc[1][0]=1.;
  Fc[1][1]=0.;
  Fc[1][2]=1.;

  L[0]=0.114467314039127;
  L[1]=0.456961257389444;
  L[2]=0.314104114532301;

}

/*multiply partial likelihoods by stationary probabilities. F is a matrix of partial likelihoods (rows are sites, cols are states). pi_i is a vector of stationary probabilities.*/
void WeightByPi(MrBFlt **F,MrBFlt *pi_i, int nsite, int nstate)
{
  int i,j;

  for(i=0;i<nsite;i++)
    for(j=0;j<nstate;j++)
      F[i][j]*=pi_i[j];
}

/*write the results to an output file*/
void WriteResults(FILE *outfile,MrBFlt *****partials, int nbranch, int nproc, int nsite, MrBFlt ***condE, MrBFlt **priorE, MrBFlt **priorV, int *multiplicities, int *sitemap, int ncols, MrBFlt *tbranch)
{
  int i,j,k,m,b;
  MrBFlt *totalb,*totalbp,*totals,*totalsp,*vbp,zij,*zs,z;
  MrBFlt tol=1e-15;

  totalb=AllocateDoubleVector(nbranch);
  totalbp=AllocateDoubleVector(nbranch);
  vbp=AllocateDoubleVector(nbranch);
  totals=AllocateDoubleVector(ncols);
  totalsp=AllocateDoubleVector(ncols);
  zs=AllocateDoubleVector(ncols);

  fprintf(outfile,"\n");
  fprintf(outfile,"Branch\tProcess\tSite\tConditional_expectation\tPrior_expectation\tPrior_variance\tzscore\n");
  for(i=0;i<ncols;i++){
    totals[i]=0.;
    totalsp[i]=0.;
  }
  b=0;
  for(i=0;i<nbranch;i++){
    totalb[i]=0.;
    vbp[i]=0.;
    if(tbranch[i]>tol) b++;/*number of nonzero edges*/
    for(j=0;j<nproc;j++){/*SUMMING OVER PROCESSES WON'T BE RIGHT: WE DON'T YET HAVE PRIOR PROBABILITIES SUPPLIED*/
      for(k=0;k<ncols;k++){
	m=sitemap[k];/*we need to look up the pattern corresponding to this alignment column*/
	zij=computeZscore(condE[i][j][m],priorE[i][j],priorV[i][j]);
	fprintf(outfile,"%d\t%d\t%d\t%.16e\t%.16e\t%.16e\t%.16e\n",i,j,k,condE[i][j][m],priorE[i][j],priorV[i][j],zij);
	if(tbranch[i]>tol) zs[k]+=zij;/*add to site score, only if branch has length and therefore prior variance greater than zero*/
	totalb[i]+=condE[i][j][m];/*branch conditional expectatation*/
	totals[k]+=condE[i][j][m];/*site conditional expectation*/
	totalsp[k]+=priorE[i][j];/*site prior expectation*/
      }
      totalbp[i]=priorE[i][j]*ncols;/*branch prior expectation*/
      vbp[i]=priorV[i][j]*ncols;/*branch prior variance*/
    }
  }
  fprintf(outfile,"\nBranch\tTotal_conditional_expectation\tTotal_prior_expectation\tzscore\tlength\n");
  for(i=0;i<nbranch;i++){
    fprintf(outfile,"%d\t%.16e\t%.16e\t%.16e\t%.16e\n",i,totalb[i],totalbp[i],computeZscore(totalb[i],totalbp[i],vbp[i]),tbranch[i]);
  }
  
  fprintf(outfile,"\nSite\tTotal_conditional_expectation\tTotal_prior_expectation\tzscore\n");
  z=0.;
  for(i=0;i<ncols;i++){
    fprintf(outfile,"%d\t%.16e\t%.16e\t%.16e\n",i,totals[i],totalsp[i],zs[i]/sqrt(b));
    z+=zs[i]/sqrt(b);
  }
  fprintf(outfile,"\nOverall_zscore\t%.16e\n",z/sqrt(ncols));
  free(totalb);
  free(totalbp);
  free(vbp);
  free(totals);
  free(totalsp);
  free(zs);
}

void CalculateAndWrite(int nsite, int nstate, int nbranch, int nproc, int ncols, int ****scalefact, int **L, int *multiplicities, int *sitemap,
        MrBFlt *****partials,
        MrBFlt ***Qset, 
        MrBFlt **sitelikes, MrBFlt **pi_i,
        MrBFlt *tbranch, MrBFlt *mixprobs,
        FILE *outfile
                ){
    MrBFlt **ENLt, **ENLtD, **Pt,t; 
    MrBFlt ***condE;
    MrBFlt **EigenValues,  ***EigVecs, ***inverseEigVecs;
    MrBFlt ***QLset;
    MrBFlt **priorE, **priorV;
    int i,j;
    priorE=AllocateDoubleMatrix(nbranch,nproc);/*prior expectation for each branch and process*/
    priorV=AllocateDoubleMatrix(nbranch,nproc);/*prior variance for each branch and process*/

    EigVecs=AllocateQset(nproc,nstate);
    inverseEigVecs=AllocateQset(nproc,nstate);
    
    /*allocate space for calculations*/
    ENLt=AllocateSquareDoubleMatrix(nstate);
    ENLtD=AllocateSquareDoubleMatrix(nstate);
    Pt=AllocateSquareDoubleMatrix(nstate);
    condE=AllocatecondE(nbranch,nproc,nsite);/*allocate memory for conditional expectations*/
    EigenValues=AllocateDoubleMatrix(nproc,nstate);/*eigenvalues for each process in rows*/
    QLset=AllocateQset(nproc,nstate);

    for(i=0;i<nproc;i++){/*eigendecomposition, stationary distribution, and QL for each process*/
      EigenDecomp(nstate,Qset[i],EigenValues[i],EigVecs[i],inverseEigVecs[i]);/* eigendecomposition of rate matrix*/
      makeQL(Qset[i],L,QLset[i],nstate);/*construct QL matrix*/
      /*computepi(nstate, inverseEigVecs[i], EigenValues[i], pi_i[i]);*//*NOT NEEDED ANY MORE: SHOULD BE IN FILE stationary distributions*/
    }
 
    /*do the computations for each branch, process, and site*/
    for(i=0;i<nbranch;i++){
      t=tbranch[i];/*get branch length*/
      for(j=0;j<nproc;j++){
	computeENLt(nstate, QLset[j], t, ENLt, EigVecs[j], inverseEigVecs[j], EigenValues[j]);/*joint mean expectation*/
	computePt(nstate, t, EigVecs[j], inverseEigVecs[j], EigenValues[j], Pt);/*transition probabilities*/
 	computeENLtD(nstate,ENLt,Pt,ENLtD);/*conditional mean matrix: we don't need this explicitly. See eq. 3.8 in Minin and Suchard*/
	priorE[i][j]=computepriorE(pi_i[j],QLset[j],t,nstate);/*prior mean*/
	/*priorV[i][j]=computepriorV(priorE[i][j],pi_i[j],QLset[j],t,nstate,EigVecs[j],inverseEigVecs[j],EigenValues[j]);*/ /*ORIGINAL VERSION: THIS IS THE WRONG VARIANCE (VARIANCE OF NUMBER OF EVENTS, NOT EXPECTED NUMBER OF EVENTS*/
	WeightByPi(partials[i][j][0],pi_i[j],nsite,nstate);/*stationary probabilities need to go in on one end or the other (doesn't matter which provided the substitution model is reversible and the counting set is symmetrical)*/
	/*computeEHD(partials[i][j][0], partials[i][j][1], Pt, sitelikes[j], condE[i][j], nstate, nsite,scalefact[i][j]);*//*TEST: should get a vector of 1s. Real code is below*/
	computeEHD(partials[i][j][0], partials[i][j][1], ENLt, sitelikes[j], condE[i][j], nstate, nsite,scalefact[i][j]);/*conditional expected numbers of labelled changes*/
	priorV[i][j]=computepriorVE(ENLtD,priorE[i][j],pi_i[j],Pt,nstate);/*prior variance of MEAN number of events: ADDED 25/6/10*/
      }
    }
    WriteResults(outfile,partials,nbranch,nproc,nsite,condE,priorE,priorV,multiplicities,sitemap,ncols,tbranch);
    FreeSquareDoubleMatrix(ENLt);
    FreeSquareDoubleMatrix(ENLtD);
    FreeSquareDoubleMatrix(Pt);
    FreecondE(condE,nbranch);
    FreeQset(EigVecs,nproc);
    FreeQset(inverseEigVecs,nproc);
    FreeQset(QLset,nproc);
    FreeDoubleMatrix(EigenValues);
    FreeDoubleMatrix(priorE);
    FreeDoubleMatrix(priorV);
 
}


#ifndef BUILDLIBRARY
int main(int argc, char * argv[])
{
  int nsite,nstate,nbranch,nproc,ncols;
  int ****scalefact;
  int **L;
  int *multiplicities, *sitemap;
  int i;
  MrBFlt *****partials;
  MrBFlt ***Qset;
  MrBFlt **sitelikes, **pi_i;
  MrBFlt *tbranch, *mixprobs;
  settings      *sets;
  char          *infilename, *outfilename, *lfilename;
  FILE          *outfile;

  sets=(settings *)calloc(1,sizeof(settings));
  infilename=(char *)calloc(LEN,sizeof(char));
  outfilename=(char *)calloc(LEN,sizeof(char));
  lfilename=(char *)calloc(LEN,sizeof(char));
   
  ReadParams(argc,argv,sets,infilename,outfilename,lfilename);/*read parameters*/

  if(sets->printHelp){/*display help and exit*/
    printHelp();
    free(infilename);
    free(outfilename);
    free(lfilename);
    free(sets);
  }
  else{

    /*open output file*/
    outfile=fopen(outfilename,"w");
    if(outfile==NULL){
      MrBayesPrint ("%s   Error: Problem writing output file %s.\n", spacer,outfilename);
      exit(1);
    }
    ncols=CountCols(infilename);/*number of columns in original alignment*/
    nstate=CountStates(infilename);/*number of states in Q matrix*/
    nbranch=CountBranches(infilename);/*number of branches in tree*/
    nproc=CountProc(infilename);/*number of processes in model*/
    nsite=CountSites(infilename);/*number of sites in file*/
    printf("%d states in partial file\n",nstate);
    printf("%d branches in partial file\n",nbranch);
    printf("%d processes in partial file\n",nproc);
    printf("%d sites read in partial file\n",nsite);
    printf("%d sites in original alignment\n",ncols);

    Qset=AllocateQset(nproc,nstate);/*allocate memory for set of Q matrices*/
    pi_i=AllocateDoubleMatrix(nproc,nstate);/*stationary probabilities for each process in rows*/
    partials=AllocatePartials(nbranch,nproc,nsite,nstate);/*allocate memory for partial likelihoods*/
    sitelikes=AllocateDoubleMatrix(nproc,nsite);/*allocate memory for site likelihoods*/
    tbranch=AllocateDoubleVector(nbranch);/*allocate memory for branch lengths*/
    scalefact=Allocatescalefact(nbranch,nproc,nsite);/*allocate memory for scale factors*/
    mixprobs=AllocateDoubleVector(nproc);/*allocate memory for mixing probabilities*/
    multiplicities=AllocateIntegerVector(nsite);/*allocate memory for site multiplicities*/
    sitemap=AllocateIntegerVector(ncols);/*allocate memory for site map*/
    L=AllocateSquareIntegerMatrix(nstate);/*matrix of labelled transitions*/
    ReadLMatrix(lfilename,L,nstate);/*1 for transitions we want to count, 0 for others*/
    ReadPartials(infilename,Qset,partials,sitelikes,tbranch,mixprobs,nstate,nproc,nsite,multiplicities,scalefact,sitemap,ncols,pi_i);/*read the partial likelihood file*/
  
    fprintf(outfile,"%s version %s\n%s\nCommand:\n",PROGRAM_NAME,PROGRAM_VERSION,COPYRIGHT);
    for(i=0;i<argc;i++) fprintf(outfile,"%s%s",argv[i],(i<argc-1) ? " " : "");/*write arguments to output file*/
    CalculateAndWrite( nsite,  nstate,  nbranch,  nproc,  ncols,
         scalefact, L,  multiplicities, sitemap,
         partials,
         Qset, 
         sitelikes, pi_i,  
         tbranch,  mixprobs, outfile);
 


    /* free memory */
    FreeQset(Qset,nproc);
    FreeDoubleMatrix(pi_i);
    FreePartials(partials,nbranch,nproc);
    FreeDoubleMatrix(sitelikes);
    free(tbranch);
    Freescalefact(scalefact,nbranch,nproc,nsite);
    free(mixprobs);
    free(multiplicities);
    free(sitemap);
   FreeSquareIntegerMatrix(L);
    free(infilename);
    free(outfilename);
    free(lfilename);
    free(sets);
    fclose(outfile);
  }
  return (0);
}
#endif
