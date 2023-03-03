#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define k2       (0.01720209895*0.01720209895)
#define MS       1.988705522e30
#define mJ       1.8988e27
#define mZ       5.9742E24
#define MAXN     20
#define MAXn     20
#define year     365.242190419

struct params_str
{
 int dim;
 double *X, *Y, *Z, *U, *V, *W, *R, *S, *D, *P;
 int npl, order;
 double *m, *mu, *beta, m0;
};

struct params_str params;

#include "functions.h"

int main()
{
 char nazwa[100];
 FILE *plik1, *plik2, *plik3;
 int i, N, n, lp, nstep, Nstep, j, ii, kierunek;
 double cart[6*MAXN], kepl_a[6*MAXN], kepl_P[6*MAXN]; 
 double E0, E, dt, Amin, t, h, tend;
 int jowisze, stopnie, liczymy;
 double maxDE = 0.0, DE;

 sprintf(nazwa, "start_taylor_Nbody.in");
 plik1 = fopen(nazwa, "r");
 ii=fscanf(plik1, "%d %d %le %d %le %d %d %d", &lp, &N, &tend, &Nstep, &dt, &n, &jowisze, &stopnie);
 tend = tend*year;
 params.npl = N;
 params.order = n;

 params.m = malloc(params.npl*sizeof(double));
 params.mu = malloc(params.npl*sizeof(double));
 params.beta = malloc(params.npl*sizeof(double));
 
 ii=fscanf(plik1, "%le", &params.m0);
 for (i=0; i<N; i++) {
     ii=fscanf(plik1, "%le", &params.m[i]);
     for (j=0; j<6; j++)  fscanf(plik1, "%le", &kepl_a[6*i+j]);
     if (jowisze == 1)      params.m[i] *= mJ/MS;
     else if (jowisze == 0) params.m[i] *= mZ/MS;
     params.mu[i] = k2*(params.m0 + params.m[i]);
     params.beta[i] = 1./(1./params.m0 + 1./params.m[i]);
     if (stopnie == 1) {
         kepl_a[6*i+2] *= M_PI/180.0;
         kepl_a[6*i+3] *= M_PI/180.0;
         kepl_a[6*i+4] *= M_PI/180.0;
         kepl_a[6*i+5] *= M_PI/180.0;
     }
 }
 fclose(plik1);

 params.X = malloc(params.npl*params.order*sizeof(double));
 params.Y = malloc(params.npl*params.order*sizeof(double));
 params.Z = malloc(params.npl*params.order*sizeof(double));
 params.U = malloc(params.npl*params.order*sizeof(double));
 params.V = malloc(params.npl*params.order*sizeof(double));
 params.W = malloc(params.npl*params.order*sizeof(double));
 params.R = malloc(params.npl*params.order*sizeof(double));
 params.S = malloc(params.npl*params.order*sizeof(double));
 params.D = malloc(params.npl*params.npl*params.order*sizeof(double));
 params.P = malloc(params.npl*params.npl*params.order*sizeof(double));
 
 kepl2cart(kepl_a, cart, params);

 t = 0.0;
 ii = 0;
 
 E0 = energy(cart, params);
 E = E0;
 
 sprintf(nazwa, "results/kepl_a_%i.out", lp);
 plik1 = fopen(nazwa, "w");
 sprintf(nazwa, "results/kepl_P_%i.out", lp);
 plik2 = fopen(nazwa, "w");
 sprintf(nazwa, "results/cart_%i.out", lp);
 plik3 = fopen(nazwa, "w");

 Amin = kepl_a[0];
 for (i=1; i<N; i++) 
     if ((kepl_a[6*i] < Amin) && (kepl_a[6*i] > 0.0))  Amin = kepl_a[6*i];
 h = 2.*M_PI*sqrt(Amin*Amin*Amin/(k2*params.m0))*dt;

 fprintf(plik1, "%le %le ", t/year, (E-E0)/E0);
 for (i=0; i<6*N; i++)  fprintf(plik1, "%.12le ", kepl_a[i]);
 fprintf(plik1, "\n");

 cart2kepl_P (cart, kepl_P, params);
 fprintf(plik2, "%le %le ", t/year, (E-E0)/E0);
 for (i=0; i<6*N; i++)  fprintf(plik2, "%.12le ", kepl_P[i]);
 fprintf(plik2, "\n");

 fprintf(plik3, "%le %le ", t/year, (E-E0)/E0);
 for (i=0; i<6*N; i++)  fprintf(plik3, "%.12le ", cart[i]);
 fprintf(plik3, "\n");

 kierunek = 1;
 if (tend < 0.0) {
     kierunek = -1;
     h = -h;
 }
 
 liczymy = 1;

 clock_t begin = clock();

 while ((kierunek*t < kierunek*tend) && (liczymy == 1)) {
     if (kierunek*(t+h) > kierunek*tend) {
         h = tend - t;
         liczymy = 0;
     }
     taylor3d(cart, &t, h, &params);

     ii++;
     if (ii == Nstep) {
         cart2kepl_a(cart, kepl_a, params);
         cart2kepl_P(cart, kepl_P, params);
         E = energy(cart, params);
	 DE = fabs((E-E0)/E0);
	 if (maxDE < DE)  maxDE = DE;

         fprintf(plik1, "%le %le ", t/year, (E-E0)/E0);
         for (i=0; i<6*N; i++)  fprintf(plik1, "%.12le ", kepl_a[i]);
         fprintf(plik1, "\n");

         fprintf(plik2, "%le %le ", t/year, (E-E0)/E0);
         for (i=0; i<6*N; i++)  fprintf(plik2, "%.12le ", kepl_P[i]);
         fprintf(plik2, "\n");

         fprintf(plik3, "%le %le ", t/year, (E-E0)/E0);
         for (i=0; i<6*N; i++)  fprintf(plik3, "%.12le ", cart[i]);
         fprintf(plik3, "\n");

	 ii=0;
     }
 }
 
 clock_t end = clock();
 printf("t = %lf s\n", (double)(end-begin)/CLOCKS_PER_SEC);
 
 fclose(plik1);
 fclose(plik2);
 fclose(plik3);

 free(params.X);
 free(params.Y);
 free(params.Z);
 free(params.U);
 free(params.V);
 free(params.W);
 free(params.R);
 free(params.S);
 free(params.D);
 free(params.P);
 free(params.m);
 free(params.mu);
 free(params.beta);

 return 0;
}
