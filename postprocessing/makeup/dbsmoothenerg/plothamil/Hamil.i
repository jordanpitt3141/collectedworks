 /* Hamil.i */
 %module Hamil
 %{
 /* Put header files here or function declarations like below */
 extern double *mallocPy(int n);
 extern void writetomem(double *x, int i , double f);
 extern double readfrommem(double *x,int i);
 extern void deallocPy(double *x);
 extern void HankEnergyall(double *x,double *h,double *u,double g,int n, int nBC,double dx, double *Ham);
 extern void HankEnergyallPT(double *x,double *h,double *u,double g,int n, int nBC,double dx, double *HamFT, double *HamST, double *HamTT);
 %} 
 extern double *mallocPy(int n);
 extern void writetomem(double *x, int i , double f);
 extern double readfrommem(double *x,int i);
 extern void deallocPy(double *x);
 extern void HankEnergyall(double *x,double *h,double *u,double g,int n, int nBC,double dx, double *Ham);
 extern void HankEnergyallPT(double *x,double *h,double *u,double g,int n, int nBC,double dx, double *HamFT, double *HamST, double *HamTT);
