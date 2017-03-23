 /* Hamil.i */
 %module Hamil
 %{
 /* Put header files here or function declarations like below */
 extern double *mallocPy(int n);
 extern void writetomem(double *x, int i , double f);
 extern double readfrommem(double *x,int i);
 extern void deallocPy(double *x);
 extern double HankEnergyall(double *x,double *h,double *u,double g,int n, int nBC,double dx);
 extern double uhall(double *x,double *h,double *u,int n, int nBC,double dx);
 extern double hall(double *x,double *h,int n, int nBC,double dx);
 %} 
 extern double *mallocPy(int n);
 extern void writetomem(double *x, int i , double f);
 extern double readfrommem(double *x,int i);
 extern void deallocPy(double *x);
 extern double HankEnergyall(double *x,double *h,double *u,double g,int n, int nBC,double dx); 
 extern double uhall(double *x,double *h,double *u,int n, int nBC,double dx);
 extern double hall(double *x,double *h,int n, int nBC,double dx);
