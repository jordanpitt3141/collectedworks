 /* Serre2GR.i */
 %module Serre2GR
 %{
 /* Put header files here or function declarations like below */
 extern double *mallocPy(int n);
 extern void writetomem(double *x, int i , double f);
 extern double readfrommem(double *x,int i);
 extern void deallocPy(double *x);
 extern void evolvewrap(double *u, double *h, double *pubc, double *h0, double *h1, double *u0, double *u1, double g, double dx, double dt, int nBC, int n, int nBCs);
 extern double HankEnergyall(double *x,double *h,double *u,double g,int n, int nBC,double dx); 
 extern void conc(double *a , double *b, double *c,int n,int m ,int k, double *d);
 %}
 extern double *mallocPy(int n);
 extern void writetomem(double *x, int i , double f);
 extern double readfrommem(double *x,int i);
 extern void deallocPy(double *x);
 extern void evolvewrap(double *u, double *h, double *pubc, double *h0, double *h1, double *u0, double *u1, double g, double dx, double dt, int nBC, int n, int nBCs);
 extern double HankEnergyall(double *x,double *h,double *u,double g,int n, int nBC,double dx); 
 extern void conc(double *a , double *b, double *c,int n,int m ,int k, double *d);
