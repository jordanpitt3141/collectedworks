 /* Serre2dc.i */
 %module Serre2dc
 %{
 /* Put header files here or function declarations like below */
 extern double *mallocPy(int n);
 extern void conc(double *a , double *b, double *c,int n,int m ,int k, double *d);
 extern void writetomem(double *x, int i , double f);
 extern double readfrommem(double *x,int i);
 extern void deallocPy(double *x);
 extern void evolvewrap(double *G, double *h, double *h0, double *h1, double *u0, double *u1, double g, double dx, double dt, int nBC, int n, int nBCs, double theta);
 extern void getufromG(double *h, double *G, double u0, double u1, double h0, double h1, double dx , int n, double *u);
 extern double minmod(double a, double b, double c);
 extern double HankEnergyall(double *x,double *h,double *u,double g,int n, int nBC,double dx); 
 %}
 extern double *mallocPy(int n);
 extern void conc(double *a , double *b, double *c,int n,int m ,int k, double *d);
 extern void writetomem(double *x, int i , double f);
 extern double readfrommem(double *x,int i);
 extern void deallocPy(double *x);
 extern void evolvewrap(double *G, double *h, double *h0, double *h1, double *u0, double *u1, double g, double dx, double dt, int nBC, int n, int nBCs, double theta);
 extern void getufromG(double *h, double *G, double u0, double u1, double h0, double h1, double dx , int n, double *u);
 extern double minmod(double a, double b, double c);
 extern double HankEnergyall(double *x,double *h,double *u,double g,int n, int nBC,double dx); 
