 /* Serre2.i */
 %module Serre2
 %{
 /* Put header files here or function declarations like below */
 extern double *mallocPy(int n);
 extern void writetomem(double *x, int i , double f);
 extern double readfrommem(double *x,int i);
 extern void deallocPy(double *x);
 extern double minmod(double a, double b, double c);
 extern double GNall(double *x,double *h,double *u, double *b,double g,int n, int nBC,double dx);
 extern void conc(double *a , double *b, double *c,int n,int m ,int k, double *d);
 extern void getufromG(double *h, double *G,double *bed, double u0, double u1, double h0, double h1, double b0, double b1,double dx, int n,double *ublank);
 extern void getGfromu(double *h, double *u,double *bed, double u0, double u1, double h0, double h1, double b0, double b1,double dx,int n, double *Gblank );
 extern void evolvewrap(double *G,double *h,double *bed,double *h0,double *h1,double *u0,double *u1,double *b0,double *b1, double g,double dx,double dt, int n, int nBC,double theta);
 %}
 extern double *mallocPy(int n);
 extern void writetomem(double *x, int i , double f);
 extern double readfrommem(double *x,int i);
 extern void deallocPy(double *x);
 extern double minmod(double a, double b, double c);
 extern double GNall(double *x,double *h,double *u, double *b,double g,int n, int nBC,double dx);
 extern void conc(double *a , double *b, double *c,int n,int m ,int k, double *d);
 extern void getufromG(double *h, double *G,double *bed, double u0, double u1, double h0, double h1, double b0, double b1,double dx, int n,double *ublank);
 extern void getGfromu(double *h, double *u,double *bed, double u0, double u1, double h0, double h1, double b0, double b1,double dx,int n, double *Gblank );
 extern void evolvewrap(double *G,double *h,double *bed,double *h0,double *h1,double *u0,double *u1,double *b0,double *b1, double g,double dx,double dt, int n, int nBC,double theta);
