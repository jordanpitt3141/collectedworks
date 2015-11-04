 /* Serre2.i */
 %module Serre2
 %{
 /* Put header files here or function declarations like below */
 extern double *mallocPy(int n);
 extern void writetomem(double *x, int i , double f);
 extern double readfrommem(double *x,int i);
 extern void deallocPy(double *x);
 extern void evolvewrap(double *G, double *h, double *bed, double *h0, double *h1, double *u0, double *u1, double *b0, double *b1, double g, double dx, double dt, int nBC, int n, int nBCs, double theta);
 extern void getufromG(double *h, double *G, double *bed, double u0, double u1, double h0, double h1, double b0, double b1, double dx , int n, double *u);
 extern double minmod(double a, double b, double c);
 %}
 extern double *mallocPy(int n);
 extern void writetomem(double *x, int i , double f);
 extern double readfrommem(double *x,int i);
 extern void deallocPy(double *x);
 extern void evolvewrap(double *G, double *h, double *bed, double *h0, double *h1, double *u0, double *u1, double *b0, double *b1, double g, double dx, double dt, int nBC, int n, int nBCs, double theta);
 extern void getufromG(double *h, double *G, double *bed, double u0, double u1, double h0, double h1, double b0, double b1, double dx , int n, double *u);
 extern double minmod(double a, double b, double c);
