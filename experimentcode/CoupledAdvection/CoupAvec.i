 /* CoupAvec.i */
 %module CoupAvec
 %{
 /* Put header files here or function declarations like below */
 extern double *mallocPy(int n);
 extern void conc(double *a , double *b, double *c,int n,int m ,int k, double *d);
 extern void writetomem(double *x, int i , double f);
 extern double readfrommem(double *x,int i);
 extern void deallocPy(double *x);
 extern double minmod(double a, double b, double c);
 extern void evolvewrap(double *u, double *v, double *u0, double *u1, double *v0, double *v1, double a, double b, double dx, double dt, int nBC, int n, int nBCs, double theta);
 extern double HankEnergyall(double *x,double *u,int n, int nBC,double dx);
 extern void evolvewrapFD(double *u, double *v, double *pubc, double *pvbc, double *u0, double *u1, double *v0, double *v1, double a, double b, double dx, double dt, int nBC, int n, int nBCs);
 %}
 extern double *mallocPy(int n);
 extern void conc(double *a , double *b, double *c,int n,int m ,int k, double *d);
 extern void writetomem(double *x, int i , double f);
 extern double readfrommem(double *x,int i);
 extern void deallocPy(double *x);
 extern double minmod(double a, double b, double c);
 extern void evolvewrap(double *u, double *v, double *u0, double *u1, double *v0, double *v1, double a, double b, double dx, double dt, int nBC, int n, int nBCs, double theta);
 extern double HankEnergyall(double *x,double *u,int n, int nBC,double dx);
 extern void evolvewrapFD(double *u, double *v, double *pubc, double *pvbc, double *u0, double *u1, double *v0, double *v1, double a, double b, double dx, double dt, int nBC, int n, int nBCs);
