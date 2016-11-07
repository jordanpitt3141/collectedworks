 /* Int.i */
 %module Int
 %{
 /* Put header files here or function declarations like below */
 extern double *mallocPy(int n);
 extern void writetomem(double *x, int i , double f);
 extern double readfrommem(double *x,int i);
 extern void deallocPy(double *x);
 extern double Intall(double *x,double *h,int n, int nBC,double dx);
 %} 
 extern double *mallocPy(int n);
 extern void writetomem(double *x, int i , double f);
 extern double readfrommem(double *x,int i);
 extern void deallocPy(double *x);
 extern double Intall(double *x,double *h,int n, int nBC,double dx);
