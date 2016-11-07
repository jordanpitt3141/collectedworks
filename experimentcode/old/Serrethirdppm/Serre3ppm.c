#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>



const double i24 = 1.0/24.0;
const double i12 = 1.0/12.0;
const double i3 = 1.0/3.0;

void conc(double *a , double *b, double *c,int n,int m ,int k, double *d)
{
    //replace with memcpy for performance?
    memcpy(d,a,n*sizeof(double));
    memcpy(d+n,b,m*sizeof(double));
    memcpy(d+n+m,c,k*sizeof(double));
}

void TDMA(double *a, double *b, double *c, double *d, int n, double *x)
{
	double *alpha = malloc(n*sizeof(double));
	double *beta = malloc(n*sizeof(double));

    alpha[0] = c[0] / b[0];
    beta[0] = d[0] / b[0];


    int i;
    double m;
    for(i=1; i<n-1; i++)
    {
        m = 1.0 / (b[i] - a[i-1]*alpha[i-1]);
        alpha[i] = c[i]*m;
        beta[i] = (d[i] - a[i-1]*beta[i-1])*m;   
        
    }

    m = 1.0 / (b[n-1] - a[n-2]*alpha[n-2]);
    beta[n-1] = (d[n-1] - a[n-2]*beta[n-2])*m;

    x[n-1] = beta[n-1];

    for (i=n-2; i > -1; i--)
    {
        x[i] = beta[i] - alpha[i]*x[i+1];  
   
    }

    free(alpha);
    free(beta);
}

void PENT(double *e, double *a, double *d, double *c,double *f, double *B, int n, double *x)
{
    //works but is destructive on inputs
    int i;
    double m;
    for(i=1; i<n-1; i++)
    {
        m = a[i-1] /d[i-1];
        
        d[i] = d[i] - m*c[i-1];
        c[i] = c[i] - m*f[i-1];
        B[i] = B[i] - m*B[i-1];

        m = e[i-1] /d[i-1];
        a[i] = a[i] - m*c[i-1];
        d[i+1] = d[i+1] - m*f[i-1];
        B[i+1] = B[i+1] - m*B[i-1];
    }

    m = a[n-2] / d[n-2];
    d[n-1] = d[n-1] - m*c[n-2];
    x[n-1] = (B[n-1] - m*B[n-2]) / d[n-1];
    x[n-2] = (B[n-2] - c[n-2]*x[n-1]) / d[n-2];

    for (i=n-3; i > -1; i--)
    {
        x[i] = (B[i] - f[i]*x[i+2] - c[i]*x[i+1]) / d[i];  
   
    }
}

void midpt2ca(double *qm , double dx , int n, double *qa)
{
    double *a = malloc((n-1)*sizeof(double));
    double *b = malloc(n*sizeof(double));
    double *c = malloc((n-1)*sizeof(double));

    int i;

    for (i=1;i<n-1;i++)
    {
        a[i-1] = -i24;
        b[i] = 26*i24;
        c[i] = -i24;


    }

    //i =0
    i = 0;
    b[i] = 1.0;
    c[i] = 0.0;

    //i=n-1
    i = n-1;
    a[i-1] = 0.0;
    b[i] = 1.0;    

    TDMA(a,b,c,qm,n,qa);
    free(a);
    free(b);
    free(c);
}

void ca2midpt(double *qa, double dx, int n,double *qm)
{
    //double *qm = malloc(n*sizeof(double));
    int i;

    for (i=1;i<n-1;i++)
    {
        qm[i] = i24*(-qa[i+1] + 26*qa[i] - qa[i-1]);
    }

    i=0;
    qm[i] = qa[i];


    i = n-1;
    qm[i] = qa[i];

}

void ufromGh(double *G, double *h,double *hbeg,double *hend,double *ubeg,double *uend, double dx , int n, int nBC, double *u)
{
    double idx = 1.0 / dx;
    double *a = malloc((n-2)*sizeof(double));
    double *b = malloc((n-1)*sizeof(double));
    double *c = malloc(n*sizeof(double));
    double *d = malloc((n-1)*sizeof(double));
    double *e = malloc((n-2)*sizeof(double));
    double *f = malloc(n*sizeof(double));


    int i,j;
    double thx,tmp1,tmp2;
    for(i =2; i < n-2 ; i++)
    {
        thx = i12*idx*(-h[i+2] + 8*h[i+1] - 8*h[i-1] + h[i-2]);
        a[i-2] = -(i12*idx)*(h[i]*h[i]*thx) + (i12*idx*idx)*(i3*h[i]*h[i]*h[i]);
        b[i-1] = (8*i12*idx)*(h[i]*h[i]*thx) - (16*i12*idx*idx)*(i3*h[i]*h[i]*h[i]);
        c[i] = h[i] + (30*i12*idx*idx)*(i3*h[i]*h[i]*h[i]); 
        d[i] = -(8*i12*idx)*(h[i]*h[i]*thx) - (16*i12*idx*idx)*(i3*h[i]*h[i]*h[i]);
        e[i] = (i12*idx)*(h[i]*h[i]*thx) + (i12*idx*idx)*(i3*h[i]*h[i]*h[i]);
        
        f[i] = G[i];        

    }
    
    //Boundaries

    i = 0;
    j = nBC;
    thx = i12*idx*(-h[i+2] + 8*h[i+1] - 8*hbeg[j-1] + hbeg[j-2]);
    tmp1 = -(i12*idx)*(h[i]*h[i]*thx) + (i12*idx*idx)*(i3*h[i]*h[i]*h[i]);
    tmp2 = (8*i12*idx)*(h[i]*h[i]*thx) - (16*i12*idx*idx)*(i3*h[i]*h[i]*h[i]);
    c[i] = h[i] + (30*i12*idx*idx)*(i3*h[i]*h[i]*h[i]); 
    d[i] = -(8*i12*idx)*(h[i]*h[i]*thx) - (16*i12*idx*idx)*(i3*h[i]*h[i]*h[i]);
    e[i] = (i12*idx)*(h[i]*h[i]*thx) + (i12*idx*idx)*(i3*h[i]*h[i]*h[i]);
        
    f[i] = G[i] - ubeg[j-1]*tmp2 - ubeg[j-2]*tmp1; 
 

    i = 1;
    j = nBC + 1;
    thx = i12*idx*(-h[i+2] + 8*h[i+1] - 8*h[i-1] + hbeg[j-2]);
    tmp1 = -(i12*idx)*(h[i]*h[i]*thx) + (i12*idx*idx)*(i3*h[i]*h[i]*h[i]);
    b[i-1] = (8*i12*idx)*(h[i]*h[i]*thx) - (16*i12*idx*idx)*(i3*h[i]*h[i]*h[i]);
    c[i] = h[i] + (30*i12*idx*idx)*(i3*h[i]*h[i]*h[i]); 
    d[i] = -(8*i12*idx)*(h[i]*h[i]*thx) - (16*i12*idx*idx)*(i3*h[i]*h[i]*h[i]);
    e[i] = (i12*idx)*(h[i]*h[i]*thx) + (i12*idx*idx)*(i3*h[i]*h[i]*h[i]);
     
    f[i] = G[i] - ubeg[j-2]*tmp1; 

    i = n-2;
    j = -2;
    thx = i12*idx*(-hend[j+2] + 8*h[i+1] - 8*h[i-1] + h[i-2]);
    a[i-2] = -(i12*idx)*(h[i]*h[i]*thx) + (i12*idx*idx)*(i3*h[i]*h[i]*h[i]);
    b[i-1] = (8*i12*idx)*(h[i]*h[i]*thx) - (16*i12*idx*idx)*(i3*h[i]*h[i]*h[i]);
    c[i] = h[i] + (30*i12*idx*idx)*(i3*h[i]*h[i]*h[i]); 
    d[i] = -(8*i12*idx)*(h[i]*h[i]*thx) - (16*i12*idx*idx)*(i3*h[i]*h[i]*h[i]);
    tmp1 = (i12*idx)*(h[i]*h[i]*thx) + (i12*idx*idx)*(i3*h[i]*h[i]*h[i]);
        
    f[i] = G[i] -uend[j+2]*tmp1;

    i = n-1;
    j = -1;
    thx = i12*idx*(-hend[j+2] + 8*hend[j+1] - 8*h[i-1] + h[i-2]);
    a[i-2] = -(i12*idx)*(h[i]*h[i]*thx) + (i12*idx*idx)*(i3*h[i]*h[i]*h[i]);
    b[i-1] = (8*i12*idx)*(h[i]*h[i]*thx) - (16*i12*idx*idx)*(i3*h[i]*h[i]*h[i]);
    c[i] = h[i] + (30*i12*idx*idx)*(i3*h[i]*h[i]*h[i]); 
    tmp1 = -(8*i12*idx)*(h[i]*h[i]*thx) - (16*i12*idx*idx)*(i3*h[i]*h[i]*h[i]);
    tmp2 = (i12*idx)*(h[i]*h[i]*thx) + (i12*idx*idx)*(i3*h[i]*h[i]*h[i]);
        
    f[i] = G[i] - uend[j+1]*tmp1 - uend[j+2]*tmp2; 

    PENT(a,b,c,d,e,f,n,u);  

    free(a);
    free(b);
    free(c);
    free(d);
    free(e);
    free(f);
}

void Gfromuh(double *u, double *h,double *hbeg,double *hend,double *ubeg,double *uend, double dx , int n, int nBC, double *G)
{
    double idx = 1.0/dx;

    int i,j;
    double ai,bi,ci,di,ei,thx;

    for(i =2;i<n-2;i++)
    {
        thx = i12*idx*(-h[i+2] + 8*h[i+1] - 8*h[i-1] + h[i-2]);
        ai = -(i12*idx)*(h[i]*h[i]*thx) + (i12*idx*idx)*(i3*h[i]*h[i]*h[i]);
        bi = (8*i12*idx)*(h[i]*h[i]*thx) - (16*i12*idx*idx)*(i3*h[i]*h[i]*h[i]);
        ci = h[i] + (30*i12*idx*idx)*(i3*h[i]*h[i]*h[i]); 
        di = -(8*i12*idx)*(h[i]*h[i]*thx) - (16*i12*idx*idx)*(i3*h[i]*h[i]*h[i]);
        ei = (i12*idx)*(h[i]*h[i]*thx) + (i12*idx*idx)*(i3*h[i]*h[i]*h[i]);

        G[i] = ai*u[i-2] + bi*u[i-1] + ci*u[i] + di*u[i+1] + ei*u[i+2];

    }

    i = 0;
    j = nBC;
    thx = i12*idx*(-h[i+2] + 8*h[i+1] - 8*hbeg[j-1] + hbeg[j-2]);
    ai = -(i12*idx)*(h[i]*h[i]*thx) + (i12*idx*idx)*(i3*h[i]*h[i]*h[i]);
    bi = (8*i12*idx)*(h[i]*h[i]*thx) - (16*i12*idx*idx)*(i3*h[i]*h[i]*h[i]);
    ci = h[i] + (30*i12*idx*idx)*(i3*h[i]*h[i]*h[i]); 
    di = -(8*i12*idx)*(h[i]*h[i]*thx) - (16*i12*idx*idx)*(i3*h[i]*h[i]*h[i]);
    ei = (i12*idx)*(h[i]*h[i]*thx) + (i12*idx*idx)*(i3*h[i]*h[i]*h[i]);

    G[i] = ai*ubeg[j-2] + bi*ubeg[j-1] + ci*u[i] + di*u[i+1] + ei*u[i+2];

    i = 1;
    j = nBC+1;
    thx = i12*idx*(-h[i+2] + 8*h[i+1] - 8*h[i-1] + hbeg[j-2]);
    ai = -(i12*idx)*(h[i]*h[i]*thx) + (i12*idx*idx)*(i3*h[i]*h[i]*h[i]);
    bi = (8*i12*idx)*(h[i]*h[i]*thx) - (16*i12*idx*idx)*(i3*h[i]*h[i]*h[i]);
    ci = h[i] + (30*i12*idx*idx)*(i3*h[i]*h[i]*h[i]); 
    di = -(8*i12*idx)*(h[i]*h[i]*thx) - (16*i12*idx*idx)*(i3*h[i]*h[i]*h[i]);
    ei = (i12*idx)*(h[i]*h[i]*thx) + (i12*idx*idx)*(i3*h[i]*h[i]*h[i]);

    G[i] = ai*ubeg[j-2] + bi*u[i-1] + ci*u[i] + di*u[i+1] + ei*u[i+2];  

    i=n-2;
    j = -2;
    thx = i12*idx*(-hend[j+2] + 8*h[i+1] - 8*h[i-1] + h[i-2]);
    ai = -(i12*idx)*(h[i]*h[i]*thx) + (i12*idx*idx)*(i3*h[i]*h[i]*h[i]);
    bi = (8*i12*idx)*(h[i]*h[i]*thx) - (16*i12*idx*idx)*(i3*h[i]*h[i]*h[i]);
    ci = h[i] + (30*i12*idx*idx)*(i3*h[i]*h[i]*h[i]); 
    di = -(8*i12*idx)*(h[i]*h[i]*thx) - (16*i12*idx*idx)*(i3*h[i]*h[i]*h[i]);
    ei = (i12*idx)*(h[i]*h[i]*thx) + (i12*idx*idx)*(i3*h[i]*h[i]*h[i]);

    G[i] = ai*u[i-2] + bi*u[i-1] + ci*u[i] + di*u[i+1] + ei*uend[j+2];  


    i = n-1;
    j = -1;
    thx = i12*idx*(-hend[j+2] + 8*hend[j+1] - 8*h[i-1] + h[i-2]);
    ai = -(i12*idx)*(h[i]*h[i]*thx) + (i12*idx*idx)*(i3*h[i]*h[i]*h[i]);
    bi = (8*i12*idx)*(h[i]*h[i]*thx) - (16*i12*idx*idx)*(i3*h[i]*h[i]*h[i]);
    ci = h[i] + (30*i12*idx*idx)*(i3*h[i]*h[i]*h[i]); 
    di = -(8*i12*idx)*(h[i]*h[i]*thx) - (16*i12*idx*idx)*(i3*h[i]*h[i]*h[i]);
    ei = (i12*idx)*(h[i]*h[i]*thx) + (i12*idx*idx)*(i3*h[i]*h[i]*h[i]);

    G[i] = ai*u[i-2] + bi*u[i-1] + ci*u[i] + di*uend[j+1] + ei*uend[j+2];

}

void reconstructppm(double *u, int i, double dx, double* uilr, double* uirr)
{
//it works yay!
    double uil, uir;
    double daip1 = 0.5*(u[i+2] - u[i]);
    double dip1 = ((u[i+2] - u[i+1])*(u[i+1] - u[i]) >0);
    dip1 = dip1*copysign(fmin(fabs(daip1),fmin(2*fabs(u[i+2] - u[i+1]), 2*fabs(u[i+1] - u[i]))) , daip1);

    double dai = 0.5*(u[i+1] - u[i-1]);
    double di = ((u[i+1] - u[i])*(u[i] - u[i-1]) >0);
    di = di*copysign(fmin(fabs(dai),fmin(2*fabs(u[i+1] - u[i]), 2*fabs(u[i] - u[i-1]))) , dai);

    uir = u[i] + 0.5*(u[i+1] - u[i]) + (0.5*i3)*(di - dip1);

    double daim1 = 0.5*(u[i] - u[i-2]);
    double dim1 = ((u[i] - u[i-1])*(u[i-1] - u[i-2]) >0);
    dim1 = dim1*copysign(fmin(fabs(daim1),fmin(2*fabs(u[i] - u[i-1]), 2*fabs(u[i-1] - u[i-2]))) , daim1);

    uil = u[i-1] + 0.5*(u[i] - u[i-1]) + (0.5*i3)*(dim1 - di);


    //local extrema
    double lce = (uir - u[i])*(u[i] - uil);

    uir = u[i]*(lce <= 0) + uir*(lce > 0);
    uil = u[i]*(lce <= 0) + uil*(lce > 0);

    //monotonicity
    double toclosellhs = (uir - uil)*(u[i] - 0.5*(uil + uir));
    double tocloselrhs = (uir - uil)*(uir -uil)*0.5*i3; 

    double tocloserlhs = (uir - uil)*(u[i] - 0.5*(uil + uir));
    double tocloserrhs = -(uir - uil)*(uir -uil)*0.5*i3;

    uil = (3*u[i] - 2*uir)*(toclosellhs > tocloselrhs) + uil*(toclosellhs <= tocloselrhs);
    uir = (3*u[i] - 2*uir)*(tocloserlhs < tocloserrhs) + uir*(tocloserlhs >= tocloserrhs);

    *uilr = uil;
    *uirr = uir;

}


void weightsum(double a,double *x, double b, double *y, int n, double *z)
{
    int i;
    for(i =0 ; i < n ;i++)
    {
        z[i] = a*x[i] + b*y[i];
    }

}


void evolve(double *G, double *h, double *u, double g, double dx, double dt,int n, int nBC, double *nh, double *nG)
{
    //Dodgy down at machine precision?
    double idx = 1.0 / dx;

    int i;

    i = nBC-1;


    //off bu 10**-15

    //printf("C , ui2mr : %f, ui1mr : %f",uim2r, uim1r);

    //i +1
    double hip1l=0,hip1r=0;
    reconstructppm(h,i+1,dx,&hip1l,&hip1r);
    double Gip1l=0,Gip1r=0;
    reconstructppm(G,i+1,dx,&Gip1l,&Gip1r);
    double uip1l=0,uip1r=0;
    reconstructppm(u,i+1,dx,&uip1l,&uip1r);

    //i
    double hil=0,hir=0;
    reconstructppm(h,i,dx,&hil,&hir);
    double Gil=0,Gir=0;
    reconstructppm(G,i,dx,&Gil,&Gir);
    double uil=0,uir=0;
    reconstructppm(u,i,dx,&uil,&uir);

    //parabola coefficients
    double uai = 3*idx*idx*(uil + uir -2*u[i]);
    double ubi = idx*(uir - uil);

    double uaip1 = 3*idx*idx*(uip1l + uip1r -2*u[i+1]);
    double ubip1 = idx*(uip1r - uip1l);


    double duel = dx*uai + ubi;
    double duer = -dx*uaip1 + ubip1;

    //printf("duer : %.8f, duel : %.8f\n",duer,duel);
    //printf("C : %f %f %f \n",uir,uim1r,uim2r);

    double sqrtghel = sqrt(g* hir);
    double sqrtgher = sqrt(g* hip1l);

    double sl = fmin(0,fmin(uir - sqrtghel, uip1l - sqrtgher));
    double sr = fmax(0,fmax(uir + sqrtghel, uip1l + sqrtgher));

    double felh = uir*hir;
    double felG = Gir*uir + 0.5*g*hir*hir - 2*i3*hir*hir*hir*duel*duel;
    double ferh = uip1l*hip1l;
    double ferG = Gip1l*uip1l + 0.5*g*hip1l*hip1l -2*i3*hip1l*hip1l*hip1l*duer*duer;

    double isrmsl = 0.0;

    if(sr != sl) isrmsl = 1.0 / (sr - sl);
   
    double foh = isrmsl*(sr*felh - sl*ferh + sl*sr*(hip1l - hir));
    double foG = isrmsl*(sr*felG - sl*ferG + sl*sr*(Gip1l - Gir));

    double fih = foh;
    double fiG = foG;
    hil = hip1l;
    hir = hip1r;
    Gil = Gip1l;
    Gir = Gip1r;
    uil = uip1l;
    uir = uip1r;
    uai = uaip1;
    ubi = ubip1;

    for (i = nBC; i < n + nBC;i++)
    {
        //i +1
        reconstructppm(h,i+1,dx,&hip1l,&hip1r);
        reconstructppm(G,i+1,dx,&Gip1l,&Gip1r);
        reconstructppm(u,i+1,dx,&uip1l,&uip1r);

        uaip1 = 3*idx*idx*(uip1l + uip1r -2*u[i+1]);
        ubip1 = idx*(uip1r - uip1l);


        duel = dx*uai + ubi;
        duer = -dx*uaip1 + ubip1;


        sqrtghel = sqrt(g* hir);
        sqrtgher = sqrt(g* hip1l);

        sl = fmin(0,fmin(uir - sqrtghel, uip1l - sqrtgher));
        sr = fmax(0,fmax(uir + sqrtghel, uip1l + sqrtgher));

        felh = uir*hir;
        felG = Gir*uir + 0.5*g*hir*hir - 2*i3*hir*hir*hir*duel*duel;
        ferh = uip1l*hip1l;
        ferG = Gip1l*uip1l + 0.5*g*hip1l*hip1l -2*i3*hip1l*hip1l*hip1l*duer*duer;

        isrmsl = 0.0;

        if(sr != sl) isrmsl = 1.0 / (sr - sl);
       
        foh = isrmsl*(sr*felh - sl*ferh + sl*sr*(hip1l - hir));
        foG = isrmsl*(sr*felG - sl*ferG + sl*sr*(Gip1l - Gir));

        nh[i -nBC] = h[i] -dt*idx*(foh - fih);
        nG[i -nBC] = G[i] -dt*idx*(foG -fiG);

        fih = foh;
        fiG = foG;
        hil = hip1l;
        hir = hip1r;
        Gil = Gip1l;
        Gir = Gip1r;
        uil = uip1l;
        uir = uip1r;
        uai = uaip1;
        ubi = ubip1;

    }
}      

void evolvewrap(double *Ga, double *ha, double *Gabeg, double *Gaend, double *habeg, double *haend, double *hmbeg, double *hmend, double *uabeg, double *uaend, double *umbeg, double *umend, int nfcBC, int nGsBC, double g, double dx, double dt, int n, int nBCa, int nBCm)
{
//again errors at machine precision, result of the division handling?
//############################### FIRST ITERATION #######################################
    double *Gm = malloc((n)*sizeof(double));
    double *hm = malloc((n)*sizeof(double));
    double *um = malloc((n)*sizeof(double));
    double *ua = malloc((n)*sizeof(double));

    ca2midpt(Ga,dx,n,Gm);
    ca2midpt(ha,dx,n,hm);

    int cnBC = nGsBC;
    
    //Boundaries might not be so good
    ufromGh(Gm,hm,hmbeg+(nBCm -cnBC),hmend,umbeg+(nBCm -cnBC),umend,dx,n, cnBC,um);
    
    midpt2ca(um ,dx ,n,ua); 

    cnBC = nfcBC;

    double *Gabc = malloc((n + 2*cnBC)*sizeof(double));
    double *habc = malloc((n + 2*cnBC)*sizeof(double));
    double *uabc = malloc((n + 2*cnBC)*sizeof(double));

    conc(Gabeg+(nBCa - cnBC), Ga, Gaend,cnBC,n,cnBC,Gabc);
    conc(habeg+(nBCa - cnBC), ha, haend,cnBC,n,cnBC,habc);
    conc(uabeg+(nBCa - cnBC), ua, uaend,cnBC,n,cnBC,uabc);

    double *nGa = malloc(n*sizeof(double));
    double *nha = malloc(n*sizeof(double));

    evolve(Gabc,habc,uabc,g,dx,dt,n,cnBC,nha,nGa);

//######################################### SECOND ITERATION #############################
    ca2midpt(nGa,dx,n,Gm);
    ca2midpt(nha,dx,n,hm);

    cnBC = nGsBC;
    
    //Boundaries might not be so good
    ufromGh(Gm,hm,hmbeg+(nBCm -cnBC),hmend,umbeg+(nBCm -cnBC),umend,dx,n, cnBC,um);
    
    midpt2ca(um ,dx ,n,ua); 

    cnBC = nfcBC;

    conc(Gabeg+(nBCa - cnBC), nGa, Gaend,cnBC,n,cnBC,Gabc);
    conc(habeg+(nBCa - cnBC), nha, haend,cnBC,n,cnBC,habc);
    conc(uabeg+(nBCa - cnBC), ua, uaend,cnBC,n,cnBC,uabc);

    double *nGap = malloc(n*sizeof(double));
    double *nhap = malloc(n*sizeof(double));

    evolve(Gabc,habc,uabc,g,dx,dt,n,cnBC,nhap,nGap);


// ################################### RK BUILD ###############################
    double *nGapp = malloc(n*sizeof(double));
    double *nhapp = malloc(n*sizeof(double));

    weightsum(0.75,ha, 0.25,nhap,n,nhapp);
    weightsum(0.75,Ga, 0.25,nGap,n,nGapp);

//######################################### THIRD ITERATION #############################

    ca2midpt(nGapp,dx,n,Gm);
    ca2midpt(nhapp,dx,n,hm);

    cnBC = nGsBC;
    
    //Boundaries might not be so good
    ufromGh(Gm,hm,hmbeg+(nBCm -cnBC),hmend,umbeg+(nBCm -cnBC),umend,dx,n, cnBC,um);
    
    midpt2ca(um ,dx ,n,ua); 

    cnBC = nfcBC;

    conc(Gabeg+(nBCa - cnBC), nGapp, Gaend,cnBC,n,cnBC,Gabc);
    conc(habeg+(nBCa - cnBC), nhapp, haend,cnBC,n,cnBC,habc);
    conc(uabeg+(nBCa - cnBC), ua, uaend,cnBC,n,cnBC,uabc);

    double *nGappp = malloc(n*sizeof(double));
    double *nhappp = malloc(n*sizeof(double));

    evolve(Gabc,habc,uabc,g,dx,dt,n,cnBC,nhappp,nGappp);

// ################################### RK BUILD ###############################
    weightsum(i3,ha,2*i3,nhappp,n,ha);
    weightsum(i3,Ga,2*i3,nGappp,n,Ga);

    free(Gm);
    free(hm);
    free(um);
    free(ua);
    free(Gabc);
    free(habc);
    free(uabc);
    free(nGa);
    free(nha);
    free(nGap);
    free(nhap);
    free(nGapp);
    free(nhapp);
    free(nGappp);
    free(nhappp);
}

double *mallocPy(int n)
{
    double *x = malloc(n*sizeof(double));
    return x;
}

void writetomem(double *x, int i , double f)
{
    x[i] = f;

}

double readfrommem(double*x,int i)
{
    return x[i];
}

void deallocPy(double *x)
{
    free(x);
}

int main()
{
    printf("h");
    return 1;
}
