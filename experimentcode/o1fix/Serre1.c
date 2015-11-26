#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

//just do a simple first order Scheme, with no bed variations

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

void getufromG(double *h, double *G, double u0, double u1, double h0, double h1, double dx , int n, double *u)
{
    double idx = 1.0 / dx;
    double ithree = 1.0 / 3.0;
    double *a = malloc((n-1)*sizeof(double));
    double *b = malloc(n*sizeof(double));
    double *c = malloc((n-1)*sizeof(double));

    int i;
    double thx;


    for (i =1;i < n-1 ; i++)
    {
        thx = 0.5*idx*(h[i+1] - h[i-1]);
        
        a[i-1] = -ithree*idx*idx*h[i]*h[i]*h[i] + 0.5*idx*h[i]*h[i]*thx;
        b[i] = h[i] + 2.0*ithree*idx*idx*h[i]*h[i]*h[i];
        c[i] = -ithree*idx*idx*h[i]*h[i]*h[i] - 0.5*idx*h[i]*h[i]*thx;

    }

    //Boundaries
    i = 0;

    thx = 0.5*idx*(h[i+1] - h0);
        
    //a[i-1] = -ithree*idx*idx*h[i]*h[i]*h[i] + 0.5*idx*th*th*thx
    b[i] = h[i] + 2.0*ithree*idx*idx*h[i]*h[i]*h[i];
    c[i] = -ithree*idx*idx*h[i]*h[i]*h[i] - 0.5*idx*h[i]*h[i]*thx;

    double tmpG1 = G[i];

    G[i] = tmpG1 - u0*(-ithree*idx*idx*h[i]*h[i]*h[i] + 0.5*idx*h[i]*h[i]*thx);

    i = n-1;

    thx = 0.5*idx*(h1 - h[i-1]);
        
    a[i-1] = -ithree*idx*idx*h[i]*h[i]*h[i] + 0.5*idx*h[i]*h[i]*thx;
    b[i] = h[i] + 2.0*ithree*idx*idx*h[i]*h[i]*h[i];
    //c[i] = -ithree*idx*idx*h[i]*h[i]*h[i] - 0.5*idx*h[i]*h[i]*thx;

    double tmpG2 = G[i];
    G[i] = tmpG2 - u1*(-ithree*idx*idx*h[i]*h[i]*h[i] - 0.5*idx*h[i]*h[i]*thx);

    TDMA(a,b,c,G,n,u);

    G[0] = tmpG1;
    G[n-1] = tmpG2;


    free(a);
    free(b);
    free(c);

}

void getGfromu(double *h, double *u, double u0, double u1, double h0, double h1, double dx , int n, double *G)
{
    double idx = 1.0 / dx;
    double ithree = 1.0 / 3.0;

    int i;
    double thx;


    for (i =1;i < n-1 ; i++)
    {
        thx = 0.5*idx*(h[i+1] - h[i-1]);

        G[i] = (-ithree*idx*idx*h[i]*h[i]*h[i] + 0.5*idx*h[i]*h[i]*thx)*u[i-1] +
                (h[i] + 2.0*ithree*idx*idx*h[i]*h[i]*h[i])*u[i] +
                (-ithree*idx*idx*h[i]*h[i]*h[i] - 0.5*idx*h[i]*h[i]*thx)*u[i+1];
    }

    //Boundaries
    i = 0;

    thx = 0.5*idx*(h[i+1] - h0);

    G[i] = (-ithree*idx*idx*h[i]*h[i]*h[i] + 0.5*idx*h[i]*h[i]*thx)*u0 +
            (h[i] + 2.0*ithree*idx*idx*h[i]*h[i]*h[i])*u[i] +
            (-ithree*idx*idx*h[i]*h[i]*h[i] - 0.5*idx*h[i]*h[i]*thx)*u[i+1];

    i = n-1;

    thx = 0.5*idx*(h1 - h[i-1]);

    G[i] = (-ithree*idx*idx*h[i]*h[i]*h[i] + 0.5*idx*h[i]*h[i]*thx)*u[i-1] +
            (h[i] + 2.0*ithree*idx*idx*h[i]*h[i]*h[i])*u[i] +
            (-ithree*idx*idx*h[i]*h[i]*h[i] - 0.5*idx*h[i]*h[i]*thx)*u1;
}

void conc(double *a , double *b, double *c,int n,int m ,int k, double *d)
{
    //replace with memcpy for performance?
    memcpy(d,a,n*sizeof(double));
    memcpy(d+n,b,m*sizeof(double));
    memcpy(d+n+m,c,k*sizeof(double));
}

void evolveBC(double *G, double *h, double *h0, double *h1, double *u0, double *u1, double g, double dx, double dt, int nBC, int n, int nBCs,double *nG, double *nh, double *nu)
{ 
    //maybe an error in calculating G
    double idx = 1.0 / dx;
    double ithree = 1.0 / 3.0;
    int j = nBCs -1;


    double *Gb = malloc(nBC*sizeof(double));
    double *Ge = malloc(nBC*sizeof(double));
    double *ub = malloc(nBC*sizeof(double));
    double *ue = malloc(nBC*sizeof(double));
    double *hb = malloc(nBC*sizeof(double));
    double *he = malloc(nBC*sizeof(double));
    double *u = malloc(n*sizeof(double));
    //ADD BC we get h,G need to solve for u then add in boundaries for G using u and h
    getufromG(h,G, u0[j], u1[0], h0[j], h1[0],dx,n,u);

    //front end
    //i keeps track of big array
    //j keeps track of small bc arrays
    //k keeps track of small new bc arrays
    int i = -1;
    int k = nBC -1;
    j = nBCs -1;


    double hx = 0.5*idx*(h[i+1] - h0[j-1]);

    double ai =0.5*idx*h0[j]*h0[j]*hx - idx*idx*ithree*h0[j]*h0[j]*h0[j];
    double bi = h0[j] + 2*ithree*idx*idx*h0[j]*h0[j]*h0[j];
    double ci = -0.5*idx*h0[j]*h0[j]*hx - idx*idx*ithree*h0[j]*h0[j]*h0[j];
    Gb[k] = ai*u0[j-1]
                 + bi*u0[j]
                 + ci*u[i+1];
    ub[k] = u0[j];
    hb[k] = h0[j];


    for(k = k-1; k > -1 ; k--)
    {
        j--;
        hx = 0.5*idx*(h0[j+1] - h0[j-1]);

        Gb[k] = (0.5*idx*h0[j]*h0[j]*hx - idx*idx*ithree*h0[j]*h0[j]*h0[j])*u0[j-1]
                    + (h0[j] + 2*ithree*idx*idx*h0[j]*h0[j]*h0[j])*u0[j]
                    + (-0.5*idx*h0[j]*h0[j]*hx - idx*idx*ithree*h0[j]*h0[j]*h0[j])*u0[j+1];

        ub[k] = u0[j];
        hb[k] = h0[j];
    }

    //back end
    i = n;
    k = 0;
    j = 0;

    hx = 0.5*idx*(h1[j+1] - h[i-1]);
    Ge[k] = (0.5*idx*h1[j]*h1[j]*hx - idx*idx*ithree*h1[j]*h1[j]*h1[j])*u[i-1]
            + (h1[j] + 2*ithree*idx*idx*h1[j]*h1[j]*h1[j])*u1[j]
            + (-0.5*idx*h1[j]*h1[j]*hx - idx*idx*ithree*h1[j]*h1[j]*h1[j])*u1[j+1];

    ue[k] = u1[j];
    he[k] = h1[j];

    for(k = k+1; k < nBC ; k++)
    {
        j++;
        hx = 0.5*idx*(h1[j+1] - h1[j-1]);
        Ge[k] = (0.5*idx*h1[j]*h1[j]*hx - idx*idx*ithree*h1[j]*h1[j]*h1[j])*u1[j-1]
                    + (h1[j] + 2*ithree*idx*idx*h1[j]*h1[j]*h1[j])*u1[j]
                    + (-0.5*idx*h1[j]*h1[j]*hx - idx*idx*ithree*h1[j]*h1[j]*h1[j])*u1[j+1];
        ue[k] = u1[j];
        he[k] = h1[j];

    }

    //bring them all together

    conc(Gb,G,Ge,nBC,n,nBC,nG);
    conc(hb,h,he,nBC,n,nBC,nh);
    conc(ub,u,ue,nBC,n,nBC,nu);

    free(Gb);
    free(Ge);
    free(hb);
    free(he);
    free(ub);
    free(ue);
    free(u);    

}

void evolve(double *G, double *h, double *u, double g, double dx, double dt, int nBC, int n,double *nG, double *nh)
{
    //modifies nh and nG to give the new values of h and G after a single time step
    double idx = 1.0 / dx;
    double ithree = 1.0 / 3.0;
    int i = nBC - 1;

    double ue = 0.5*(u[i+1]+u[i]);

    //calculate values at right of i cell
    double hir = h[i];
    double Gir = G[i];
    double uir = ue;

    //calculate values at left of i+1 cell
    double hip1l = h[i+1];
    double Gip1l = G[i+1];
    double uip1l = ue;

    //right force

    double duer = idx*(u[i+1] - u[i]);
    double duel = idx*(u[i+1] - u[i]);

    double sqrtghel = sqrt(g* hir);
    double sqrtgher = sqrt(g* hip1l);

    double sl = fmin(0,fmin(uir - sqrtghel, uip1l - sqrtgher));
    double sr = fmax(0,fmax(uir + sqrtghel, uip1l + sqrtgher));

    double felh = uir*hir;
    double felG = Gir*uir + 0.5*g*hir*hir - 2*ithree*hir*hir*hir*duel*duel;
    double ferh = uip1l*hip1l;
    double ferG = Gip1l*uip1l + 0.5*g*hip1l*hip1l -2*ithree*hip1l*hip1l*hip1l*duer*duer;

    double isrmsl = 0.0;

    if(sr != sl) isrmsl = 1.0 / (sr - sl);
   
    double foh = isrmsl*(sr*felh - sl*ferh + sl*sr*(hip1l - hir));

    double foG = isrmsl*(sr*felG - sl*ferG + sl*sr*(Gip1l - Gir));

    double fih = foh;
    double fiG = foG;
    for (i = nBC ; i < n +nBC;i++)
    {
        //i right

        ue = 0.5*(u[i+1]+u[i]);
        hir = h[i];
        Gir = G[i];
        uir = ue;

        //i+1 left

        hip1l = h[i+1];
        Gip1l = G[i+1];
        uip1l = ue;


        duer = idx*(u[i+1] - u[i]);
        duel = idx*(u[i+1] - u[i]);

        sqrtghel = sqrt(g*hir);
        sqrtgher = sqrt(g*hip1l);

        sl = fmin(0,fmin(uir - sqrtghel, uip1l - sqrtgher));
        sr = fmax(0,fmax(uir + sqrtghel, uip1l + sqrtgher));

        felh = uir*hir;
        felG = Gir*uir + 0.5*g*hir*hir - 2*ithree*hir*hir*hir*duel*duel;
        ferh = uip1l*hip1l;
        ferG = Gip1l*uip1l + 0.5*g*hip1l*hip1l -2*ithree*hip1l*hip1l*hip1l*duer*duer;

        isrmsl = 0.0;

        if(sr != sl) isrmsl = 1.0 / (sr - sl);
   
        foh = isrmsl*(sr*felh - sl*ferh + sl*sr*(hip1l - hir));
        foG = isrmsl*(sr*felG - sl*ferG + sl*sr*(Gip1l - Gir));

        //source term
        nh[i -nBC] = h[i] -dt*idx*(foh - fih);
        nG[i -nBC] = G[i] -dt*idx*(foG -fiG);

        fih = foh;
        fiG = foG;  
 
    }
    
   
}

void evolvewrap(double *G, double *h, double *h0, double *h1, double *u0, double *u1, double g, double dx, double dt, int nBC, int n, int nBCs)
{
    //first allocate memory for BC variables
    double *Gbc = malloc((n + 2*nBC)*sizeof(double));
    double *hbc = malloc((n + 2*nBC)*sizeof(double));
    double *ubc = malloc((n + 2*nBC)*sizeof(double));

    evolveBC(G,h,h0,h1,u0,u1,g,dx,dt,nBC,n,nBCs,Gbc,hbc,ubc);
    evolve(Gbc,hbc,ubc,g,dx,dt,nBC,n,G,h);


    free(Gbc);
    free(hbc);
    free(ubc);

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
