#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

//put in slope limiting.

const double i24 = 1.0/24.0;
const double i12 = 1.0/12.0;
const double i3 = 1.0/3.0;
const double i8 = 1.0/8.0;
const double i48 = 1.0/48.0;
double minmod(double a, double b, double c)
{
    if((a > 0) && (b>0) && (c>0))
    {
        return fmin(a,fmin(b,c));
    }
    else if((a < 0) && (b<0) && (c<0))
    {
        return fmax(a,fmax(b,c));
    }
        return 0.0;
}


    
double interpquarticval(double *coeff,double xj,double x)
{    
    return coeff[0]*(x -xj)*(x -xj)*(x -xj)*(x -xj) + coeff[1]*(x -xj)*(x -xj)*(x -xj)
    + coeff[2]*(x -xj)*(x -xj) + coeff[3]*(x -xj)+ coeff[4];
}  
  
double interpquarticgrad(double *coeff,double xj,double x)
{
    
    return 4*coeff[0]*(x -xj)*(x -xj)*(x -xj) + 3*coeff[1]*(x -xj)*(x -xj)
    + 2*coeff[2]*(x -xj) + coeff[3];
}    
void interpquartcoeff(double *q,double *coeff,int j,double dx)
{
    double idx = 1.0/dx;

    coeff[0] = i24*idx*idx*idx*idx*(q[j+2] - 4*q[j+1] + 6*q[j] - 4*q[j-1] + q[j-2]);
    coeff[1] = i12*idx*idx*idx*(q[j+2] - 2*q[j+1] + 2*q[j-1] - q[j-2]);
    coeff[2] = i24*idx*idx*(-q[j+2] + 16*q[j+1] - 30*q[j] + 16*q[j-1] - q[j-2]);
    coeff[3] = i12*idx*(-q[j+2] + 8*q[j+1] - 8*q[j-1] + q[j-2]);
    coeff[4] = q[j];
}
    
double HankEnergyacrosscell(double *x,double *h,double *u,double g,int j,double dx)
{
    //so we have h,u at midpoints
    //epsilon and sigma are everywhere

	double *ucoeff = malloc(5*sizeof(double));
	double *hcoeff = malloc(5*sizeof(double));
	

    //jth cell
    interpquartcoeff(u,ucoeff,j,dx);
    interpquartcoeff(h,hcoeff,j,dx);
    
    //first gauss point
    double fgp = 0.5*dx*sqrt(3.0/5.0) + x[j];
    double fgph = interpquarticval(hcoeff,x[j],fgp);
    double fgpu = interpquarticval(ucoeff,x[j],fgp);
    double fgpux = interpquarticgrad(ucoeff,x[j],fgp);
    
    double fgpe = fgph*fgpu*fgpu + g*fgph*fgph + i3*(fgph*fgph*fgph)*fgpux*fgpux;
        
    //second gauss point
    double sgp = x[j];
    double sgph = interpquarticval(hcoeff,x[j],sgp);
    double sgpu = interpquarticval(ucoeff,x[j],sgp);
    double sgpux = interpquarticgrad(ucoeff,x[j],sgp);
    
    double sgpe = sgph*sgpu*sgpu + g*sgph*sgph + i3*(sgph*sgph*sgph)*sgpux*sgpux;

    //third gauss point
    double tgp = -0.5*dx*sqrt(3.0/5.0) + x[j];
    double tgph = interpquarticval(hcoeff,x[j],tgp);
    double tgpu = interpquarticval(ucoeff,x[j],tgp);
    double tgpux = interpquarticgrad(ucoeff,x[j],tgp);
    
    double tgpe = tgph*tgpu*tgpu + g*tgph*tgph + i3*(tgph*tgph*tgph)*tgpux*tgpux;

	free(ucoeff);
	free(hcoeff);
    
    return 0.5*dx*( (5.0/9.0)*fgpe + (8.0/9.0)*sgpe + (5.0/9.0)*tgpe);
}
    
double HankEnergyall(double *x,double *h,double *u,double g,int n, int nBC,double dx)
{
    double sum1 = 0.0;
	int i;
	for(i = nBC; i < n - nBC;i++)
	{
       sum1 = sum1 + HankEnergyacrosscell(x,h,u,g,i,dx);
	}
    return 0.5*sum1; 

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

    Gb[k] = (0.5*idx*h0[j]*h0[j]*hx - idx*idx*ithree*h0[j]*h0[j]*h0[j])*u0[j-1]
                    + (h0[j] + 2*ithree*idx*idx*h0[j]*h0[j]*h0[j])*u0[j]
                    + (-0.5*idx*h0[j]*h0[j]*hx - idx*idx*ithree*h0[j]*h0[j]*h0[j])*u[i+1];
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

void evolve(double *G, double *h, double *u, double g, double dx, double dt, int nBC, int n,double *nG, double *nh, double theta)
{
    //modifies nh and nG to give the new values of h and G after a single time step
    double idx = 1.0 / dx;
    double ithree = 1.0 / 3.0;
    int i = nBC - 1;

    //calculate gradients
    double dhib = (h[i] - h[i-1]);
    double dhim = 0.5*(h[i+1] - h[i-1]);
    double dhif = (h[i+1] - h[i]);

    //double duib = (u[i] - u[i-1]);
    //double dui = 0.5*(u[i+1] - u[i-1]);
    //double duif = (u[i+1] - u[i]);

    double dGib = (G[i] - G[i-1]);
    double dGim = 0.5*(G[i+1] - G[i-1]);
    double dGif = (G[i+1] - G[i]);


    //calculate values at right of i cell
    double dhi = minmod(theta*dhib, dhim, theta*dhif);
    //double dui = minmod(theta*duib, duim, theta*duif);
    double dGi = minmod(theta*dGib, dGim, theta*dGif);

    double ue = 0.5*(u[i+1]+u[i]);

    double hir = h[i] + 0.5*dhi;
    double Gir = G[i] + 0.5*dGi;
    double uir = ue;

    //calculate values at left of i+1 cell

    //calculate gradients
    double dhip1b = (h[i+1] - h[i]);
    double dhip1m = 0.5*(h[i+2] - h[i]);
    double dhip1f = (h[i+2] - h[i+1]);

    //double duip1b = (u[i+1] - u[i]);
    //double duip1m = 0.5*(u[i+2] - u[i]);
    //double duip1f = (u[i+2] - u[i+1]);

    double dGip1b = (G[i+1] - G[i]);
    double dGip1m = 0.5*(G[i+2] - G[i]);
    double dGip1f = (G[i+2] - G[i+1]);

    double dhip1 = minmod(theta*dhip1b, dhip1m, theta*dhip1f);
    //double duip1 = minmod(theta*duip1b, duip1m, theta*duip1f);
    double dGip1 = minmod(theta*dGip1b, dGip1m, theta*dGip1f);

    double hip1l = h[i+1] - 0.5*dhip1;
    double Gip1l = G[i+1] - 0.5*dGip1;
    double uip1l = ue;

    //right force
    double duer = idx*(u[i+1] - u[i]);

    double sqrtghel = sqrt(g*hir);
    double sqrtgher = sqrt(g*hip1l);

    double sl = fmin(0,fmin(uir - sqrtghel, uip1l - sqrtgher));
    double sr = fmax(0,fmax(uir + sqrtghel, uip1l + sqrtgher));

    double felh = uir*hir;
    double felG = Gir*uir + 0.5*g*hir*hir - 2*ithree*hir*hir*hir*duer*duer;
    double ferh = uip1l*hip1l;
    double ferG = Gip1l*uip1l + 0.5*g*hip1l*hip1l -2*ithree*hip1l*hip1l*hip1l*duer*duer;

    double isrmsl = 0.0;

    if(sr != sl) isrmsl = 1.0 / (sr - sl);
   
    double foh = isrmsl*(sr*felh - sl*ferh + sl*sr*(hip1l - hir));

    double foG = isrmsl*(sr*felG - sl*ferG + sl*sr*(Gip1l - Gir));

    double fih = foh;
    double fiG = foG;

    dhi = dhip1;
    //dui = duip1;
    dGi = dGip1;

    for (i = nBC ; i < n +nBC;i++)
    {
        ue = 0.5*(u[i+1]+u[i]);
        //i right

        hir = h[i] + 0.5*dhi;
        Gir = G[i] + 0.5*dGi;
        uir = ue;

        //calculate gradients
        dhip1b = (h[i+1] - h[i]);
        dhip1m = 0.5*(h[i+2] - h[i]);
        dhip1f = (h[i+2] - h[i+1]);

        //duip1b = (u[i+1] - u[i]);
        //duip1m = 0.5*(u[i+2] - u[i]);
        //duip1f = (u[i+2] - u[i+1]);

        dGip1b = (G[i+1] - G[i]);
        dGip1m = 0.5*(G[i+2] - G[i]);
        dGip1f = (G[i+2] - G[i+1]);

        dhip1 = minmod(theta*dhip1b, dhip1m, theta*dhip1f);
        //duip1 = minmod(theta*duip1b, duip1m, theta*duip1f);
        dGip1 = minmod(theta*dGip1b, dGip1m, theta*dGip1f);

        hip1l = h[i+1] - 0.5*dhip1;
        Gip1l = G[i+1] - 0.5*dGip1;
        uip1l = ue;

        duer = idx*(u[i+1] - u[i]);

        sqrtghel = sqrt(g*hir);
        sqrtgher = sqrt(g*hip1l);

        sl = fmin(0,fmin(uir - sqrtghel, uip1l - sqrtgher));
        sr = fmax(0,fmax(uir + sqrtghel, uip1l + sqrtgher));

        felh = uir*hir;
        felG = Gir*uir + 0.5*g*hir*hir - 2*ithree*hir*hir*hir*duer*duer;
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

        dhi = dhip1;
        //dui = duip1;
        dGi = dGip1;   
 
    }
    
   
}

void evolvewrap(double *G, double *h, double *h0, double *h1, double *u0, double *u1, double g, double dx, double dt, int nBC, int n, int nBCs, double theta)
{
    //first allocate memory for BC variables
    double *Gbc = malloc((n + 2*nBC)*sizeof(double));
    double *hbc = malloc((n + 2*nBC)*sizeof(double));
    double *ubc = malloc((n + 2*nBC)*sizeof(double));

    double *Gp = malloc(n*sizeof(double));
    double *hp = malloc(n*sizeof(double));

    evolveBC(G,h,h0,h1,u0,u1,g,dx,dt,nBC,n,nBCs,Gbc,hbc,ubc);
    evolve(Gbc,hbc,ubc,g,dx,dt,nBC,n,Gp,hp,theta);

    evolveBC(Gp,hp,h0,h1,u0,u1,g,dx,dt,nBC,n,nBCs,Gbc,hbc,ubc);
    evolve(Gbc,hbc,ubc,g,dx,dt,nBC,n,Gp,hp,theta);

    int i;
    for(i=0;i<n;i++)
    {
        G[i] = 0.5*(G[i] + Gp[i]);
        h[i] = 0.5*(h[i] + hp[i]);
    }


    free(Gbc);
    free(hbc);
    free(ubc);
    free(Gp);
    free(hp);

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
