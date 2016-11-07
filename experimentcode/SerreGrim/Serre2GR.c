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


void conc(double *a , double *b, double *c,int n,int m ,int k, double *d)
{
    //replace with memcpy for performance?
    memcpy(d,a,n*sizeof(double));
    memcpy(d+n,b,m*sizeof(double));
    memcpy(d+n+m,c,k*sizeof(double));
}

void addBC(double *h, double *h0, double *h1, int nBC, int n, int nBCs, double *nh)
{ 
    //maybe an error in calculating G
    int j,k;


    double *hb = malloc(nBC*sizeof(double));
    double *he = malloc(nBC*sizeof(double));
    //ADD BC we get h,G need to solve for u then add in boundaries for G using u and h

    //front end
    //i keeps track of big array
    //j keeps track of small bc arrays
    //k keeps track of small new bc arrays
    j = nBCs-1;

    for(k = nBC -1; k > -1 ; k--)
    {
        hb[k] = h0[j];
        j--;
    }

    //back end
    j = 0;
    for(k = 0; k < nBC ; k++)
    {
        he[k] = h1[j];
        j++;

    }

    //bring them all together

    conc(hb,h,he,nBC,n,nBC,nh);

    free(hb);
    free(he);  

}
/*
void evolveh(double *h, double *u, double g, double dx, double dt, int nBC, int n,double *nh)
{
    //modifies nh and nG to give the new values of h and G after a single time step
    double idx = 1.0 / dx;
    int i = nBC - 1;
    double aiph,aimh;

    //calculate gradients

    for (i = nBC ; i < n +nBC;i++)
    {
        aiph = 0.5*(u[i+1] + u[i]);
        aimh = 0.5*(u[i] + u[i-1]);
        nh[i -nBC] = h[i] - 0.5*dt*idx*(u[i+1]*h[i+1] - u[i-1]*h[i-1]) 
            + 0.5*dt*dt*idx*idx*(aiph*(u[i+1]*h[i+1] - u[i]*h[i]) - aimh*(u[i]*h[i] - u[i-1]*h[i-1])); 
 
    }    


   
}
*/

void evolveh(double *h, double *u, double *pu, double g, double dx, double dt, int nBC, int n,double *nh)
{
    //modifies nh and nG to give the new values of h and G after a single time step
    double idx = 1.0 / dx;
    int i = nBC - 1;
    double hiph,himh, uiph,uimh;

    //calculate gradients

    for (i = nBC ; i < n +nBC;i++)
    {
        hiph = 0.5*(h[i+1] + h[i]) - 0.5*(dt*idx)*(u[i+1]*h[i+1] - u[i]*h[i]);
        himh = 0.5*(h[i] + h[i-1]) - 0.5*(dt*idx)*(u[i]*h[i] - u[i-1]*h[i-1]);

		uiph = 0.5*(0.5*(u[i+1] + pu[i+1]) + 0.5*(u[i] + pu[i]));
		uimh = 0.5*(0.5*(u[i] + pu[i]) + 0.5*(u[i-1] + pu[i-1]));

        nh[i -nBC] = h[i] - dt*idx*(hiph*uiph -  himh*uimh); 
 
    }    


   
}
void evolveu(double *h, double *u, double *pu, double g, double dx, double dt, int nBC, int n,double *fu)
{
    //modifies nh and nG to give the new values of h and G after a single time step
    double idx = 1.0 / dx;
    double *a = malloc((n-1)*sizeof(double));
    double *b = malloc(n*sizeof(double));
    double *c = malloc((n-1)*sizeof(double));
    double *f = malloc(n*sizeof(double));

    int i,j;
    double nu,nh,nux,nuxx,nuxxx,nhx,pux,puxx,F,S;


    for (i =1;i < n-1 ; i++)
    {
        //i for tridiagonals, j for the BC u,h vectors 
        j = i + nBC;
        nhx = 0.5*idx*(h[j+1] - h[j-1]);
        
        a[i-1] = 0.5*idx*h[j]*h[j]*nhx - i3*idx*idx*h[j]*h[j]*h[j];
        b[i] = h[j] + 2*i3*idx*idx*h[j]*h[j]*h[j];
        c[i] = -0.5*idx*h[j]*h[j]*nhx - i3*idx*idx*h[j]*h[j]*h[j];

        nu = u[j];
        nh = h[j];
        nux = 0.5*idx*(u[j+1] - u[j-1]);
        nuxx = idx*idx*(u[j+1] - 2*u[j] + u[j-1]);
        nuxxx = 0.5*idx*idx*idx*(u[j+2] - 2*u[j+1]  + 2*u[j-1] - u[j-2]);
        pux = 0.5*idx*(pu[j+1] - pu[j-1]);
        puxx = idx*idx*(pu[j+1] - 2*pu[j] + pu[j-1]);
        F = nu*nh*nux + g*nh*nhx + nh*nh*nhx*nux*nux + i3*nh*nh*nh*nux*nuxx 
            - nu*nh*nh*nhx*nuxx -i3*nh*nh*nh*nu*nuxxx;

        S = 2*dt*F - pu[j]*nh + nh*nh*nhx*pux + i3*nh*nh*nh*puxx;
        f[i] = -S;

    }

    //Boundaries
    i = 0;
    j = i + nBC;

    nhx = 0.5*idx*(h[j+1] - h[j-1]);
    //a[i-1] = 0.5*idx*h[j]*h[j]*nhx - i3*idx*idx*h[j]*h[j]*h[j];
    b[i] = h[j] + 2*i3*idx*idx*h[j]*h[j]*h[j];
    c[i] = -0.5*idx*h[j]*h[j]*nhx - i3*idx*idx*h[j]*h[j]*h[j];

    nu = u[j];
    nh = h[j];
    nux = 0.5*idx*(u[j+1] - u[j-1]);
    nuxx = idx*idx*(u[j+1] - 2*u[j] + u[j-1]);
    nuxxx = 0.5*idx*idx*idx*(u[j+2] - 2*u[j+1]  + 2*u[j-1] - u[j-2]);
    nhx = 0.5*idx*(h[j+1] - h[j-1]);
    pux = 0.5*idx*(pu[j+1] - pu[j-1]);
    puxx = idx*idx*(pu[j+1] - 2*pu[j] + pu[j-1]);
    F = nu*nh*nux + g*nh*nhx + nh*nh*nhx*nux*nux + i3*nh*nh*nh*nux*nuxx 
        - nu*nh*nh*nhx*nuxx -i3*nh*nh*nh*nu*nuxxx;

    S = 2*dt*F - pu[j]*nh + nh*nh*nhx*pux + i3*nh*nh*nh*puxx;
    f[i] = -S - u[j-1]*(0.5*idx*h[j]*h[j]*nhx - i3*idx*idx*h[j]*h[j]*h[j]);

    i = n-1;
    j = i + nBC;

    nhx = 0.5*idx*(h[j+1] - h[j-1]);
    a[i-1] = 0.5*idx*h[j]*h[j]*nhx - i3*idx*idx*h[j]*h[j]*h[j];
    b[i] = h[j] + 2*i3*idx*idx*h[j]*h[j]*h[j];
    //c[i] = -0.5*idx*h[j]*h[j]*nhx - i3*idx*idx*h[j]*h[j]*h[j];

    nu = u[j];
    nh = h[j];
    nux = 0.5*idx*(u[j+1] - u[j-1]);
    nuxx = idx*idx*(u[j+1] - 2*u[j] + u[j-1]);
    nuxxx = 0.5*idx*idx*idx*(u[j+2] - 2*u[j+1]  + 2*u[j-1] - u[j-2]);
    nhx = 0.5*idx*(h[j+1] - h[j-1]);
    pux = 0.5*idx*(pu[j+1] - pu[j-1]);
    puxx = idx*idx*(pu[j+1] - 2*pu[j] + pu[j-1]);
    F = nu*nh*nux + g*nh*nhx + nh*nh*nhx*nux*nux + i3*nh*nh*nh*nux*nuxx 
        - nu*nh*nh*nhx*nuxx -i3*nh*nh*nh*nu*nuxxx;

    S = 2*dt*F - pu[j]*h[j] + nh*nh*nhx*pux + i3*nh*nh*nh*puxx;
    f[i] = -S - u[j+1]*(-0.5*idx*h[j]*h[j]*nhx - i3*idx*idx*h[j]*h[j]*h[j]);

    TDMA(a,b,c,f,n,fu);

    free(a);
    free(b);
    free(c); 
    free(f);  
   
}


void evolvewrap(double *u, double *h, double *pubc, double *h0, double *h1, double *u0, double *u1, double g, double dx, double dt, int nBC, int n, int nBCs)
{

    double *ubc = malloc((n + 2*nBC)*sizeof(double));
    double *hbc = malloc((n + 2*nBC)*sizeof(double));
    addBC(h, h0, h1,nBC, n,nBCs,hbc);
    addBC(u, u0, u1,nBC, n,nBCs,ubc);
    evolveu(hbc, ubc,pubc,g,dx,dt,nBC,n,u); 
    memcpy(pubc,ubc,(n + 2*nBC)*sizeof(double));


	addBC(u, u0, u1,nBC, n,nBCs,ubc);

	evolveh(hbc, ubc,pubc,g,dx, dt,nBC, n,h); 
    free(ubc);
    free(hbc);

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
