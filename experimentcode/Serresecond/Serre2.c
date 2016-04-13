#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

//put in slope limiting.

const double i24 = 1.0/24.0;
const double i48 = 1.0/48.0;
const double i12 = 1.0/12.0;
const double i3 = 1.0/3.0;
const double i8 = 1.0/8.0;

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

double interpquarticHess(double *coeff,double xj,double x)
{
    
    return 12*coeff[0]*(x -xj)*(x -xj)+ 6*coeff[1]*(x -xj)
    + 2*coeff[2];
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
    
double MyEnergyacrosscell(double *x,double *h,double *u, double *b,double g,int j,double dx)
{
    //so we have h,u at midpoints
    //epsilon and sigma are everywhere

	double *ucoeff = malloc(5*sizeof(double));
	double *hcoeff = malloc(5*sizeof(double));
	double *bcoeff = malloc(5*sizeof(double));
	

    //jth cell
    interpquartcoeff(u,ucoeff,j,dx);
    interpquartcoeff(h,hcoeff,j,dx);
    interpquartcoeff(b,bcoeff,j,dx);
    
    //first gauss point
    double fgp = 0.5*dx*sqrt(3.0/5.0) + x[j];
    double fgph = interpquarticval(hcoeff,x[j],fgp);
    double fgpu = interpquarticval(ucoeff,x[j],fgp);
    double fgpux = interpquarticgrad(ucoeff,x[j],fgp);
	double fgpb = interpquarticval(bcoeff,x[j],fgp);
    
    double fgpe = fgph*fgpu*fgpu + 2*g*fgph*(0.5*fgph + fgpb) + i12*(fgph*fgph*fgph)*fgpux*fgpux + fgph*(0.5*fgph + fgpb);
        
    //second gauss point
    double sgp = x[j];
    double sgph = interpquarticval(hcoeff,x[j],sgp);
    double sgpu = interpquarticval(ucoeff,x[j],sgp);
    double sgpux = interpquarticgrad(ucoeff,x[j],sgp);
	double sgpb = interpquarticval(bcoeff,x[j],sgp);
    
    double sgpe = sgph*sgpu*sgpu + 2*g*sgph*(0.5*sgph + sgpb) + i12*(sgph*sgph*sgph)*sgpux*sgpux + sgph*(0.5*sgph + sgpb);

    //third gauss point
    double tgp = -0.5*dx*sqrt(3.0/5.0) + x[j];
    double tgph = interpquarticval(hcoeff,x[j],tgp);
    double tgpu = interpquarticval(ucoeff,x[j],tgp);
    double tgpux = interpquarticgrad(ucoeff,x[j],tgp);
	double tgpb = interpquarticval(bcoeff,x[j],tgp);
    
    double tgpe = tgph*tgpu*tgpu + 2*g*tgph*(0.5*tgph + tgpb) + i12*(tgph*tgph*tgph)*tgpux*tgpux + tgph*(0.5*tgph + tgpb);

	free(ucoeff);
	free(hcoeff);
	free(bcoeff);
    
    return 0.5*dx*( (5.0/9.0)*fgpe + (8.0/9.0)*sgpe + (5.0/9.0)*tgpe);
}
    
double MyEnergyall(double *x,double *h,double *u, double *b,double g,int n, int nBC,double dx)
{
	//include approximations to H(a)u(a) + p(a)u(a) - H(b)u(b) - p(b)u(b) (a end, b beg)
    double sum1 = 0.0;
	int i;
	for(i = nBC; i < n - nBC;i++)
	{
       sum1 = sum1 + MyEnergyacrosscell(x,h,u,b,g,i,dx);
	}
	sum1 = 0.5*sum1;

    return sum1; 

}

double MyCorrectionacrosstimestep(double *t,double *Corrintx,int j,double dt)
{
    //so we have h,u at midpoints
    //epsilon and sigma are everywhere

	double *Corrintxcoeff = malloc(5*sizeof(double));
	

    //jth cell
    interpquartcoeff(Corrintx,Corrintxcoeff,j,dt);
    
    //first gauss point
    double fgp = 0.5*dt*sqrt(3.0/5.0) + t[j];
    double fgpe = interpquarticval(Corrintxcoeff,t[j],fgp);
        
    //second gauss point
    double sgp = t[j];
  
    double sgpe = interpquarticval(Corrintxcoeff,t[j],sgp);

    //third gauss point
    double tgp = -0.5*dt*sqrt(3.0/5.0) + t[j];
	  
    double tgpe = interpquarticval(Corrintxcoeff,t[j],tgp);

	free(Corrintxcoeff);
    
    return 0.5*dt*( (5.0/9.0)*fgpe + (8.0/9.0)*sgpe + (5.0/9.0)*tgpe);
}

double MyCorrectionalltimes(double *t,double *Corrintx, int startj, int endj,double dt)
{
	//include approximations to H(a)u(a) + p(a)u(a) - H(b)u(b) - p(b)u(b) (a end, b beg)
    double sum1 = 0.0;
	int i;
	for(i = startj; i < endj;i++)
	{
       sum1 = sum1 + MyCorrectionacrosstimestep(t,Corrintx,i,dt);
	}
    return sum1; 

}

double MyCorrectacrosscell(double *x,double *h,double *u, double *ut, double *b,double g,int j,double dx)
{
    //so we have h,u at midpoints
    //epsilon and sigma are everywhere

	double *ucoeff = malloc(5*sizeof(double));
	double *utcoeff = malloc(5*sizeof(double));
	double *hcoeff = malloc(5*sizeof(double));
	double *bcoeff = malloc(5*sizeof(double));
	

    //jth cell
    interpquartcoeff(u,ucoeff,j,dx);
	interpquartcoeff(ut,utcoeff,j,dx);	
    interpquartcoeff(h,hcoeff,j,dx);
    interpquartcoeff(b,bcoeff,j,dx);
    
    //first gauss point
    double fgp = 0.5*dx*sqrt(3.0/5.0) + x[j];
    double fgph = interpquarticval(hcoeff,x[j],fgp);
    double fgpu = interpquarticval(ucoeff,x[j],fgp);
	double fgput = interpquarticval(utcoeff,x[j],fgp);
	double fgputx = interpquarticgrad(utcoeff,x[j],fgp);
    double fgpux = interpquarticgrad(ucoeff,x[j],fgp);
	double fgpuxx = interpquarticHess(ucoeff,x[j],fgp);
	double fgpb = interpquarticval(bcoeff,x[j],fgp);
	double fgpbx = interpquarticgrad(bcoeff,x[j],fgp);
	double fgpbxx = interpquarticHess(bcoeff,x[j],fgp);
    
	double fgpGam = fgpux*fgpux - fgpu*fgpuxx - fgputx;
	double fgpPhi = fgpbx*(fgput + fgpu*fgpux) + fgpu*fgpu*fgpbxx;
	double fgppbar = g*fgph + 0.5*fgph*fgph*fgpGam + fgph*fgpPhi;
    double fgpe = fgppbar*(-fgpbx*fgpu + 0.5*fgph + fgpb + 0.5*fgph*fgpux);
        
    //second gauss point
    fgp = x[j];
    fgph = interpquarticval(hcoeff,x[j],fgp);
    fgpu = interpquarticval(ucoeff,x[j],fgp);
	fgput = interpquarticval(utcoeff,x[j],fgp);
	fgputx = interpquarticgrad(utcoeff,x[j],fgp);
    fgpux = interpquarticgrad(ucoeff,x[j],fgp);
	fgpuxx = interpquarticHess(ucoeff,x[j],fgp);
	fgpb = interpquarticval(bcoeff,x[j],fgp);
	fgpbx = interpquarticgrad(bcoeff,x[j],fgp);
	fgpbxx = interpquarticHess(bcoeff,x[j],fgp);
    
	fgpGam = fgpux*fgpux - fgpu*fgpuxx - fgputx;
	fgpPhi = fgpbx*(fgput + fgpu*fgpux) + fgpu*fgpu*fgpbxx;
	fgppbar = g*fgph + 0.5*fgph*fgph*fgpGam + fgph*fgpPhi;
    double sgpe = fgppbar*(-fgpbx*fgpu + 0.5*fgph + fgpb + 0.5*fgph*fgpux);

    //third gauss point
    fgp = -0.5*dx*sqrt(3.0/5.0) + x[j];
    fgph = interpquarticval(hcoeff,x[j],fgp);
    fgpu = interpquarticval(ucoeff,x[j],fgp);
	fgput = interpquarticval(utcoeff,x[j],fgp);
	fgputx = interpquarticgrad(utcoeff,x[j],fgp);
    fgpux = interpquarticgrad(ucoeff,x[j],fgp);
	fgpuxx = interpquarticHess(ucoeff,x[j],fgp);
	fgpb = interpquarticval(bcoeff,x[j],fgp);
	fgpbx = interpquarticgrad(bcoeff,x[j],fgp);
	fgpbxx = interpquarticHess(bcoeff,x[j],fgp);
    
	fgpGam = fgpux*fgpux - fgpu*fgpuxx - fgputx;
	fgpPhi = fgpbx*(fgput + fgpu*fgpux) + fgpu*fgpu*fgpbxx;
	fgppbar = g*fgph + 0.5*fgph*fgph*fgpGam + fgph*fgpPhi;
    double tgpe = fgppbar*(-fgpbx*fgpu + 0.5*fgph + fgpb + 0.5*fgph*fgpux);

	free(ucoeff);
	free(hcoeff);
	free(bcoeff);
    
    return 0.5*dx*( (5.0/9.0)*fgpe + (8.0/9.0)*sgpe + (5.0/9.0)*tgpe);
}


double MyCorrectionallcells(double *x,double *h,double *u, double *ut, double *b,double g,int n, int nBC,double dx)
{
	//include approximations to H(a)u(a) + p(a)u(a) - H(b)u(b) - p(b)u(b) (a end, b beg)
    double sum1 = 0.0;
	int i;
	for(i = nBC; i < n - nBC;i++)
	{
       sum1 = sum1 + MyCorrectacrosscell(x,h,u,ut,b,g,i,dx);
	}
	sum1 = sum1;

    return sum1; 

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

void getufromG(double *h, double *G, double *bed, double u0, double u1, double h0, double h1, double b0, double b1, double dx , int n, double *u)
{
    double idx = 1.0 / dx;
    double ithree = 1.0 / 3.0;
    double *a = malloc((n-1)*sizeof(double));
    double *b = malloc(n*sizeof(double));
    double *c = malloc((n-1)*sizeof(double));

    int i;
    double thx,tbx,tbxx,D;


    for (i =1;i < n-1 ; i++)
    {
        thx = 0.5*idx*(h[i+1] - h[i-1]);
        tbx = 0.5*idx*(bed[i+1] - bed[i-1]);
        tbxx = idx*idx*(bed[i+1] - 2*bed[i] + bed[i-1]);
        D = h[i] + h[i]*thx*tbx + 0.5*h[i]*h[i]*tbxx + h[i]*tbx*tbx;
        
        a[i-1] = -ithree*idx*idx*h[i]*h[i]*h[i] + 0.5*idx*h[i]*h[i]*thx;
        b[i] = D + 2.0*ithree*idx*idx*h[i]*h[i]*h[i];
        c[i] = -ithree*idx*idx*h[i]*h[i]*h[i] - 0.5*idx*h[i]*h[i]*thx;

    }

    //Boundaries
    i = 0;

    thx = 0.5*idx*(h[i+1] - h0);
    tbx = 0.5*idx*(bed[i+1] - b0);
    tbxx = idx*idx*(bed[i+1] - 2*bed[i] + b0);
    D = h[i] + h[i]*thx*tbx + 0.5*h[i]*h[i]*tbxx + h[i]*tbx*tbx;
        
    //a[i-1] = -ithree*idx*idx*h[i]*h[i]*h[i] + 0.5*idx*th*th*thx
    b[i] = D + 2.0*ithree*idx*idx*h[i]*h[i]*h[i];
    c[i] = -ithree*idx*idx*h[i]*h[i]*h[i] - 0.5*idx*h[i]*h[i]*thx;

    double tmpG1 = G[i];

    G[i] = tmpG1 - u0*(-ithree*idx*idx*h[i]*h[i]*h[i] + 0.5*idx*h[i]*h[i]*thx);

    i = n-1;

    thx = 0.5*idx*(h1 - h[i-1]);
    tbx = 0.5*idx*(b1 - bed[i-1]);
    tbxx = idx*idx*(b1 - 2*bed[i] + bed[i-1]);
    D = h[i] + h[i]*thx*tbx + 0.5*h[i]*h[i]*tbxx + h[i]*tbx*tbx;
        
    a[i-1] = -ithree*idx*idx*h[i]*h[i]*h[i] + 0.5*idx*h[i]*h[i]*thx;
    b[i] = D + 2.0*ithree*idx*idx*h[i]*h[i]*h[i];
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

void getGfromu(double *h, double *u, double *bed, double u0, double u1, double h0, double h1, double b0, double b1, double dx , int n, double *G)
{
    double idx = 1.0 / dx;
    double ithree = 1.0 / 3.0;

    int i;
    double thx,tbx,tbxx,D;


    for (i =1;i < n-1 ; i++)
    {
        thx = 0.5*idx*(h[i+1] - h[i-1]);
        tbx = 0.5*idx*(bed[i+1] - bed[i-1]);
        tbxx = idx*idx*(bed[i+1] - 2*bed[i] + bed[i-1]);
        D = h[i] + h[i]*thx*tbx + 0.5*h[i]*h[i]*tbxx + h[i]*tbx*tbx;

        G[i] = (-ithree*idx*idx*h[i]*h[i]*h[i] + 0.5*idx*h[i]*h[i]*thx)*u[i-1] +
                (D + 2.0*ithree*idx*idx*h[i]*h[i]*h[i])*u[i] +
                (-ithree*idx*idx*h[i]*h[i]*h[i] - 0.5*idx*h[i]*h[i]*thx)*u[i+1];
    }

    //Boundaries
    i = 0;

    thx = 0.5*idx*(h[i+1] - h0);
    tbx = 0.5*idx*(bed[i+1] - b0);
    tbxx = idx*idx*(bed[i+1] - 2*bed[i] + b0);
    D = h[i] + h[i]*thx*tbx + 0.5*h[i]*h[i]*tbxx + h[i]*tbx*tbx;

    G[i] = (-ithree*idx*idx*h[i]*h[i]*h[i] + 0.5*idx*h[i]*h[i]*thx)*u0 +
            (D + 2.0*ithree*idx*idx*h[i]*h[i]*h[i])*u[i] +
            (-ithree*idx*idx*h[i]*h[i]*h[i] - 0.5*idx*h[i]*h[i]*thx)*u[i+1];

    i = n-1;

    thx = 0.5*idx*(h1 - h[i-1]);
    tbx = 0.5*idx*(b1 - bed[i-1]);
    tbxx = idx*idx*(b1 - 2*bed[i] + bed[i-1]);
    D = h[i] + h[i]*thx*tbx + 0.5*h[i]*h[i]*tbxx + h[i]*tbx*tbx;

    G[i] = (-ithree*idx*idx*h[i]*h[i]*h[i] + 0.5*idx*h[i]*h[i]*thx)*u[i-1] +
            (D + 2.0*ithree*idx*idx*h[i]*h[i]*h[i])*u[i] +
            (-ithree*idx*idx*h[i]*h[i]*h[i] - 0.5*idx*h[i]*h[i]*thx)*u1;
}

void conc(double *a , double *b, double *c,int n,int m ,int k, double *d)
{
    //replace with memcpy for performance?
    memcpy(d,a,n*sizeof(double));
    memcpy(d+n,b,m*sizeof(double));
    memcpy(d+n+m,c,k*sizeof(double));
}

void evolveBC(double *G, double *h, double *bed, double *h0, double *h1, double *u0, double *u1, double *b0, double *b1, double g, double dx, double dt, int nBC, int n, int nBCs,double *nG, double *nh, double *nu, double *nbed)
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
    double *bedb = malloc(nBC*sizeof(double));
    double *bede = malloc(nBC*sizeof(double));
    double *u = malloc(n*sizeof(double));
    //ADD BC we get h,G need to solve for u then add in boundaries for G using u and h
    getufromG(h,G,bed, u0[j], u1[0], h0[j], h1[0], b0[j], b1[0],dx,n,u);

    //front end
    //i keeps track of big array
    //j keeps track of small bc arrays
    //k keeps track of small new bc arrays
    int i = -1;
    int k = nBC -1;
    j = nBCs -1;


    double hx = 0.5*idx*(h[i+1] - h0[j-1]);
    double bx = 0.5*idx*(bed[i+1] - b0[j-1]);
    double bxx = idx*idx*(bed[i+1] -2*b0[j] + b0[j-1]);

    double D = h0[j] + h0[j]*hx*bx + 0.5*h0[j]*h0[j]*bxx + h0[j]*bx*bx;

    double ai =0.5*idx*h0[j]*h0[j]*hx - idx*idx*ithree*h0[j]*h0[j]*h0[j];
    double bi = D + 2*ithree*idx*idx*h0[j]*h0[j]*h0[j];
    double ci = -0.5*idx*h0[j]*h0[j]*hx - idx*idx*ithree*h0[j]*h0[j]*h0[j];
    Gb[k] = ai*u0[j-1]
                 + bi*u0[j]
                 + ci*u[i+1];
    ub[k] = u0[j];
    hb[k] = h0[j];
    bedb[k] = b0[j];


    for(k = k-1; k > -1 ; k--)
    {
        j--;
        hx = 0.5*idx*(h0[j+1] - h0[j-1]);
        bx = 0.5*idx*(b0[j+1] - b0[j-1]);
        bxx = idx*idx*(b0[j+1] -2*b0[j] + b0[j-1]);

        D = h0[j] + h0[j]*hx*bx + 0.5*h0[j]*h0[j]*bxx + h0[j]*bx*bx;
        Gb[k] = (0.5*idx*h0[j]*h0[j]*hx - idx*idx*ithree*h0[j]*h0[j]*h0[j])*u0[j-1]
                    + (D + 2*ithree*idx*idx*h0[j]*h0[j]*h0[j])*u0[j]
                    + (-0.5*idx*h0[j]*h0[j]*hx - idx*idx*ithree*h0[j]*h0[j]*h0[j])*u0[j+1];

        ub[k] = u0[j];
        hb[k] = h0[j];
        bedb[k] = b0[j];
    }

    //back end
    i = n;
    k = 0;
    j = 0;

    hx = 0.5*idx*(h1[j+1] - h[i-1]);
    bx = 0.5*idx*(b1[j+1] - bed[i-1]);
    bxx = idx*idx*(b1[j+1] -2*b1[j] + bed[i-1]);

    D = h1[j] + h1[j]*hx*bx + 0.5*h1[j]*h1[j]*bxx + h1[j]*bx*bx;
    Ge[k] = (0.5*idx*h1[j]*h1[j]*hx - idx*idx*ithree*h1[j]*h1[j]*h1[j])*u[i-1]
            + (D + 2*ithree*idx*idx*h1[j]*h1[j]*h1[j])*u1[j]
            + (-0.5*idx*h1[j]*h1[j]*hx - idx*idx*ithree*h1[j]*h1[j]*h1[j])*u1[j+1];

    ue[k] = u1[j];
    he[k] = h1[j];
    bede[k] = b1[j];

    for(k = k+1; k < nBC ; k++)
    {
        j++;
        hx = 0.5*idx*(h1[j+1] - h1[j-1]);
        bx = 0.5*idx*(b1[j+1] - b1[j-1]);
        bxx = idx*idx*(b1[j+1] -2*b1[j] + b1[j-1]);

        D = h1[j] + h1[j]*hx*bx + 0.5*h1[j]*h1[j]*bxx + h1[j]*bx*bx;
        Ge[k] = (0.5*idx*h1[j]*h1[j]*hx - idx*idx*ithree*h1[j]*h1[j]*h1[j])*u1[j-1]
                    + (D + 2*ithree*idx*idx*h1[j]*h1[j]*h1[j])*u1[j]
                    + (-0.5*idx*h1[j]*h1[j]*hx - idx*idx*ithree*h1[j]*h1[j]*h1[j])*u1[j+1];
        ue[k] = u1[j];
        he[k] = h1[j];
        bede[k] = b1[j];

    }

    //bring them all together

    conc(Gb,G,Ge,nBC,n,nBC,nG);
    conc(hb,h,he,nBC,n,nBC,nh);
    conc(ub,u,ue,nBC,n,nBC,nu);
    conc(bedb,bed,bede,nBC,n,nBC,nbed);

    free(Gb);
    free(Ge);
    free(hb);
    free(he);
    free(ub);
    free(ue);
    free(bedb);
    free(bede);
    free(u);    

}

void evolve(double *G, double *h, double *u, double *bed, double g, double dx, double dt, int nBC, int n,double *nG, double *nh, double theta)
{

    double tux,tbx,tbxx, srcr,srcc,srcl;
    //modifies nh and nG to give the new values of h and G after a single time step
    double idx = 1.0 / dx;
    double ithree = 1.0 / 3.0;
    int i = nBC - 1;

	double wim2 = h[i-2] + bed[i-2];
    double wim1 = h[i-1] + bed[i-1];
    double wi = h[i] + bed[i];
    double wip1 = h[i+1] + bed[i+1];
    double wip2 = h[i+2] + bed[i+2];
	double wip3 = h[i+3] + bed[i+3];

    //calculate gradients
    double dwib = (wi - wim1);
    double dwim = 0.5*(wip1 - wim1);
    double dwif = (wip1 - wi);

    double dhib = (h[i] - h[i-1]);
    double dhim = 0.5*(h[i+1] - h[i-1]);
    double dhif = (h[i+1] - h[i]);

    double duib = (u[i] - u[i-1]);
    double duim = 0.5*(u[i+1] - u[i-1]);
    double duif = (u[i+1] - u[i]);

    double dGib = (G[i] - G[i-1]);
    double dGim = 0.5*(G[i+1] - G[i-1]);
    double dGif = (G[i+1] - G[i]);


    //calculate values at right of i cell
    double dwi = minmod(theta*dwib, dwim, theta*dwif);
    double dhi = minmod(theta*dhib, dhim, theta*dhif);
    double dui = minmod(theta*duib, duim, theta*duif);
    double dGi = minmod(theta*dGib, dGim, theta*dGif);

    double hir = h[i] + 0.5*dhi;
    double wir = wi + 0.5*dwi;
    double Gir = G[i] + 0.5*dGi;
    double uir = u[i] + 0.5*dui;
    double bir = wir - hir;

    double hil = h[i] - 0.5*dhi;
    double wil = wi - 0.5*dwi;
    double Gil = G[i] - 0.5*dGi;
    double uil = u[i] - 0.5*dui;
    double bil = wil - hil;

    //calculate values at left of i+1 cell

    //calculate gradients
    double dwip1b = (wip1 - wi);
    double dwip1m = 0.5*(wip2 - wi);
    double dwip1f = (wip2 - wip1);

    double dhip1b = (h[i+1] - h[i]);
    double dhip1m = 0.5*(h[i+2] - h[i]);
    double dhip1f = (h[i+2] - h[i+1]);

    double duip1b = (u[i+1] - u[i]);
    double duip1m = 0.5*(u[i+2] - u[i]);
    double duip1f = (u[i+2] - u[i+1]);

    double dGip1b = (G[i+1] - G[i]);
    double dGip1m = 0.5*(G[i+2] - G[i]);
    double dGip1f = (G[i+2] - G[i+1]);

    double dwip1 = minmod(theta*dwip1b, dwip1m, theta*dwip1f);
    double dhip1 = minmod(theta*dhip1b, dhip1m, theta*dhip1f);
    double duip1 = minmod(theta*duip1b, duip1m, theta*duip1f);
    double dGip1 = minmod(theta*dGip1b, dGip1m, theta*dGip1f);

	//LEFT
    double hip1l = h[i+1] - 0.5*dhip1;
    double wip1l = wip1 - 0.5*dwip1;
    double Gip1l = G[i+1] - 0.5*dGip1;
    double uip1l = u[i+1] - 0.5*duip1;
    double bip1l = wip1l - hip1l;

	//RIGHT
    double hip1r = h[i+1] + 0.5*dhip1;
    double wip1r = wip1 + 0.5*dwip1;
    double Gip1r = G[i+1] + 0.5*dGip1;
    double uip1r = u[i+1] + 0.5*duip1;
    double bip1r = wip1r - hip1r;

	// want bim1r/uim1r and bip1l/uip1l

    //calculate gradients
    double dwim1b = (wim1 - wim2);
    double dwim1m = 0.5*(wi - wim2);
    double dwim1f = (wi - wim1);

    double dhim1b = (h[i-1] - h[i-2]);
    double dhim1m = 0.5*(h[i] - h[i-2]);
    double dhim1f = (h[i] - h[i-1]);

    double duim1b = (u[i-1] - u[i-2]);
    double duim1m = 0.5*(u[i] - u[i-2]);
    double duim1f = (u[i] - u[i-1]);

    double dGim1b = (G[i-1] - G[i-2]);
    double dGim1m = 0.5*(G[i] - G[i-2]);
    double dGim1f = (G[i] - G[i-1]);

    //calculate values at right of i cell
    double dwim1 = minmod(theta*dwim1b, dwim1m, theta*dwim1f);
    double dhim1 = minmod(theta*dhim1b, dhim1m, theta*dhim1f);
    double duim1 = minmod(theta*duim1b, duim1m, theta*duim1f);
    double dGim1 = minmod(theta*dGim1b, dGim1m, theta*dGim1f);

    double him1r = h[i-1] + 0.5*dhim1;
    double wim1r = wim1 + 0.5*dwim1;
	double Gim1r = G[i-1] + 0.5*dGim1;
    double uim1r = u[i-1] + 0.5*duim1;
    double bim1r = wim1r - him1r;

    //calculate gradients
    double dwip2b = (wip2 - wip1);
    double dwip2m = 0.5*(wip3 - wip1);
    double dwip2f = (wip3 - wip2);

    double dhip2b = (h[i+2] - h[i+1]);
    double dhip2m = 0.5*(h[i+3] - h[i+1]);
    double dhip2f = (h[i+3] - h[i+2]);

    double duip2b = (u[i+2] - u[i+1]);
    double duip2m = 0.5*(u[i+3] - u[i+1]);
    double duip2f = (u[i+3] - u[i+2]);

    double dGip2b = (G[i+2] - G[i+1]);
    double dGip2m = 0.5*(G[i+3] - G[i+1]);
    double dGip2f = (G[i+3] - G[i+2]);

    double dwip2 = minmod(theta*dwip2b, dwip2m, theta*dwip2f);
    double dhip2 = minmod(theta*dhip2b, dhip2m, theta*dhip2f);
    double duip2 = minmod(theta*duip2b, duip2m, theta*duip2f);
    double dGip2 = minmod(theta*dGip2b, dGip2m, theta*dGip2f);

    double hip2l = h[i+2] - 0.5*dhip2;
    double wip2l = wip2 - 0.5*dwip2;
    double Gip2l = G[i+2] - 0.5*dGip2;
    double uip2l = u[i+2] - 0.5*duip2;
    double bip2l = wip2l - hip2l;

    //right force
    double nbi  = fmax(bip1l,bir);
    double hihm = fmax(0, wir - nbi);
    double hihp = fmax(0, wip1l -nbi);

	//This is the only difference, change it
	double duer = idx*(uip2l - uip1l);
	double dber = idx*(bip2l - bip1l);
		        
	double duel = idx*(uir - uim1r);
	double dbel = idx*(bir - bim1r);

    double sqrtghel = sqrt(g* hihm);
    double sqrtgher = sqrt(g* hihp);

    double sl = fmin(0,fmin(uir - sqrtghel, uip1l - sqrtgher));
    double sr = fmax(0,fmax(uir + sqrtghel, uip1l + sqrtgher));

    double felh = uir*hihm;
    double felG = Gir*uir + 0.5*g*hihm*hihm - 2*ithree*hihm*hihm*hihm*duel*duel + hihm*hihm*uir*duel*dbel;
    double ferh = uip1l*hihp;
    double ferG = Gip1l*uip1l + 0.5*g*hihp*hihp -2*ithree*hihp*hihp*hihp*duer*duer + hihp*hihp*uip1l*duer*dber;

    double isrmsl = 0.0;

    if(sr != sl) isrmsl = 1.0 / (sr - sl);
   
    double foh = isrmsl*(sr*felh - sl*ferh + sl*sr*(hihp - hihm));

    double foG = isrmsl*(sr*felG - sl*ferG + sl*sr*(Gip1l - Gir));

    double fih = foh;
    double fiG = foG;
    double himhp = hihp;

    him1r = hir;
    wim1r = wir;
    uim1r = uir;
    bim1r = bir;

    hil = hip1l;
    wil = wip1l;
    bil = bip1l;
    Gil = Gip1l;
    uil = uip1l;

    hir = hip1r;
    wir = wip1r;
    bir = bip1r;
    Gir = Gip1r;
    uir = uip1r;

    hip1l = hip2l;
    wip1l = wip2l;
    uip1l = uip2l;
	Gip1l = Gip2l;
    bip1l = bip2l;

    dwip1 = dwip2;
    dhip1 = dhip2;
    duip1 = duip2;
    dGip1 = dGip2;
	

    for (i = nBC ; i < n +nBC;i++)
    {

		wim2 = h[i-2] + bed[i-2];
		wim1 = h[i-1] + bed[i-1];
		wi = h[i] + bed[i];
		wip1 = h[i+1] + bed[i+1];
		wip2 = h[i+2] + bed[i+2];
		wip3 = h[i+3] + bed[i+3];

        //i right

        hip1r = h[i+1] + 0.5*dhip1;
        wip1r = wip1 + 0.5*dwip1;
        Gip1r = G[i+1] + 0.5*dGip1;
        uip1r = u[i+1] + 0.5*duip1;
        bip1r = wip1r - hip1r;

        //calculate gradients
		/*
        dwip1b = (wip1 - wi);
        dwip1m = 0.5*(wip2 - wi);
        dwip1f = (wip2 - wip1);

        dhip1b = (h[i+1] - h[i]);
        dhip1m = 0.5*(h[i+2] - h[i]);
        dhip1f = (h[i+2] - h[i+1]);

        duip1b = (u[i+1] - u[i]);
        duip1m = 0.5*(u[i+2] - u[i]);
        duip1f = (u[i+2] - u[i+1]);

        dGip1b = (G[i+1] - G[i]);
        dGip1m = 0.5*(G[i+2] - G[i]);
        dGip1f = (G[i+2] - G[i+1]);

        dwip1 = minmod(theta*dwip1b, dwip1m, theta*dwip1f);
        dhip1 = minmod(theta*dhip1b, dhip1m, theta*dhip1f);
        duip1 = minmod(theta*duip1b, duip1m, theta*duip1f);
        dGip1 = minmod(theta*dGip1b, dGip1m, theta*dGip1f);

        hip1l = h[i+1] - 0.5*dhip1;
        wip1l = wip1 - 0.5*dwip1;
        Gip1l = G[i+1] - 0.5*dGip1;
        uip1l = u[i+1] - 0.5*duip1;
        bip1l = wip1l - hip1l;
		*/

		dwip2b = (wip2 - wip1);
		dwip2m = 0.5*(wip3 - wip1);
		dwip2f = (wip3 - wip2);

		dhip2b = (h[i+2] - h[i+1]);
		dhip2m = 0.5*(h[i+3] - h[i+1]);
		dhip2f = (h[i+3] - h[i+2]);

		duip2b = (u[i+2] - u[i+1]);
		duip2m = 0.5*(u[i+3] - u[i+1]);
		duip2f = (u[i+3] - u[i+2]);

		dGip2b = (G[i+2] - G[i+1]);
		dGip2m = 0.5*(G[i+3] - G[i+1]);
		dGip2f = (G[i+3] - G[i+2]);

		dwip2 = minmod(theta*dwip2b, dwip2m, theta*dwip2f);
		dhip2 = minmod(theta*dhip2b, dhip2m, theta*dhip2f);
		duip2 = minmod(theta*duip2b, duip2m, theta*duip2f);
		dGip2 = minmod(theta*dGip2b, dGip2m, theta*dGip2f);

		hip2l = h[i+2] - 0.5*dhip2;
		wip2l = wip2 - 0.5*dwip2;
		uip2l = u[i+2] - 0.5*duip2;
		Gip2l = G[i+2] - 0.5*dGip2;
		bip2l = wip2l - hip2l;

		//dGip1b = (G[i+1] - G[i]);
        //dGip1m = 0.5*(G[i+2] - G[i]);
        //dGip1f = (G[i+2] - G[i+1]);

        //dGip1 = minmod(theta*dGip1b, dGip1m, theta*dGip1f);

        //Gip1l = G[i+1] - 0.5*dGip1;


        nbi  = fmax(bip1l,bir);
        hihm = fmax(0, wir - nbi);
        hihp = fmax(0, wip1l -nbi);

		duer = idx*(uip2l - uip1l);
		dber = idx*(bip2l - bip1l);
		        
		duel = idx*(uir - uim1r);
		dbel = idx*(bir - bim1r);

        sqrtghel = sqrt(g*hihm);
        sqrtgher = sqrt(g*hihp);

        sl = fmin(0,fmin(uir - sqrtghel, uip1l - sqrtgher));
        sr = fmax(0,fmax(uir + sqrtghel, uip1l + sqrtgher));

        felh = uir*hihm;
        felG = Gir*uir + 0.5*g*hihm*hihm - 2*ithree*hihm*hihm*hihm*duel*duel + hihm*hihm*uir*duel*dbel;
        ferh = uip1l*hihp;
        ferG = Gip1l*uip1l + 0.5*g*hihp*hihp -2*ithree*hihp*hihp*hihp*duer*duer + hihp*hihp*uip1l*duer*dber;

        isrmsl = 0.0;

        if(sr != sl) isrmsl = 1.0 / (sr - sl);
   
        foh = isrmsl*(sr*felh - sl*ferh + sl*sr*(hihp - hihm));
        foG = isrmsl*(sr*felG - sl*ferG + sl*sr*(Gip1l - Gir));

        //source term
        //tbx = 0.5*idx*(bed[i+1] - bed[i-1]);
		tux = (uil - uir);
		tbx = (bil - bir);
        tbxx = idx*idx*(bed[i+1] -2*bed[i] + bed[i-1]);
        srcr = g*0.5*(hihm*hihm - hir*hir);
        srcc = g*h[i]*tbx + 0.5*h[i]*h[i]*u[i]*tux*tbxx -h[i]*u[i]*u[i]*tbx*tbxx;
        srcl =g*0.5*(hil*hil -himhp*himhp);

		//printf("G[%d] : %f ||", i, G[i]);

        nh[i -nBC] = h[i] -dt*idx*(foh - fih);
        nG[i -nBC] = G[i] -dt*idx*(foG -fiG) +dt*idx*(srcr + srcc + srcl);
		//printf("nG[%d] : %f \n", i, nG[i - nBC]);

        fih = foh;
        fiG = foG;
        himhp = hihp;

		him1r = hir;
		wim1r = wir;
		uim1r = uir;
		bim1r = bir;

		hil = hip1l;
		wil = wip1l;
		bil = bip1l;
		Gil = Gip1l;
		uil = uip1l;

		hir = hip1r;
		wir = wip1r;
		bir = bip1r;
		Gir = Gip1r;
		uir = uip1r;

		hip1l = hip2l;
		wip1l = wip2l;
		uip1l = uip2l;
		Gip1l = Gip2l;
		bip1l = bip2l;

		dwip1 = dwip2;
		dhip1 = dhip2;
		duip1 = duip2;
		dGip1 = dGip2;
 
    }
    
   
}

void evolvewrap(double *G, double *h, double *bed, double *h0, double *h1, double *u0, double *u1, double *b0, double *b1, double g, double dx, double dt, int nBC, int n, int nBCs, double theta)
{
    //first allocate memory for BC variables
    double *Gbc = malloc((n + 2*nBC)*sizeof(double));
    double *hbc = malloc((n + 2*nBC)*sizeof(double));
    double *ubc = malloc((n + 2*nBC)*sizeof(double));
    double *bedbc = malloc((n + 2*nBC)*sizeof(double));

    double *Gp = malloc(n*sizeof(double));
    double *hp = malloc(n*sizeof(double));

    evolveBC(G,h,bed,h0,h1,u0,u1,b0,b1,g,dx,dt,nBC,n,nBCs,Gbc,hbc,ubc,bedbc);
    evolve(Gbc,hbc,ubc,bedbc,g,dx,dt,nBC,n,Gp,hp,theta);

    evolveBC(Gp,hp,bed,h0,h1,u0,u1,b0,b1,g,dx,dt,nBC,n,nBCs,Gbc,hbc,ubc,bedbc);
    evolve(Gbc,hbc,ubc,bedbc,g,dx,dt,nBC,n,Gp,hp,theta);

    int i;
    for(i=0;i<n;i++)
    {
        G[i] = 0.5*(G[i] + Gp[i]);
        h[i] = 0.5*(h[i] + hp[i]);
    }

    free(Gbc);
    free(hbc);
    free(ubc);
    free(bedbc);
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
