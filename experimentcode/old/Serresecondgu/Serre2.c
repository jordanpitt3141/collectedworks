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

//GENERAL STUFF
void conc(double *a , double *b, double *c,int n,int m ,int k, double *d)
{
    //replace with memcpy for performance?
    memcpy(d,a,n*sizeof(double));
    memcpy(d+n,b,m*sizeof(double));
    memcpy(d+n+m,c,k*sizeof(double));
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


// APPLICATION SPECIFIC
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
	else
	{
        return 0.0;
	}
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
    
double GNacrosscell(double *x,double *h,double *u, double *b,double g,int j,double dx)
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
	double fgpbx = interpquarticval(bcoeff,x[j],fgp);
    
    double fgpe = fgph*fgpu*fgpu + 2*g*fgph*(0.5*fgph + fgpb) + i12*(fgph*fgph*fgph)*fgpux*fgpux 
				+ fgph*(fgpu*fgpbx - 0.5*fgpux*fgph)*(fgpu*fgpbx - 0.5*fgpux*fgph);
        
    //second gauss point
    double sgp = x[j];
    double sgph = interpquarticval(hcoeff,x[j],sgp);
    double sgpu = interpquarticval(ucoeff,x[j],sgp);
    double sgpux = interpquarticgrad(ucoeff,x[j],sgp);
	double sgpb = interpquarticval(bcoeff,x[j],sgp);
	double sgpbx = interpquarticval(bcoeff,x[j],sgp);
    
    double sgpe = sgph*sgpu*sgpu + 2*g*sgph*(0.5*sgph + sgpb) + i12*(sgph*sgph*sgph)*sgpux*sgpux 
				+ sgph*(sgpu*sgpbx - 0.5*sgpux*sgph)*(sgpu*sgpbx - 0.5*sgpux*sgph);

    //third gauss point
    double tgp = -0.5*dx*sqrt(3.0/5.0) + x[j];
    double tgph = interpquarticval(hcoeff,x[j],tgp);
    double tgpu = interpquarticval(ucoeff,x[j],tgp);
    double tgpux = interpquarticgrad(ucoeff,x[j],tgp);
	double tgpb = interpquarticval(bcoeff,x[j],tgp);
	double tgpbx = interpquarticval(bcoeff,x[j],tgp);
    
    double tgpe = tgph*tgpu*tgpu + 2*g*tgph*(0.5*tgph + tgpb) + i12*(tgph*tgph*tgph)*tgpux*tgpux 
				+ tgph*(tgpu*tgpbx - 0.5*tgpux*tgph)*(tgpu*tgpbx - 0.5*tgpux*tgph);

	free(ucoeff);
	free(hcoeff);
	free(bcoeff);
    
    return 0.5*dx*( (5.0/9.0)*fgpe + (8.0/9.0)*sgpe + (5.0/9.0)*tgpe);
}
    
double GNall(double *x,double *h,double *u, double *b,double g,int n, int nBC,double dx)
{
	//include approximations to H(a)u(a) + p(a)u(a) - H(b)u(b) - p(b)u(b) (a end, b beg)
    double sum1 = 0.0;
	int i;
	for(i = nBC; i < n - nBC;i++)
	{
       sum1 = sum1 + GNacrosscell(x,h,u,b,g,i,dx);
	}
	sum1 = 0.5*sum1;

    return sum1; 

}

void TDMA(double *a, double *b, double *c,double *d, int n, double *x)
{
	//tridiag matrix solver for [0abc0][x] = [d]
	double *alpha = malloc(n*sizeof(double));
	double *beta = malloc(n*sizeof(double));
	
	//int n will be size of b,d,x

	alpha[0] = c[0] / b[0];
	beta[0] = d[0] / b[0];
	int i;
	double m;
	for(i = 1; i < n -1;i++)
	{
		m = 1.0 / (b[i] - a[i-1]*alpha[i-1]);
		alpha[i] = c[i]*m;
		beta[i] = (d[i] - a[i-1]*beta[i-1])*m;
       
	}

	m = 1.0 / (b[n-1] - a[n-2]*alpha[n-2]);
	beta[n-1] = (d[n-1] - a[n-2]*beta[n-2])*m;

	x[n-1] = beta[n-1];

	for(i = n-2; i > -1; i--)
	{
		x[i] = beta[i] - alpha[i]*x[i+1];
       
	}

    free(alpha);
    free(beta);		
}

void getufromG(double *h, double *G,double *bed, double u0, double u1, double h0, double h1, double b0, double b1,double dx, int n,double *ublank)
{
	double idx = 1.0 / dx;
	double *a = malloc((n-1)*sizeof(double));
	double *b = malloc(n*sizeof(double));
	double *c = malloc((n-1)*sizeof(double));

	double th,thx,tbx,tbxx,D,ai,bi,ci;
	int i;
	for(i = 1; i < n-1; i++)
	{
		th = h[i];
		thx = 0.5*idx*(h[i+1] - h[i-1]);
		tbx = 0.5*idx*(bed[i+1] - bed[i-1]);
		tbxx = idx*idx*(bed[i+1] -2*bed[i]+ bed[i-1]);

		D = th + th*th*tbx + 0.5*th*th*tbxx + th*tbx*tbx;
		
		ai = -i3*idx*idx*th*th*th + 0.5*idx*th*th*thx;
		bi = D + 2.0*i3*idx*idx*th*th*th;
		ci = -i3*idx*idx*th*th*th - 0.5*idx*th*th*thx;

		a[i-1] = ai;
		b[i] = bi;
		c[i] = ci;
	}

	//Boundary 
	//i = 0
	i = 0;
	th = h[i];
	thx = 0.5*idx*(h[i+1] - h0);
	tbx = 0.5*idx*(bed[i+1] - b0);
	tbxx = idx*idx*(bed[i+1] -2*bed[i]+ b0);

	D = th + th*th*tbx + 0.5*th*th*tbxx + th*tbx*tbx;
	
	ai = -i3*idx*idx*th*th*th + 0.5*idx*th*th*thx;
	bi = D + 2.0*i3*idx*idx*th*th*th;
	ci = -i3*idx*idx*th*th*th - 0.5*idx*th*th*thx;

	b[i] = bi;
	c[i] = ci;	

	double tmpG1 = G[i];
	G[i] = tmpG1 - u0*ai;

	i = n-1;
	
	th = h[i];
	thx = 0.5*idx*(h1- h[i-1]);
	tbx = 0.5*idx*(b1 - bed[i-1]);
	tbxx = idx*idx*(b1 -2*bed[i]+ bed[i-1]);

	D = th + th*th*tbx + 0.5*th*th*tbxx + th*tbx*tbx;
	
	ai = -i3*idx*idx*th*th*th + 0.5*idx*th*th*thx;
	bi = D + 2.0*i3*idx*idx*th*th*th;
	ci = -i3*idx*idx*th*th*th - 0.5*idx*th*th*thx;

	a[i-1] = ai;
	b[i] = bi;	

	double tmpG2 = G[i];
	G[i] = tmpG2 - u1*ci;

	TDMA(a,b,c,G,n,ublank);
	
	G[0] = tmpG1;
	G[n-1] = tmpG2;


    free(a);
    free(b);
    free(c);
}

void getGfromu(double *h, double *u,double *bed, double u0, double u1, double h0, double h1, double b0, double b1,double dx,int n, double *Gblank )
{
	double idx = 1.0 / dx;

	double th,thx,tbx,tbxx,D,ai,bi,ci;
	int i;
	for(i = 1; i < n-1; i++)
	{
		th = h[i];
		thx = 0.5*idx*(h[i+1] - h[i-1]);
		tbx = 0.5*idx*(bed[i+1] - bed[i-1]);
		tbxx = idx*idx*(bed[i+1] -2*bed[i]+ bed[i-1]);

		D = th + th*th*tbx + 0.5*th*th*tbxx + th*tbx*tbx;
		
		ai = -i3*idx*idx*th*th*th + 0.5*idx*th*th*thx;
		bi = D + 2.0*i3*idx*idx*th*th*th;
		ci = -i3*idx*idx*th*th*th - 0.5*idx*th*th*thx;

		Gblank[i] = ai*u[i-1] + bi*u[i] + ci*u[i+1];
	}

	//Boundary 
	//i = 0
	i = 0;
	th = h[i];
	thx = 0.5*idx*(h[i+1] - h0);
	tbx = 0.5*idx*(bed[i+1] - b0);
	tbxx = idx*idx*(bed[i+1] -2*bed[i]+ b0);

	D = th + th*th*tbx + 0.5*th*th*tbxx + th*tbx*tbx;
	
	ai = -i3*idx*idx*th*th*th + 0.5*idx*th*th*thx;
	bi = D + 2.0*i3*idx*idx*th*th*th;
	ci = -i3*idx*idx*th*th*th - 0.5*idx*th*th*thx;

	Gblank[i] = ai*u0 + bi*u[i] + ci*u[i+1];

	//i = n-1
	i = n-1;
	
	th = h[i];
	thx = 0.5*idx*(h1- h[i-1]);
	tbx = 0.5*idx*(b1 - bed[i-1]);
	tbxx = idx*idx*(b1 -2*bed[i]+ bed[i-1]);

	D = th + th*th*tbx + 0.5*th*th*tbxx + th*tbx*tbx;
	
	ai = -i3*idx*idx*th*th*th + 0.5*idx*th*th*thx;
	bi = D + 2.0*i3*idx*idx*th*th*th;
	ci = -i3*idx*idx*th*th*th - 0.5*idx*th*th*thx;

	Gblank[i] = ai*u[i-1] + bi*u[i] + ci*u1;
}

void evolve(double *h,double *G, double *hblank,double *Gblank,double *bed,double g,double *u0,double *u1,double *h0,double *h1,double *b0,double *b1,double theta, int n, int nBC,double dx,double dt)
{
    
    //get averages
    double idx = 1.0 / dx;  

	double *u = malloc(n*sizeof(double));
	int i,j;
	int nBCn = 3;
    double th,hx,bx,bxx,D,ai,bi,ci,gb1,gb2,gb3,ge1,ge2,ge3;
    
    //calculate u at time step
	//getufromG(h,G,bed, u, u0[nBC - 1], u1[0], h0[nBC - 1], h1[0], b0[nBC - 1], b1[0], n,dx);
	// u seems to be the problem here
	getufromG(h,G,bed,u0[nBC - 1], u1[0], h0[nBC - 1], h1[0],b0[nBC - 1], b1[0],dx,n,u);
    
    //boundaries
    //beginning
    //i=-1
    i = -1;
    j = nBC-1;
    th = h0[j];
    hx = 0.5*idx*(h[i+1] - h0[j-1]);
    bx = 0.5*idx*(bed[i+1] - b0[j-1]);
    bxx =idx*idx*(bed[i+1] -2*b0[j] + b0[j-1]);
    
    D = th + th*hx*bx + 0.5*th*th*bxx + th*bx*bx;
    
    ai = th*th*hx*(0.5*idx) - i3*th*th*th*(idx*idx);
    bi = D + 2.0*i3*th*th*th*(idx*idx);
    ci = -1.0*th*th*hx*(0.5*idx) - i3*th*th*th*(idx*idx);
    
    gb3 = ai*u0[j-1] + bi*u0[j] +ci*u[i+1];
	//printf("uim1 : %f | ui : %f | uip1 : %f\n",u0[j-1],u0[j],u[i+1]);
	//printf("uip2 : %f | uip3 : %f | uip4 : %f\n",u[i+2],u[i+3], u[i+4]);
    
    //i=-2
    i = -2;
    j = nBC -2;
    th = h0[j];
    hx = 0.5*idx*(h0[j+1] - h0[j-1]);
    bx = 0.5*idx*(b0[j+1] - b0[j-1]);
    bxx = idx*idx*(b0[j+1] - 2*b0[j] + b0[j+1]);
    
    D = th + th*hx*bx + 0.5*th*th*bxx + th*bx*bx;
    
    ai = th*th*hx*(0.5*idx) - i3*th*th*th*(idx*idx);
    bi = D + 2.0*i3*th*th*th*(idx*idx);
    ci = -1.0*th*th*hx*(0.5*idx) - i3*th*th*th*(idx*idx);
    
    gb2 = ai*u0[j-1] + bi*u0[j] +ci*u0[j+1];
    
    //i=-3
    i = -3;
    j = nBC -3;
    th = h0[j];
    hx = 0.5*idx*(h0[j+1] - h0[j-1]);
    bx = 0.5*idx*(b0[j+1] - b0[j-1]);
    bxx = idx*idx*(b0[j+1] - 2*b0[j] + b0[j+1]);
    
    D = th + th*hx*bx + 0.5*th*th*bxx + th*bx*bx;
    
    ai = th*th*hx*(0.5*idx) - i3*th*th*th*(idx*idx);
    bi = D + 2.0*i3*th*th*th*(idx*idx);
    ci = -1.0*th*th*hx*(0.5*idx) - i3*th*th*th*(idx*idx);
    
    gb1 = ai*u0[j-1] + bi*u0[j] +ci*u0[j+1];
    
    //i = n
    i = n;
    j = 0;
    th = h1[j];
    hx = 0.5*idx*(h1[j+1] - h[i-1]);
    bx = 0.5*idx*(b1[j+1] - bed[i-1]);
    bxx = idx*idx*(b1[j+1] - 2*b1[j] + bed[i-1]);
    
    D = th + th*hx*bx + 0.5*th*th*bxx + th*bx*bx;
    
    ai = th*th*hx*(0.5*idx) - i3*th*th*th*(idx*idx);
    bi = D + 2.0*i3*th*th*th*(idx*idx);
    ci = -1.0*th*th*hx*(0.5*idx) - i3*th*th*th*(idx*idx);
    
    ge1 = ai*u[i-1] + bi*u1[j] + ci*u1[j+1];
    
    //i = n+1
    j = 1;
    th = h1[j];
    hx = 0.5*idx*(h1[j+1] - h1[j-1]);
    bx = 0.5*idx*(b1[j+1] - b1[j-1]);
    bxx = idx*idx*(b1[j+1] - 2*b1[j] + b1[j-1]);
    
    D = th + th*hx*bx + 0.5*th*th*bxx + th*bx*bx;
    
    ai = th*th*hx*(0.5*idx) - i3*th*th*th*(idx*idx);
    bi = D + 2.0*i3*th*th*th*(idx*idx);
    ci = -1.0*th*th*hx*(0.5*idx) - i3*th*th*th*(idx*idx);
    
    ge2 = ai*u1[j-1] + bi*u1[j] +ci*u1[j+1];
    
    //i = n+2    
    j = 2;
    th = h1[j];
    hx = 0.5*idx*(h1[j+1] - h1[j-1]);
    bx = 0.5*idx*(b1[j+1] - b1[j-1]);
    bxx = idx*idx*(b1[j+1] - 2*b1[j] + b1[j-1]);
    
    D = th + th*hx*bx + 0.5*th*th*bxx + th*bx*bx;
    
    ai = th*th*hx*(0.5*idx) - i3*th*th*th*(idx*idx);
    bi = D + 2.0*i3*th*th*th*(idx*idx);
    ci = -1.0*th*th*hx*(0.5*idx) - i3*th*th*th*(idx*idx);
    
    ge3 = ai*u1[j-1] + bi*u1[j] +ci*u1[j+1];
    
    

    double *ubeg = malloc(nBCn*sizeof(double));
    double *uend = malloc(nBCn*sizeof(double));
    double *bbeg = malloc(nBCn*sizeof(double));
    double *bend = malloc(nBCn*sizeof(double));
    double *hbeg = malloc(nBCn*sizeof(double));
    double *hend = malloc(nBCn*sizeof(double));
    double *gbeg = malloc(nBCn*sizeof(double));
    double *gend = malloc(nBCn*sizeof(double));


	ubeg[0] = u0[1];
	ubeg[1] = u0[2];
	ubeg[2] = u0[3];

	hbeg[0] = h0[1];
	hbeg[1] = h0[2];
	hbeg[2] = h0[3];

	bbeg[0] = b0[1];
	bbeg[1] = b0[2];
	bbeg[2] = b0[3];

	uend[0] = u1[0];
	uend[1] = u1[1];
	uend[2] = u1[2];

	hend[0] = h1[0];
	hend[1] = h1[1];
	hend[2] = h1[2];

	bend[0] = b1[0];
	bend[1] = b1[1];
	bend[2] = b1[2];

    gbeg[0] = gb1;
    gbeg[1] = gb2;
    gbeg[2] = gb3;
    
    gend[0] = ge1;
    gend[1] = ge2;
    gend[2] = ge3;

	//printf("gb1 : %f | gb2 : %f | gb3 : %f \n",gb1,gb2,gb3);
	//printf("ge1 : %f | ge2 : %f | ge3 : %f \n",ge1,ge2,ge3);
    
    double *Gbc = malloc((n + 2*nBCn)*sizeof(double));
    double *hbc = malloc((n + 2*nBCn)*sizeof(double));
    double *ubc = malloc((n + 2*nBCn)*sizeof(double));
    double *bedbc = malloc((n + 2*nBCn)*sizeof(double));

    conc(gbeg,G,gend,nBCn,n,nBCn,Gbc);
    conc(hbeg,h,hend,nBCn,n,nBCn,hbc);
    conc(bbeg,bed,bend,nBCn,n,nBCn,bedbc);
    conc(ubeg,u,uend,nBCn,n,nBCn,ubc);
        
    //do normal stuff 

    double wi,wip1,wip2,wip3,wim1,wim2,dwib,dwif,dwim,dhib,dhif,dhim,dGib,dGif,dGim;
    double duib,duif,duim, dwi,dhi, dGi, dui;
    double hir,wir,Gir,uir,bir,hil,wil,Gil,uil,bil;  
        
    //i = 2
    i = nBCn -1;
    //define the stage
    wi = hbc[i] + bedbc[i];
    wip1 = hbc[i+1] + bedbc[i+1];
    wip2 = hbc[i+2] + bedbc[i+2];
    wip3 = hbc[i+3] + bedbc[i+3];
    wim1 = hbc[i-1] + bedbc[i-1];
    wim2 = hbc[i-2] + bedbc[i-2];
        
    //reconstruct common values first
        
    //i left and right
        
    //gradients
    dwib = (wi - wim1);
    dwif = (wip1 - wi);
    dwim = 0.5*(wip1 - wim1);
    dhib = (hbc[i] - hbc[i-1]);
    dhif = (hbc[i+1] - hbc[i]);
    dhim = 0.5*(hbc[i+1] - hbc[i-1]);
    dGib = (Gbc[i] - Gbc[i-1]);
    dGif = (Gbc[i+1] - Gbc[i]);
    dGim = 0.5*(Gbc[i+1] - Gbc[i-1]);
    duib = (ubc[i] - ubc[i-1]);
    duif = (ubc[i+1] - ubc[i]);
    duim = 0.5*(ubc[i+1] - ubc[i-1]);
        
    //limiting
    dwi = minmod(theta*dwib,theta*dwif,dwim);
    dhi = minmod(theta*dhib,theta*dhif,dhim);
    dGi = minmod(theta*dGib,theta*dGif,dGim);
    dui = minmod(theta*duib,theta*duif,duim);
      
    //reconstruct right
    hir = hbc[i]+ 0.5*dhi;
    wir = wi + 0.5*dwi;
    Gir = Gbc[i] + 0.5*dGi;
    uir = ubc[i] + 0.5*dui;
    bir = wir - hir;
        
    //reconstruct left
    hil = hbc[i] - 0.5*dhi;
    wil = wi - 0.5*dwi;
    Gil = Gbc[i] - 0.5*dGi;
    uil = ubc[i] - 0.5*dui;
    bil = wil - hil;
        
    //only left of i+1 common but do both

    double dwip1b,dwip1f,dwip1m,dhip1b,dhip1f,dhip1m,dGip1b,dGip1f,dGip1m;
    double duip1b,duip1f,duip1m, dwip1,dhip1, dGip1, duip1;
    double hip1r,wip1r,Gip1r,uip1r,bip1r,hip1l,wip1l,Gip1l,uip1l,bip1l; 
        
    //gradients
    dwip1b = (wip1 - wi);
    dwip1f = (wip2 - wip1);
    dwip1m = 0.5*(wip2 - wi);
    dhip1b = (hbc[i+1] - hbc[i]);
    dhip1f = (hbc[i+2] - hbc[i+1]);
    dhip1m = 0.5*(hbc[i+2] - hbc[i]);
    dGip1b = (Gbc[i+1] - Gbc[i]);
    dGip1f = (Gbc[i+2] - Gbc[i+1]);
    dGip1m = 0.5*(Gbc[i+2] - Gbc[i]);
    duip1b = (ubc[i+1] - ubc[i]);
    duip1f = (ubc[i+2] - ubc[i+1]);
    duip1m = 0.5*(ubc[i+2] - ubc[i]);
        
    //limiting
    dwip1 = minmod(theta*dwip1b,theta*dwip1f,dwip1m);
    dhip1 = minmod(theta*dhip1b,theta*dhip1f,dhip1m);
    dGip1 = minmod(theta*dGip1b,theta*dGip1f,dGip1m);
    duip1 = minmod(theta*duip1b,theta*duip1f,duip1m);
        
    //reconstruct right
    hip1r = hbc[i+1] + 0.5*dhip1;
    wip1r = wip1 + 0.5*dwip1;
    Gip1r = Gbc[i+1] + 0.5*dGip1;
    uip1r = ubc[i+1] + 0.5*duip1;
    bip1r = wip1r - hip1r;
        
    //reconstruct left
    hip1l = hbc[i+1] - 0.5*dhip1;
    wip1l = wip1 - 0.5*dwip1;
    Gip1l = Gbc[i+1] - 0.5*dGip1;
    uip1l = ubc[i+1] - 0.5*duip1;
    bip1l = wip1l - hip1l;
        
        
    //only right of i-1
    //i-1  right

    double dwim1b,dwim1f,dwim1m,dhim1b,dhim1f,dhim1m,dGim1b,dGim1f,dGim1m;
    double duim1b,duim1f,duim1m, dwim1,dhim1, dGim1, duim1;
    double him1r,wim1r,Gim1r,uim1r,bim1r,him1l,wim1l,Gim1l,uim1l,bim1l; 
        
    //gradients
    dwim1b = (wim1 - wim2);
    dwim1f = (wi - wim1);
    dwim1m = 0.5*(wi - wim2);
    dhim1b = (hbc[i-1] - hbc[i-2]);
    dhim1f = (hbc[i] - hbc[i-1]);
    dhim1m = 0.5*(hbc[i] - hbc[i-2]);
    dGim1b = (Gbc[i-1] - Gbc[i-2]);
    dGim1f = (Gbc[i] - Gbc[i-1]);
    dGim1m = 0.5*(Gbc[i] - Gbc[i-2]);
    duim1b = (ubc[i-1] - ubc[i-2]);
    duim1f = (ubc[i] - ubc[i-1]);
    duim1m = 0.5*(ubc[i] - ubc[i-2]);
        
    //limiting
    dwim1 = minmod(theta*dwim1b,theta*dwim1f,dwim1m);
    dhim1 = minmod(theta*dhim1b,theta*dhim1f,dhim1m);
    dGim1 = minmod(theta*dGim1b,theta*dGim1f,dGim1m);
    duim1 = minmod(theta*duim1b,theta*duim1f,duim1m);
        
    //reconstruct right
    him1r = hbc[i-1] + 0.5*dhim1;
    wim1r = wim1 + 0.5*dwim1;
    Gim1r = Gbc[i-1] + 0.5*dGim1;
    uim1r = ubc[i-1] + 0.5*duim1;
    bim1r = wim1r - him1r;
        
    //reconstruct i+2 left
    double dwip2b,dwip2f,dwip2m,dhip2b,dhip2f,dhip2m,dGip2b,dGip2f,dGip2m;
    double duip2b,duip2f,duip2m,dwip2,dhip2, dGip2, duip2;
    double hip2r,wip2r,Gip2r,uip2r,bip2r,hip2l,wip2l,Gip2l,uip2l,bip2l; 
    
    //gradients
    dwip2b = (wip2 - wip1);
    dwip2f = (wip3 - wip2);
    dwip2m = 0.5*(wip3 - wip1);
    dhip2b = (hbc[i+2] - hbc[i+1]);
    dhip2f = (hbc[i+3] - hbc[i+2]);
    dhip2m = 0.5*(hbc[i+3] - hbc[i+1]);
    dGip2b = (Gbc[i+2] - Gbc[i+1]);
    dGip2f = (Gbc[i+3] - Gbc[i+2]);
    dGip2m = 0.5*(Gbc[i+3] - Gbc[i+1]);
    duip2b = (ubc[i+2] - ubc[i+1]);
    duip2f = (ubc[i+3] - ubc[i+2]);
    duip2m = 0.5*(ubc[i+3] - ubc[i+1]);
        
    //limiting
    dwip2 = minmod(theta*dwip2b,theta*dwip2f,dwip2m);
    dhip2 = minmod(theta*dhip2b,theta*dhip2f,dhip2m);
    dGip2 = minmod(theta*dGip2b,theta*dGip2f,dGip2m);
    duip2 = minmod(theta*duip2b,theta*duip2f,duip2m);
                
    //reconstruct left
    hip2l = hbc[i+2] - 0.5*dhip2;
    wip2l = wip2 - 0.5*dwip2;
    Gip2l = Gbc[i+2] - 0.5*dGip2;
    uip2l = ubc[i+2] - 0.5*duip2;
    bip2l = wip2l - hip2l;
    
    //calculate forces          
    double nbi,hihm,hihp,her,Ger,uer,hel,Gel,uel, himhp; 
    double duer,duel,dber,dbel,sqrtghel,sqrtgher,sl,sr,felh,felG,ferh,ferG,isrmsl;   
    double foh,foG,fih,fiG;
        
    //right force i
    nbi = fmax(bip1l,bir);
    hihm = fmax(0,wir-nbi);
    hihp = fmax(0,wip1l-nbi);

    her = hihp;
    Ger = Gip1l;
    uer = uip1l;
        
    hel = hihm;
    Gel = Gir;
    uel = uir;

	//do u like o2fix
	double ue = 0.5*(ubc[i+1]+ubc[i]);
	double due = idx*(ubc[i+1] - ubc[i]);

	double dbe = idx*(bedbc[i+1] - bedbc[i]);

	uer = ue;
	uel = ue;
	duer = due;
	duel = due;
	dber = dbe;
	dbel = dbe;

        
    //duer = idx*(uip2l - uip1l);
    dber = idx*(bip2l - bip1l);
            
    //duel = idx*(uir - uim1r);
    dbel = idx*(bir - bim1r);
        
    sqrtghel = sqrt(g*hel);
    sqrtgher = sqrt(g*her);
    sl = fmin(0, fmin(uel - sqrtghel, uer - sqrtgher));
    sr = fmax(0,fmax(uel + sqrtghel, uer + sqrtgher));
        
    felh = uel*hel;
    felG = Gel*uel + 0.5*g*hel*hel - 2*i3*hel*hel*hel*duel*duel + hel*hel*uel*duel*dbel;
    ferh = uer*her;
    ferG = Ger*uer + 0.5*g*her*her - 2*i3*her*her*her*duer*duer + her*her*uer*duer*dber;

	//printf("%d || hihm - hir : %e || hihp - wip1l : %e\n", i, hihm - hir,hihp - wip1l);
  
    isrmsl = 0.0;

    if(sr != sl) isrmsl = 1.0 / (sr - sl);

    foh = isrmsl*(sr*felh - sl*ferh + sl*sr*(her - hel ));
    foG = isrmsl*(sr*felG - sl*ferG + sl*sr*(Ger - Gel ));
    
    fih = foh;
    fiG = foG;
    himhp = hihp;
        
    him1r = hir;
    wim1r = wir;
    bim1r = bir;
    Gim1r = Gir;
    uim1r = uir;
        
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
    bip1l = bip2l;
    Gip1l = Gip2l;
    uip1l = uip2l;      
    
    dhip1 = dhip2;
    dwip1 = dwip2;
    duip1 = duip2;
    dGip1 = dGip2;
    
	for(i = nBCn; i < n + nBCn; i++)
    {
        //update both forces at same time

        //define the stage
        wi = hbc[i] + bedbc[i];
        wip1 = hbc[i+1]  + bedbc[i+1];
        wip2 = hbc[i+2]  + bedbc[i+2];
        wip3 = hbc[i+3] + bedbc[i+3];
        wim1 = hbc[i-1]  + bedbc[i-1];
        wim2 = hbc[i-2]  + bedbc[i-2];
        
        //reconstruct common values first
                
        //only left of i+1 common but do both   
        
        //reconstruct right
        hip1r = hbc[i+1] + 0.5*dhip1;
        wip1r = wip1 + 0.5*dwip1;
        Gip1r = Gbc[i+1] + 0.5*dGip1;
        uip1r = ubc[i+1] + 0.5*duip1;
        bip1r = wip1r - hip1r;
        
        
        //reconstruct i+2 left
        
        //gradients
        dwip2b = (wip2 - wip1);
        dwip2f = (wip3 - wip2);
        dwip2m = 0.5*(wip3 - wip1);
        dhip2b = (hbc[i+2] - hbc[i+1]);
        dhip2f = (hbc[i+3] - hbc[i+2]);
        dhip2m = 0.5*(hbc[i+3] - hbc[i+1]);
        dGip2b = (Gbc[i+2] - Gbc[i+1]);
        dGip2f = (Gbc[i+3] - Gbc[i+2]);
        dGip2m = 0.5*(Gbc[i+3] - Gbc[i+1]);
        duip2b = (ubc[i+2] - ubc[i+1]);
        duip2f = (ubc[i+3] - ubc[i+2]);
        duip2m = 0.5*(ubc[i+3] - ubc[i+1]);
        
        //limiting
        dwip2 = minmod(theta*dwip2b,theta*dwip2f,dwip2m);
        dhip2 = minmod(theta*dhip2b,theta*dhip2f,dhip2m);
        dGip2 = minmod(theta*dGip2b,theta*dGip2f,dGip2m);
        duip2 = minmod(theta*duip2b,theta*duip2f,duip2m);
                
        //reconstruct left
        hip2l = hbc[i+2] - 0.5*dhip2;
        wip2l = wip2 - 0.5*dwip2;
        Gip2l = Gbc[i+2] - 0.5*dGip2;
        uip2l = ubc[i+2] - 0.5*duip2;
        bip2l = wip2l - hip2l;


                
        //calculate forces              
        
        //right force i
        nbi = fmax(bip1l,bir);
        hihm = fmax(0.0,wir-nbi);
        hihp = fmax(0.0,wip1l-nbi);

		//printf("%d || hihm - hir : %e || hihp - wip1l : %e\n", i, hihm - hir,hihp - wip1l);

        her = hihp;
        Ger = Gip1l;
        uer = uip1l;
        
        hel = hihm;
        Gel = Gir;
        uel = uir;
        
        duer = idx*(uip2l - uip1l);
        dber = idx*(bip2l - bip1l);

            
        duel = idx*(uir - uim1r);
        dbel = idx*(bir - bim1r);

		ue = 0.5*(ubc[i+1]+ubc[i]);
		due = idx*(ubc[i+1] - ubc[i]);

		dbe = idx*(bedbc[i+1] - bedbc[i]);

		uer = ue;
		uel = ue;
		duer = due;
		duel = due;
		dber = dbe;
		dbel = dbe;
        
        sqrtghel = sqrt(g*hel);
        sqrtgher = sqrt(g*her);
        sl = fmin(0,fmin(uel - sqrtghel, uer - sqrtgher));
        sr = fmax(0,fmax(uel + sqrtghel, uer + sqrtgher));
        
        felh = uel*hel;
        felG = Gel*uel + 0.5*g*hel*hel - 2*i3*hel*hel*hel*duel*duel + hel*hel*uel*duel*dbel;
        ferh = uer*her;
        ferG = Ger*uer + 0.5*g*her*her - 2*i3*her*her*her*duer*duer + her*her*uer*duer*dber;

        isrmsl = 0.0;

        if(sr != sl) isrmsl = 1.0 / (sr - sl);
        foh = isrmsl*(sr*felh - sl*ferh + sl*sr*(her - hel ));
        foG = isrmsl*(sr*felG - sl*ferG + sl*sr*(Ger - Gel ));
        
        double th,tu,tux,tbx,tbxx, sourcer, sourcec, sourcel;
        //calculate the source term
        th = hbc[i];
        tu = ubc[i];
        tux = 0.5*idx*(ubc[i+1] - ubc[i-1]);
        tbx = 0.5*idx*(bedbc[i+1] - bedbc[i-1]);
        tbxx = idx*idx*(bedbc[i+1] - 2*bedbc[i] + bedbc[i-1]);
        
        sourcer = g*0.5*(hihm*hihm - hir*hir);
        sourcec = g*th*tbx +  0.5*th*th*tu*tux*tbxx - th*tu*tu*tbx*tbxx ;
        sourcel = g*0.5*(hil*hil - himhp*himhp);

		//printf("%d | Source: %e\n",i, dt*idx*(sourcer+sourcel + sourcec));
        
        hblank[i-nBCn] = hbc[i] - dt*idx*(foh - fih);
        Gblank[i-nBCn] = Gbc[i] - dt*idx*(foG - fiG) + dt*idx*(sourcer+sourcel + sourcec);
        
        fih = foh;
        fiG = foG;
        himhp = hihp;
        
        him1r = hir;
        wim1r = wir;
        bim1r = bir;
        Gim1r = Gir;
        uim1r = uir;
        
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
        bip1l = bip2l;
        Gip1l = Gip2l;
        uip1l = uip2l;      
        
    
        dhip1 = dhip2;
        dwip1 = dwip2;
        duip1 = duip2;
        dGip1 = dGip2;
    }
        
    free(u);
    free(ubeg);
    free(uend);
    free(bbeg);
    free(bend);
    free(hbeg);
    free(hend);
    free(gbeg);
    free(gend);
    free(Gbc);
    free(hbc);
    free(ubc);
    free(bedbc);

}

void evolvewrap(double *G,double *h,double *bed,double *h0,double *h1,double *u0,double *u1,double *b0,double *b1, double g,double dx,double dt, int n, int nBC,double theta)
{ 

//UNFINISHED~~~
    double *hp = malloc(n*sizeof(double));
    double *Gp = malloc(n*sizeof(double));
    double *hpp = malloc(n*sizeof(double));
    double *Gpp = malloc(n*sizeof(double));
    //update h' and G'    
    evolve(h,G, hp,Gp, bed,g,u0,u1,h0,h1,b0,b1,theta,n,nBC,dx,dt);
        
    //update h'' and G''
    evolve(hp,Gp, hpp,Gpp, bed,g,u0,u1,h0,h1,b0,b1,theta,n,nBC,dx,dt);

    int i;
    for(i=0;i<n;i++)
    {
		//printf("Before G[%d] : %f || h[%d] : %f \n",i, G[i],i, h[i]);
		//printf("Before Gp[%d] : %f || hp[%d] : %f \n",i, Gp[i],i, hp[i]);
		//printf("Before Gpp[%d] : %f || hpp[%d] : %f \n",i, Gpp[i],i, hpp[i]);
        G[i] = 0.5*(G[i] + Gpp[i]);
        h[i] = 0.5*(h[i] + hpp[i]);
		//printf("After G[%d] : %f || h[%d] : %f \n",i, G[i],i, h[i]);
    }
    
    free(hp);
    free(Gp);
    free(Gpp);
    free(hpp);
}

