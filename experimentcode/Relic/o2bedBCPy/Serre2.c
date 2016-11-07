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
    int i;
    for(i = 0; i < n;i++)
	{
       d[i] = a[i];
	}
    for(i = n; i < n + m;i++)
	{
       d[i] = b[i-n];
	}
    for(i = n + m; i < n + m + k;i++)
	{
       d[i] = c[i-n - m];
	}
    //memcpy(d,a,n*sizeof(double));
    //memcpy(d+n,b,m*sizeof(double));
    //memcpy(d+n+m,c,k*sizeof(double));
}

/*
void conc(double *a , double *b, double *c,int n,int m ,int k, double *d)
{
    //replace with memcpy for performance?
    memcpy(d,a,n*sizeof(double));
    memcpy(d+n,b,m*sizeof(double));
    memcpy(d+n+m,c,k*sizeof(double));
}
*/
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

		D = th + th*thx*tbx + 0.5*th*th*tbxx + th*tbx*tbx;
		
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

	D = th + th*thx*tbx + 0.5*th*th*tbxx + th*tbx*tbx;
	
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

	D = th + th*thx*tbx + 0.5*th*th*tbxx + th*tbx*tbx;
	
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

void getufromGtankBC(double *h, double *G,double *bed,double u0, double u1, double h0, double h1, double b0, double b1,double g,double dx, int n,double *ublank)
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

		D = th + th*thx*tbx + 0.5*th*th*tbxx + th*tbx*tbx;
		
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

	D = th + th*thx*tbx + 0.5*th*th*tbxx + th*tbx*tbx;
	
	ai = -i3*idx*idx*th*th*th + 0.5*idx*th*th*thx;
	bi = D + 2.0*i3*idx*idx*th*th*th;
	ci = -i3*idx*idx*th*th*th - 0.5*idx*th*th*thx;

	b[i] = bi;
	c[i] = ci;	
	double tmpG1 = G[i];
	G[i] = tmpG1 - ai*u0;
	//printf("ai : %f | g : %f | h0 : %f | th : %f  \n",ai,g,h0,th);
	//printf("2*ai*(sqrt(g*h0) - sqrt(g*th)) :  %f \n",2*ai*(sqrt(g*h0) - sqrt(g*th)));

	i = n-1;
	
	th = h[i];
	thx = 0.5*idx*(h1- h[i-1]);
	tbx = 0.5*idx*(b1 - bed[i-1]);
	tbxx = idx*idx*(b1 -2*bed[i]+ bed[i-1]);

	D = th + th*thx*tbx + 0.5*th*th*tbxx + th*tbx*tbx;
	
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

void getufromGperiodic(double *h, double *G,double *bed,double dx, int n,double *ublank)
{

	double idx = 1.0 / dx;
	double *a = malloc((n-1)*sizeof(double));
	double *b = malloc((n)*sizeof(double));
	double *c = malloc((n-1)*sizeof(double));

	double th,thx,tbx,tbxx,D,ai,bi,ci,a1,cn,b1;
	int i,j;
	for(i = 1; i < n-1; i++)
	{
		th = h[i];
		thx = 0.5*idx*(h[i+1] - h[i-1]);
		tbx = 0.5*idx*(bed[i+1] - bed[i-1]);
		tbxx = idx*idx*(bed[i+1] -2*bed[i]+ bed[i-1]);

		D = th + th*thx*tbx + 0.5*th*th*tbxx + th*tbx*tbx;
		
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
	j = n-1;
	th = h[i];
	thx = 0.5*idx*(h[i+1] - h[j]);
	tbx = 0.5*idx*(bed[i+1] - bed[j]);
	tbxx = idx*idx*(bed[i+1] -2*bed[i]+ bed[j]);

	D = th + th*thx*tbx + 0.5*th*th*tbxx + th*tbx*tbx;
	
	ai = -i3*idx*idx*th*th*th + 0.5*idx*th*th*thx;
	bi = D + 2.0*i3*idx*idx*th*th*th;
	ci = -i3*idx*idx*th*th*th - 0.5*idx*th*th*thx;

	b[i] = 2*bi;
	c[i] = ci;
	b1 = bi;
	a1 = ai;

	//i = 0
	i = n-1;
	j = 0;
	th = h[i];
	thx = 0.5*idx*(h[j] - h[i-1]);
	tbx = 0.5*idx*(bed[j] - bed[i-1]);
	tbxx = idx*idx*(bed[j] -2*bed[i]+ bed[i-1]);

	D = th + th*thx*tbx + 0.5*th*th*tbxx + th*tbx*tbx;
	
	ai = -i3*idx*idx*th*th*th + 0.5*idx*th*th*thx;
	bi = D + 2.0*i3*idx*idx*th*th*th;
	ci = -i3*idx*idx*th*th*th - 0.5*idx*th*th*thx;

	a[i-1] = ai;
	b[i] = bi + ci*(a1/b1);
	cn = ci;
//good to here

double *u1 = malloc((n)*sizeof(double));
double *v1 = malloc((n)*sizeof(double));

	for(i = 0; i < n; i++)
	{
		u1[i] = 0;
		v1[i] = 0;
		
		
	}

	u1[0] = - b1;
	u1[n-1] =  cn;
	v1[0] = 1.0;
	v1[n-1] = - a1/b1;

	double *ap = malloc((n-1)*sizeof(double));
	double *bp = malloc((n)*sizeof(double));
	double *cp = malloc((n-1)*sizeof(double));	

	double *y1 = malloc((n)*sizeof(double));	
	double *q1 = malloc((n)*sizeof(double));

	memcpy(ap,a,(n-1)*sizeof(double));
	memcpy(bp,b,(n)*sizeof(double));
	memcpy(cp,c,(n-1)*sizeof(double));

	double *Gp = malloc((n)*sizeof(double));
	memcpy(Gp,G,(n)*sizeof(double));

	TDMA(a,b,c,Gp,n,y1);

	TDMA(ap,bp,cp,u1,n,q1);

	double dotprodvy, dotprodvq ;

	dotprodvy = 0;
	dotprodvq = 0;
	for(i = 0; i < n; i++)
	{
		dotprodvy = dotprodvy + v1[i]*y1[i];
		dotprodvq = dotprodvq + v1[i]*q1[i];
		
	}


	for(i = 0; i < n; i++)
	{
		ublank[i] = y1[i] - (dotprodvy / (1 + dotprodvq))*q1[i];
		
	}



    free(a);
    free(b);
    free(c);
    free(ap);
    free(bp);
    free(cp);
    free(y1);
    free(q1);
    free(v1);
	free(u1);
	free(Gp);
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

		D = th + th*thx*tbx + 0.5*th*th*tbxx + th*tbx*tbx;
		
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

	D = th + th*thx*tbx + 0.5*th*th*tbxx + th*tbx*tbx;
	
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

	D = th + th*thx*tbx + 0.5*th*th*tbxx + th*tbx*tbx;
	
	ai = -i3*idx*idx*th*th*th + 0.5*idx*th*th*thx;
	bi = D + 2.0*i3*idx*idx*th*th*th;
	ci = -i3*idx*idx*th*th*th - 0.5*idx*th*th*thx;

	Gblank[i] = ai*u[i-1] + bi*u[i] + ci*u1;
}


void evolveBC(double *h,double *G,double *bed,double *u0,double *u1,double *h0,double *h1 , double *G0, double *G1,double *b0,double *b1, int n, int nBC, int nBCn,double dx, double *nhbc, double *nGbc, double *nubc, double *nbbc)
{

    //get averages
	double *u = malloc(n*sizeof(double));
	int i,j;
    
	getufromG(h,G,bed,u0[nBC - 1], u1[0], h0[nBC - 1], h1[0],b0[nBC - 1], b1[0],dx,n,u);   
    

    double *ubeg = malloc(nBCn*sizeof(double));
    double *uend = malloc(nBCn*sizeof(double));
    double *bbeg = malloc(nBCn*sizeof(double));
    double *bend = malloc(nBCn*sizeof(double));
    double *hbeg = malloc(nBCn*sizeof(double));
    double *hend = malloc(nBCn*sizeof(double));
    double *gbeg = malloc(nBCn*sizeof(double));
    double *gend = malloc(nBCn*sizeof(double));

	for(i = 0 ; i < nBCn;i++)
	{

		ubeg[i] = u0[nBC - nBCn + i];
		hbeg[i] = h0[nBC - nBCn + i];
		gbeg[i] = G0[nBC - nBCn + i];
		bbeg[i] = b0[nBC - nBCn + i];

		uend[i] = u1[i];
		hend[i] = h1[i];
		gend[i] = G1[i];
		bend[i] = b1[i];

	}

	conc(hbeg,h,hend,nBCn,n,nBCn,nhbc);
	conc(ubeg,u,uend,nBCn,n,nBCn,nubc);
	conc(gbeg,G,gend,nBCn,n,nBCn,nGbc);
	conc(bbeg,bed,bend,nBCn,n,nBCn,nbbc);
    
	free(ubeg);
	free(hbeg);
	free(gbeg);
	free(bbeg);

	free(uend);
	free(hend);
	free(gend);
	free(bend);

	free(u);
}

void evolvewavetankBC(double *h,double *G,double *bed,double *h0,double *b0,double *u1,double *h1, double *G1,double *b1, double g, int n, int nBC, int nBCn,double dx, double *nhbc, double *nGbc, double *nubc, double *nbbc)
{

	double u0D = sqrt(g*h0[nBC - 1]);

    //get averages
	double *u = malloc(n*sizeof(double));
	int i,j;
    
	getufromGtankBC(h,G,bed,u0D, u1[0], h0[nBC - 1], h1[0],b0[nBC - 1], b1[0],g,dx,n,u);   
    

    double *ubeg = malloc(nBCn*sizeof(double));
    double *uend = malloc(nBCn*sizeof(double));
    double *bbeg = malloc(nBCn*sizeof(double));
    double *bend = malloc(nBCn*sizeof(double));
    double *hbeg = malloc(nBCn*sizeof(double));
    double *hend = malloc(nBCn*sizeof(double));
    double *gbeg = malloc(nBCn*sizeof(double));
    double *gend = malloc(nBCn*sizeof(double));

	double th, thx,tbx,tbxx,D,ai,bi,ci;
	double idx = 1.0 / dx;

	for(i = 0 ; i < nBCn;i++)
	{
		uend[i] = u1[i];
		hend[i] = h1[i];
		gend[i] = G1[i];
		bend[i] = b1[i];

		hbeg[i] = h0[nBC - nBCn + i];
		bbeg[i] = b0[nBC - nBCn + i];

	}
	for(i = nBCn - 1 ; i > -1;i--)
	{

		ubeg[i] = sqrt(g*hbeg[i]);

	}

	for(i = 1 ; i < nBCn - 1;i++)
	{

		th = hbeg[i];
		thx = 0.5*idx*(hbeg[i+1] - hbeg[i-1]);
		tbx = 0.5*idx*(bbeg[i+1] - bbeg[i-1]);
		tbxx = idx*idx*(bbeg[i+1] -2*bbeg[i]+ bbeg[i-1]);

		D = th + th*thx*tbx + 0.5*th*th*tbxx + th*tbx*tbx;
	
		ai = -i3*idx*idx*th*th*th + 0.5*idx*th*th*thx;
		bi = D + 2.0*i3*idx*idx*th*th*th;
		ci = -i3*idx*idx*th*th*th - 0.5*idx*th*th*thx;

		gbeg[i] = ai*ubeg[i-1] + bi*ubeg[i] + ci*ubeg[i+1];

	}

	//i = 0
	i = 0 ;

	th = hbeg[i];
	thx = 0.5*idx*(hbeg[i+1] - h0[nBC - nBCn -1]);
	tbx = 0.5*idx*(bbeg[i+1] - b0[nBC - nBCn -1]);
	tbxx = idx*idx*(bbeg[i+1] -2*bbeg[i]+ b0[nBC - nBCn -1]);

	D = th + th*thx*tbx + 0.5*th*th*tbxx + th*tbx*tbx;

	ai = -i3*idx*idx*th*th*th + 0.5*idx*th*th*thx;
	bi = D + 2.0*i3*idx*idx*th*th*th;
	ci = -i3*idx*idx*th*th*th - 0.5*idx*th*th*thx;

	gbeg[i] = ai*(sqrt(g*h0[nBC - nBCn -1])) + bi*ubeg[i] + ci*ubeg[i+1];

	//i = nBCn - 1 
	i = nBCn - 1;
	th = hbeg[i];
	thx = 0.5*idx*(h[0] - hbeg[i-1]);
	tbx = 0.5*idx*(bed[0] - bbeg[i-1]);
	tbxx = idx*idx*(bed[0] -2*bbeg[i]+ bbeg[i-1]);

	D = th + th*thx*tbx + 0.5*th*th*tbxx + th*tbx*tbx;

	ai = -i3*idx*idx*th*th*th + 0.5*idx*th*th*thx;
	bi = D + 2.0*i3*idx*idx*th*th*th;
	ci = -i3*idx*idx*th*th*th - 0.5*idx*th*th*thx;

	gbeg[i] = ai*ubeg[i-1] + bi*ubeg[i] + ci*u[0];

	conc(hbeg,h,hend,nBCn,n,nBCn,nhbc);
	conc(ubeg,u,uend,nBCn,n,nBCn,nubc);
	conc(gbeg,G,gend,nBCn,n,nBCn,nGbc);
	conc(bbeg,bed,bend,nBCn,n,nBCn,nbbc);
    
	free(ubeg);
	free(hbeg);
	free(gbeg);
	free(bbeg);

	free(uend);
	free(hend);
	free(gend);
	free(bend);

	free(u);
}


//Try this
void evolveperiodicBC(double *h,double *G,double *bed, int n, int nBCn,double dx, double *nhbc, double *nGbc, double *nubc, double *nbbc)
{

    //get averages
	double *u = malloc(n*sizeof(double));
	int i,j;
    
	getufromGperiodic(h,G,bed,dx,n,u);   
    

    double *ubeg = malloc(nBCn*sizeof(double));
    double *uend = malloc(nBCn*sizeof(double));
    double *bbeg = malloc(nBCn*sizeof(double));
    double *bend = malloc(nBCn*sizeof(double));
    double *hbeg = malloc(nBCn*sizeof(double));
    double *hend = malloc(nBCn*sizeof(double));
    double *gbeg = malloc(nBCn*sizeof(double));
    double *gend = malloc(nBCn*sizeof(double));

	for(i = 0 ; i < nBCn;i++)
	{

		ubeg[i] = u[n - nBCn + i];
		hbeg[i] = h[n - nBCn + i];
		gbeg[i] = G[n - nBCn + i];
		bbeg[i] = bed[n - nBCn + i];

		uend[i] = u[i];
		hend[i] = h[i];
		gend[i] = G[i];
		bend[i] = bed[i];

	}

	conc(hbeg,h,hend,nBCn,n,nBCn,nhbc);
	conc(ubeg,u,uend,nBCn,n,nBCn,nubc);
	conc(gbeg,G,gend,nBCn,n,nBCn,nGbc);
	conc(bbeg,bed,bend,nBCn,n,nBCn,nbbc);
    
	free(ubeg);
	free(hbeg);
	free(gbeg);
	free(bbeg);

	free(uend);
	free(hend);
	free(gend);
	free(bend);

	free(u);
}


void evolve(double *hbc,double *Gbc, double *hblank,double *Gblank,double *bedbc, double *ubc ,double g,double theta, int n, int nBCn,double dx,double dt)
{
    
    //get averages
    double idx = 1.0 / dx;  
	int i,j;

	//reconstruction at the edges
	double dhif,dhib,dhim,dhilim, hil,hir;
	double dGif,dGib,dGim,dGilim, Gil,Gir;
	double dwif,dwib,dwim,dwilim, wil,wir;
	double bil, bir;
	double wi, wip1,wim1;

	i = nBCn-1;

	// hil and hir
	dhif = idx*(hbc[i+1] - hbc[i]) ;
	dhib = idx*(hbc[i] - hbc[i-1]) ;
	dhim = 0.5*idx*(hbc[i+1] - hbc[i-1]) ;
	dhilim = minmod(theta*dhif,dhim,theta*dhib);
	hil = hbc[i] - 0.5*dx*dhilim;
	hir = hbc[i] + 0.5*dx*dhilim;

	// Gil and Gir
	dGif = idx*(Gbc[i+1] - Gbc[i]) ;
	dGib = idx*(Gbc[i] - Gbc[i-1]) ;
	dGim = 0.5*idx*(Gbc[i+1] - Gbc[i-1]) ;
	dGilim = minmod(theta*dGif,dGim,theta*dGib);
	Gil = Gbc[i] - 0.5*dx*dGilim;
	Gir = Gbc[i] + 0.5*dx*dGilim;

	// wil and wir
	wi = hbc[i] + bedbc[i];
	wip1 = hbc[i+1] + bedbc[i+1];
	wim1 = hbc[i-1] + bedbc[i-1];
	dwif = idx*(wip1 - wi) ;
	dwib = idx*(wi - wim1) ;
	dwim = 0.5*idx*(wip1 - wim1) ;
	dwilim = minmod(theta*dwif,dwim,theta*dwib);
	wil = wi - 0.5*dx*dwilim;
	wir = wi + 0.5*dx*dwilim;

	// bil and bir
	bil = wil - hil;
	bir = wir - hir;

	double dhip1f,dhip1b,dhip1m,dhip1lim, hip1l,hip1r;
	double dGip1f,dGip1b,dGip1m,dGip1lim, Gip1l,Gip1r;
	double dwip1f,dwip1b,dwip1m,dwip1lim, wip1l,wip1r;
	double bip1l, bip1r;
	double wip2;

	// hip1l and hip1r
	dhip1f = idx*(hbc[i+2] - hbc[i+1]) ;
	dhip1b = idx*(hbc[i+1] - hbc[i]) ;
	dhip1m = 0.5*idx*(hbc[i+2] - hbc[i]) ;
	dhip1lim = minmod(theta*dhip1f,dhip1m,theta*dhip1b);
	hip1l = hbc[i+1] - 0.5*dx*dhip1lim;
	hip1r = hbc[i+1] + 0.5*dx*dhip1lim;

	// Gip1l and Gip1r
	dGip1f = idx*(Gbc[i+2] - Gbc[i+1]) ;
	dGip1b = idx*(Gbc[i+1] - Gbc[i]) ;
	dGip1m = 0.5*idx*(Gbc[i+2] - Gbc[i]) ;
	dGip1lim = minmod(theta*dGip1f,dGip1m,theta*dGip1b);
	Gip1l = Gbc[i+1] - 0.5*dx*dGip1lim;
	Gip1r = Gbc[i+1] + 0.5*dx*dGip1lim;

	// wip1l and wip1r
	wip1 = hbc[i+1] + bedbc[i+1];
	wip2 = hbc[i+2] + bedbc[i+2];
	wi = hbc[i] + bedbc[i];
	dwip1f = idx*(wip2 - wip1) ;
	dwip1b = idx*(wip1 - wi) ;
	dwip1m = 0.5*idx*(wip2 - wi) ;
	dwip1lim = minmod(theta*dwip1f,dwip1m,theta*dwip1b);
	wip1l = wip1 - 0.5*dx*dwip1lim;
	wip1r = wip1 + 0.5*dx*dwip1lim;

	// bip1l and bip1r
	bip1l = wip1l - hip1l;
	bip1r = wip1r - hip1r;

	double dhip2f,dhip2b,dhip2m,dhip2lim, hip2l,hip2r;
	double dGip2f,dGip2b,dGip2m,dGip2lim, Gip2l,Gip2r;
	double dwip2f,dwip2b,dwip2m,dwip2lim, wip2l,wip2r;
	double bip2l, bip2r;
	double wip3;

	// hip2l and hip2r
	dhip2f = idx*(hbc[i+3] - hbc[i+2]) ;
	dhip2b = idx*(hbc[i+2] - hbc[i+1]) ;
	dhip2m = 0.5*idx*(hbc[i+3] - hbc[i+1]) ;
	dhip2lim = minmod(theta*dhip2f,dhip2m,theta*dhip2b);
	hip2l = hbc[i+2] - 0.5*dx*dhip2lim;
	hip2r = hbc[i+2] + 0.5*dx*dhip2lim;

	// Gip2l and Gip2r
	dGip2f = idx*(Gbc[i+3] - Gbc[i+2]) ;
	dGip2b = idx*(Gbc[i+2] - Gbc[i+1]) ;
	dGip2m = 0.5*idx*(Gbc[i+3] - Gbc[i+1]) ;
	dGip2lim = minmod(theta*dGip2f,dGip2m,theta*dGip2b);
	Gip2l = Gbc[i+2] - 0.5*dx*dGip2lim;
	Gip2r = Gbc[i+2] + 0.5*dx*dGip2lim;

	// wip2l and wip2r
	wip2 = hbc[i+2] + bedbc[i+2];
	wip3 = hbc[i+3] + bedbc[i+3];
	wip1 = hbc[i+1] + bedbc[i+1];
	dwip2f = idx*(wip3 - wip2) ;
	dwip2b = idx*(wip2 - wip1) ;
	dwip2m = 0.5*idx*(wip3 - wip1) ;
	dwip2lim = minmod(theta*dwip2f,dwip2m,theta*dwip2b);
	wip2l = wip2 - 0.5*dx*dwip2lim;
	wip2r = wip2 + 0.5*dx*dwip2lim;

	// bip2l and bip2r
	bip2l = wip2l - hip2l;
	bip2r = wip2r - hip2r;


	double dhim1f,dhim1b,dhim1m,dhim1lim, him1l,him1r;
	double dGim1f,dGim1b,dGim1m,dGim1lim, Gim1l,Gim1r;
	double dwim1f,dwim1b,dwim1m,dwim1lim, wim1l,wim1r;
	double bim1l, bim1r;
	double wim2;

	// him1l and him1r
	dhim1f = idx*(hbc[i] - hbc[i-1]) ;
	dhim1b = idx*(hbc[i-1] - hbc[i-2]) ;
	dhim1m = 0.5*idx*(hbc[i] - hbc[i-2]) ;
	dhim1lim = minmod(theta*dhim1f,dhim1m,theta*dhim1b);
	him1l = hbc[i-1] - 0.5*dx*dhim1lim;
	him1r = hbc[i-1] + 0.5*dx*dhim1lim;

	// Gim1l and Gim1r
	dGim1f = idx*(Gbc[i] - Gbc[i-1]) ;
	dGim1b = idx*(Gbc[i-1] - Gbc[i-2]) ;
	dGim1m = 0.5*idx*(Gbc[i] - Gbc[i-2]) ;
	dGim1lim = minmod(theta*dGim1f,dGim1m,theta*dGim1b);
	Gim1l = Gbc[i-1] - 0.5*dx*dGim1lim;
	Gim1r = Gbc[i-1] + 0.5*dx*dGim1lim;

	// wim1l and wim1r
	wim1 = hbc[i-1] + bedbc[i-1];
	wi = hbc[i] + bedbc[i];
	wim2 = hbc[i-2] + bedbc[i-2];
	dwim1f = idx*(wi - wim1) ;
	dwim1b = idx*(wim1 - wim2) ;
	dwim1m = 0.5*idx*(wi - wim2) ;
	dwim1lim = minmod(theta*dwim1f,dwim1m,theta*dwim1b);
	wim1l = wim1 - 0.5*dx*dwim1lim;
	wim1r = wim1 + 0.5*dx*dwim1lim;

	// bim1l and bim1r
	bim1l = wim1l - him1l;
	bim1r = wim1r - him1r;

	double nbi,hihm,hihp;
	double her, Ger, ber,dber, uer, duer, hel, Gel, bel,dbel, uel, duel;
	double fhel,fher,fGer, fGel,sqrtghel, sqrtgher, sl, sr,isrmsl,foh,foG;
	nbi = fmax(bip1l, bir);
	hihm = fmax(0, wir - nbi);
	hihp = fmax(0, wip1l - nbi);

	her = hihp;
	Ger = Gip1l;
	ber = bip1l;
	dber = idx*(bip2l - bip1l);
	uer  = 0.5*(ubc[i+1] + ubc[i]);
	duer = idx*(ubc[i+1] - ubc[i]);

	hel = hihm;
	Gel = Gir;
	bel = bir;
	dbel = idx*(bir - bim1r);
	uel  = 0.5*(ubc[i+1] + ubc[i]);
	duel = idx*(ubc[i+1] - ubc[i]);

	fhel = uel*hel;
	fher = uer*her;
	
	fGel = Gel*uel + 0.5*g*hel*hel - 2*i3*hel*hel*hel*duel*duel + hel*hel*uel*duel*dbel;
	fGer = Ger*uer + 0.5*g*her*her - 2*i3*her*her*her*duer*duer + her*her*uer*duer*dber;

    sqrtghel = sqrt(g* hel);
    sqrtgher = sqrt(g* her);

    sl = fmin(0,fmin(uel - sqrtghel, uer - sqrtgher));
    sr = fmax(0,fmax(uel + sqrtghel, uer + sqrtgher));

    isrmsl = 0.0;

    if(sr != sl) isrmsl = 1.0 / (sr - sl);	

	foh =isrmsl*( sr*fhel - sl*fher + sl*sr*(her - hel));
	foG =isrmsl*( sr*fGel - sl*fGer + sl*sr*(Ger - Gel));

	double fih,fiG,himhp;
	double th,tu,tux,tbx,tbxx,sourcer,sourcel,sourcec;
	fih = foh;
	fiG = foG;
	himhp = hihp;

	for(i = nBCn; i < n + nBCn; i++)
    {

		// hil and hir
		dhif = idx*(hbc[i+1] - hbc[i]) ;
		dhib = idx*(hbc[i] - hbc[i-1]) ;
		dhim = 0.5*idx*(hbc[i+1] - hbc[i-1]) ;
		dhilim = minmod(theta*dhif,dhim,theta*dhib);
		hil = hbc[i] - 0.5*dx*dhilim;
		hir = hbc[i] + 0.5*dx*dhilim;

		// Gil and Gir
		dGif = idx*(Gbc[i+1] - Gbc[i]) ;
		dGib = idx*(Gbc[i] - Gbc[i-1]) ;
		dGim = 0.5*idx*(Gbc[i+1] - Gbc[i-1]) ;
		dGilim = minmod(theta*dGif,dGim,theta*dGib);
		Gil = Gbc[i] - 0.5*dx*dGilim;
		Gir = Gbc[i] + 0.5*dx*dGilim;

		// wil and wir
		wi = hbc[i] + bedbc[i];
		wip1 = hbc[i+1] + bedbc[i+1];
		wim1 = hbc[i-1] + bedbc[i-1];
		dwif = idx*(wip1 - wi) ;
		dwib = idx*(wi - wim1) ;
		dwim = 0.5*idx*(wip1 - wim1) ;
		dwilim = minmod(theta*dwif,dwim,theta*dwib);
		wil = wi - 0.5*dx*dwilim;
		wir = wi + 0.5*dx*dwilim;

		// bil and bir
		bil = wil - hil;
		bir = wir - hir;


		// hip1l and hip1r
		dhip1f = idx*(hbc[i+2] - hbc[i+1]) ;
		dhip1b = idx*(hbc[i+1] - hbc[i]) ;
		dhip1m = 0.5*idx*(hbc[i+2] - hbc[i]) ;
		dhip1lim = minmod(theta*dhip1f,dhip1m,theta*dhip1b);
		hip1l = hbc[i+1] - 0.5*dx*dhip1lim;
		hip1r = hbc[i+1] + 0.5*dx*dhip1lim;

		// Gip1l and Gip1r
		dGip1f = idx*(Gbc[i+2] - Gbc[i+1]) ;
		dGip1b = idx*(Gbc[i+1] - Gbc[i]) ;
		dGip1m = 0.5*idx*(Gbc[i+2] - Gbc[i]) ;
		dGip1lim = minmod(theta*dGip1f,dGip1m,theta*dGip1b);
		Gip1l = Gbc[i+1] - 0.5*dx*dGip1lim;
		Gip1r = Gbc[i+1] + 0.5*dx*dGip1lim;

		// wip1l and wip1r
		wip1 = hbc[i+1] + bedbc[i+1];
		wip2 = hbc[i+2] + bedbc[i+2];
		wi = hbc[i] + bedbc[i];
		dwip1f = idx*(wip2 - wip1) ;
		dwip1b = idx*(wip1 - wi) ;
		dwip1m = 0.5*idx*(wip2 - wi) ;
		dwip1lim = minmod(theta*dwip1f,dwip1m,theta*dwip1b);
		wip1l = wip1 - 0.5*dx*dwip1lim;
		wip1r = wip1 + 0.5*dx*dwip1lim;

		// bip1l and bip1r
		bip1l = wip1l - hip1l;
		bip1r = wip1r - hip1r;


		// hip2l and hip2r
		dhip2f = idx*(hbc[i+3] - hbc[i+2]) ;
		dhip2b = idx*(hbc[i+2] - hbc[i+1]) ;
		dhip2m = 0.5*idx*(hbc[i+3] - hbc[i+1]) ;
		dhip2lim = minmod(theta*dhip2f,dhip2m,theta*dhip2b);
		hip2l = hbc[i+2] - 0.5*dx*dhip2lim;
		hip2r = hbc[i+2] + 0.5*dx*dhip2lim;

		// Gip2l and Gip2r
		dGip2f = idx*(Gbc[i+3] - Gbc[i+2]) ;
		dGip2b = idx*(Gbc[i+2] - Gbc[i+1]) ;
		dGip2m = 0.5*idx*(Gbc[i+3] - Gbc[i+1]) ;
		dGip2lim = minmod(theta*dGip2f,dGip2m,theta*dGip2b);
		Gip2l = Gbc[i+2] - 0.5*dx*dGip2lim;
		Gip2r = Gbc[i+2] + 0.5*dx*dGip2lim;

		// wip2l and wip2r
		wip2 = hbc[i+2] + bedbc[i+2];
		wip3 = hbc[i+3] + bedbc[i+3];
		wip1 = hbc[i+1] + bedbc[i+1];
		dwip2f = idx*(wip3 - wip2) ;
		dwip2b = idx*(wip2 - wip1) ;
		dwip2m = 0.5*idx*(wip3 - wip1) ;
		dwip2lim = minmod(theta*dwip2f,dwip2m,theta*dwip2b);
		wip2l = wip2 - 0.5*dx*dwip2lim;
		wip2r = wip2 + 0.5*dx*dwip2lim;

		// bip2l and bip2r
		bip2l = wip2l - hip2l;
		bip2r = wip2r - hip2r;

		// him1l and him1r
		dhim1f = idx*(hbc[i] - hbc[i-1]) ;
		dhim1b = idx*(hbc[i-1] - hbc[i-2]) ;
		dhim1m = 0.5*idx*(hbc[i] - hbc[i-2]) ;
		dhim1lim = minmod(theta*dhim1f,dhim1m,theta*dhim1b);
		him1l = hbc[i-1] - 0.5*dx*dhim1lim;
		him1r = hbc[i-1] + 0.5*dx*dhim1lim;

		// Gim1l and Gim1r
		dGim1f = idx*(Gbc[i] - Gbc[i-1]) ;
		dGim1b = idx*(Gbc[i-1] - Gbc[i-2]) ;
		dGim1m = 0.5*idx*(Gbc[i] - Gbc[i-2]) ;
		dGim1lim = minmod(theta*dGim1f,dGim1m,theta*dGim1b);
		Gim1l = Gbc[i-1] - 0.5*dx*dGim1lim;
		Gim1r = Gbc[i-1] + 0.5*dx*dGim1lim;

		// wim1l and wim1r
		wim1 = hbc[i-1] + bedbc[i-1];
		wi = hbc[i] + bedbc[i];
		wim2 = hbc[i-2] + bedbc[i-2];
		dwim1f = idx*(wi - wim1) ;
		dwim1b = idx*(wim1 - wim2) ;
		dwim1m = 0.5*idx*(wi - wim2) ;
		dwim1lim = minmod(theta*dwim1f,dwim1m,theta*dwim1b);
		wim1l = wim1 - 0.5*dx*dwim1lim;
		wim1r = wim1 + 0.5*dx*dwim1lim;

		// bim1l and bim1r
		bim1l = wim1l - him1l;
		bim1r = wim1r - him1r;

		nbi = fmax(bip1l, bir);
		hihm = fmax(0, wir - nbi);
		hihp = fmax(0, wip1l - nbi);

		her = hihp;
		Ger = Gip1l;
		ber = bip1l;
		dber = idx*(bip2l - bip1l);
		uer  = 0.5*(ubc[i+1] + ubc[i]);
		duer = idx*(ubc[i+1] - ubc[i]);

		hel = hihm;
		Gel = Gir;
		bel = bir;
		dbel = idx*(bir - bim1r);
		uel  = 0.5*(ubc[i+1] + ubc[i]);
		duel = idx*(ubc[i+1] - ubc[i]);

		fhel = uel*hel;
		fher = uer*her;
	
		fGel = Gel*uel + 0.5*g*hel*hel - 2*i3*hel*hel*hel*duel*duel + hel*hel*uel*duel*dbel;
		fGer = Ger*uer + 0.5*g*her*her - 2*i3*her*her*her*duer*duer + her*her*uer*duer*dber;

		sqrtghel = sqrt(g* hel);
		sqrtgher = sqrt(g* her);

		sl = fmin(0,fmin(uel - sqrtghel, uer - sqrtgher));
		sr = fmax(0,fmax(uel + sqrtghel, uer + sqrtgher));

		isrmsl = 0.0;

		if(sr != sl) isrmsl = 1.0 / (sr - sl);	

		foh =isrmsl*( sr*fhel - sl*fher + sl*sr*(her - hel));
		foG =isrmsl*( sr*fGel - sl*fGer + sl*sr*(Ger - Gel));

		th = hbc[i];
		tu = ubc[i];
		tux = -0.5*(ubc[i+1] - ubc[i-1]);
		tbx = (bil - bir);
		tbxx = idx*idx*(bedbc[i+1] - 2*bedbc[i] + bedbc[i-1]);
		
		sourcer = g*0.5*(hihm*hihm - hir*hir);
		sourcec = g*th*tbx +  0.5*th*th*tu*tux*tbxx - th*tu*tu*tbx*tbxx ;
		sourcel = g*0.5*(hil*hil - himhp*himhp);


	
		hblank[i - nBCn] = hbc[i] - dt*idx*(foh - fih);
		Gblank[i - nBCn] = Gbc[i] - dt*idx*(foG -fiG) + dt*idx*(sourcer+sourcel + sourcec);

		//printf("nhv : %2.12f  | nGv : %2.12f \n",nhv,nGv);

		fih = foh;
		fiG = foG;
		himhp = hihp;
	}
}


void evolvewrapperiodic(double *G,double *h,double *bed, double g,double dx,double dt, int n, int nBCn,double theta, double *hbc, double *Gbc, double *ubc)
{ 

//UNFINISHED~~~
    double *hp = malloc((n)*sizeof(double));
    double *Gp = malloc((n)*sizeof(double));
    double *hpp = malloc((n)*sizeof(double));
    double *Gpp = malloc((n)*sizeof(double));

	double *ihbc = malloc((n + 2*nBCn)*sizeof(double));
	double *iubc = malloc((n + 2*nBCn)*sizeof(double));
	double *iGbc = malloc((n + 2*nBCn)*sizeof(double));
	double *ibbc = malloc((n + 2*nBCn)*sizeof(double));

    //update h' and G'  

	evolveperiodicBC(h,G,bed,n,nBCn,dx,ihbc, iGbc, iubc, ibbc);
 
    evolve(ihbc,iGbc, hp,Gp,ibbc,iubc,g,theta,n,nBCn,dx,dt);

	evolveperiodicBC(hp,Gp,bed,n,nBCn,dx,ihbc, iGbc, iubc, ibbc);

	evolve(ihbc,iGbc, hpp,Gpp,ibbc,iubc,g,theta,n,nBCn,dx,dt);

	int i;
	for(i = 0; i < n ; i++)
	{
		h[i] =0.5*(h[i] +hpp[i]);
		G[i] =0.5*(G[i] +Gpp[i]);
	}
	
    
    free(hp);
    free(Gp);
    free(hpp);
    free(Gpp);
    free(ihbc);
    free(iubc);
    free(iGbc);
    free(ibbc);
}

void evolvewrapBC(double *G,double *h,double *bed,double *h0,double *h1,double *u0,double *u1, double *G0, double *G1,double *h0h,double *h1h,double *u0h,double *u1h, double *G0h, double *G1h,double *b0,double *b1, double g,double dx,double dt, int n, int nBC, int nBCn,double theta, double *hbc, double *Gbc, double *ubc)
{ 

//UNFINISHED~~~
    double *hp = malloc((n)*sizeof(double));
    double *Gp = malloc((n)*sizeof(double));
    double *hpp = malloc((n)*sizeof(double));
    double *Gpp = malloc((n)*sizeof(double));

	double *ihbc = malloc((n + 2*nBCn)*sizeof(double));
	double *iubc = malloc((n + 2*nBCn)*sizeof(double));
	double *iGbc = malloc((n + 2*nBCn)*sizeof(double));
	double *ibbc = malloc((n + 2*nBCn)*sizeof(double));

    //update h' and G' 
	evolveBC(h,G,bed,u0,u1,h0,h1,G0,G1,b0,b1,n,nBC,nBCn,dx,ihbc, iGbc, iubc, ibbc);
	evolve(ihbc,iGbc, hp,Gp,ibbc,iubc,g,theta,n,nBCn,dx,dt);


	evolveBC(hp,Gp,bed,u0h,u1h,h0h,h1h,G0h,G1h,b0,b1,n,nBC,nBCn,dx,ihbc, iGbc, iubc, ibbc);
	evolve(ihbc,iGbc, hpp,Gpp,ibbc,iubc,g,theta,n,nBCn,dx,dt);



	int i;
	for(i = 0; i < n ; i++)
	{
		h[i] =0.5*(h[i] +hpp[i]);
		G[i] =0.5*(G[i] +Gpp[i]);
	}
	
    
    free(hp);
    free(Gp);
    free(hpp);
    free(Gpp);
    free(ihbc);
    free(iubc);
    free(iGbc);
    free(ibbc);
}

void evolvewrapBCwavetank(double *G,double *h,double *bed,double *h0,double *h1,double *u1, double *G1,double *h0h,double *h1h,double *u1h, double *G1h,double *b0,double *b1, double g,double dx,double dt, int n, int nBC, int nBCn,double theta, double *hbc, double *Gbc, double *ubc)
{ 

//UNFINISHED~~~
    double *hp = malloc((n)*sizeof(double));
    double *Gp = malloc((n)*sizeof(double));
    double *hpp = malloc((n)*sizeof(double));
    double *Gpp = malloc((n)*sizeof(double));

	double *ihbc = malloc((n + 2*nBCn)*sizeof(double));
	double *iubc = malloc((n + 2*nBCn)*sizeof(double));
	double *iGbc = malloc((n + 2*nBCn)*sizeof(double));
	double *ibbc = malloc((n + 2*nBCn)*sizeof(double));

    //update h' and G' 
//evolvewavetankBC(double *h,double *G,double *bed,double *h0,double *b0,double *u1,double *h1, double *G1,double *b1, double g, int n, int nBC, int nBCn,double dx, double *nhbc, double *nGbc, double *nubc, double *nbbc)
	evolvewavetankBC(h,G,bed,h0,b0,u1,h1,G1,b1,g,n,nBC,nBCn,dx,ihbc, iGbc, iubc, ibbc);

	memcpy(hbc,ihbc,(n + 2*nBCn)*sizeof(double));
	memcpy(Gbc,iGbc,(n + 2*nBCn)*sizeof(double));
	memcpy(ubc,iubc,(n + 2*nBCn)*sizeof(double));
	evolve(ihbc,iGbc, hp,Gp,ibbc,iubc,g,theta,n,nBCn,dx,dt);


	evolvewavetankBC(hp,Gp,bed,h0h,b0,u1h,h1h,G1h,b1,g,n,nBC,nBCn,dx,ihbc, iGbc, iubc, ibbc);

	evolve(ihbc,iGbc, hpp,Gpp,ibbc,iubc,g,theta,n,nBCn,dx,dt);



	int i;
	for(i = 0; i < n ; i++)
	{
		h[i] =0.5*(h[i] +hpp[i]);
		G[i] =0.5*(G[i] +Gpp[i]);
	}
	
    
    free(hp);
    free(Gp);
    free(hpp);
    free(Gpp);
    free(ihbc);
    free(iubc);
    free(iGbc);
    free(ibbc);
}




