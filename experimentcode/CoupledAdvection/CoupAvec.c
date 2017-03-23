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
      
void interpquartcoeff(double *q,double *coeff,int j,double dx)
{
    double idx = 1.0/dx;

    coeff[0] = i24*idx*idx*idx*idx*(q[j+2] - 4*q[j+1] + 6*q[j] - 4*q[j-1] + q[j-2]);
    coeff[1] = i12*idx*idx*idx*(q[j+2] - 2*q[j+1] + 2*q[j-1] - q[j-2]);
    coeff[2] = i24*idx*idx*(-q[j+2] + 16*q[j+1] - 30*q[j] + 16*q[j-1] - q[j-2]);
    coeff[3] = i12*idx*(-q[j+2] + 8*q[j+1] - 8*q[j-1] + q[j-2]);
    coeff[4] = q[j];
}
    
double HankEnergyacrosscell(double *x,double *u,int j,double dx)
{
    //so we have h,u at midpoints
    //epsilon and sigma are everywhere

	double *ucoeff = malloc(5*sizeof(double));

    //jth cell
    interpquartcoeff(u,ucoeff,j,dx);
    
    //first gauss point
    double fgp = 0.5*dx*sqrt(3.0/5.0) + x[j];
    double fgpu = interpquarticval(ucoeff,x[j],fgp);
    
    double fgpe = fgpu;
        
    //second gauss point
    double sgp = x[j];
    double sgpu = interpquarticval(ucoeff,x[j],sgp);
    
    double sgpe = sgpu;

    //third gauss point
    double tgp = -0.5*dx*sqrt(3.0/5.0) + x[j];
    double tgpu = interpquarticval(ucoeff,x[j],tgp);
    
    double tgpe = tgpu;

	free(ucoeff);
    
    return 0.5*dx*( (5.0/9.0)*fgpe + (8.0/9.0)*sgpe + (5.0/9.0)*tgpe);
}
    
double HankEnergyall(double *x,double *u,int n, int nBC,double dx)
{
    double sum1 = 0.0;
	int i;
	for(i = nBC; i < n - nBC;i++)
	{
       sum1 = sum1 + HankEnergyacrosscell(x,u,i,dx);
	}
    return sum1; 

}


void conc(double *a , double *b, double *c,int n,int m ,int k, double *d)
{
    //replace with memcpy for performance?
    memcpy(d,a,n*sizeof(double));
    memcpy(d+n,b,m*sizeof(double));
    memcpy(d+n+m,c,k*sizeof(double));
}

void evolveBC(double *u, double *v, double *u0, double *u1, double *v0, double *v1, double dx, double dt, int nBC, int n, int nBCs,double *nu, double *nv)
{ 
    //maybe an error in calculating G
    int j = nBCs -1;


    double *ub = malloc(nBC*sizeof(double));
    double *ue = malloc(nBC*sizeof(double));
    double *vb = malloc(nBC*sizeof(double));
    double *ve = malloc(nBC*sizeof(double));

    //front end
    //i keeps track of big array
    //j keeps track of small bc arrays
    //k keeps track of small new bc arrays
    int k = nBC -1;
    j = nBCs -1;


    ub[k] = u0[j];
    vb[k] = v0[j];


    for(k = k-1; k > -1 ; k--)
    {
        j--;

        ub[k] = u0[j];
        vb[k] = v0[j];
    }

    //back end
    k = 0;
    j = 0;


    ue[k] = u1[j];
    ve[k] = v1[j];

    for(k = k+1; k < nBC ; k++)
    {
        j++;
        ue[k] = u1[j];
        ve[k] = v1[j];

    }

    //bring them all together

    conc(vb,v,ve,nBC,n,nBC,nv);
    conc(ub,u,ue,nBC,n,nBC,nu);

    free(vb);
    free(ve);
    free(ub);
    free(ue);  

}

void evolve(double *u, double *v, double a, double b, double dx, double dt, int nBC, int n,double *nu, double *nv, double theta)
{
    //modifies nh and nG to give the new values of h and G after a single time step
    double idx = 1.0 / dx;
    int i = nBC - 1;

    //calculate gradients
    double duib = (u[i] - u[i-1]);
    double duim = 0.5*(u[i+1] - u[i-1]);
    double duif = (u[i+1] - u[i]);

    double dvib = (v[i] - v[i-1]);
    double dvim = 0.5*(v[i+1] - v[i-1]);
    double dvif = (v[i+1] - v[i]);


    //calculate values at right of i cell
    double dui = minmod(theta*duib, duim, theta*duif);
    double dvi = minmod(theta*dvib, dvim, theta*dvif);


    double uir = u[i] + 0.5*dui;
    double vir = v[i] + 0.5*dvi;

    //calculate values at left of i+1 cell

    //calculate gradients
    double duip1b = (u[i+1] - u[i]);
    double duip1m = 0.5*(u[i+2] -u[i]);
    double duip1f = (u[i+2] - u[i+1]);

    double dvip1b = (v[i+1] - v[i]);
    double dvip1m = 0.5*(v[i+2] - v[i]);
    double dvip1f = (v[i+2] - v[i+1]);

    double duip1 = minmod(theta*duip1b, duip1m, theta*duip1f);
    double dvip1 = minmod(theta*dvip1b, dvip1m, theta*dvip1f);

    double uip1l = u[i+1] - 0.5*duip1;
    double vip1l = v[i+1] - 0.5*dvip1;

    double sl = fmin(0,fmin(0, 0));
    double sr = fmax(0,fmax(b*uir + a*vir,b*uip1l + a*vip1l));

    double felu = a*uir*vir;
    double felv = b*uir*vir;
    double feru = a*uip1l*vip1l;
    double ferv = b*uip1l*vip1l;

    double isrmsl = 0.0;

    if(sr != sl) isrmsl = 1.0 / (sr - sl);
   
    double fou = isrmsl*(sr*felu - sl*feru + sl*sr*(uip1l - uir));

    double fov = isrmsl*(sr*felv - sl*ferv + sl*sr*(vip1l - vir));

    double fiu = fou;
    double fiv = fov;

    dui = duip1;
    //dui = duip1;
    dvi = dvip1;

    for (i = nBC ; i < n +nBC;i++)
    {

        uir = u[i] + 0.5*dui;
        vir = v[i] + 0.5*dvi;

        //calculate gradients
        duip1b = (u[i+1] - u[i]);
        duip1m = 0.5*(u[i+2] - u[i]);
        duip1f = (u[i+2] - u[i+1]);

        dvip1b = (v[i+1] - v[i]);
        dvip1m = 0.5*(v[i+2] - v[i]);
        dvip1f = (v[i+2] - v[i+1]);

        duip1 = minmod(theta*duip1b, duip1m, theta*duip1f);
        dvip1 = minmod(theta*dvip1b, dvip1m, theta*dvip1f);

        uip1l = u[i+1] - 0.5*duip1;
        vip1l = v[i+1] - 0.5*dvip1;

    	sl = fmin(0,fmin(0, 0));
    	sr = fmax(0,fmax(b*uir + a*vir,b*uip1l + a*vip1l));

    	felu = a*uir*vir;
    	felv = b*uir*vir;
    	feru = a*uip1l*vip1l;
    	ferv = b*uip1l*vip1l;

   	isrmsl = 0.0;

    	if(sr != sl) isrmsl = 1.0 / (sr - sl);
   
    	fou = isrmsl*(sr*felu - sl*feru + sl*sr*(uip1l - uir));

    	fov = isrmsl*(sr*felv - sl*ferv + sl*sr*(vip1l - vir));
        //source term

        nu[i -nBC] = u[i] -dt*idx*(fou - fiu);
        nv[i -nBC] = v[i] -dt*idx*(fov -fiv);

        fiu = fou;
        fiv = fov;

        dui = duip1;
        //dui = duip1;
        dvi = dvip1;   
 
    }
    
   
}

void evolveu(double *u, double *v, double *pu, double a, double dx, double dt, int nBC, int n,double *nu)
{
    //modifies nh and nG to give the new values of h and G after a single time step
    double idx = 1.0 / dx;
    int i = nBC - 1;
    double cux,cvx;

    //calculate gradients

    for (i = nBC ; i < n +nBC;i++)
    {
        cux = 0.5*idx*(u[i+1] - u[i-1]);
        cvx = 0.5*idx*(v[i+1] - v[i-1]);
        nu[i -nBC] = pu[i] - 2*dt*a*(v[i]*cux + u[i]*cvx); 
 
    }    
   
}


void evolvewrap(double *u, double *v, double *u0, double *u1, double *v0, double *v1, double a, double b, double dx, double dt, int nBC, int n, int nBCs, double theta)
{

    //first allocate memory for BC variables
    double *ubc = malloc((n + 2*nBC)*sizeof(double));
    double *vbc = malloc((n + 2*nBC)*sizeof(double));

    double *up = malloc(n*sizeof(double));
    double *vp = malloc(n*sizeof(double));


    evolveBC(u,v,u0,u1,v0,v1,dx,dt,nBC,n,nBCs,ubc,vbc);
    evolve(ubc,vbc,a,b,dx,dt,nBC,n,up,vp,theta);

    evolveBC(up,vp,u0,u1,v0,v1,dx,dt,nBC,n,nBCs,ubc,vbc);
    evolve(ubc,vbc,a,b,dx,dt,nBC,n,up,vp,theta);


    int i;
    for(i=0;i<n;i++)
    {
        u[i] = 0.5*(u[i] + up[i]);
        v[i] = 0.5*(v[i] + vp[i]);
    }

    free(ubc);
    free(vbc);
    free(up); 
    free(vp);
}

void evolvewrapFD(double *u, double *v, double *pubc, double *pvbc, double *u0, double *u1, double *v0, double *v1, double a, double b, double dx, double dt, int nBC, int n, int nBCs)
{

    //first allocate memory for BC variables
    double *ubc = malloc((n + 2*nBC)*sizeof(double));
    double *vbc = malloc((n + 2*nBC)*sizeof(double));



    evolveBC(u,v,u0,u1,v0,v1,dx,dt,nBC,n,nBCs,ubc,vbc);

    evolveu(ubc, vbc, pubc, a, dx,dt,nBC,n,u);
    evolveu(vbc, ubc, pvbc, b, dx,dt,nBC,n,v);

    memcpy(pubc,ubc,(n + 2*nBC)*sizeof(double));
    memcpy(pvbc,vbc,(n + 2*nBC)*sizeof(double));

    free(ubc);
    free(vbc);
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
