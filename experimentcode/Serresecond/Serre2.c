#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

//put in slope limiting.

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

    double tbx,tbxx, srcr,srcc,srcl;
    //modifies nh and nG to give the new values of h and G after a single time step
    double idx = 1.0 / dx;
    double ithree = 1.0 / 3.0;
    int i = nBC - 1;

    double wim1 = h[i-1] + bed[i-1];
    double wi = h[i] + bed[i];
    double wip1 = h[i+1] + bed[i+1];
    double wip2 = h[i+2] + bed[i+2];

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

    double hip1l = h[i+1] - 0.5*dhip1;
    double wip1l = wip1 - 0.5*dwip1;
    double Gip1l = G[i+1] - 0.5*dGip1;
    double uip1l = u[i+1] - 0.5*duip1;
    double bip1l = wip1l - hip1l;

    //right force
    double nbi  = fmax(bip1l,bir);
    double hihm = fmax(0, wir - nbi);
    double hihp = fmax(0, wip1l -nbi);

    double duer = idx*(u[i+1] - u[i]);
    double dber = idx*(bed[i+1] - bed[i]);

    double sqrtghel = sqrt(g* hihm);
    double sqrtgher = sqrt(g* hihp);

    double sl = fmin(0,fmin(uir - sqrtghel, uip1l - sqrtgher));
    double sr = fmax(0,fmax(uir + sqrtghel, uip1l + sqrtgher));

    double felh = uir*hihm;
    double felG = Gir*uir + 0.5*g*hihm*hihm - 2*ithree*hihm*hihm*hihm*duer*duer + hihm*hihm*uir*duer*dber;
    double ferh = uip1l*hihp;
    double ferG = Gip1l*uip1l + 0.5*g*hihp*hihp -2*ithree*hihp*hihp*hihp*duer*duer + hihp*hihp*uip1l*duer*dber;

    double isrmsl = 0.0;

    if(sr != sl) isrmsl = 1.0 / (sr - sl);
   
    double foh = isrmsl*(sr*felh - sl*ferh + sl*sr*(hihp - hihm));

    double foG = isrmsl*(sr*felG - sl*ferG + sl*sr*(Gip1l - Gir));

    double fih = foh;
    double fiG = foG;
    double himhp = hihp;

    double hil = hip1l;

    dwi = dwip1;
    dhi = dhip1;
    dui = duip1;
    dGi = dGip1;

    for (i = nBC ; i < n +nBC;i++)
    {

        wim1 = h[i-1] + bed[i-1];
        wi = h[i] + bed[i];
        wip1 = h[i+1] + bed[i+1];
        wip2 = h[i+2] + bed[i+2];

        //i right

        hir = h[i] + 0.5*dhi;
        wir = wi + 0.5*dwi;
        Gir = G[i] + 0.5*dGi;
        uir = u[i] + 0.5*dui;
        bir = wir - hir;

        //calculate gradients
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


        nbi  = fmax(bip1l,bir);
        hihm = fmax(0, wir - nbi);
        hihp = fmax(0, wip1l -nbi);

        duer = idx*(u[i+1] - u[i]);
        dber = idx*(bed[i+1] - bed[i]);

        sqrtghel = sqrt(g*hihm);
        sqrtgher = sqrt(g*hihp);

        sl = fmin(0,fmin(uir - sqrtghel, uip1l - sqrtgher));
        sr = fmax(0,fmax(uir + sqrtghel, uip1l + sqrtgher));

        felh = uir*hihm;
        felG = Gir*uir + 0.5*g*hihm*hihm - 2*ithree*hihm*hihm*hihm*duer*duer + hihm*hihm*uir*duer*dber;
        ferh = uip1l*hihp;
        ferG = Gip1l*uip1l + 0.5*g*hihp*hihp -2*ithree*hihp*hihp*hihp*duer*duer + hihp*hihp*uip1l*duer*dber;

        isrmsl = 0.0;

        if(sr != sl) isrmsl = 1.0 / (sr - sl);
   
        foh = isrmsl*(sr*felh - sl*ferh + sl*sr*(hihp - hihm));
        foG = isrmsl*(sr*felG - sl*ferG + sl*sr*(Gip1l - Gir));

        //source term
        tbx = 0.5*idx*(bed[i+1] - bed[i-1]);
        tbxx = idx*idx*(bed[i+1] -2*bed[i] + bed[i-1]);
        srcr = g*0.5*(hihm*hihm - hir*hir);
        srcc = g*h[i]*tbx + 0.5*h[i]*h[i]*u[i]*dui*tbxx -h[i]*u[i]*u[i]*tbx*tbxx;
        srcl =g*0.5*(hil*hil -himhp*himhp);

        nh[i -nBC] = h[i] -dt*idx*(foh - fih);
        nG[i -nBC] = G[i] -dt*idx*(foG -fiG) +dt*idx*(srcr + srcc + srcl);

        fih = foh;
        fiG = foG;
        himhp = hihp;

        hil = hip1l;

        dwi = dwip1;
        dhi = dhip1;
        dui = duip1;
        dGi = dGip1;   
 
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
