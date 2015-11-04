#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

//'Integrate' between the 3 cell values to get u bar

const double i3780 = 1.0/3780.0;
const double i420 = 1.0/420.0;
const double i30 = 1.0/30.0;
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

double phikm(double r)
{
    return fmax(0.,fmin(fmin(2.*r,i3*(1.+2.*r)),2.0));
}
double phikp(double r)
{
    return fmax(0.,fmin(fmin(2.*r,i3*(2.+r)),2.));
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
    int m = 2*n + 1;
    double idx = 1.0 / dx;

    double *a = malloc((m-2)*sizeof(double));
    double *b = malloc((m-1)*sizeof(double));
    double *c = malloc(m*sizeof(double));
    double *d = malloc((m-1)*sizeof(double));
    double *e = malloc((m-2)*sizeof(double));
    double *f = malloc(m*sizeof(double));
    double *ue = malloc(m*sizeof(double));

    double pja;


    int i,j;

    //have to clear out all so they are all zero
    for(i = 0; i < m-2 ; i++)
    {
        a[i] = 0.0;
        b[i] = 0.0;
        c[i] = 0.0;
        d[i] = 0.0;
        e[i] = 0.0;
        f[i] = 0.0;
    }

    //boundary zeros
    b[m-2] = 0.0;
    c[m-2] = 0.0;
    c[m-1] = 0.0;
    d[m-2] = 0.0;     
    f[m-2] = 0.0;
    f[m-1] = 0.0;

    j = 3;
    double Gci,hci,Gri,hri,Gir,hir,Gil,hil,Gim,him,p11,p12,p13,p21,p22,p23,p31,p32,p33;
    double q11,q12,q13,q21,q22,q23,q31,q32,q33,r1,r2,r3;
    for(i = nBC + 1; i < n + nBC - 1 ; i++)
    {
        // Get h and G over an element
        //increment
        
        //value that fits over cell and its neighbours
        Gci = G[i] - i24*(G[i+1] - 2*G[i] + G[i-1]);
        hci = h[i] - i24*(h[i+1] - 2*h[i] + h[i-1]);
        
        Gri = (G[i+1] - G[i]) / (G[i] - G[i-1]); 
        hri = (h[i+1] - h[i]) / (h[i] - h[i-1]); 

        Gir = G[i] + 0.5*phikm(Gri)*(G[i] - G[i-1]);
        hir = h[i] + 0.5*phikm(hri)*(h[i] - h[i-1]);
                
        //i - 1/2 +
        Gil = G[i] - 0.5*phikp(Gri)*(G[i] - G[i-1]);
        hil = h[i] - 0.5*phikp(hri)*(h[i] - h[i-1]);
        
        //i
        Gim = Gci;
        him = hci;

        
        //can now calculate Qe, Pe, Re for e = i
        p11 = (dx*i420)*(39*hil + 20*him - 3*hir);
        p12 = (dx*i420)*(20*hil + 16*him - 8*hir);
        p13 = (dx*i420)*(-3*hil - 8*him - 3*hir);
        p21 = (dx*i420)*(20*hil + 16*him - 8*hir);
        p22 = (dx*i420)*(16*hil + 192*him + 16*hir);
        p23 = (dx*i420)*(-8*hil + 16*him + 20*hir);   
        p31 = (dx*i420)*(-3*hil - 8*him - 3*hir);
        p32 = (dx*i420)*(-8*hil + 16*him + 20*hir);
        p33 = (dx*i420)*(-3*hil + 20*him + 39*hir);
        q11 = (i3780*idx)*(912*(hil*him*him) + 948*(hil*hil*him) + 61*(hir*hir*hir) 
                + 21*(hil*hir*hir) - 195*(hil*hil*hir) - 240*(him*him*hir) -336*(hil*him*hir) 
                + 832*(him*him*him) + 853*(hil*hil*hil) + 84*(him*hir*hir));
        q12 = (i3780 *idx)*(-240*(him*hir*hir) -284*(hir*hir*hir) -1104*(hil*hil*him) + 384*(hil*him*hir) 
                -1076*(hil*hil*hil) - 512*(him*him*him) - 960*(hil*him*him) + 228*(hil*hil*hir) 
                +12*(hil*hir*hir) + 192*(him*him*hir));
        q13 = (i3780 *idx)*(156*(hil*hil*him) - 33*(hil*hil*hir) - 33*(hil*hir*hir) + 156*(him*hir*hir)
                + 48*(him*him*hir) - 48*(hil*him*hir) + 223*(hir*hir*hir) - 320*(him*him*him) + 223*(hil*hil*hil) 
                + 48*(hil*him*him));
        q21 = (i3780 *idx)*(-240*(hir*him*hir) -284*(hir*hir*hir) -1104*(hil*hil*him) + 384*(hil*him*hir) 
                -1076*(hil*hil*hil) - 512*(him*him*him) - 960*(hil*him*him) + 228*(hil*hil*hir) 
                +12*(hil*hir*hir) + 192*(him*him*hir));
        q22 = (i3780 *idx)*(1344*(hil*hil*him) - 240*(hil*hil*hir) + 768*(hil*him*him) - 240*(hil*hir*hir) 
                + 768*(him*him*hir) + 1344*(him*hir*hir) + 1360*(hil*hil*hil) + 1024*(him*him*him) 
                + 1360*(hir*hir*hir) - 768*(hil*him*hir)) ;
        q23 = (i3780 *idx)*(12*(hil*hil*hir) - 240*(hil*hil*him) - 1104*(him*hir*hir) + 192*(hil*him*him) 
                - 960*(him*him*hir) -512*(him*him*him) - 1076*(hir*hir*hir) + 228*(hil*hir*hir) 
                - 284*(hil*hil*hil) + 384*(hil*him*hir));   
        q31 = (i3780 *idx)*(-48*(hil*him*hir) - 320*(him*him*him) - 33*(hil*hir*hir) + 48*(him*him*hir) 
                + 156*(him*hir*hir)+ 48*(hil*him*him) + 156*(hil*hil*him) - 33*(hil*hil*hir) + 223*(hil*hil*hil) 
                + 223*(hir*hir*hir));
        q32 = (i3780 *idx)*(12*(hil*hil*hir) - 240*(hil*hil*him) - 1104*(him*hir*hir) + 192*(hil*him*him) 
                - 960*(him*him*hir)-512*(him*him*him) - 1076*(hir*hir*hir) + 228*(hil*hir*hir) 
                - 284*(hil*hil*hil) + 384*(hil*him*hir));
        q33 = (i3780 *idx)*(-336*(hil*him*hir) - 195*(hil*hir*hir) + 912*(him*him*hir) + 948*(him*hir*hir) 
                + 832*(him*him*him) + 853*(hir*hir*hir) + 61*(hil*hil*hil) + 84*(hil*hil*him) 
                + 21*(hil*hil*hir) - 240*(hil*him*him));
        r1 = (dx*i30)*(4*Gil + 2*Gim - Gir);
        r2 = (dx*i30)*(2*Gil + 16*Gim + 2*Gir);
        r3 = (dx*i30)*(-Gil + 2*Gim + 4*Gir);
        
        //assembly
        f[j-1] = f[j-1] + r1; 
        f[j] = r2;
        f[j+1] = f[j+1] + r3;
        
        c[j-1] = c[j-1] + q11 + p11;
        c[j] =  q22 + p22;
        c[j+1] = c[j+1] + q33 + p33;
        
        //only change, so dont really need such an update
        d[j-1] = p12 + q12;
        d[j] = p23 + q23;
        
        e[j-1] = p13 + q13;
        
        //need to push back a and b because they are 'behind' by the counting system
        b[j-1] = p21 + q21;
        b[j] = p32 + q32;
        
        a[j-1] = p31 + q31;
        
        j = j + 2;      

    }

    i = nBC;
    j = 1;

    // Get h and G over an element
    //increment
        
    //value that fits over cell and its neighbours
    Gci = G[i] - i24*(G[i+1] - 2*G[i] + G[i-1]);
    hci = h[i] - i24*(h[i+1] - 2*h[i] + h[i-1]);
        
    Gri = (G[i+1] - G[i]) / (G[i] - G[i-1]); 
    hri = (h[i+1] - h[i]) / (h[i] - h[i-1]); 

    Gir = G[i] + 0.5*phikm(Gri)*(G[i] - G[i-1]);
    hir = h[i] + 0.5*phikm(hri)*(h[i] - h[i-1]);
                
    //i - 1/2 +
    Gil = G[i] - 0.5*phikp(Gri)*(G[i] - G[i-1]);
    hil = h[i] - 0.5*phikp(hri)*(h[i] - h[i-1]);
        
    //i
    Gim = Gci;
    him = hci;

        
    //can now calculate Qe, Pe, Re for e = i
    p11 = (dx*i420)*(39*hil + 20*him - 3*hir);
    p12 = (dx*i420)*(20*hil + 16*him - 8*hir);
    p13 = (dx*i420)*(-3*hil - 8*him - 3*hir);
    p21 = (dx*i420)*(20*hil + 16*him - 8*hir);
    p22 = (dx*i420)*(16*hil + 192*him + 16*hir);
    p23 = (dx*i420)*(-8*hil + 16*him + 20*hir) ;  
    p31 = (dx*i420)*(-3*hil - 8*him - 3*hir);
    p32 = (dx*i420)*(-8*hil + 16*him + 20*hir);
    p33 = (dx*i420)*(-3*hil + 20*him + 39*hir);
    q11 = (i3780*idx)*(912*(hil*him*him) + 948*(hil*hil*him) + 61*(hir*hir*hir) 
            + 21*(hil*hir*hir) - 195*(hil*hil*hir) - 240*(him*him*hir) -336*(hil*him*hir) 
            + 832*(him*him*him) + 853*(hil*hil*hil) + 84*(him*hir*hir));
    q12 = (i3780 *idx)*(-240*(him*hir*hir) -284*(hir*hir*hir) -1104*(hil*hil*him) + 384*(hil*him*hir) 
            - 1076*(hil*hil*hil) - 512*(him*him*him) - 960*(hil*him*him) + 228*(hil*hil*hir) 
            + 12*(hil*hir*hir) + 192*(him*him*hir));
    q13 = (i3780 *idx)*(156*(hil*hil*him) - 33*(hil*hil*hir) - 33*(hil*hir*hir) + 156*(him*hir*hir)
            + 48*(him*him*hir) - 48*(hil*him*hir) + 223*(hir*hir*hir) - 320*(him*him*him) + 223*(hil*hil*hil) 
            + 48*(hil*him*him));
    q21 = (i3780 *idx)*(-240*(hir*him*hir) -284*(hir*hir*hir) -1104*(hil*hil*him) + 384*(hil*him*hir) 
            -1076*(hil*hil*hil) - 512*(him*him*him) - 960*(hil*him*him) + 228*(hil*hil*hir) 
            +12*(hil*hir*hir) + 192*(him*him*hir));
    q22 = (i3780 *idx)*(1344*(hil*hil*him) - 240*(hil*hil*hir) + 768*(hil*him*him) - 240*(hil*hir*hir) 
            + 768*(him*him*hir) + 1344*(him*hir*hir) + 1360*(hil*hil*hil) + 1024*(him*him*him) 
            + 1360*(hir*hir*hir) - 768*(hil*him*hir)) ;
    q23 = (i3780 *idx)*(12*(hil*hil*hir) - 240*(hil*hil*him) - 1104*(him*hir*hir) + 192*(hil*him*him) 
            - 960*(him*him*hir) -512*(him*him*him) - 1076*(hir*hir*hir) + 228*(hil*hir*hir) 
            - 284*(hil*hil*hil) + 384*(hil*him*hir));   
    q31 = (i3780 *idx)*(-48*(hil*him*hir) - 320*(him*him*him) - 33*(hil*hir*hir) + 48*(him*him*hir) 
            + 156*(him*hir*hir)+ 48*(hil*him*him) + 156*(hil*hil*him) - 33*(hil*hil*hir) + 223*(hil*hil*hil) 
            + 223*(hir*hir*hir));
    q32 = (i3780 *idx)*(12*(hil*hil*hir) - 240*(hil*hil*him) - 1104*(him*hir*hir) + 192*(hil*him*him) 
            - 960*(him*him*hir)-512*(him*him*him) - 1076*(hir*hir*hir) + 228*(hil*hir*hir) 
            - 284*(hil*hil*hil) + 384*(hil*him*hir));
    q33 = (i3780 *idx)*(-336*(hil*him*hir) - 195*(hil*hir*hir) + 912*(him*him*hir) + 948*(him*hir*hir) 
            + 832*(him*him*him) + 853*(hir*hir*hir) + 61*(hil*hil*hil) + 84*(hil*hil*him) 
            + 21*(hil*hil*hir) - 240*(hil*him*him));
    r1 = (dx*i30)*(4*Gil + 2*Gim - Gir);
    r2 = (dx*i30)*(2*Gil + 16*Gim + 2*Gir);
    r3 = (dx*i30)*(-Gil + 2*Gim + 4*Gir);
        
    //assembly
    f[j-1] = ubeg[nBC-1]; 
    f[j] = r2;
    f[j+1] = f[j+1] + r3;
        
    c[j-1] = 1.0;
    c[j] =  q22 + p22;
    c[j+1] = c[j+1] + q33 + p33;
        
    //only change, so dont really need such an update
    d[j] = p23 + q23;
        
    //need to push back a and b because they are 'behind' by the counting system
    b[j-1] = p21 + q21;
    b[j] = p32 + q32;
        
    a[j-1] = p31 + q31;




    //END ############################
    i = n + nBC - 1;
    j = 2*n -1;

    // Get h and G over an element
    //increment
        
    //value that fits over cell and its neighbours
    Gci = G[i] - i24*(G[i+1] - 2*G[i] + G[i-1]);
    hci = h[i] - i24*(h[i+1] - 2*h[i] + h[i-1]);
        
    Gri = (G[i+1] - G[i]) / (G[i] - G[i-1]); 
    hri = (h[i+1] - h[i]) / (h[i] - h[i-1]); 

    Gir = G[i] + 0.5*phikm(Gri)*(G[i] - G[i-1]);
    hir = h[i] + 0.5*phikm(hri)*(h[i] - h[i-1]);
                
    //i - 1/2 +
    Gil = G[i] - 0.5*phikp(Gri)*(G[i] - G[i-1]);
    hil = h[i] - 0.5*phikp(hri)*(h[i] - h[i-1]);
        
    //i
    Gim = Gci;
    him = hci;

        
    //can now calculate Qe, Pe, Re for e = i
    p11 = (dx*i420)*(39*hil + 20*him - 3*hir);
    p12 = (dx*i420)*(20*hil + 16*him - 8*hir);
    p13 = (dx*i420)*(-3*hil - 8*him - 3*hir);
    p21 = (dx*i420)*(20*hil + 16*him - 8*hir);
    p22 = (dx*i420)*(16*hil + 192*him + 16*hir);
    p23 = (dx*i420)*(-8*hil + 16*him + 20*hir) ;  
    p31 = (dx*i420)*(-3*hil - 8*him - 3*hir);
    p32 = (dx*i420)*(-8*hil + 16*him + 20*hir);
    p33 = (dx*i420)*(-3*hil + 20*him + 39*hir);
    q11 = (i3780*idx)*(912*(hil*him*him) + 948*(hil*hil*him) + 61*(hir*hir*hir) 
            + 21*(hil*hir*hir) - 195*(hil*hil*hir) - 240*(him*him*hir) -336*(hil*him*hir) 
            + 832*(him*him*him) + 853*(hil*hil*hil) + 84*(him*hir*hir));
    q12 = (i3780 *idx)*(-240*(him*hir*hir) -284*(hir*hir*hir) -1104*(hil*hil*him) + 384*(hil*him*hir) 
            - 1076*(hil*hil*hil) - 512*(him*him*him) - 960*(hil*him*him) + 228*(hil*hil*hir) 
            + 12*(hil*hir*hir) + 192*(him*him*hir));
    q13 = (i3780 *idx)*(156*(hil*hil*him) - 33*(hil*hil*hir) - 33*(hil*hir*hir) + 156*(him*hir*hir)
            + 48*(him*him*hir) - 48*(hil*him*hir) + 223*(hir*hir*hir) - 320*(him*him*him) + 223*(hil*hil*hil) 
            + 48*(hil*him*him));
    q21 = (i3780 *idx)*(-240*(hir*him*hir) -284*(hir*hir*hir) -1104*(hil*hil*him) + 384*(hil*him*hir) 
            -1076*(hil*hil*hil) - 512*(him*him*him) - 960*(hil*him*him) + 228*(hil*hil*hir) 
            +12*(hil*hir*hir) + 192*(him*him*hir));
    q22 = (i3780 *idx)*(1344*(hil*hil*him) - 240*(hil*hil*hir) + 768*(hil*him*him) - 240*(hil*hir*hir) 
            + 768*(him*him*hir) + 1344*(him*hir*hir) + 1360*(hil*hil*hil) + 1024*(him*him*him) 
            + 1360*(hir*hir*hir) - 768*(hil*him*hir)) ;
    q23 = (i3780 *idx)*(12*(hil*hil*hir) - 240*(hil*hil*him) - 1104*(him*hir*hir) + 192*(hil*him*him) 
            - 960*(him*him*hir) -512*(him*him*him) - 1076*(hir*hir*hir) + 228*(hil*hir*hir) 
            - 284*(hil*hil*hil) + 384*(hil*him*hir));   
    q31 = (i3780 *idx)*(-48*(hil*him*hir) - 320*(him*him*him) - 33*(hil*hir*hir) + 48*(him*him*hir) 
            + 156*(him*hir*hir)+ 48*(hil*him*him) + 156*(hil*hil*him) - 33*(hil*hil*hir) + 223*(hil*hil*hil) 
            + 223*(hir*hir*hir));
    q32 = (i3780 *idx)*(12*(hil*hil*hir) - 240*(hil*hil*him) - 1104*(him*hir*hir) + 192*(hil*him*him) 
            - 960*(him*him*hir)-512*(him*him*him) - 1076*(hir*hir*hir) + 228*(hil*hir*hir) 
            - 284*(hil*hil*hil) + 384*(hil*him*hir));
    q33 = (i3780 *idx)*(-336*(hil*him*hir) - 195*(hil*hir*hir) + 912*(him*him*hir) + 948*(him*hir*hir) 
            + 832*(him*him*him) + 853*(hir*hir*hir) + 61*(hil*hil*hil) + 84*(hil*hil*him) 
            + 21*(hil*hil*hir) - 240*(hil*him*him));
    r1 = (dx*i30)*(4*Gil + 2*Gim - Gir);
    r2 = (dx*i30)*(2*Gil + 16*Gim + 2*Gir);
    r3 = (dx*i30)*(-Gil + 2*Gim + 4*Gir);
        
    //assembly
    f[j-1] = f[j-1] + r1 ; 
    f[j] = r2;
    f[j+1] = uend[0];
        
    c[j-1] = c[j-1] + q11 + p11;
    c[j] =  q22 + p22;
    c[j+1] = 1.0;
        
    //only change, so dont really need such an update
    d[j] = p23 + q23;
    d[j-1] = p12 + q12;
        
    //need to push back a and b because they are 'behind' by the counting system
    e[j-1] = p13 + q13;    
    

    PENT(a,b,c,d,e,f,m,ue); 
    j  = 0;
    for(i=1 ; i < m ; i=i+2)
    {
        pja = 2*(ue[i+1] -2*ue[i] + ue[i-1]); 
        u[j] = i12*pja*dx + ue[i]*dx;
        
        j++;
        //printf("nums = i : %d | n : %d | j : %d | m : %d \n",i,n,j,m);
        //printf("us = u[%d] : %e | ue[%d] : %e \n",j,u[j],i,ue[i]);

    }    

    free(a);
    free(b);
    free(c);
    free(d);
    free(e);
    free(f);
    free(ue);
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

    //need these to reconstruct the derivatives at cell edges
    double urim2 = (u[i-1] - u[i-2]) / (u[i-2] - u[i-3]);
    double uim2r = u[i-2] + 0.5*phikm(urim2)*(u[i-2] - u[i-3]);

    double urim1 = (u[i] - u[i-1]) / (u[i-1] - u[i-2]);
    double uim1r = u[i-1] + 0.5*phikm(urim1)*(u[i-1] - u[i-2]);

    //printf("C , ui2mr : %f, ui1mr : %f",uim2r, uim1r);

    //i right
    double Gri = (G[i+1] - G[i]) / (G[i] - G[i-1]); 
    double hri = (h[i+1] - h[i]) / (h[i] - h[i-1]); 
    double uri = (u[i+1] - u[i]) / (u[i] - u[i-1]); 

    double Gir = G[i] + 0.5*phikm(Gri)*(G[i] - G[i-1]);
    double hir = h[i] + 0.5*phikm(hri)*(h[i] - h[i-1]);
    double uir = u[i] + 0.5*phikm(uri)*(u[i] - u[i-1]);
    //printf("C , Gir : %f, hir : %f , uir : %f",Gir, hir,uir);

    // i+1 left
    double Grip1 = (G[i+2] - G[i+1]) / (G[i+1] - G[i]); 
    double hrip1 = (h[i+2] - h[i+1]) / (h[i+1] - h[i]); 
    double urip1 = (u[i+2] - u[i+1]) / (u[i+1] - u[i]);

    double Gip1l = G[i+1] - 0.5*phikp(Grip1)*(G[i+1] - G[i]);
    double hip1l = h[i+1] - 0.5*phikp(hrip1)*(h[i+1] - h[i]);
    double uip1l = u[i+1] - 0.5*phikp(urip1)*(u[i+1] - u[i]); 

    //for derivatives
    double urip2 = (u[i+3] - u[i+2]) / (u[i+2] - u[i+1]);
    double uip2l = u[i+2] - 0.5*phikp(urip2)*(u[i+2] - u[i+1]); 

    double urip3 = (u[i+4] - u[i+3]) / (u[i+3] - u[i+2]);
    double uip3l = u[i+3] - 0.5*phikp(urip3)*(u[i+3] - u[i+2]);
    //printf("C , uip3l : %f, uip2l : %f",uip3l, uip2l); 



    double duer = 0.5*idx*(-uip3l + 4*uip2l - 3*uip1l);
    double duel = 0.5*idx*(3*uir - 4*uim1r + uim2r);

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

    for (i = nBC; i < n + nBC;i++)
    {
        //need these to reconstruct the derivatives at cell edges
        urim2 = (u[i-1] - u[i-2]) / (u[i-2] - u[i-3]);
        uim2r = u[i-2] + 0.5*phikm(urim2)*(u[i-2] - u[i-3]);

        urim1 = (u[i] - u[i-1]) / (u[i-1] - u[i-2]);
        uim1r = u[i-1] + 0.5*phikm(urim1)*(u[i-1] - u[i-2]);

        //printf("C , ui2mr : %f, ui1mr : %f",uim2r, uim1r);

        //i right
        Gri = (G[i+1] - G[i]) / (G[i] - G[i-1]); 
        hri = (h[i+1] - h[i]) / (h[i] - h[i-1]); 
        uri = (u[i+1] - u[i]) / (u[i] - u[i-1]); 

        Gir = G[i] + 0.5*phikm(Gri)*(G[i] - G[i-1]);
        hir = h[i] + 0.5*phikm(hri)*(h[i] - h[i-1]);
        uir = u[i] + 0.5*phikm(uri)*(u[i] - u[i-1]);
        //printf("C , Gir : %f, hir : %f , uir : %f",Gir, hir,uir);

        // i+1 left
        Grip1 = (G[i+2] - G[i+1]) / (G[i+1] - G[i]); 
        hrip1 = (h[i+2] - h[i+1]) / (h[i+1] - h[i]); 
        urip1 = (u[i+2] - u[i+1]) / (u[i+1] - u[i]);

        Gip1l = G[i+1] - 0.5*phikp(Grip1)*(G[i+1] - G[i]);
        hip1l = h[i+1] - 0.5*phikp(hrip1)*(h[i+1] - h[i]);
        uip1l = u[i+1] - 0.5*phikp(urip1)*(u[i+1] - u[i]); 

        //for derivatives
        urip2 = (u[i+3] - u[i+2]) / (u[i+2] - u[i+1]);
        uip2l = u[i+2] - 0.5*phikp(urip2)*(u[i+2] - u[i+1]); 

        urip3 = (u[i+4] - u[i+3]) / (u[i+3] - u[i+2]);
        uip3l = u[i+3] - 0.5*phikp(urip3)*(u[i+3] - u[i+2]);
        //printf("C , uip3l : %f, uip2l : %f",uip3l, uip2l); 

        duer = 0.5*idx*(-uip3l + 4*uip2l - 3*uip1l);
        duel = 0.5*idx*(3*uir - 4*uim1r + uim2r);

        //printf("duer : %.8f, duel : %.8f\n",duer,duel);
        //printf("C : %f %f %f \n",uir,uim1r,uim2r);

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

    }
}      

void evolvewrap(double *Ga, double *ha, double *Gabeg, double *Gaend, double *habeg, double *haend, double *hmbeg, double *hmend, double *uabeg, double *uaend, double *umbeg, double *umend, int nfcBC, int nGsBC, double g, double dx, double dt, int n, int nBCa, int nBCm)
{
//again errors at machine precision, result of the division handling?
//############################### FIRST ITERATION #######################################
    double *um = malloc((n)*sizeof(double));
    double *ua = malloc((n)*sizeof(double));
    int cnBC;

    cnBC = nfcBC;

    double *Gabc = malloc((n + 2*cnBC)*sizeof(double));
    double *habc = malloc((n + 2*cnBC)*sizeof(double));
    double *uabc = malloc((n + 2*cnBC)*sizeof(double));

    conc(Gabeg+(nBCa - cnBC), Ga, Gaend,cnBC,n,cnBC,Gabc);
    conc(habeg+(nBCa - cnBC), ha, haend,cnBC,n,cnBC,habc);

    ufromGh(Gabc,habc,hmbeg+(nBCm -cnBC),hmend,umbeg+(nBCm -cnBC),umend,dx,n,cnBC,ua);

    conc(uabeg+(nBCa - cnBC), ua, uaend,cnBC,n,cnBC,uabc);

    double *nGa = malloc(n*sizeof(double));
    double *nha = malloc(n*sizeof(double));

    evolve(Gabc,habc,uabc,g,dx,dt,n,cnBC,nha,nGa);

//######################################### SECOND ITERATION #############################

    cnBC = nfcBC;

    conc(Gabeg+(nBCa - cnBC), nGa, Gaend,cnBC,n,cnBC,Gabc);
    conc(habeg+(nBCa - cnBC), nha, haend,cnBC,n,cnBC,habc);

    ufromGh(Gabc,habc,hmbeg+(nBCm -cnBC),hmend,umbeg+(nBCm -cnBC),umend,dx,n,cnBC,ua);

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

    cnBC = nfcBC;

    conc(Gabeg+(nBCa - cnBC), nGapp, Gaend,cnBC,n,cnBC,Gabc);
    conc(habeg+(nBCa - cnBC), nhapp, haend,cnBC,n,cnBC,habc);

    ufromGh(Gabc,habc,hmbeg+(nBCm -cnBC),hmend,umbeg+(nBCm -cnBC),umend,dx,n,cnBC,ua);

    conc(uabeg+(nBCa - cnBC), ua, uaend,cnBC,n,cnBC,uabc);

    double *nGappp = malloc(n*sizeof(double));
    double *nhappp = malloc(n*sizeof(double));

    evolve(Gabc,habc,uabc,g,dx,dt,n,cnBC,nhappp,nGappp);

// ################################### RK BUILD ###############################
    weightsum(i3,ha,2*i3,nhappp,n,ha);
    weightsum(i3,Ga,2*i3,nGappp,n,Ga);

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
