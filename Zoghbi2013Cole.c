/* Implements procedure from section 3.1 of Zoghbi, Reynolds, and Cackett 2013, ApJ, 777, 24 (which was taken from Bond, Jaffe, and Knox 1998, PRD, 57, 2117 and Lance Miller et al. 2010) for computing a correlation matrix and then computing a likelihood for a candidate underlying power spectrum. */

/* Code of Zoghbi et al. 2013 algorithm written by Coleman Miller. University of Maryland, College Park. 2015. */

/* In this version we go away from the recursive form for the determinant, because it appears to be causing problems when there are too many (>10?) data points. */

/* In this version we integrate over the log of the frequency.  This appears to be necessary because, at least for some types of red noise, the lowest frequencies contribute considerably.  Also, we subtract the mean of the light curve before doing analysis; this makes a profound difference to the results (e.g., without doing so, low-frequency modes completely dominate). */

/* In this version the hard code is changed to accept arguments. Changes performed by Peter Teuben. University of Maryland, College Park. 2016. */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI      3.14159265359
#define TWO_PI  6.28318530718
#define Nparams 2                /* Number of parameters in power density model */
#define ftol    1.0e-6           /* Parameter for integration */
#define Ninteg  10000            /* Number of steps in integral */

double logdeterminant(double **m, int n);
void inverse(double **m, double **inv, int n);
void exchange(double **m, int n, int k, int l);
void divide(double **m, int n, int k, double r);
void addrow(double **m, int n, int k, int l, double r);
void multiply(double **a, double **b, double **c, int l, int m, int n);
double power_slow(double freq, double *params);

inline static double sqr(double x) {
    return x*x;
}


inline static double power(double freq, double tau, double sigma)
{
    /* Given a model for the power density and values of the model's
     parameters, compute the power density at frequency freq.
     */
    
    double pow0;
    pow0=2.0*sqr(sigma*tau)/(1.0+sqr(TWO_PI*tau*freq));
    return pow0;
}


int main(int argc, char **argv)
{
    int i,j,k,kmax,Ndat,flag;
    double tmp;
    double *params;
    double *t,*X,*err,avg;
    double **Cx,**Cs,**Cn;
    double **invCx,**xvec,**xtransvec,**Cmultx,**product,exponent,logdet;
    double *pow0;
    double mindiff,logf,dlogf,logfmin,fmin,fmax,freq,integrand,powmax,logl;
    FILE *data,*wc;
    char *data_file;
    
    if (argc < 4) {
        printf("Usage: %s DATAFILE TAU SIGMA\n",argv[0]);
        printf("DATAFILE should contain time,flux,flux_error\n");
        printf("TAU and SIGMA are model parameters, eg. 56.2341 0.00728546\n");
        exit(1);
    }
    params=calloc(Nparams,sizeof(double));
    
    /* grab command line arguments */
    data_file = argv[1];           /* data file name with 3 columns */
    params[0] = atof(argv[2]);     /* tau */
    params[1] = atof(argv[3]);     /* sigma */
    
    /* open file once and count lines, so we know Ndat */
    data=fopen(data_file,"r");
    Ndat = 0;
    while (fscanf(data,"%lg %lg %lg", &tmp, &tmp, &tmp) > 0)
        Ndat++;
    fclose(data);
    
    /* allocate all other spaces */
    t=calloc(Ndat,sizeof(double));
    X=calloc(Ndat,sizeof(double));
    err=calloc(Ndat,sizeof(double));
    Cx=calloc(Ndat,sizeof(double *));
    invCx=calloc(Ndat,sizeof(double *));
    xvec=calloc(Ndat,sizeof(double *));
    Cmultx=calloc(Ndat,sizeof(double *));
    product=calloc(1,sizeof(double *));
    product[0]=calloc(1,sizeof(double));
    xtransvec=calloc(1,sizeof(double *));
    xtransvec[0]=calloc(Ndat,sizeof(double));
    Cs=calloc(Ndat,sizeof(double *));
    Cn=calloc(Ndat,sizeof(double *));
    for (i=0; i<Ndat; i++)
    {
        Cx[i]=calloc(Ndat,sizeof(double));
        invCx[i]=calloc(Ndat,sizeof(double));
        xvec[i]=calloc(1,sizeof(double));
        Cmultx[i]=calloc(1,sizeof(double));
        Cs[i]=calloc(Ndat,sizeof(double));
        Cn[i]=calloc(Ndat,sizeof(double));
    }
    pow0=calloc(Ninteg+1,sizeof(double));
    mindiff=1.0e20;
    avg=0.0;
    /* open data for real and read the data */
    data=fopen(data_file,"r");
    for (i=0; i<Ndat; i++)
    {
        fscanf(data,"%lg %lg %lg",&t[i],&X[i],&err[i]);
        avg+=X[i];
        if (i>0)
        {
            if (t[i]-t[i-1]<mindiff) mindiff=t[i]-t[i-1];
        }
    }
    fclose(data);
    /* continueing on... */
    avg/=(1.0*Ndat);
    for (i=0; i<Ndat; i++)
    {
        X[i]-=avg;
        xvec[i][0]=X[i];
        xtransvec[0][i]=X[i];
    }
    fmax=1.0e10/mindiff;
    fmin=1.0e-10/(t[Ndat-1]-t[0]);
    logfmin=log(fmin);
    dlogf=log(fmax/fmin)/(1.0*Ninteg);
    /* Now we have read in the data and declared the arrays.
     We need to compute the correlation matrix for the noise and the
     candidate signal and add them together to get the correlation
     matrix for the model.  We then compute and output the log likelihood.
     */
    for (i=0; i<Ndat; i++)
        Cn[i][i]=err[i]*err[i];   /* Noise correlation matrix */
    flag=0;
    powmax=0.0;
    for (k=0; k<=Ninteg && !flag; k++) /* Integral over frequency */
    {
        logf=logfmin+k*dlogf;
        freq=exp(logf);
        //pow0[k]=power_slow(freq,params);
        pow0[k]=power(freq,params[0],params[1]);
        if (freq*fabs(pow0[k])>powmax) powmax=freq*fabs(pow0[k]);
        if (freq*fabs(pow0[k])<ftol*powmax)
        {
            flag=1;
            kmax=k;
        }
    }
    if (!flag) kmax=Ninteg;
    
    /* Now compute signal correlation matrix for candidate model */
    flag=0;
    for (i=0; i<Ndat; i++)
    {
        for (j=i; j<Ndat; j++)
        {
            Cs[i][j]=0.0;  /* Might need to compute for many param combinations */
            Cs[j][i]=0.0;
            flag=0;
            for (k=0; k<=kmax; k++) /* Integral over frequency */
            {
                logf=logfmin+k*dlogf;
                freq=exp(logf);
                if (freq*fabs(pow0[k])>ftol*powmax)
                {
                    integrand=2.0*freq*pow0[k]*cos(2.0*PI*freq*(t[j]-t[i]));
                    if (k==0 || k==Ninteg)
                    {
                        Cs[i][j]+=integrand;
                        if (i != j) Cs[j][i]+=integrand;
                    }
                    else if (k%2==0)
                    {
                        Cs[i][j]+=2.0*integrand;
                        if (i != j) Cs[j][i]+=2.0*integrand;
                    }
                    else
                    {
                        Cs[i][j]+=4.0*integrand;
                        if (i != j) Cs[j][i]+=4.0*integrand;
                    }
                }
                else
                    flag=1;
            }
            Cs[i][j]*=dlogf/3.0;
            if (i != j) Cs[j][i]*=dlogf/3.0;
            Cx[i][j]=Cs[i][j]+Cn[i][j];
            if (i != j) Cx[j][i]=Cs[j][i]+Cn[j][i];
            // printf("%3d %3d %lg\n",i,j,Cx[i][j]);
        }
    }
    logdet=logdeterminant(Cx,Ndat);
    inverse(Cx,invCx,Ndat);
    multiply(invCx,xvec,Cmultx,Ndat,Ndat,1);
    multiply(xtransvec,Cmultx,product,1,Ndat,1);
    exponent=product[0][0];
    logl=-0.5*Ndat*log(2.0*PI)-0.5*logdet-0.5*exponent;
    printf("%lg\n",logl);
    
    for (i=0; i<Ndat; i++)
    {
        free(Cx[i]);
        free(invCx[i]);
        free(xvec[i]);
        free(Cmultx[i]);
        free(Cs[i]);
        free(Cn[i]);
    }
    free(Cx);
    free(invCx);
    free(xvec);
    free(Cmultx);
    free(xtransvec[0]);
    free(xtransvec);
    free(product[0]);
    free(product);
    free(Cs);
    free(Cn);
    free(t);
    free(X);
    free(err);
    free(params);
    free(pow0);
}

double power_slow(double freq, double *params)
{
    /* Given a model for the power density and values of the model's
     parameters, compute the power density at frequency freq.
     */
    
    double tau,sigma,pow0;
    
    tau=params[0];
    sigma=params[1];
    pow0=2.0*sigma*sigma*tau*tau/(1.0+(2.0*PI*tau*freq)*(2.0*PI*tau*freq));
    
    /*   double amp,slope,pow0;
     
     amp=params[0];
     slope=params[1];
     pow0=amp*pow(freq,slope); */
    
    return pow0;
}

double logdeterminant(double **a, int n)
{
    /* Calculates the log of the determinant of an nxn matrix */
    
    int i,j,k;
    double logdet,**detmat;
    
    detmat=calloc(n,sizeof(double *));
    for (i=0; i<n; i++)
        detmat[i]=calloc(n,sizeof(double));
    for (i=0; i<n; i++)
        for (j=0; j<n; j++)
            detmat[i][j]=a[i][j];
    for (i=0; i<n; i++)
    {
        if (detmat[i][i]==0.0)
        {
            k=i+1;
            do
            {
                exchange(detmat,n,i,k);
                k++;
            } while (detmat[i][i]==0.0);
        }
        for (j=0; j<n; j++)
            if (j!=i)
                addrow(detmat,n,i,j,-detmat[j][i]/detmat[i][i]);
    }
    /* Exchanges of rows multiply the determinant by -1, but in this
     specific case the determinant had better be positive.  Thus we
     ignore the sign flips.
     */
    
    logdet=0.0;
    for (i=0; i<n; i++)
        logdet+=log(fabs(detmat[i][i]));
    
    for (i=0; i<n; i++)
        free(detmat[i]);
    free(detmat);
    
    return logdet;
}

void inverse(double **m, double **inv, int n)
{
    /* Calculates the inverse of an nxn matrix using Gaussian elimination */
    
    int i,j,k;
    double r,**mtemp;
    
    mtemp=calloc(n,sizeof(double *));
    for (i=0; i<n; i++)
        mtemp[i]=calloc(n,sizeof(double));
    for (i=0; i<n; i++)
    {
        for (j=0; j<n; j++)
        {
            if (i==j)
                inv[i][j]=1.0;
            else
                inv[i][j]=0.0;
            mtemp[i][j]=m[i][j];
        }
    }
    for (i=0; i<n; i++)
    {
        if (mtemp[i][i]==0.0)
        {
            k=i+1;
            do
            {
                exchange(mtemp,n,i,k);
                exchange(inv,n,i,k);
                k++;
            } while (mtemp[i][i]==0.0);
        }
        r=mtemp[i][i];
        divide(inv,n,i,r);
        divide(mtemp,n,i,r);
        for (j=0; j<n; j++)
        {
            if (j!=i)
            {
                addrow(inv,n,i,j,-mtemp[j][i]);
                addrow(mtemp,n,i,j,-mtemp[j][i]);
            }
        }
    }
    for (i=0; i<n; i++)
        free(mtemp[i]);
    free(mtemp);
}

void exchange(double **m, int n, int k, int l)
{
    /* Exchanges rows k and l of an nxn matrix */
    
    int j;
    double temp;
    
    for (j=0; j<n; j++)
    {
        temp=m[k][j];
        m[k][j]=m[l][j];
        m[l][j]=temp;
    }
}

void divide(double **m, int n, int k, double r)
{
    /* Divides row k of an nxn matrix by r */
    
    int j;
    
    for (j=0; j<n; j++)
    {
        m[k][j]/=r;
    }
}

void addrow(double **m, int n, int k, int l, double r)
{
    /* Add r times row k of an nxn matrix to row l */
    
    int j;
    
    for (j=0; j<n; j++)
        m[l][j]+=r*m[k][j];
}

void multiply(double **a, double **b, double **c, int l, int m, int n)
{
    /* Matrix multiplication of lxm matrix a by mxn matrix b to get lxn matrix c */
    
    int i,j,k;
    
    for (i=0; i<l; i++)
    {
        for (j=0; j<n; j++)
        {
            c[i][j]=0.0;
            for (k=0; k<m; k++)
            {
                if (k==0) c[i][j]=0.0;
                c[i][j]+=a[i][k]*b[k][j];
            }
        }
    }
}
