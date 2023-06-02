/********************************************************************/
/*                                                                  */
/* MODULE  : kt11.c                                                  */
/* MODIFIED:    Interpolates in LISA noise curve                    */
 /*                   (using esa-LISA 2018 curve)                   */
/*              SNR part completed (May 2020)                       */
/*              corrected chirping SNR calculation (May 2020)       */
/*              corrected Shane's factor of 4 in ho (June 2020)     */
/*              modified for COSMIC output ***                      */
/*              added tracking of type of WDs in binaries           */
/*              outputs NS and BH binaries to separate file         */
/*              outputs SNR between 0 and 1 to separate file        */
/*                                                                  */
/* AUTHORS : Patricia Purdue (Colorado College)                     */
/*                                                                  */
/********************************************************************/


/* ============ INCLUDE LIBRARIES ============ */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* ============ CONSTANTS DEFINED ============ */

#define PI              3.14159265358979323846      /* Why ask pi?              */
#define G				6.67384e-11                 /* N m^2 kg^-2 (2010 CODATA value)  */
#define Msun			1.9891e30					/* mass of sun in kilograms */
#define C               299792458                   /* meters per second        */
#define mkpc            3.08567758e19               /* meters in a kiloparsec- maybe off by factor of 1000 */
#define yrsec           31557600                    /* seconds in a year */
#define daysec          86400                       /* seconds in a day */

#define MAXLINELENGTH   255                         /* maximum string length */

/* ============== constants for LISA data and manipulation ------------- */

double  Lf[900], Lh[900];                        /* freq and strain arrays for LISA noise curve */
                                                 /* (postioning this command here makes it accessible */
                                                 /* inside LISA interpolation function) */

/* ========================================================= */
/*   MAIN  ROUTINE                                           */
/* ========================================================= */

int main()
{
    /* ============= VARIABLE DECLARATIONS ============= */
    
    char    hdrRead[MAXLINELENGTH+1];
    
    FILE *ifp, *ofp1, *ofp2, *ofp3, *ofp4, *LISAdat;
    
    ifp = fopen("gx_dat_full.csv","r");      /* specify input file */
    ofp1 = fopen("full_DWD_1.txt","w");           /* specify output files */
    ofp2 = fopen("full_DWD_0.txt","w");
    ofp3 = fopen("full_nonDWD_1.txt","w");
    ofp4 = fopen("full_nonDWD_0.txt","w");
    
                                /* NOTE: When changing input file, check the following: */
                                /*        - number, type, & contents of column variables matches */
                                /*        - number of lines of data is correct */
    
    int      binnum;
    float    tphys, m1, m2, kstar1, kstar2, sep, ecc, Porbday, fGWpeak, xGx, yGx, zGx, dist;
                                                            /* variables from input file */
                                                            /* m1, m2 = masses in Msun */
                                                            /* Porbday = orbital period in DAYS */
                                                            /* kstar1/2 indicate star type */
                                                            /*  10 = He WD, 11 = C-O WD, 12 = O-Ne WD */
                                                            /*  13 = NS, 14 = BH, 15 = massless remnant */
                                                            /* rest not used by this program */
    
    float   snrcutoff = 1.0;                                /* cutoff value for SNR */
    int     inputfilelines = 39388147;                      /* number of lines in input data file */
                                                            /* (assumes 1 header line) */
                                                            /* 121 for gx_dat_test.csv   */

   /* ------------------ iteration variables --------------------- */
    
    int     i = 1;                                          /* iteration variable for binaries */
    int     j = 1, k = 0;                                   /* iteration variables for LISA curve */
    int     ii = 0;                                         /* iteration variable for chirping SNR */
    
   /* ------------------ useful fractions ----------------------- */
    
     double  threeelevenths = 3.0/11.0;                      /* exponent for sfreq calculation */
     double  eightthirds = 8.0/3.0;
     double  threeeighths = 3.0/8.0;
     double  negeleventhirds = -11.0/3.0;
     double  twothirds = 2.0/3.0;
    
    /* ---------- constants used in h and f calculation ------------- */
    
    double  Jconst, hconst, fconst, ms5, G2;                /* constants to make calculations faster */
    
       ms5 = pow(Msun,5);                              /* Msun terms convert units to kg */
       G2 = pow(G,2);
       Jconst = cbrt(G2*ms5/(2*PI));                   /* constants appearing in Jorb calculation */
       hconst = 4.0 * pow(G,3) * ms5/pow(C,4);         /* constants appearing in hnorm calculation */
       fconst = G2 * ms5/PI;                           /* constants appearing in f calculation */
       
    /* ----- variables for determining type of WDs in binary ------ */
    
    int     WDbintype;                                      /* category of source DWD binary: */
                                                            /* 0 = non-WDs, which are being excluded */
                                                            /* 1 = He-He */
                                                            /* 2 = He-C/O */
                                                            /* 3 = He-O/Ne */
                                                            /* 4 = C/O-C/O */
                                                            /* 5 = C/O-O/Ne */
                                                            /* 6 = O/Ne-O/Ne */
    
    int     DWDtype(float,float);                     /* define function to characterize type of binary */
    
    /* ---------- variables for h and f calculation ------------- */
    
    double  ma, md, q, Q, Q3, mt, mt5, dhnorm, f;           /* variables calculated in program */
    double  Porb, Jorb, t3, t5;
    
    /* ---------- constants used in chirping calculations ------------- */
       
    double  kapcoeftop, kapcoefbot, kapcoef;             /* values to calculate kappa */
    double  Tobs, Tobssq;                                    /* observation time  */
      
        Tobs = 1.0 * yrsec;                            /* observation time in years converted to seconds */
        Tobssq = pow(Tobs,2);
        kapcoeftop = 5.0/(256.0 * cbrt(pow(PI,8)));
        kapcoefbot = pow(cbrt(G/pow(C,3)),5);
        kapcoef = kapcoeftop/kapcoefbot;              /* kappa appears in stationary frequency calculation */
    
    
    /* ---------- variables for chirping calculations ------------- */
    
    double  mchirp, sfreq;                                  /* chirp mass, stationary frequency for SNR */
    double  kappa;
    
    double  fbinlo;                                         /* lower boundary frequency for bins */
    double  df, fi, ff, nbin;                               /* freq increment, starting & ending freq, 
                                                                number of bins */
    double  nbininit, hcoef;                                /* rounded off number of bins, coef needed for h*/
    
    
    double  ho;                                             /* scaling amplitude */
    double  hsrc;                                           /* amplitude of source */
    
    double tbin, hobin, hchirpbin, hLISAbin, snrbin;        /* for chirping sources */
    
    /* ---------- variables for LISA data and manipulation ------------- */

    double LISAinter(double);                               /* define function to interpolate LISA */
                                                            /* sensitivity for specific frequency */
    
    double hLISA;                                           /* sensitivity of LISA for freq f */
                                                            /* (output of LISAinter() function) */
    
    double LISAfmin, LISAfmax;                              /* min and max values of LISA freq range */
    
    /* -------- variables for SNR calculation -------------------------- */
    
    int     sigtype;                                        /* category of source signal: */
                                                            /* 0 = out of range */
                                                            /* 1 = monochromatic */
                                                            /* 2 = chirping */
    
    double  snr, snrsq;                                     /* signal to noise ratio*/
    
    /* ----------------- counting variables ------------------------- */
    
    int outcnt = 0;                                         /* count number of sources outside LISA band */
    int monocnt = 0;                                        /* count number of monochromatic sources */
    int chirpcnt = 0;                                       /* count number of chirping sources */
    int totDWDcnt = 0;                                      /* count number of DWD sources checked */
    int detectcnt = 0;                                      /* count number of detectable sources */
                                                            /*  (above SNR cutoff) */
    int totcnt = 0;                                         /* count total number of sources checked */
    int nonDWDcnt = 0;                                      /* count number of non-DWD binary sources */
    int othercnt = 0;                                       /* count number of non-binaries */
    
    /* =============== Read in LISA Noise Curve =============== */
    
    LISAdat = fopen("../LISA2018_esaSTDp.txt","r");            /* input LISA noise curve */
                                                                 /* Note: Using updated LISA curve, */
                                                                 /* had to make file plain csv for code to read */

    while (j <= 26) {                                       /* 26 header lines in LISA file     */
        fgets(hdrRead,MAXLINELENGTH,LISAdat);
        ++j;
    }

    while (k < 900) {                                        /* 900 lines of data in LISA file     */
        fscanf(LISAdat, "%lf,%lf", &Lf[k], &Lh[k]);          /* reading in LISA noise curve (freq & h) */
        ++k;
    };
    
    /*printf("Lf[0]: %g   Lh[0]: %g \n",Lf[0],Lh[0]);
    printf("Lf[899]: %g   Lh[899]: %g \n",Lf[899],Lh[899]);*/

    fclose(LISAdat);
    
    /* LISAfmin = Lf[0]; */
    LISAfmin = 1.0e-5;                                      /* set low-frequency cutoff */
    LISAfmax = Lf[899];                                     /* set high-frequency cutoff */
    /* printf("LISA freq range: %g  to %g\n",LISAfmin,LISAfmax); */
    
    /* =============== START ROUTINE HERE =============== */
    
    /* ------------------ preliminaries ----------------- */
    
    
    /*printf("ms5: %g\n",ms5);
    printf("Jconst: %g\n",Jconst);
    printf("hconst: %g\n",hconst);
    printf("fconst: %g\n",fconst);
    printf("kapcoeftop: %g\n",kapcoeftop);
    printf("kapcoefbot: %g\n",kapcoefbot);
    printf("kapcoef: %g\n",kapcoef);*/

    
    /* ----------------- put header lines in output files --------------------- */
    
        /* --- high SNR DWD sources --- */
    
    fprintf(ofp1,"# DWD sources above SNR cutoff: %f\n",snrcutoff);
    
    fputs("# Signal Type 0 = source out of detector range \n",ofp1);
    fputs("# Signal Type 1 = monochromatic source \n",ofp1);
    fputs("# Signal Type 2 = chirping source \n",ofp1);
    
    fputs("# DWD Binary Type 1 = He-He \n",ofp1);
    fputs("# DWD Binary Type 2 = He-C/O \n",ofp1);
    fputs("# DWD Binary Type 3 = He-O/Ne \n",ofp1);
    fputs("# DWD Binary Type 4 = C/O-C/O \n",ofp1);
    fputs("# DWD Binary Type 5 = C/O-O/Ne \n",ofp1);
    fputs("# DWD Binary Type 6 = O/Ne-O/Ne \n",ofp1);

    fputs("# BinNum SigType WDbintype Freq(Hz) dhnorm(m) md(solar) ma(solar) mchirp Porb(s) dist SNR\n",ofp1);
    fputs("BinNum SigType WDbintype Freq dhnorm md ma mchirp Porb dist SNR\n",ofp1);     /* this line is for pandas dataframe headers */
    
        /* --- low SNR DWD sources --- */
    
    fprintf(ofp2,"# DWD sources below SNR cutoff: %f\n",snrcutoff);
    
    fputs("# Signal Type 0 = source out of detector range \n",ofp2);
    fputs("# Signal Type 1 = monochromatic source \n",ofp2);
    fputs("# Signal Type 2 = chirping source \n",ofp2);
    
    fputs("# DWD Binary Type 1 = He-He \n",ofp2);
    fputs("# DWD Binary Type 2 = He-C/O \n",ofp2);
    fputs("# DWD Binary Type 3 = He-O/Ne \n",ofp2);
    fputs("# DWD Binary Type 4 = C/O-C/O \n",ofp2);
    fputs("# DWD Binary Type 5 = C/O-O/Ne \n",ofp2);
    fputs("# DWD Binary Type 6 = O/Ne-O/Ne \n",ofp2);

    fputs("# BinNum SigType WDbintype Freq(Hz) dhnorm(m) md(solar) ma(solar) mchirp Porb(s) dist SNR\n",ofp2);
    fputs("BinNum SigType WDbintype Freq dhnorm md ma mchirp Porb dist SNR\n",ofp2);   /* this line is for pandas dataframe headers */
    
        /* --- non-DWD sources above cutoff --- */
    
    fprintf(ofp3,"# non-DWD sources, above SNR cutoff: %f\n",snrcutoff);
    
    fputs("# Signal Type 0 = source out of detector range \n",ofp3);
    fputs("# Signal Type 1 = monochromatic source \n",ofp3);
    fputs("# Signal Type 2 = chirping source \n",ofp3);
    
    fputs("# Binary Type 7 = He WD-NS \n",ofp3);
    fputs("# Binary Type 8 = C/O WD-NS \n",ofp3);
    fputs("# Binary Type 9 = O/Ne WD-NS \n",ofp3);
    fputs("# Binary Type 10 = He WD-BH \n",ofp3);
    fputs("# Binary Type 11 = C/O WD-BH \n",ofp3);
    fputs("# Binary Type 12 = O/Ne WD-BH \n",ofp3);
    fputs("# Binary Type 13 = NS-NS \n",ofp3);
    fputs("# Binary Type 14 = NS-BH \n",ofp3);
    fputs("# Binary Type 15 = BH-BH \n",ofp3);
   

    fputs("# BinNum SigType WDbintype Freq(Hz) dhnorm(m) md(solar) ma(solar) mchirp(solar) Porb(s) dist SNR\n",ofp3);
    fputs("BinNum SigType WDbintype Freq dhnorm md ma mchirp Porb dist SNR\n",ofp3);   /* this line is for pandas dataframe headers */
    
    /* --- non-DWD sources below cutoff --- */
     
     fprintf(ofp4,"# non-DWD sources, SNR cutoff: %f\n",snrcutoff);
     
     fputs("# Signal Type 0 = source out of detector range \n",ofp4);
     fputs("# Signal Type 1 = monochromatic source \n",ofp4);
     fputs("# Signal Type 2 = chirping source \n",ofp4);
     
     fputs("# Binary Type 7 = He WD-NS \n",ofp4);
     fputs("# Binary Type 8 = C/O WD-NS \n",ofp4);
     fputs("# Binary Type 9 = O/Ne WD-NS \n",ofp4);
     fputs("# Binary Type 10 = He WD-BH \n",ofp4);
     fputs("# Binary Type 11 = C/O WD-BH \n",ofp4);
     fputs("# Binary Type 12 = O/Ne WD-BH \n",ofp4);
     fputs("# Binary Type 13 = NS-NS \n",ofp4);
     fputs("# Binary Type 14 = NS-BH \n",ofp4);
     fputs("# Binary Type 15 = BH-BH \n",ofp4);
    

     fputs("# BinNum SigType WDbintype Freq(Hz) dhnorm(m) md(solar) ma(solar) mchirp(solar) Porb(s) dist SNR\n",ofp4);
     fputs("BinNum SigType WDbintype Freq dhnorm md ma mchirp Porb dist SNR\n",ofp4);   /* this line is for pandas dataframe headers */
    
    /* ---------------- start calculations on data  ----------------- */

    while (i <= 1) {                                /* skip the header lines of input data file */
           fgets(hdrRead,MAXLINELENGTH,ifp);
           ++i;
       }
       
/*  printf("i: %i\n",i); */
    
    while (i <= inputfilelines) {                 /* iteration maximum specified above   */
                                                   
        
        fscanf(ifp,"%i,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f",
               &binnum, &tphys, &m1, &m2, &kstar1, &kstar2, &sep, &ecc, &Porbday, &fGWpeak,
               &xGx, &yGx, &zGx, &dist);                /* read in data point */
        
        WDbintype = DWDtype(kstar1,kstar2);
  
     if (WDbintype != 0){
        if (m1 < m2)                                /* make sure that md < ma */
        {md = m1, ma = m2;}
            else
            {md = m2, ma = m1;}
        q = md/ma;                                  /* ratio of the masses, note: 0 < q <= 1 */
        Q = q/pow(1+q,2);                           /* ratio of reduced mass to total mass */
        Q3 = pow(Q,3);
        mt = ma + md;                               /* total mass */
        mt5 = pow(mt,5);
      
        Porb = Porbday * daysec;                    /* convert period to seconds */
        Jorb = cbrt(Porb/mt)*ma*md*Jconst;          /* orbital angular momentum of binary */
        
        t3 = mt5 * Q3 * hconst;
        dhnorm = t3 / pow(Jorb,2);                  /* (distance)*(GW amplitude) (vertical axis value) */
        
        t5 = mt5 * Q3 * fconst;
        f = t5/pow(Jorb,3);                         /* GW frequency (horizontal axis value) */
           
    if (f < LISAfmax && f > LISAfmin) {         /* only calculate SNR if in LISA's band */
            
        /* -------------- find SNR and make cutoff ----------------- */
        
        mchirp = pow(m1*m2,0.6)/pow(mt,0.2);        /* calculate chirp mass */

        kappa = kapcoef/pow(cbrt(mchirp * Msun),5);  /* find kappa for stationary freqency */
        
        sfreq = pow(eightthirds*kappa/Tobssq,threeelevenths);   /* stationary frequency */
        
        ho = dhnorm/8/mkpc;                         /*  scaling amplitude assuming source at 8kpc */
        
        if (f < sfreq) {                            /* if freq low enough, source is monochromatic */
            sigtype = 1;
            hsrc = ho * sqrt(Tobs);                     /* estimate of source strenth for mono source */
            /* printf("mono source with hsrc: %g\n",hsrc); */
            /*printf("monochromatic source\n");*/
            hLISA = LISAinter(f);
            snr = hsrc/hLISA;                           /* SNR estimate for mono source */
            ++monocnt;                                  /* count number of monochromatic sources */
        } else {
            sigtype = 2;
            fi = f;                                         /* initial frequency is given frequency */
            ff = 1/pow(1/pow(fi,eightthirds) - Tobs/kappa, threeeighths);
                                                            /* find final frequency during observation time */
            nbin = Tobs*(ff-fi);
            nbininit = ceil(nbin);                          /* calculate number of bins (rounding up) */
            /*printf("number of bins: %g\n", nbininit);*/
            hcoef = sqrt(eightthirds*kappa);
            snrsq = 0.;
            for (ii = 0; ii < nbininit; ++ii) {
                fbinlo = fi + ii * df;
                tbin = eightthirds * kappa * pow(fbinlo,negeleventhirds)/Tobs;
                hobin = ho*pow(fbinlo/f, twothirds);
                hchirpbin = pow(ho,2) * tbin;                  /* actually square of numerator for snr sum */
                hLISAbin = LISAinter(fbinlo);                  /* use interpolation rountine to find h_LISA_bin */
                snrbin = hchirpbin/pow(hLISAbin,2);             /* actually the square of snr for bin */
                snrsq = snrsq + snrbin;                         /* sum of squares for snr */
            };                                                  /* end of for (bins) */
            snr = sqrt(snrsq);                                  /* corrected chirping SNR (May 2020) to be */
                                                                /* square root of sum of the squares */
            ++chirpcnt;                                         /* count number of chirping sources */
            };                                                      /* end of else for SNR */
            
          /* end SNR part */
        
    } else {                    /* end of if for range check */
        sigtype = 0;               /* sigtype 0 is out of range */
        ++outcnt;               /* count number of sources outside LISA band */
    }                       /* end of if-else for sources in LISA frequency range */
         
     };                      /* end of if for source type */
       
    /* ------------- printing results to files --------------- */
         
    
        if (sigtype !=0) {                  /* only print points that are within the LISA band */
        if (WDbintype > 0 && WDbintype < 7) {
            if (snr > snrcutoff){             /* only keeping points with SNR above cutoff  */
                                                        /* and within LISA band */
                fprintf(ofp1, "%i %i %i %g %g %g %g %g %g %g %g\n",
                    binnum, sigtype, WDbintype, f, dhnorm, md, ma, mchirp, Porb, dist, snr);
                                                    /* write results to file */
            ++detectcnt;                             /* count number of detectable sources */
            } else {
                fprintf(ofp2, "%i %i %i %g %g %g %g %g %g %g %g\n",
                        binnum, sigtype, WDbintype, f, dhnorm, md, ma, mchirp, Porb, dist, snr);
            };
        ++totDWDcnt;                                   /* count total number of sources checked */
         } else {
         if (WDbintype > 6 && WDbintype < 16) {
             if (snr > snrcutoff){             /* only keeping points with SNR above cutoff  */
                                                         /* and within LISA band */
                 fprintf(ofp3, "%i %i %i %g %g %g %g %g %g %g %g\n",
                     binnum, sigtype, WDbintype, f, dhnorm, md, ma, mchirp, Porb, dist, snr);
                                                     /* write results to file */
             ++detectcnt;                             /* count number of detectable sources */
             } else {
                 fprintf(ofp4, "%i %i %i %g %g %g %g %g %g %g %g\n",
                         binnum, sigtype, WDbintype, f, dhnorm, md, ma, mchirp, Porb, dist, snr);
             };
             /*fprintf(ofp3, "%i %i %i %g %g %g %g %g %g %g %g\n",
                     binnum, sigtype, WDbintype, f, dhnorm, md, ma, mchirp, Porb, dist, snr); */
             ++nonDWDcnt;
         } else {++othercnt;}
         }};      /* end of if-else statement for printing results to file */
        
         ++i;
     }      /* end of while statement */
    
    totcnt = totDWDcnt + nonDWDcnt;
    
    printf("\n");
    printf("total number of sources checked: %i\n",totcnt);
    printf("number of DWD sources checked: %i\n",totDWDcnt);
    printf("number of non-DWD binary sources checked: %i\n",totDWDcnt);
    printf("number of non-binary sources checked: %i\n",othercnt);
    printf("number of out-of-band sources found: %i\n",outcnt);
    printf("number of monochromatic sources found: %i\n",monocnt);
    printf("number of chirping sources found: %i\n",chirpcnt);
    printf("SNR cutoff: %f\n",snrcutoff);
    printf("number of detectable sources found: %i\n",detectcnt);
    
    fclose(ifp);
    fclose(ofp1);
    fclose(ofp2);
    fclose(ofp3);
    fclose(ofp4);
    
	return 0;
}

/* ========================================================= */
/*   LISA INTERPOLATION ROUTINE                              */
/* ========================================================= */

double LISAinter(double f)
{
    /* ---------- variables for LISA data and manipulation ------------- */

    extern double  Lf[900], Lh[900];                        /* freq and strain arrays for LISA noise curve */
                                                            /* extern command should make it accessible in */
                                                            /* LISA interpolation function */
    
    int     k = 0;                                          /* iteration variable for LISA curve */
    double  flo, fhi, hlo, hhi;                             /* points on either side of target frequency */
    double  logflo, logfhi, loghlo, loghhi;                 /* logarithms of above points */
    double  logslope, logint;                               /* logs of slope and intercept for interpolation */
    double  loghLISA, hLISA;                                /* h of LISA for given freqency */
    
    /* ------- find and interpolate relevant LISA points ------- */
 
    k = 0;
    while (Lf[k]-f < 0) {                                   /* locate relevant part of LISA noise curve */
        flo = Lf[k];
        hlo = Lh[k];                                        /* outputs point just below freq of interest */
        ++k;
    }
 
    fhi = Lf[k]; hhi = Lh[k];                               /* assign point just above freq of interest */
 
    /* printf("k: %i    flo: %g   hlo: %g \n",k,flo,hlo);*/
    /* printf("k: %i    fhi: %g   hhi: %g \n",k,fhi,hhi);*/
 
    logflo = log10(flo);
    logfhi = log10(fhi);
    loghlo = log10(hlo);
    loghhi = log10(hhi);
 
    logslope = (loghhi - loghlo)/(logfhi - logflo);         /* find slope and intercept in logspace */
    logint = loghlo - logslope * logflo;
    /* printf("logslope: %g    logint: %g\n",logslope,logint);*/
 
    loghLISA = logint + logslope * log10(f);                /* find sensitivity at freq of interest */
    hLISA = pow(10, loghLISA);
 
    /* printf("loghLISA: %g    hLISA: %g\n",loghLISA,hLISA); */
    return hLISA;                                           /* returns LISA sensitivity at input freq */
}

/* ========================================================= */
/*   DWD BINARY TYPE SORTING ROUTINE                         */
/* ========================================================= */

int DWDtype(float k1, float k2)
{
    float   klg,ksm;
    int     bintype;
    
    if (k1 < k2)                                /* make sure that ksm < klg */
    {klg = k2, ksm = k1;}
        else
        {klg = k1, ksm = k2;}
    
    if (ksm < 10.0) {bintype = 0;}
    if (klg > 14.0) {bintype = 0;}
    if (ksm == 10.0 && klg == 10.0) {bintype = 1;}      /* He-He DWD */
    if (ksm == 10.0 && klg == 11.0) {bintype = 2;}      /* He-C/O DWD */
    if (ksm == 10.0 && klg == 12.0) {bintype = 3;}      /* He-O/Ne DWD */
    if (ksm == 11.0 && klg == 11.0) {bintype = 4;}      /* C/O-C/O DWD */
    if (ksm == 11.0 && klg == 12.0) {bintype = 5;}      /* C/O-O/Ne DWD*/
    if (ksm == 12.0 && klg == 12.0) {bintype = 6;}      /* O/Ne-O/Ne DWD */
    if (ksm == 10.0 && klg == 13.0) {bintype = 7;}      /* He WD-NS */
    if (ksm == 11.0 && klg == 13.0) {bintype = 8;}      /* C/O WD-NS */
    if (ksm == 12.0 && klg == 13.0) {bintype = 9;}      /* O/Ne WD-NS */
    if (ksm == 10.0 && klg == 14.0) {bintype = 10;}     /* He WD-BH */
    if (ksm == 11.0 && klg == 14.0) {bintype = 11;}     /* C/O WD-BH */
    if (ksm == 12.0 && klg == 14.0) {bintype = 12;}     /* O/Ne WD-BH */
    if (ksm == 13.0 && klg == 13.0) {bintype = 13;}     /* NS-NS */
    if (ksm == 13.0 && klg == 14.0) {bintype = 14;}     /* NS-BH */
    if (ksm == 14.0 && klg == 14.0) {bintype = 15;}     /* BH-BH */

    
    /*printf("kstar1: %f      kstar2: %f      bintype: %i\n",k1,k2,bintype);*/
    return bintype;
}
