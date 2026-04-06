/* ts2phima.c */

/*
++++++++++++++++++++++++++++++++++++++++++++++++++
 
.IDENTIFICATION
  program TS2PHIMA			version 1.00	940328
  A. Kaufer                  		LSW - Heidelberg
  1.00	   940328	
  1.10     940415
 
.KEYWORDS
  bulk data frames, time series
 
.PURPOSE
  convert a time series of 1-dim spectra to 2-dim frame 
 
.ALGORITHM
  extract input frames from catalogue,
  renorm extracted frame 
  filter median 3x1 in t direction 
  resample to phase bins
 
.INPUT/OUTPUT
  the following keys are used:
 
  OUT_A/C/1/60			result frame
  P3/C/1/80			catalog of frames in time series
  INPUTD/D/1/1                  lower limit of x interval
  INPUTD/D/2/1                  upper limit of x interval
  INPUTI/I/1/2                  catalog interval
  INPUTD/D/4/1                  period
  INPUTD/D/5/1                  epoch
  INPUTD/D/6/1                  no. of phase-bins
 
  P7/C/1/1                      renorm flag
  P8/C/1/1                      filter flag  
    

.VERSIONS
  1.00		first version, created from ts2ima.c
  1.10          averaging in phase-bins
--------------------------------------------------
*/

#include <midas_def.h>
#include <math.h>

#define NRANSI
#include "nr.h"
#include "nrutil.h"

#define MAXCAT 1000
#define MAXWIN 25
#define HMAXMED 1
#define BLACK -1.0

void median(float x[],int n,float *xmed);
void statistics(float data[], int n, float *ave, float *adev, float *sdev,
	float *var, float *min, float *max);


main() {

int      renorm=0, filter=0, first=1;
int      t_limits[2];  
int      frmcnt, iav, iseq, x, xstart, t, tt, tmax, r, i, n;
int      naxis, npix[2], tnaxis, tnpix[2], rnaxis, rnpix[2];
int      imno, imnot, imnor, null, dum, noent, noelem, bytelem;
int      tid, xocol, tocol, idocol, entry;
int      nbin, phase_bin, *yy;
int      unit, uni;
 
char     flag[1], dtype[1];
char     line[84];
char     msg[80];
char     frame[64], rframe[64], catfil[64];
char     cunit[16*(2+1)], ident[72], output[64], history[64];

float     *inima, *tempima, *resima, *tcol, *data;
float     yp1, ypn, out;
float     mjd[MAXCAT], medwin[2*HMAXMED+1], leftwin[MAXWIN], rightwin[MAXWIN];
float     med, medrwin, medlwin, a, b, phase[MAXCAT], catnr[MAXCAT];
float     lhcuts[4], min=1.e99, max=-1.e99, ave, adev, sdev, var;

double   start[2], step[2], tstart[2], tstep[2], rstart[2], rstep[2];
double   x_limits[3], lam_c, v_lo, v_hi, period, epoch;
double   time[MAXCAT][6], dphase, dbin;
double   xo1, xo2, mjdo;
 
 
/* set up MIDAS environment + enable automatic error abort */

SCSPRO("ts2phima");

/* get result frame, history, input frame list, flags  */

SCKGETC("OUT_A",1L,60L,&iav,rframe);
SCKGETC("HISTORY",1L,64L,&iav,history);
SCKGETC("P3",1L,80L,&iav,line);
SCKGETC("P7",1L,1L,&iav,flag);
if (strcmp(flag,"y")==0 || strcmp(flag,"Y")==0) 
 renorm=1; 
SCKGETC("P8",1L,1L,&iav,flag);
if (strcmp(flag,"y")==0 || strcmp(flag,"Y")==0) 
 filter=1; 



/*  get the data intervals, period, epoch, phase_bin  */

SCKRDD("INPUTD",1L,2L,&iav,x_limits,&unit,&null);
SCKRDD("INPUTD",4L,1L,&iav,&period,&unit,&null);
SCKRDD("INPUTD",5L,1L,&iav,&epoch,&unit,&null);
SCKRDD("INPUTD",6L,1L,&iav,&dbin,&unit,&null);
SCKRDI("INPUTI",1L,2L,&iav,t_limits,&unit,&null);
phase_bin = (long int)dbin;

if (x_limits[1] < x_limits[0]){
  dum = x_limits[1];
  x_limits[1] = x_limits[0];
  x_limits[0] = dum;
  sprintf(msg," correction: using x-interval from %6.1f to %6.1f  \n",x_limits[0],x_limits[1]);
  SCTPUT (msg);
}

/*   handle input from catalog - get name of first image file in catalog  */

strcpy(catfil,line);

SCCSHO(catfil,&noent,&dum);

if (t_limits[1] < t_limits[0]){
  dum = t_limits[1];
  t_limits[1] = t_limits[0];
  t_limits[0] = dum;
  sprintf (msg, " correction: using catalog entries %d - %d \n",t_limits[0],t_limits[1]);
  SCTPUT (msg);
}

if (noent < t_limits[1]) {
  t_limits[1] = noent;
  sprintf (msg, " correction: using catalog entries %d - %d \n",t_limits[0],t_limits[1]);
  SCTPUT (msg);
}

/* default: use all active entries    */
if (t_limits[0]==0 || t_limits[1]==0) {
  t_limits[0] = 1;
  t_limits[1] = noent;
}

/* loop on catalogue */

for (frmcnt = 1; frmcnt < t_limits[0]; frmcnt++) {
   SCCGET(catfil,0L,frame,output,&iseq);
   if (frame[0] == ' ') break;		
}


t = 0; 
for (frmcnt = t_limits[0]; frmcnt <= t_limits[1]; frmcnt++) {

   SCCGET(catfil,0L,frame,output,&iseq);
   if (frame[0] == ' ') break;		

   catnr[t] = (float)t;
   SCIGET(frame,D_R4_FORMAT,F_I_MODE,F_IMA_TYPE,1L,&naxis,npix,
        	 start,step,ident,cunit,(char **)&inima,&imno);

   SCDFND (imno,"DTIME", dtype, &noelem, &bytelem);
   if (strcmp (dtype, " ")!=NULL) {
     SCDRDD(imno,"DTIME",1L,6L,&iav,&time[t][0],&unit,&null);   
   } else {
     sprintf(msg," Warning: descriptor DTIME not present in %s, using O_TIME"
	     , frame);
     SCTPUT (msg);
     SCDRDD(imno,"O_TIME",1L,6L,&iav,&time[t][0],&unit,&null);
   }

/*    take the step in x from first frame and create temporary frame  
      with 2 dimensions (x,trow)                                           */
   if (t==0) {
     SCDRDC(imno,"IDENT",1L,1L,72L,&iav,ident,&unit,&null);
     tnaxis = 2L;
     tstart[0] = x_limits[0];
     tstep[0] = step[0];
     tnpix[0] = (long int)((x_limits[1]-x_limits[0]) / tstep[0]) + 1;
     tstart[1] = t_limits[0];
     tstep[1] =  1.0;
     tnpix[1] = (long int)(t_limits[1]-t_limits[0]) + 1;
     SCIPUT("middumma",D_R4_FORMAT,F_O_MODE,F_IMA_TYPE,tnaxis,tnpix,
             tstart,tstep,ident,cunit,(char **)&tempima,&imnot);
   }  
   t++;       
   SCFCLO(imno); 
}


/* compute phase and sort catalog numbers by phase */

for (t = 0; t < tnpix[1]; t++) {
  mjd[t] = time[t][3];
  phase[t] = (float)((time[t][3]-epoch)/period);
  phase[t] = (float)(phase[t]-(int)phase[t]);
  if (phase[t]<0.) phase[t] += 1.;
}
sort2 (tnpix[1],(phase-1),(catnr-1));

/* read images from catalogue again sorted by phase,
   extract selected wavelength range and 
   insert in temp image                                           */   

t = 0;
for (frmcnt = t_limits[0]; frmcnt <= t_limits[1]; frmcnt++) {

   entry = (long)catnr[t]+1L;
   SCCFND(catfil,entry,frame);
   SCIGET(frame,D_R4_FORMAT,F_I_MODE,F_IMA_TYPE,1L,&naxis,npix,
        	 start,step,ident,cunit,(char **)&inima,&imno);

   /* check range */  
   if (start[0]>x_limits[0] || (start[0]+npix[0]*step[0])<x_limits[1]){
       sprintf (msg, " error: selected region out of wavelength range \n");
       SCETER(1,msg);
   }

   xstart = (long int)((x_limits[0] - start[0]) / step[0]);
   for (x = 0; x < tnpix[0]; x++) {
     *(tempima+x+t*tnpix[0]) = *(inima+xstart+x);
   }
   t++;       
   SCFCLO(imno); 
}


if (renorm) {

/* find local continuum and renorm temp_ima  */

for (t=0; t<tnpix[1]; t++){
  for (i = 0; i < MAXWIN; i++){
    leftwin[i]  = *(tempima+i+t*tnpix[0]); 
    rightwin[i] = *(tempima-(i+1)+(t+1)*tnpix[0]);
  }  
  median(leftwin-1,MAXWIN,&medlwin);
  median(rightwin-1,MAXWIN,&medrwin);
  a = (medrwin-medlwin)/(tnpix[0]-0-MAXWIN);
  b = medlwin - a * MAXWIN/2;
  for (x = 0; x < tnpix[0]; x++) {
   *(tempima+x+t*tnpix[0]) = *(tempima+x+t*tnpix[0]) / (a * x + b);      
  }
}
}

if (filter) {

/* filter median 1x3 in t direction         */


tcol = vector(1,tnpix[1]);

for (x = 0; x < tnpix[0]; x++){
  for (t = 0; t < tnpix[1]; t++) {
    *(tcol+t) = *(tempima+t*tnpix[0]+x);
  }

  for (t = HMAXMED; t < tnpix[1]-HMAXMED; t++) {
    for (i = 0; i < 2*HMAXMED+1; i++){
      medwin[i] = *(tcol+t+i-HMAXMED);
    }
    median(medwin-1,2*HMAXMED+1,&med);
    *(tempima+t*tnpix[0]+x) = med;
  }
}

free_vector(tcol,1,tnpix[1]);
}


/* average spectra to phase-bins and write to result image  */

/* create and open result image             */
rnaxis = 2L;
rstart[0] = tstart[0];
rstep[0] =  tstep[0];
rnpix[0] =  tnpix[0];
rstart[1] = 0.0;
rstep[1] =  (double)(1./phase_bin);
rnpix[1] =  phase_bin;
strcpy(cunit,"rel.flux        wavelength      phase           ");
SCIPUT(rframe,D_R4_FORMAT,F_O_MODE,F_IMA_TYPE,rnaxis,rnpix,
       rstart,rstep,ident,cunit,(char **)&resima,&imnor);

/* init result images with BLACK pixels */
for (x = 0; x < tnpix[0]; x++){
 for (r = 0; r < rnpix[1]; r++){
  *(resima+r*rnpix[0]+x) = BLACK;
 }
}

data = vector (1,rnpix[1]*rnpix[0]);
n = 1;

tcol = vector(1,rnpix[1]);
  yy = ivector(1,rnpix[1]);

for (x = 0; x < tnpix[0]; x++){

  for (nbin = 0; nbin < phase_bin; nbin++) {
    *(tcol+nbin) = 0.0;
    *(yy+nbin) = 0L;
  }

  for (t = 0; t < tnpix[1]; t++){
    nbin = (long int)((phase[t]-rstart[1])/rstep[1]);
    *(tcol+nbin) = *(tcol+nbin)+*(tempima+t*tnpix[0]+x);
    *(yy+nbin) = *(yy+nbin) +1L;
  }

  for (r = 0; r < rnpix[1]; r++){
    if (*(yy+r) > 0){
     out = *(tcol+r)/(float)*(yy+r);
     *(resima+r*rnpix[0]+x) = out;
     data[n] = (float)out; n++;
    } else {
     if (first) {
      sprintf (msg, "Warning: empty bin at phase = %f", rstart[1]+r*rstep[1]); 
      SCTPUT (msg);
    }
   }
  }
 first = 0;
}
free_ivector(yy,1,rnpix[1]);
free_vector(tcol,1,rnpix[1]);

/* compute and write additional descriptors of image header */
SCDWRC(imnor,"HISTORY",1L,history,1L,64L,&uni);
SCDWRI(imnor,"N_SPEC",tnpix+1,1L,1L,&uni);
SCDWRI(imnor,"CATINTV",t_limits,1L,2L,&uni);
SCDWRD (imnor,"XINTV",x_limits,1L,2L,&uni);
SCDWRD(imnor,"PERIOD",&period,1L,1L,&uni);
SCDWRD(imnor,"EPOCH",&epoch,1L,1L,&uni);
statistics (data, n-1, &ave, &adev, &sdev, &var, &min, &max);
lhcuts[0] = ave-3*sdev;
lhcuts[1] = ave+3*sdev;
lhcuts[2] = min;
lhcuts[3] = max;
if (lhcuts[0] < min) lhcuts[0] = min;
if (lhcuts[1] > max) lhcuts[1] = max;
SCDWRR (imnor,"LHCUTS",lhcuts,1L,4L,&uni);


/* init table for overlay */
TCTINI (rframe, F_TRANS, F_IO_MODE, 3L, 200L, &tid);
TCCINI (tid, D_R4_FORMAT, 1L, "F8.3", "  ", "x", &xocol);
TCCINI (tid, D_R4_FORMAT, 1L, "F8.3", "  ", "t", &tocol);
TCCINI (tid, D_R4_FORMAT, 1L, "I4", "  ", "ident", &idocol);

/* write x and t for overlay to table */
xo2 = rstart[0] + rnpix[0]*rstep[0];
xo1 = xo2 + rnpix[0]*rstep[0]/50.;

tt = 1;
for (t = 0; t < tnpix[1]; t++) {
  dphase = (double)phase[t];
  entry = catnr[t];
  TCEWRD(tid,tt,tocol,&dphase);
  TCEWRD(tid,tt,xocol,&xo1);
  TCEWRI(tid,tt,idocol,&entry);
  TCEWRD(tid,tt+1,tocol,&dphase);
  TCEWRD(tid,tt+1,xocol,&xo2);
  TCEWRI(tid,tt+1,idocol,&entry);
  tt += 2;
}

SCSEPI();

} /* end of main */




/*  subroutines   */          

void statistics(float data[], int n, float *ave, float *adev, float *sdev,
	float *var, float *min, float *max)
{
	void nrerror(char error_text[]);
	int j;
	float ep=0.0,s,p;
	char msg[20];

	if (n <= 1) nrerror("n must be at least 2 in moment");
	s=0.0;
	for (j=1;j<=n;j++){
          if (data[j] > *max) *max = data[j];
          if (data[j] < *min) *min = data[j];
	  s += data[j];
        }
	*ave=s/n;
	*adev=(*var)=0.0;
	for (j=1;j<=n;j++) {
		*adev += fabs(s=data[j]-(*ave));
		*var += (p=s*s);
	}
	*adev /= n;
	*var=(*var-ep*ep/n)/(n-1);
	*sdev=(float)sqrt((double)*var);
}



#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;
#define M 7
#define NSTACK 50

void sort2(unsigned long n, float arr[], float brr[])
{
	unsigned long i,ir=n,j,k,l=1;
	int *istack,jstack=0;
	float a,b,temp;

	istack=ivector(1,NSTACK);
	for (;;) {
		if (ir-l < M) {
			for (j=l+1;j<=ir;j++) {
				a=arr[j];
				b=brr[j];
				for (i=j-1;i>=1;i--) {
					if (arr[i] <= a) break;
					arr[i+1]=arr[i];
					brr[i+1]=brr[i];
				}
				arr[i+1]=a;
				brr[i+1]=b;
			}
			if (!jstack) {
				free_ivector(istack,1,NSTACK);
				return;
			}
			ir=istack[jstack];
			l=istack[jstack-1];
			jstack -= 2;
		} else {
			k=(l+ir) >> 1;
			SWAP(arr[k],arr[l+1])
			SWAP(brr[k],brr[l+1])
			if (arr[l+1] > arr[ir]) {
				SWAP(arr[l+1],arr[ir])
				SWAP(brr[l+1],brr[ir])
			}
			if (arr[l] > arr[ir]) {
				SWAP(arr[l],arr[ir])
				SWAP(brr[l],brr[ir])
			}
			if (arr[l+1] > arr[l]) {
				SWAP(arr[l+1],arr[l])
				SWAP(brr[l+1],brr[l])
			}
			i=l+1;
			j=ir;
			a=arr[l];
			b=brr[l];
			for (;;) {
				do i++; while (arr[i] < a);
				do j--; while (arr[j] > a);
				if (j < i) break;
				SWAP(arr[i],arr[j])
				SWAP(brr[i],brr[j])
			}
			arr[l]=arr[j];
			arr[j]=a;
			brr[l]=brr[j];
			brr[j]=b;
			jstack += 2;
			if (jstack > NSTACK) nrerror("NSTACK too small in sort2.");
			if (ir-i+1 >= j-l) {
				istack[jstack]=ir;
				istack[jstack-1]=i;
				ir=j-1;
			} else {
				istack[jstack]=j-1;
				istack[jstack-1]=l;
				l=i;
			}
		}
	}
}
#undef M
#undef NSTACK
#undef SWAP


void median(float x[],int n,float *xmed)
{
	int n2,n2p;

	sort(n,x);
	n2p=(n2=n/2)+1;
	*xmed=(n % 2 ? x[n2p] : 0.5*(x[n2]+x[n2p]));
}


#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;
#define M 7
#define NSTACK 50

void sort(unsigned long n, float arr[])
{
	unsigned long i,ir=n,j,k,l=1;
	int jstack=0,*istack;
	float a,temp;

	istack=ivector(1,NSTACK);
	for (;;) {
		if (ir-l < M) {
			for (j=l+1;j<=ir;j++) {
				a=arr[j];
				for (i=j-1;i>=1;i--) {
					if (arr[i] <= a) break;
					arr[i+1]=arr[i];
				}
				arr[i+1]=a;
			}
			if (jstack == 0) break;
			ir=istack[jstack--];
			l=istack[jstack--];
		} else {
			k=(l+ir) >> 1;
			SWAP(arr[k],arr[l+1])
			if (arr[l+1] > arr[ir]) {
				SWAP(arr[l+1],arr[ir])
			}
			if (arr[l] > arr[ir]) {
				SWAP(arr[l],arr[ir])
			}
			if (arr[l+1] > arr[l]) {
				SWAP(arr[l+1],arr[l])
			}
			i=l+1;
			j=ir;
			a=arr[l];
			for (;;) {
				do i++; while (arr[i] < a);
				do j--; while (arr[j] > a);
				if (j < i) break;
				SWAP(arr[i],arr[j]);
			}
			arr[l]=arr[j];
			arr[j]=a;
			jstack += 2;
			if (jstack > NSTACK) nrerror("NSTACK too small in sort.");
			if (ir-i+1 >= j-l) {
				istack[jstack]=ir;
				istack[jstack-1]=i;
				ir=j-1;
			} else {
				istack[jstack]=j-1;
				istack[jstack-1]=l;
				l=i;
			}
		}
	}
	free_ivector(istack,1,NSTACK);
}
#undef M
#undef NSTACK


/* (C) Copr. 1986-92 Numerical Recipes Software *1#R.+,V+%. */


#include "nrutil.c"


#undef NRANSI




