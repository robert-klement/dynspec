
/* ts2tab.c */

/*
++++++++++++++++++++++++++++++++++++++++++++++++++
 
.IDENTIFICATION
  program TS2TAB			version 1.00	930318
  A. Kaufer                  		LSW - Heidelberg
  1.00	   940318	
 
.KEYWORDS
  bulk data frames, time series
 
.PURPOSE
  convert a time series of 1-dim spectra to table with time and colums
    for each wavelength step
 
.ALGORITHM
  extract input frames from catalogue,
  renorm extracted frame 
  filter median 3x1 in t direction 
  extract columns and write to table
 
.INPUT/OUTPUT
  the following keys are used:
 
  OUT_A/C/1/60			result table
  P3/C/1/80			catalog of frames in time series
  INPUTD/D/1/1                  lower limit of x interval
  INPUTD/D/2/1                  upper limit of x interval
  INPUTI/I/1/2                  catalog interval
  P6/C/1/1                      renorm flag
  P7/C/1/1                      filter flag  
    

.VERSIONS
  1.00		first version; transcription of ts2ima.c 
----------------------------------------------------------------
*/

#include <midas_def.h>

#define NRANSI
#include "nr.h"
#include "nrutil.h"

#define MAXCAT 1000
#define MAXWIN 25
#define HMAXMED 1

void median(float x[],int n,float *xmed);


main() {

int      renorm=0, filter=0;
int      t_limits[2];  
int      frmcnt, iav, iseq, x, xstart, t, tt, tmax, r, i;
int      naxis, npix[2], tnaxis, tnpix[2], rnaxis, rnpix[2];
int      imno, imnot, imnor, null, dum, noent, noelem, bytelem;
int      tid, xocol, tocol, idocol;
int      unit, uni;
 
char     flag[1], dtype[1];
char     msg[80];
char     line[84];
char     frame[64], rframe[64], catfil[64];
char     cunit[16*(2+1)], ident[6], output[64], history[64];

float     *inima, *tempima, *resima, *tcol;
float     yp1, ypn, *yy, in, out;
float     mjd[MAXCAT], medwin[2*HMAXMED+1], leftwin[MAXWIN], rightwin[MAXWIN];
float     med, medrwin, medlwin, a, b;

double   start[2], step[2], tstart[2], tstep[2], rstart[2], rstep[2];
double   x_limits[3], lam_c, v_lo, v_hi, t_step;
double   time[MAXCAT][6];
double   xo1, xo2, mjdo;
 
 
/* set up MIDAS environment + enable automatic error abort */

SCSPRO("ts2tab");

/* get result frame, history, input frame list, flags  */

SCKGETC("OUT_A",1L,60L,&iav,rframe);
SCKGETC("HISTORY",1L,64L,&iav,history);
SCKGETC("P3",1L,80L,&iav,line);
SCKGETC("P6",1L,1L,&iav,flag);
if (strcmp(flag,"y")==0 || strcmp(flag,"Y")==0) 
 renorm=1; 
SCKGETC("P7",1L,1L,&iav,flag);
if (strcmp(flag,"y")==0 || strcmp(flag,"Y")==0) 
 filter=1; 


/*  get the data intervals  */

SCKRDD("INPUTD",1L,2L,&iav,x_limits,&unit,&null);
SCKRDI("INPUTI",1L,2L,&iav,t_limits,&unit,&null);

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
  printf(" correction: using catalog entries %d - %d \n",t_limits[0],t_limits[1]);
}

if (noent < t_limits[1]) {
  t_limits[1] = noent;
  printf(" correction: using catalog entries %d - %d \n",t_limits[0],t_limits[1]);
}

/* default: use all active entries    */
if (t_limits[0]==0 || t_limits[1]==0) {
  t_limits[0] = 1;
  t_limits[1] = noent;
}

for (frmcnt = 1; frmcnt < t_limits[0]; frmcnt++) {
   SCCGET(catfil,0L,frame,output,&iseq);
   if (frame[0] == ' ') break;		
}


t = 0; 
for (frmcnt = t_limits[0]; frmcnt <= t_limits[1]; frmcnt++) {

   SCCGET(catfil,0L,frame,output,&iseq);
   if (frame[0] == ' ') break;		

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

   
/*   extract selected wavelength range and insert in temp image           */   

   xstart = (long int)((x_limits[0] - start[0]) / step[0]);
   for (x = 0; x < tnpix[0]; x++) {
     *(tempima+x+t*tnpix[0]) = *(inima+xstart+x);
   }
   t++;       
   SCFCLO(imno); 
}


for (t = 0; t < tnpix[1]; t++) {
  mjd[t] = time[t][3];
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


/* extract columns and write to result table  */

/* create and open result table             */

TCTINI (rframe, F_TRANS, F_IO_MODE, tnpix[0]+1, tnpix[1], &tid);
SCDWRI (tid, "LNPIX", tnpix, 1L, 1L, &uni);
SCDWRD (tid, "LSTART",tstart,1L, 1L, &uni);
SCDWRD (tid, "LSTEP", tstep, 1L, 1L, &uni);
SCDWRD (tid, "XINTV",x_limits,1L,2L,&uni);
SCDWRI (tid, "CATNTV",t_limits,1L,2L,&uni);

/*  inserting mjd as column :jd24   */ 

TCCINI (tid, D_R8_FORMAT, 1L, "G10.4", "  ", "jd24", &tocol);
for (t = 0; t < tnpix[1]; t++) {
  mjdo = (double)mjd[t];
  TCEWRD(tid,t+1,tocol,&mjdo);
}

/*  insert columns as columns :liii  */

for (x = 0; x < tnpix[0]; x++){
  sprintf (ident,"l%04i",x+1);
  TCCINI (tid, D_R4_FORMAT, 1L, "G8.3", "  ", ident, &xocol);
  TCCMAP (tid, xocol, (char **)&tcol); 
  for (t = 0; t < tnpix[1]; t++) {
    *(tcol+t) = *(tempima+t*tnpix[0]+x);
  }
}


SCSEPI();

} /* end of main */




/*  subroutines   */          


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



