/* ts2vtab.c */

/*
++++++++++++++++++++++++++++++++++++++++++++++++++
 
.IDENTIFICATION
  program TS2TAB			version 1.00	940802
  A. Kaufer                  		LSW - Heidelberg
  1.00	   940802	

.KEYWORDS
  bulk data frames, time series
 
.PURPOSE
  convert a time series of 1-dim spectra to table with time and colums
    for each velocity step
 
.ALGORITHM
  extract input frames from catalogue,
  rebin to velocities
  renorm extracted frame 
  filter median 3x1 in t direction 
  extract columns and write to table
 
.INPUT/OUTPUT
  the following keys are used:
 
  OUT_A/C/1/60			result table
  P3/C/1/80			catalog of frames in time series
  INPUTD/D/1/1                  central wavelength 
  INPUTD/D/2/1                  lower velocity limit of interval
  INPUTD/D/3/1                  upper velocity limit of interval
  INPUTD/D/4/1                  velocity step size
  INPUTI/I/1/2                  catalog interval
  P7/C/1/1                      renorm flag
  P8/C/1/1                      filter flag  
    

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

int    t_limits[2];  
int    frmcnt, iav, iseq, x, xstart, t, tt, tmax, r, i;
int    naxis, npix[2], tnaxis, tnpix[2], rnaxis, rnpix[2];
int    imno, imnot, imnor, null, dum, noent, noelem, bytelem;
int    tid, xocol, tocol, idocol;
int    inpix;
int    unit, uni;
 
char     flag[1], dtype[1];
char     msg[80];
char     line[84];
char     frame[64], rframe[64], catfil[64];
char     cunit[16*(2+1)], ident[6], output[64], history[64];

float     *inima, *tempima, *resima, *xrow, *xint, *tcol;
float     yp1, ypn, *yy, in, out;
float     mjd[MAXCAT], medwin[2*HMAXMED+1], leftwin[MAXWIN], rightwin[MAXWIN];
float     med, medrwin, medlwin, a, b;

double   start[2], step[2], tstart[2], tstep[2], rstart[2], rstep[2];
double   x_limits[3], lam_c, v_lo, v_hi, t_step, v_step;
double   time[MAXCAT][6];
double   xo1, xo2, mjdo;
 
 
/* set up MIDAS environment + enable automatic error abort */

SCSPRO("ts2vtab");

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


/*  get the data intervals  */

SCKRDD("INPUTD",1L,1L,&iav,&lam_c,&unit,&null);
SCKRDD("INPUTD",2L,1L,&iav,&v_lo,&unit,&null);
SCKRDD("INPUTD",3L,1L,&iav,&v_hi,&unit,&null);
SCKRDD("INPUTD",4L,1L,&iav,&v_step,&unit,&null);
SCKRDI("INPUTI",1L,2L,&iav,t_limits,&unit,&null);

if (v_hi < v_lo){
  dum = v_hi;
  v_hi = v_lo;
  v_lo = dum;
  sprintf(msg, " correction: using velocity interval from %6.1f to %6.1f km/s \n",v_lo,v_hi);
  SCTPUT (msg);
}

x_limits[0] = lam_c * (1.0 + v_lo/3.0e05);
x_limits[1] = lam_c * (1.0 + v_hi/3.0e05);


/*   handle input from catalog - get name of first image file in catalog  */

strcpy(catfil,line);

SCCSHO(catfil,&noent,&dum);

if (t_limits[1] < t_limits[0]){
  dum = t_limits[1];
  t_limits[1] = t_limits[0];
  t_limits[0] = dum;
  sprintf(msg, " correction: using catalog entries %d - %d \n",t_limits[0],t_limits[1]);
  SCTPUT (msg);  
}

if (noent < t_limits[1]) {
  t_limits[1] = noent;
  sprintf(msg, " correction: using catalog entries %d - %d \n",t_limits[0],t_limits[1]);
  SCTPUT (msg);
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

/*    create temporary frame with 2 dimensions (x,trow = velocity, time)   */

tnaxis = 2L;
tstart[0] = v_lo;
tstep[0] = v_step;
tnpix[0] = (long int)((v_hi - v_lo) / v_step) + 1;
tstart[1] = t_limits[0];
tstep[1] =  1.0;
tnpix[1] = (long int)(t_limits[1]-t_limits[0]) + 1;
SCIPUT("middumma",D_R4_FORMAT,F_O_MODE,F_IMA_TYPE,tnaxis,tnpix,
        tstart,tstep,ident,cunit,(char **)&tempima,&imnot);


t = 0; 
for (frmcnt = t_limits[0]; frmcnt <= t_limits[1]; frmcnt++) {

   SCCGET(catfil,0L,frame,output,&iseq);
   if (frame[0] == ' ') break;		

   SCIGET(frame,D_R4_FORMAT,F_I_MODE,F_IMA_TYPE,1L,&naxis,npix,
        	 start,step,ident,cunit,(char **)&inima,&imno);

   /* check range */  
   if (start[0]>x_limits[0] || (start[0]+npix[0]*step[0])<x_limits[1]){
       sprintf (msg, " error: selected region out of wavelength range \n");
       SCETER(1,msg);
   }

   SCDFND (imno,"DTIME", dtype, &noelem, &bytelem);
   if (strcmp (dtype, " ")!=NULL) {
     SCDRDD(imno,"DTIME",1L,6L,&iav,&time[t][0],&unit,&null);   
   } else {
     sprintf(msg," Warning: descriptor DTIME not present in %s, using O_TIME"
	     , frame);
     SCTPUT (msg);
     SCDRDD(imno,"O_TIME",1L,6L,&iav,&time[t][0],&unit,&null);
   }
   
/*   extract selected wavelength range and insert in temp image           */   

   inpix = (long int)((x_limits[1] - x_limits[0]) / step[0]);
   xrow = vector(1,inpix);
   xint = vector(1,inpix);
   yy = vector(1,inpix);

   xstart = (long int)((x_limits[0] - start[0]) / step[0])+1;
   for (r = 0; r < inpix; r++){
    xrow[r] = 3.0e05*((start[0] + (xstart+r) * step[0])/lam_c-1.0);  
    xint[r] = *(inima+xstart+r);  
   }

   /* rebin with cubic spline interpolation */
   spline((xrow-1),(xint-1),inpix,0.0,0.0,(yy-1));
   for (x = 0; x < tnpix[0]; x++){
    in = v_lo + x * v_step;  
    splint((xrow-1),(xint-1),(yy-1),inpix,in,&out);
    *(tempima+x+t*tnpix[0]) = out;
   }
   free_vector(yy,1,inpix);
   free_vector(xrow,1,inpix);
   free_vector(xint,1,inpix);

   t++;       
   SCFCLO(imno); 
}

/* set time vector mjd                       */
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
SCDWRI (tid, "VNPIX", tnpix, 1L, 1L, &uni);
SCDWRD (tid, "VSTART",tstart,1L, 1L, &uni);
SCDWRD (tid, "VSTEP", tstep, 1L, 1L, &uni);
SCDWRD (tid, "LAMC",  &lam_c, 1L, 1L, &uni);

/*  inserting mjd as column :jd24   */ 

TCCINI (tid, D_R8_FORMAT, 1L, "G10.4", "  ", "jd24", &tocol);
for (t = 0; t < tnpix[1]; t++) {
  mjdo = (double)mjd[t];
  TCEWRD(tid,t+1,tocol,&mjdo);
}

/*  insert columns as columns :liii  */

for (x = 0; x < tnpix[0]; x++){
  sprintf (ident,"l%04i",x+1);
  TCCINI (tid, D_R4_FORMAT, 1L, "G10.5", "  ", ident, &xocol);
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


void spline(float x[], float y[], int n, float yp1, float ypn, float y2[])
{
	int i,k;
	float p,qn,sig,un,*u;
	
	u=vector(1,n-1);
	if (yp1 > 0.99e30)
		y2[1]=u[1]=0.0;
	else {
		y2[1] = -0.5;
		u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
	}
	for (i=2;i<=n-1;i++) {
		sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
		p=sig*y2[i-1]+2.0;
		y2[i]=(sig-1.0)/p;
		u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
		u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
	}
	if (ypn > 0.99e30)
		qn=un=0.0;
	else {
		qn=0.5;
		un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
	}
	y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);
	for (k=n-1;k>=1;k--)
		y2[k]=y2[k]*y2[k+1]+u[k];
	free_vector(u,1,n-1);
}

/* (C) Copr. 1986-92 Numerical Recipes Software *1#R.+,V+%. */



void splint(float xa[], float ya[], float y2a[], int n, float x, float *y)
{
	void nrerror(char error_text[]);
	int klo,khi,k,i;
	float h,b,a;

	klo=1;
	khi=n;
	while (khi-klo > 1) {
		k=(khi+klo) >> 1;
		if (xa[k] > x) khi=k;
		else klo=k;
	}
	h=xa[khi]-xa[klo];
	if (h == 0.0) nrerror("Bad xa input to routine splint");
	a=(xa[khi]-x)/h;
	b=(x-xa[klo])/h;
	*y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}



/* (C) Copr. 1986-92 Numerical Recipes Software *1#R.+,V+%. */


#include "nrutil.c"


#undef NRANSI



