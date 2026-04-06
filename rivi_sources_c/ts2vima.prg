! @@ ts2ima.prg	  
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! MIDAS procedure  ts2vima.prg  
! A. Kaufer	930817
! 
! use via
! TS2VIMA out = inframes x_interval t_interval [t_step] [renorm_flag] [fil_flag]
! with inframes:  catalog.cat
!      x_interval: lambda_center,v_low,v_high,v_step  
!      t_interval: cat_entry_start,cat_entry_stop, defaulted to all 
!      v_step,t_step : vel,time_step for sampling, defaulted to 5,0.5 [km/s,d] 
!      renorm flag y/n  : renorm spectra, defaulted to y
!      filter flag y/n  : filter median 3x1 in t, defaulted to y 
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
DEFINE/PAR P1 ? IMA "Enter result frame: "
!
DEFINE/PAR P3 ? C     "Enter catalog.cat with 1-dim frames    : "
DEFINE/PAR P4 ? N     "Enter central wavelength, vel. interval lamc,vlo,vhi: "
DEFINE/PAR P5 0,0 N   "Enter valid catalog interval tlo,thi : "
DEFINE/PAR P6 5,0.5 N "Enter velocity,time step for sampling : "
DEFINE/PAR P7 Y C     "Renorm spectra [y/n] : "
DEFINE/PAR P8 Y C     "Filter median 3x1 in time direction [y/n] : "
! 
WRITE/KEY HISTORY "TS2VIMA {P3}"
WRITE/KEY OUT_A {P1}
WRITE/KEY INPUTD/D/1/5 {P4},{P6}
WRITE/KEY INPUTI/I/1/2 {P5}
!
RUN TSEXEC:TS2VIMA
!
DELETE/TEMP
DISPLAY/TS {P1} 0 2,1 "Velocity" ?
!END 








