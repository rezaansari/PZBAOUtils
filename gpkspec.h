/* ----
   Project   LSST/BAO/PhotoZ
   Class to compute power spectrum from from mass density of 
   galaxy count number grid (array) 
   R.Ansari for A. Choyer , Feb 2015 
                                                     -------  */

#ifndef GPKSPEC_SEEN
#define SPECPK_SEEN

#include "machdefs.h"      
#include "sopnamsp.h"       
#include <math.h>
#include <iostream>
#include <vector>
#include <string> 

#include "classfunc.h"  
#include "array.h"         
#include "histats.h"       
#include "fftwserver.h"    
#include "randinterf.h"      



//--- Change this to r_8 if one needs to work with double precision arrays 
#define TF  r_4 

// -- GFour3DPk class :  3D fourier amplitudes and power spectrum 
class GFour3DPk {
public:
// Constructor
  GFour3DPk(TArray< TF > & mgrid);
  virtual ~GFour3DPk(); 

  // Set the grid cell size (in Mpc) 
  void SetGridCellSize(double dx=1., double dy=1., double dz=1., bool fgprt=false);

  // Define the print level 0,1,2 ...
  inline int SetPrtLevel(int lev=0, int prtmod=10) 
  { int olev=prtlev_; prtlev_=lev; prtmodulo_=prtmod; return olev; }

// Return the total number of modes , and modes along each axis  
  inline sa_size_t TotN_k() const { return fourAmp.Size(); }
  inline sa_size_t N_kX()   const { return fourAmp.SizeX(); }
  inline sa_size_t N_kY()   const { return fourAmp.SizeY(); }
  inline sa_size_t N_kZ()   const { return fourAmp.SizeZ(); }

// return the grid cell size 
  inline double getdX() const {  return dx_; } 
  inline double getdY() const {  return dy_; } 
  inline double getdZ() const {  return dz_; } 

// return the cell/step of Fourier modes 
  inline double getdkX() const {  return dkx_; } 
  inline double getdkY() const {  return dky_; } 
  inline double getdkZ() const {  return dkz_; } 

// Return the fourier amplitude matrix  
  TArray< complex<TF> > GetFourierAmp()  { return fourAmp; }

// Return the reconstructed power spectrum as a profile histogram   
  HProf ComputePk(int nbin=100, double kmin=0., double kmax=-1., bool fgmodcnt=false);

  // Fills a data table from the computed P(k) profile histogram and mode count 
  Histo FillPkDataTable(DataTable& dt, double rfac=1.);
  inline HProf& getPk() { return *hp_pk_p_; }

  // perform the FFT of the input array, it is called by ComputePk() if not called already
  void doFFT();
  // Changes mass/ngal grid cells < thr to thr 
  size_t CleanNegatives(TF thr=-1.);

protected:
  void  HisPkCumul();

  // member attribute
  TArray< TF > mgal_grid_;          // drho/rho or galaxy count array (grid)
  TArray< complex<TF> > fourAmp;    // complex array of fourier coefficients
  double dx_, dy_, dz_;       // mgal_grid_ cell size
  double dkx_, dky_, dkz_;    // fourier mode steps 
  int prtlev_;
  int prtmodulo_;
  // Profile histograms for power spectrum and number of modes 
  HProf* hp_pk_p_;
  Histo* hmcnt_p_;
  Histo* hmcntok_p_;
  // double s2cut_;    for later 
};


#endif 
