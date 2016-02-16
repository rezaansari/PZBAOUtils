/* ----
   Project   LSST/BAO/PhotoZ
   Class to compute power spectrum from from mass density of 
   galaxy count number grid (array) 
   R.Ansari for A. Choyer , Feb 2015 
                                                     -------  */


#include "gpkspec.h"
#include "randr48.h"      
#include "ctimer.h"      
#include "fftwserver.h"
#include "randfmt.h"

#include "corfunc.h"

#define DeuxPI 2.*M_PI 

//---------------------------------------------------------------
// -- GFour3DPk class : power spectrum computation from a 3D-grid 
//---------------------------------------------------------------
// Constructeur a partir du tableau d rho/rho ou n-gal 
GFour3DPk::GFour3DPk(TArray< TF > & mgrid)
: mgal_grid_(mgrid) 
{
  SetPrtLevel();
  // Get the d rho/rho , n_gal array cell size
  double dx = mgal_grid_.Info().GetD("DX",1.);
  double dy = mgal_grid_.Info().GetD("DY",1.);
  double dz = mgal_grid_.Info().GetD("DZ",1.);
  cout<<"GFour3DPk() - calling SetGridCellSize() with cell size from the input grid Info() object ..."<<endl;
  SetGridCellSize(dx, dy, dz, true);
  hp_pk_p_=NULL;  hmcnt_p_=NULL;  hmcntok_p_=NULL;
  hpk2_p_=h2mcnt_p_=NULL;
}

// Destructor
GFour3DPk::~GFour3DPk()
{
  if (hp_pk_p_) delete hp_pk_p_;
  if (hmcnt_p_) delete hmcnt_p_;
  if (hmcntok_p_) delete hmcntok_p_;
}

void GFour3DPk::SetGridCellSize(double dx, double dy, double dz, bool fgprt)
{
  dx_=dx; dy_=dy;  dz_=dz;
  dkx_ = DeuxPI/(mgal_grid_.SizeX()*dx);
  dky_ = DeuxPI/(mgal_grid_.SizeY()*dy);
  dkz_ = DeuxPI/(mgal_grid_.SizeZ()*dz);
  if ((prtlev_>0)||(fgprt))  {
    cout<<"GFour3DPk::SetGridCellSize() dxyz="<<dx_<<"x"<<dy_<<"x"<<dz_<<" --> dkxyz="
	<<dkx_<<"x"<<dky_<<"x"<<dkz_<<endl;
  }
  return;
}

// Generate mass field Fourier Coefficient
void GFour3DPk::doFFT()
{
  FFTWServer ffts;
  // ATTENTION : il faut qu'on soit normalize a true pour avoir la bonne normalisation (c'est le mode par defaut)
  ffts.setNormalize(true);
  ffts.FFTForward(mgal_grid_, fourAmp);
  if (prtlev_>1)
    cout << " GFour3DPk::doFFT done ..." << endl;
  return;
}

// Generate mass field from Fourier Coefficient
void GFour3DPk::doInverseFFT(bool fgNOk0)
{
  if (fgNOk0) fourAmp(0,0,0)=complex<TF>(0.,0.);
  FFTWServer ffts;
  // ATTENTION : il faut qu'on soit normalize a true pour avoir la bonne normalisation (c'est le mode par defaut)
  ffts.setNormalize(true);
  ffts.FFTBackward(fourAmp, mgal_grid_, true);
  if (prtlev_>1)
    cout << " GFour3DPk::doInverseFFT done ..." << endl;
  return;
}

// Generate mass field Fourier Coefficient
void GFour3DPk::generateFourierAmp(ClassFunc1D & pk, double sigma_z)
{
  if (! fourAmp.IsAllocated())  doFFT();
  FMTRandGen rg;
  bool fgsmg=false;  // True, gaussian smoothing 
  double dsig2=0.;
  double Asmooth=1.;
  if (sigma_z>1.e-19) {
    dsig2 = 0.5*sigma_z*sigma_z;
    Asmooth=1./sqrt(2.*M_PI);
    fgsmg=true;
  }
  // fourAmp represent 3-D fourier transform of a real input array. 
  // The second half of the array along Y and Z contain negative frequencies
  double kxx, kyy, kzz;
  // sa_size_t is large integer type
  cout << " Four3DPk::generateFourierAmp/Info : generating fourier coefficients ..."<<endl;
  ProgressBar pgb(fourAmp.SizeZ());
  // We ignore 0th term in all frequency directions ...
  for(sa_size_t kz=0; kz<fourAmp.SizeZ(); kz++) {
    kzz =  (kz > fourAmp.SizeZ()/2) ? (double)(fourAmp.SizeZ()-kz)*dkz_ : (double)kz*dkz_; 
    for(sa_size_t ky=0; ky<fourAmp.SizeY(); ky++) {
      kyy =  (ky > fourAmp.SizeY()/2) ? (double)(fourAmp.SizeY()-ky)*dky_ : (double)ky*dky_; 
      for(sa_size_t kx=0; kx<fourAmp.SizeX(); kx++) {  // ignore the 0th coefficient (constant term)
	kxx=(double)kx*dkx_;
	complex<TF> za = fourAmp(kx, ky, kz);
	//	if (za.real()>8.e19) continue;
	double wk = sqrt(kxx*kxx+kyy*kyy+kzz*kzz);
        double amp = sqrt(pk(wk)/2.);
	if (fgsmg)  amp *= (Asmooth*exp(-(kzz*kzz*dsig2)));
	fourAmp(kx,ky,kz) = complex<TF>(rg.Gaussian(amp), rg.Gaussian(amp));
      }
    }
    pgb.update(kz);
  }
  return;
}

// Compute power spectrum as a function of wave number k 
// cells with amp^2=re^2+im^2>s2cut are ignored
// Output : power spectrum (profile histogram)
HProf GFour3DPk::ComputePk(int nbin, double kmin, double kmax, bool fgmodcnt)
{
  // The second half of the array along Y (matrix rows) contain
  // negative frequencies
  //  int nbh = sqrt(fourAmp.SizeX()*fourAmp.SizeX()+fourAmp.SizeY()*fourAmp.SizeY()/4.+fourAmp.SizeZ()*fourAmp.SizeY()/4.);
  // The profile histogram will contain the mean value of FFT amplitude
  // as a function of wave-number k = sqrt((double)(kx*kx+ky*ky))
  if ((kmax<0.)||(kmax<kmin)) {
    kmin=0.;
    double maxx=fourAmp.SizeX()*dkx_;
    double maxy=fourAmp.SizeY()/2*dky_;
    double maxz=fourAmp.SizeZ()/2*dkz_;
    kmax=sqrt(maxx*maxx+maxy*maxy+maxz*maxz);
  }
  if (nbin<2) nbin=100;
  if (hp_pk_p_) delete hp_pk_p_;
  hp_pk_p_ = new HProf(kmin, kmax, nbin);
  hp_pk_p_->SetErrOpt(false);
  if (fgmodcnt) {
    if (hmcnt_p_) delete hmcnt_p_;
    hmcnt_p_ = new Histo(kmin, kmax, nbin);
    if (hmcntok_p_) delete hmcntok_p_;
    hmcntok_p_ = new Histo(kmin, kmax, nbin);
  }
  if (! fourAmp.IsAllocated())  doFFT();
  if (prtlev_>1)  
    cout << " GFour3DPk::ComputePk()  calling HisPkCumul() kmin,kmax,nbin="
	 <<kmin<<","<<kmax<<","<<nbin<<" ..."<<endl;
  HisPkCumul();
  return *hp_pk_p_;
}

Histo2D GFour3DPk::ComputePk2D(int nbin, double kmax)
{
  if (kmax<1.e-19)  {
    double maxx=fourAmp.SizeX()*dkx_;
    double maxy=fourAmp.SizeY()/2*dky_;
    double maxz=fourAmp.SizeZ()/2*dkz_;
    kmax=sqrt((maxx*maxx+maxy*maxy+maxz*maxz)*0.5);
  }
  if (hpk2_p_) delete hpk2_p_;
  hpk2_p_ = new Histo2D(-kmax,kmax,nbin,-kmax,kmax,nbin);
  if (h2mcnt_p_)  delete h2mcnt_p_;
  h2mcnt_p_ = new Histo2D(-kmax,kmax,nbin,-kmax,kmax,nbin);
  
  uint_8 nmodeok=0;
  // fourAmp represent 3-D fourier transform of a real input array. 
  // The second half of the array along Y and Z contain negative frequencies
  double kxx, kyy, kzz;
  // sa_size_t is large integer type  
  // We ignore 0th term in all frequency directions ...
  for(sa_size_t kz=0; kz<fourAmp.SizeZ(); kz++) {
    kzz =  (kz > fourAmp.SizeZ()/2) ? -(double)(fourAmp.SizeZ()-kz)*dkz_ : (double)kz*dkz_; 
    for(sa_size_t ky=0; ky<fourAmp.SizeY(); ky++) {
      kyy =  (ky > fourAmp.SizeY()/2) ? -(double)(fourAmp.SizeY()-ky)*dky_ : (double)ky*dky_; 
      for(sa_size_t kx=0; kx<fourAmp.SizeX(); kx++) {  // ignore the 0th coefficient (constant term)
	kxx=(double)kx*dkx_;
	complex<TF> za = fourAmp(kx, ky, kz);
	//	if (za.real()>8.e19) continue;
	double ktrans = sqrt(kxx*kxx+kyy*kyy);
	double klong = kzz;
	double amp2 = norm(za); 
	if (h2mcnt_p_) {
	  h2mcnt_p_->Add(ktrans, klong);
	  h2mcnt_p_->Add(-ktrans, klong);  // pour les kx negatifs qu'on n'a pas F(-kx) = conj(F(kx)) 
	}
	hpk2_p_->Add(ktrans, klong, amp2);
	hpk2_p_->Add(-ktrans, klong, amp2);
	nmodeok++;
      }
    }
  }
  (*hpk2_p_) /= (*h2mcnt_p_);
  //  if ((prtlev_>1)||((prtlev_>0)&&(s2cut_>1.e-9))) {
  if (prtlev_>0) {
    cout << " Four3DPk::ComputePk2D/Info : NModeOK=" << nmodeok << " / NMode=" << fourAmp.Size() 
	 << " -> " << 100.*(double)nmodeok/(double)fourAmp.Size() << "%" << endl;
  }
  return (*hpk2_p_);
}

// Compute power spectrum as a function of wave number k 
// Cumul dans hp - cells with amp^2=re^2+im^2>s2cut are ignored
void GFour3DPk::HisPkCumul()
{
  uint_8 nmodeok=0;
  // fourAmp represent 3-D fourier transform of a real input array. 
  // The second half of the array along Y and Z contain negative frequencies
  double kxx, kyy, kzz;
  // sa_size_t is large integer type  
  // We ignore 0th term in all frequency directions ...
  for(sa_size_t kz=0; kz<fourAmp.SizeZ(); kz++) {
    kzz =  (kz > fourAmp.SizeZ()/2) ? (double)(fourAmp.SizeZ()-kz)*dkz_ : (double)kz*dkz_; 
    for(sa_size_t ky=0; ky<fourAmp.SizeY(); ky++) {
      kyy =  (ky > fourAmp.SizeY()/2) ? (double)(fourAmp.SizeY()-ky)*dky_ : (double)ky*dky_; 
      for(sa_size_t kx=0; kx<fourAmp.SizeX(); kx++) {  // ignore the 0th coefficient (constant term)
	kxx=(double)kx*dkx_;
	complex<TF> za = fourAmp(kx, ky, kz);
	//	if (za.real()>8.e19) continue;
	double wk = sqrt(kxx*kxx+kyy*kyy+kzz*kzz);
	double amp2 = za.real()*za.real()+za.imag()*za.imag();
	if (hmcnt_p_) hmcnt_p_->Add(wk);
	//	if ((s2cut_>1.e-9)&&(amp2>s2cut_))  continue;
	if (hmcntok_p_) hmcntok_p_->Add(wk);
	hp_pk_p_->Add(wk, amp2);
	nmodeok++;
      }
    }
  }
  //  if ((prtlev_>1)||((prtlev_>0)&&(s2cut_>1.e-9))) {
  if (prtlev_>0) {
    cout << " Four3DPk::HisPkCumul/Info : NModeOK=" << nmodeok << " / NMode=" << fourAmp.Size() 
	 << " -> " << 100.*(double)nmodeok/(double)fourAmp.Size() << "%" << endl;
  }
  return;
}

size_t GFour3DPk::CleanNegatives(TF seuil)
{
  size_t nneg = 0.;
  for(sa_size_t kz=0; kz<mgal_grid_.SizeZ(); kz++) 
    for(sa_size_t ky=0; ky<mgal_grid_.SizeY(); ky++) 
      for(sa_size_t kx=0; kx<mgal_grid_.SizeX(); kx++) 
        if (mgal_grid_(kx, ky, kz) < seuil)  {
          nneg++; mgal_grid_(kx, ky, kz)=seuil;
        }
  cout << " GFour3DPk::CleanNegatives " << nneg << " cells <" << seuil << " changed to" << seuil << endl;
  return nneg;
}



// Fills a data table from the computed P(k) profile histogram and mode count 
Histo GFour3DPk::FillPkDataTable(DataTable& dt, double rfac)
{
  if (hp_pk_p_==NULL) throw ParmError("Four3DPk::FillPkDataTable P(k) NOT computed");
  const char* nomcol[5] = {"k","Pk","nmode","nmodok","fracmodok"};
  dt.Clear();
  dt.AddDoubleColumn(nomcol[0]);    
  dt.AddDoubleColumn(nomcol[1]);    
  double Vol = (dx_*mgal_grid_.SizeX())*(dy_*mgal_grid_.SizeY())*(dz_*mgal_grid_.SizeZ());
  if (prtlev_>0) 
    cout << " GFour3DPk::FillPkDataTable()/Info: normalizing P(k) by V=dx*dy*dz="<<Vol<<" Mpc^3 x fac="<<rfac<<endl;
  Vol *= rfac;

  bool fgokmodcnt=true;
  if ((hmcnt_p_==NULL)||(hmcntok_p_==NULL)) {
    cout << " GFour3DPk::FillPkDataTable()/Warning Mode count histos NOT filled, using only P(k) ..."<<endl; 
    fgokmodcnt=false;
    HProf& hp=(*hp_pk_p_);
    DataTableRow 	dtr = dt.EmptyRow();
    for(int_4 ib=0; ib<hp.NBins(); ib++) {
      dtr[0]=hp.BinCenter(ib);
      dtr[1]=hp(ib)*Vol;
      dt.AddRow(dtr);
    }
    Histo fracmodok(hp.XMin(), hp.XMax(), hp.NBins());
    return fracmodok;
  } 
  HProf& hp=(*hp_pk_p_);
  Histo& hmcnt=(*hmcnt_p_);
  Histo& hmcntok=(*hmcntok_p_);
  Histo fracmodok=hmcntok/hmcnt;
  dt.AddIntegerColumn(nomcol[2]);    
  dt.AddIntegerColumn(nomcol[3]);    
  dt.AddFloatColumn(nomcol[4]);    
  DataTableRow 	dtr = dt.EmptyRow();
  for(int_4 ib=0; ib<hp.NBins(); ib++) {
    dtr[0]=hp.BinCenter(ib);
    dtr[1]=hp(ib)*Vol;
    dtr[2]=hmcnt(ib);
    dtr[3]=hmcntok(ib);
    dtr[4]=fracmodok(ib);
    dt.AddRow(dtr);
  }
return fracmodok;
}

