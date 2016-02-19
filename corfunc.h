/* ----
   Project   LSST/BAO/PhotoZ
   Class to compute power spectrum from from mass density of 
   galaxy count number grid (array) 
   R.Ansari for A. Choyer , Feb 2015 
                                                     -------  */

#ifndef CORFUNC_SEEN
#define CORFUNC_SEEN

#include "classfunc.h"
#include "slininterp.h"

#include "myinteg2d.h"

#include "hsplfit.h"

//--------------------------------------------------------------------------------
//-----------  Fonction Xsi1 : fonction d'auto-correlation 
//--------------------------------------------------------------------------------
class Xsi1 : public ClassFunc1D , public ClassFunc2D {
public:
  Xsi1(double A_=1., double ra=0.15, double B_=0.05, double sigma_b=0.06, double r0b=1., double C_=0., double lenc=0.075)
  {
    A=A_;  rA=ra;
    B=B_;  sigmaB=sigma_b; r0B=r0b;
    C=C_;  lenC=lenc;
    a_=1.; b_=2.; 
    AA=1./(b_*sqrt(M_PI))/4.;
  }
  Xsi1(Xsi1 const & a)
  {
    Copy(a);
  }
  Xsi1& operator= (Xsi1 const & a)
  {
    Copy(a);  return(*this);
  }
  inline double MyValue(double r) const
  {
    r/=100.;  // pour avoir une echelle similaire a la cosmo - s_bao ~ 100 Mpc
    double ra=r/rA;
    double raa=(r-a_*rA)/(b_*rA);
    double rb=(r-r0B)/sigmaB;
    return A*(exp(-ra)-AA*exp(-raa*raa))+B*(exp(-rb*rb));

   //    double rbb=(r-b*r0B)/(b_*sigmaB);
   //    double rc=r/lenC;
   // return (A*(exp(-ra)-AA*exp(-raa*raa))+B*(exp(-rb*rb))+C*exp(-rc));
    //    return (A*(exp(-ra)-AA*exp(-raa*raa))+B*(exp(-rb*rb)*cos(0.5*rb))+C*exp(-rc));
    //    double va=((ra>1.e-9)?sin(ra)/ra:1.);
    //    return (A*va+B*exp(-0.5*rb*rb)+C*exp(-r/lenC));

  }
  virtual double operator()(double r)  const
  {
    return MyValue(r);
  }
  virtual double operator()(double rlong, double rtrans)  const
  {
    double r = sqrt(rlong*rlong+rtrans*rtrans);
    return MyValue(r);
  }
    
  void Copy(Xsi1 const & a)
  {
    A=a.A;  rA=a.rA;
    B=a.B;  sigmaB=a.sigmaB;   r0B=a.r0B;
    C=a.C;  lenC=a.lenC;
    a_=a.a_;  b_=a.b_;  AA=a.AA;
  }
  Vector FillVec(double rmax=300., int N=200)
  {
    Vector vxsi(N);
    double dr=rmax/(double)N;
    for(sa_size_t i=0; i<N; i++) {
      double r=(double)i*dr;
      vxsi(i)=MyValue(r);
    }
    vxsi.Info()["rmax"]=rmax;
    return vxsi;
  }
  Histo2D FillHisto2D(double rmax=300., int N=100)
  {
    Histo2D h2(-rmax,rmax,N,-rmax,rmax,N);
    cout<<"Xsi1::FillHisto() ..."<<endl;
    //    ProgressBar pgb(h2.NBinX());
    for(int i=0; i<h2.NBinX() ; i++) {
      // pgb.update(i);
      for(int j=0; j<h2.NBinY() ; j++) {
	double rtrans, rlong; 
	h2.BinCenter(i, j, rtrans, rlong);
	double valxsi=this->operator()(rlong, rtrans);
	h2.Add(rtrans, rlong, valxsi);
      }
    }
    return h2;
  }

  double A,rA, B, sigmaB, r0B, C, lenC;
  double a_,b_,AA;
};

//--------------------------------------------------------------------------------
//---- Enrobage d'une fonction 1D (P(k)/Xsi) pour en faire l'argument d'integration
//---  en klong, ktrans pour calculer la fonction d'auto-correlation 2D 
//--------------------------------------------------------------------------------

class FctPkXsi4Integ : public ClassFunc1D {
public:
  FctPkXsi4Integ(ClassFunc1D const& pkxsi, double r_ou_k)
    : pkxsi_(pkxsi), rok_(r_ou_k)
  { }
  virtual double operator()(double kor)  const
  {
    double rk=rok_*kor;
    if (rk<1.e-19) return  pkxsi_(kor)*kor*kor; 
    else return pkxsi_(kor)*kor*kor*sin(rk)/rk;
  }    
  ClassFunc1D const & pkxsi_;
  double rok_;
};

//--------------------------------------------------------------------------------
//-----------  Fonction P(k) calcule a partir d'une fonction d'auto-correlation 
//--------------------------------------------------------------------------------
class PkFrXsi : public ClassFunc1D, public ClassFunc2D {
public:
  //  PkFrXsi(ClassFunc1D const & xsi, double rmax=4., int glorder=75)
  //    : xsi_(xsi), rmax_(rmax), gli1d_(xsi_,0.,rmax)
  PkFrXsi(ClassFunc1D const & xsi, double rmax=500., double dr=0.1) : xsi_(xsi), rmax_(rmax), dr_(dr)
  {
    //    gli1d_.SetOrder(glorder);
    //    C_=2./sqrt(2.*M_PI);
    C_=2.*dr_/sqrt(2.*M_PI);
  }
  inline double MyValue(double k) const
  {
    if (k<0.)  throw ParmError("PkFrXsi::operator()(double k)  k<0 !");
    double rv=0;
    if (k<1.e-9) {
      for(double r=dr_; r<rmax_; r+=dr_)  rv+=xsi_(r)*r*r; 
    }
    else {
      for(double r=dr_; r<rmax_; r+=dr_)  rv+=xsi_(r)*(sin(k*r))*r/k; 
    }
    return (C_*rv); 
    /*
    FctPkXsi4Integ fintxsi(xsi_, k);
    gli1d_.SetFunction(fintxsi);
    return C_*gli1d_.Value();
    */
  }
  virtual double operator()(double k)  const
  {
    return MyValue(k);
  }
  virtual double operator()(double klong, double ktrans)  const
  {
    double k = sqrt(klong*klong+ktrans*ktrans);
    return MyValue(k);
  }
  Vector FillVec(double kmax=0.5, int N=200)
  {
    Vector vpk(N);
    double dk=kmax/(double)N;
    for(sa_size_t i=0; i<N; i++) {
      double k=(double)i*dk;
      vpk(i)=MyValue(k);
    }
    vpk.Info()["kmax"]=kmax;
    return vpk;
  }
  Histo2D FillHisto(double kmax=0.5, int N=100)
  {
    cout<<"PkFrXsi::FillHisto() ..."<<endl;
    Histo2D h2(-kmax,kmax,N,-kmax,kmax,N);
    ProgressBar pgb(h2.NBinX());
    for(int i=0; i<h2.NBinX() ; i++) {
      pgb.update(i);
      for(int j=0; j<h2.NBinY() ; j++) {
	double klong,ktrans; 
	h2.BinCenter(i, j, ktrans, klong);
	double valpk=this->operator()(klong, ktrans);
	h2.Add(ktrans, klong, valpk);
      }
    }
    return h2;
  }

  ClassFunc1D const & xsi_;
  double rmax_, dr_;
  //  mutable GLInteg1D gli1d_;
  double C_;
};

//-----------------------------------------------------------------------------------------
//-----------  Fonction P(k) interpole calcule a partir d'un spectre de puissance P(k) 
//-----------------------------------------------------------------------------------------
class InterpPk : public SLinInterp1D, public ClassFunc2D {
public:
  InterpPk(ClassFunc1D const & pk, double kmin=0., double kmax=2., int npt=2000, bool fgckneg=true)
  {
    double dk=(kmax-kmin)/(double)npt;
    vector<double> vpk(npt+1);
    double k=kmin;
    int cntneg=0;
    for(size_t i=0; i<=npt; i++) {
      vpk[i]=pk(k); // +1.e-4*(1.-k/100.);
      if (fgckneg&&(vpk[i]<0.)) {
	vpk[i]=0.;  cntneg++;
      }
      kmax=k;  k+=dk;
    }
    DefinePoints(kmin,kmax,vpk);
    if (fgckneg && cntneg>0)
      cout << " InterpPk() - negative P(k) set to zero , NegCount="<<cntneg<<" /NPT="<<npt+1<<endl;
  }
  InterpPk(string const & pkfilename)
  {
    cout << " InterpPk() - Creating from (k Pk) pairs from input file"<<pkfilename;
    ReadXYFromFile(pkfilename);
    //DBG    Print(1);
  }
  virtual double operator()(double k)  const
  {
    return YInterp(k);
    //    return SLinInterp1D::operator()(k);
  }
  virtual double operator()(double klong, double ktrans)  const
  {
    double k = sqrt(klong*klong+ktrans*ktrans);
    return YInterp(k);
    //    return SLinInterp1D::operator()(k);
  }
  Vector FillVec(double kmax=0.5, int N=200)
  {
    Vector vpk(N);
    double dk=kmax/(double)N;
    for(sa_size_t i=0; i<N; i++) {
      double k=(double)i*dk;
      vpk(i)=this->operator()(k);
    }
    vpk.Info()["kmax"]=kmax;
    return vpk;
  }  
  Histo2D FillHisto(double kmax=0.5, int N=100)
  {
    cout<<"InterpPk::FillHisto() ..."<<endl;
    Histo2D h2(-kmax,kmax,N,-kmax,kmax,N);
    //    ProgressBar pgb(h2.NBinX());
    for(int i=0; i<h2.NBinX() ; i++) {
      //      pgb.update(i);
      for(int j=0; j<h2.NBinY() ; j++) {
	double klong,ktrans; 
	h2.BinCenter(i, j, ktrans, klong);
	double valpk=this->operator()(klong, ktrans);
	h2.Add(ktrans, klong, valpk);
      }
    }
    return h2;
  }

};


//--------------------------------------------------------------------------------
//-----------  Fonction d'auto-correlation isotrope calculee a partir d'un spectre de puissance P(k)
//--------------------------------------------------------------------------------
class XsiFrPk : public ClassFunc1D {
public:
  XsiFrPk(ClassFunc1D const & pk, double kmax=2.0, double dk=0.0005) : pk_(pk), kmax_(kmax), dk_(dk)
  {
    C_=2.*dk_/sqrt(2.*M_PI);
  }
  virtual double operator()(double r)  const
  {
    if (r<0.)  throw ParmError("XsiFrPk::operator()(double r)  r<0 !");
    double rv=0;
    if (r<1.e-9) {
      for(double k=dk_; k<kmax_; k+=dk_)  rv+=pk_(k)*k*k; 
    }
    else {
      for(double k=dk_; k<kmax_; k+=dk_)  rv+=pk_(k)*(sin(k*r))*k/r; 
    }
    return (C_*rv);
  }
  virtual double operator()(double rlong, double rtrans)  const
  {
    double r = sqrt(rlong*rlong+rtrans*rtrans);
    return this->operator()(r);
  }

  Vector FillVec(double rmax=300., int N=200)
  {
    Vector vxsi(N);
    double dr=rmax/(double)N;
    for(sa_size_t i=0; i<N; i++) {
      double r=(double)i*dr;
      vxsi(i)=this->operator()(r);
    }
    vxsi.Info()["rmax"]=rmax;
    return vxsi;
  }
  Histo2D FillHisto(double rmax=300., int N=100)
  {
    Histo2D h2(-rmax,rmax,N,-rmax,rmax,N);
    cout<<"XsiFrPk::FillHisto() ..."<<endl;
    //    ProgressBar pgb(h2.NBinX());
    for(int i=0; i<h2.NBinX() ; i++) {
      // pgb.update(i);
      for(int j=0; j<h2.NBinY() ; j++) {
	double rtrans, rlong; 
	h2.BinCenter(i, j, rtrans, rlong);
	double valxsi=this->operator()(rlong, rtrans);
	h2.Add(rtrans, rlong, valxsi);
      }
    }
    return h2;
  }

  ClassFunc1D const & pk_;
  double kmax_, dk_;
  double C_;
};

//--------------------------------------------------------------------------------
//---- Enrobage d'une fonction 2D (P(k)) pour en faire l'argument d'integration
//---  en klong, ktrans pour calculer la fonction d'auto-correlation 2D 
//--------------------------------------------------------------------------------

class FctPk2D4Integ : public ClassFunc2D {
public:
  FctPk2D4Integ(ClassFunc2D const& pk, double rlong, double rtrans)
    : pk_(pk), rlong_(rlong), rtrans_(rtrans)
  { }
  virtual double operator()(double klong, double ktrans)  const
  {
    return pk_(klong, ktrans)*j0(ktrans*rtrans_)*ktrans*cos(klong*rlong_);
  }    
  ClassFunc2D const & pk_;
  double rlong_, rtrans_;
};


//--------------------------------------------------------------------------------
//-----------  Fonction d'auto-correlation anisotrope calculee a partir d'un spectre de puissance P(k//,k_|_)
//--------------------------------------------------------------------------------
class XsiFrPk2D : public ClassFunc2D {
public:
  XsiFrPk2D(ClassFunc2D const & pk, double kmax=2., int glorder=75) : pk_(pk), kmax_(kmax), glorder_(glorder) { }
  
  inline void setGLOrder(int glorder=75) { glorder_=glorder; }
  inline int  getGLOrder() { return glorder_; }
  inline void setKmax(double kmax=200.) { kmax_=kmax; }
  inline double getKmax() { return kmax_; }

  virtual double operator()(double rlong, double rtrans)  const
  {
    FctPk2D4Integ fun2(pk_, rlong, rtrans);
    GLInteg2D gli2d(fun2, 0., kmax_, 0., kmax_);
    gli2d.SetOrder(glorder_);
    double rv=gli2d.Value();
    /*    
    // Integrale sur klong
    for(double kl=0.; kl<kmax_; kl+=dk_) {
      double tki=0.;
      // Integrale sur ktrans
      for(double kt=0.; kt<kmax_; kt+=dk_) {
	tki += pk_(kl,kt)*j0(kt*rtrans)*kt;
      }
      rv+=(tki*dk_*cos(kl*rlong));
    }
    return (2.*rv*dk_)/(sqrt(M_PI));
    */
    return (2.*rv)/(sqrt(M_PI));
  }

  Histo2D FillHisto(double rmax=500., int N=100)
  {
    Histo2D h2(-rmax,rmax,N,-rmax,rmax,N);
    cout<<"XsiFrPk2D::FillHisto() ..."<<endl;
    ProgressBar pgb(h2.NBinX());
    for(int i=0; i<h2.NBinX() ; i++) {
      pgb.update(i);
      for(int j=0; j<h2.NBinY() ; j++) {
	double rtrans, rlong; 
	h2.BinCenter(i, j, rtrans, rlong);
	double valxsi=this->operator()(rlong, rtrans);
	h2.Add(rtrans, rlong, valxsi);
      }
    }
    return h2;
  }
  
  ClassFunc2D const & pk_;
  double kmax_;
  int glorder_;
};

/* --------------------  
class FctPk4Integ : public ClassFunc2D {
public:
  FctPk4Integ(ClassFunc2D const& pk, double rlong)
    : pk_(pk), rlong_(rlong)
  { }
  virtual double operator()(double klong, double ktrans)  const
  {
    return pk_(klong, ktrans)*ktrans*cos(klong*rlong_);
  }    
  ClassFunc2D const & pk_;
  double rlong_;
};
-------------------------- */

//-----------------------------------------------------------------------------------------
//-----------  Fonction Pk1(k) : un spectre de puissance 
//-----------------------------------------------------------------------------------------
class Pk1 : public ClassFunc1D {
public:
  //  Pk1(double A_=1.e-2, double B_=2.e-3, double db=0.15, double C_=1e-3, double fc=1., double dc=0.2)
  //  Pk1(double A_=3.e-3, double B_=1.e-2, double db=2., double C_=1e-3, double fc=0.05, double dc=0.2)
  Pk1(double A_=0., double B_=1., double db=1., double C_=0.9, double fc=40., double dc=0.)

  {
    A=A_;  B=B_;  dB=db;  C=C_;   fC=fc;  dC=dc;
  }
  Pk1(Pk1 const & a)
  {
    Copy(a);
  }
  Pk1& operator= (Pk1 const & a)
  {
    Copy(a);  return(*this);
  }
  virtual double operator()(double k)  const
  {
    // return (A+B*k/(1.+exp(dB*k))+C*exp(-dC*k)*sin(fC*k));
    // return (A+B/(1.+k/dC)+C*exp(-dC*k)*sin(fC*k));
    //    return (A+B/(1.+k/dC)+C*sin(fC*k));
    // return (A+C*sin(fC*k));
    double s=sin(fC*k);
    //    return (A*100.+B*0.1/(1.e-5+k*k*sqrt(k))+C*100.*exp(-dC*k)*s*s);
    // return (A*100.+B*0.1/(1.e-5+k*k*sqrt(k))+C*100.*exp(-dC*k)*s);
    return (A*100.+B*0.1/(1.e-5+k*k*sqrt(k)))*(1.+C*exp(-dC*k)*s);
    //    return (A*100.+B*0.1/(1.e-6+k*k*k))*(1.+C*exp(-dC*k)*s);

    // return (A+B*0.1/(1.e-6+k*k*k)+C*exp(-dC*k)*sin(fC*k));

  }
  void Copy(Pk1 const & a)
  {
    A=a.A;  B=a.B;   dB=a.dB;  
    C=a.C;   dC=a.dC;  fC=a.fC; 
  }
  double A, B, dB, C, dC, fC;
};


#endif 
