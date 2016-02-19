// -------------------------------------------------------------
//       Fit spline a un histo 2D 
//--------------------------------------------------------------
#ifndef HSPLFIT_H_SEEN
#define HSPLFIT_H_SEEN

#include "machdefs.h"
#include "cspline.h"

namespace SOPHYA {
class H2DSpline : public ClassFunc2D {
public:
  H2DSpline(Histo2D & h2d, bool fgexy=false) : fgexy_(fgexy) // : csp2_()
  {
    h2d.GetXCoor(lesx1);
    h2d.GetYCoor(lesx2); 
    h2d.GetValue(lesv);
    // On passe par un pointeur pour controler l'appel du destructeur
    // il y a en effet plantage lors de l'appel du destruteur _ pb avec des free 
    csp2_ = new CSpline2 ; 
    csp2_->SetNewTab(lesx1.Size(), lesx1.Data(), lesx2.Size(), lesx2.Data(), lesv.Data(), true, true);
    csp2_->ComputeCSpline();
  }
  ~H2DSpline()
  {
    // cout << " ----*DBG*---- Getting in  ~H2DSpline()"<<endl;
    // delete csp2_ ;  on ne fait pas le delete - ca fait planter !
  } 
  virtual double operator()(double x, double y)  const
  {
    if (fgexy_) { double tmp=x; x=y; y=tmp; }
    return csp2_->CSplineInt(x,y);
  }
  void FillHisto(Histo2D & h2)
  {
    cout<<"H2DSpline::FillHisto() ..."<<endl;
    //    ProgressBar pgb(h2.NBinX());
    for(int i=0; i<h2.NBinX() ; i++) {
      //      pgb.update(i);
      for(int j=0; j<h2.NBinY() ; j++) {
	double x, y; 
	h2.BinCenter(i, j, x, y);
	double val=this->operator()(x, y);
	h2.Add(x,y,val);
      }
    }
    return;
  }

protected:
  mutable CSpline2 * csp2_;
  mutable Vector lesx1, lesx2;
  mutable Matrix lesv;
  bool fgexy_;  // true -> exchange x,y with respect to the Histo2D

};

class H1DLinInterpol : public ClassFunc1D {
public:
  H1DLinInterpol(Histo & h)  : h_(h)
  {
  }
  ~H1DLinInterpol()
  {
  } 
  virtual double operator()(double x)  const
  {
    int ibx = (x-h_.XMin())/h_.BinWidth ();
    if (ibx<0)  ibx=0;
    if (ibx >= h_.NBins())  ibx=h_.NBins()-1;
    return h_(ibx);
  }
  void FillHisto(Histo & h)
  {
    cout<<"H1DLinInterpol::FillHisto() ..."<<endl;
    //    ProgressBar pgb(h2.NBinX());
    for(int i=0; i<h.NBins() ; i++) {
      double x; 
      x=h.BinCenter(i); 
      double val=this->operator()(x);
      h.Add(x,val);
    }
    return;
  }

protected:
  Histo & h_;
};

class H2DLinInterpol : public ClassFunc2D {
public:
  H2DLinInterpol(Histo2D & h2d, bool fgexy=false)  : h2_(h2d), fgexy_(fgexy)
  {
  }
  ~H2DLinInterpol()
  {
  } 
  virtual double operator()(double x, double y)  const
  {
    if (fgexy_) { double tmp=x; x=y; y=tmp; }
    int ibx = (x-h2_.XMin())/h2_.WBinX();
    int jby = (y-h2_.YMin())/h2_.WBinY();
    if (ibx<0)  ibx=0;
    if (jby<0)  jby=0;
    if (ibx >= h2_.NBinX())  ibx=h2_.NBinX()-1;
    if (jby >= h2_.NBinY())  jby=h2_.NBinY()-1;
    return h2_(ibx,jby);
  }
  void FillHisto(Histo2D & h2)
  {
    cout<<"H2DLinInterpol::FillHisto() ..."<<endl;
    //    ProgressBar pgb(h2.NBinX());
    for(int i=0; i<h2.NBinX() ; i++) {
      //      pgb.update(i);
      for(int j=0; j<h2.NBinY() ; j++) {
	double x, y; 
	h2.BinCenter(i, j, x, y);
	double val=this->operator()(x, y);
	h2.Add(x,y,val);
      }
    }
    return;
  }

protected:
  Histo2D & h2_;
  bool fgexy_;  // true -> exchange x,y with respect to the Histo2D
};

} // Fin du namespace

#endif
