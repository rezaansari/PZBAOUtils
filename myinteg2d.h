// -------------------------------------------------------------
//       Adaptation GLInteg a 2D , Reza , Fevrier 2016 

#ifndef MYINTEG2D_H_SEEN
#define MYINTEG2D_H_SEEN

#include "machdefs.h"
#include "classfunc.h"

namespace SOPHYA {


class GLInteg1D {
public:
  GLInteg1D(ClassFunc1D const & f, double xmin, double xmax)
    : mFunc_p(&f), mXmin(xmin), mXmax(xmax), mXPos(NULL), mWeights(NULL), mOrder(8)
  {
  }
  ~GLInteg1D()
  {
    InvalWeights();
  }

  void SetOrder(int order)
  {
    mOrder = order;
    InvalWeights();
  }
  inline  int GetOrder(void) {return mOrder;}

  inline  void SetFunction(const ClassFunc1D & f) { mFunc_p=&f; }
  inline  void SetBounds(double xmin, double xmax)
  {  mXmin=xmin; mXmax=xmax;  }

  double ComputeIntegral(double xmin, double xmax);
  double Value() 
  {
    return ComputeIntegral(mXmin, mXmax);
  }
  virtual void Print(int lp=0);
  
  static void  Global_ComputeWeights(int morder, double * mxPos, double* mweights);
protected:
  
  const ClassFunc1D * mFunc_p;
  double mXmin, mXmax; 
  double* mXPos;     // integration entre 0 et 1
  double* mWeights;  // integration entre 0 et 1
  int     mOrder;    // ordre GL
  //  set<double> mBounds; // limites des intervalles d'integration
  
  void ComputeWeights()
  {
    if (mXPos)  InvalWeights();
    mXPos    = new double[mOrder];
    mWeights = new double[mOrder];
    Global_ComputeWeights(mOrder, mXPos, mWeights);
  }
  //  void LimitsChanged();
  //  void StepsChanged();
  void InvalWeights()
  {
    if(mXPos)    delete[] mXPos;     mXPos    = NULL;
    if(mWeights) delete[] mWeights;  mWeights = NULL;
  }

  //  void ComputeBounds();
};

class GLInteg2D {
public:
  GLInteg2D(ClassFunc2D const & f, double xmin, double xmax, double ymin, double ymax)
    : mFunc_p(&f), mXmin(xmin), mXmax(xmax), mYmin(ymin), mYmax(ymax), mXPos(NULL), mWeights(NULL), mOrder(8)
  {
  }
  ~GLInteg2D()
  {
    InvalWeights();
  }
  void SetOrder(int order)
  {
    mOrder = order;
    InvalWeights();
  }
  inline  int GetOrder(void) {return mOrder;}

  inline  void SetFunction(ClassFunc2D const & f) { mFunc_p=&f; }
  inline  void SetBounds(double xmin, double xmax, double ymin, double ymax)
  {  mXmin=xmin; mXmax=xmax; mYmin=ymin; mYmax=ymax; }
  
  double ComputeIntegral(double xmin, double xmax, double ymin, double ymax);
  virtual double Value()
  {
    return ComputeIntegral(mXmin, mXmax, mYmin, mYmax);
  }
  virtual void Print(int lp=0);
  //   virtual GLInteg& AddBound(double x);
protected:
  
  const ClassFunc2D * mFunc_p;
  double mXmin, mXmax, mYmin, mYmax; 
  double* mXPos;     // integration entre 0 et 1
  double* mWeights;  // integration entre 0 et 1
  int     mOrder;    // ordre GL
  //  set<double> mBounds; // limites des intervalles d'integration
  
  void ComputeWeights()
  {
    if (mXPos)  InvalWeights();
    mXPos    = new double[mOrder];
    mWeights = new double[mOrder];
    GLInteg1D::Global_ComputeWeights(mOrder, mXPos, mWeights);
  }
  //  void LimitsChanged();
  //  void StepsChanged();
  void InvalWeights()
  {
    if(mXPos)    delete[] mXPos;     mXPos    = NULL;
    if(mWeights) delete[] mWeights;  mWeights = NULL;
  }
  //  void ComputeBounds();
};

} // Fin du namespace

#endif
