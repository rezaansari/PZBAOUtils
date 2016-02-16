#include "myinteg2d.h"
#include <iostream> 
#include <math.h> 



/*!
  \ingroup NTools
  \class SOPHYA::GLInteg2D 
  
  \brief Implementation basique d'integration 2D par la methode de Gauss-Legendre.

  Le principe de la méthode est de calculer les valeurs de la 
  fonction aux zéros des polynomes de Legendre. Avec les poids
  qui vont bien, GL d'ordre n est exacte pour des polynomes de 
  degré <= 2n+1 (monome le + haut x^(2*n-1).
  Impossible de demander une précision donnée.

  \sa SOPHYA::Integrator

  \warning statut EXPERIMENTAL , NON TESTE
*/

namespace SOPHYA {

void
GLInteg1D::Global_ComputeWeights(int order, double * xPos, double* weights)
{
   const double EPS_gauleg = 5.0e-12;

   int m=(order+1)/2;
   const double xxMin = 0;
   const double xxMax = 1;
      
   double xm=0.5*(xxMax+xxMin);
   double xl=0.5*(xxMax-xxMin);
   for (int i=1;i<=m;i++)  {
      double z=cos(3.141592654*(i-0.25)/(order+0.5));
      double z1, pp;
      do {
         double p1=1.0;
         double p2=0.0;
         for (int j=1;j<=order;j++) {
            double p3=p2;
            p2=p1;
            p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
         }
         pp=order*(z*p1-p2)/(z*z-1.0);
         z1=z;
         z=z1-p1/pp;
      } while (fabs(z-z1) > EPS_gauleg);
      xPos[i-1]         = xm-xl*z;
      xPos[order-i]    = xm+xl*z;
      weights[i-1]      = 2.0*xl/((1.0-z*z)*pp*pp);
      weights[order-i] = weights[i-1];
   }
   return;
}

  

//! Retourne la valeur de l'integrale
double
GLInteg1D::ComputeIntegral(double xmin, double xmax) 
{
  if (!mXPos) ComputeWeights();
  double s = 0;
 // set<double>::iterator i=mBounds.begin();
 // set<double>::iterator j=i; j++;
 // while (j != mBounds.end()) {
   double x1 = xmin;
   double x2 = xmax;
   for(int k=0; k<mOrder; k++)
     s += mWeights[k] * (*mFunc_p)(x1 + (x2-x1)*mXPos[k]);
   //   i++; j++;
 return s*(x2-x1);
}

void
GLInteg1D::Print(int lp)
{
  // Integrator::Print(lp);
 cout<<"GLInteg1D order="<<mOrder<<endl;
 if(lp>0 && mOrder>0) {
   for(int i=0;i<mOrder;i++)
     cout<<" ("<<mXPos[i]<<","<<mWeights[i]<<")";
   cout<<endl;
 }
}

/*
void
GLInteg::ComputeBounds()
{
  mBounds.erase(mBounds.begin(), mBounds.end());
  for (int i=0; i<=mNStep; i++) 
    mBounds.insert(mXMin + (mXMax-mXMin)*i/mNStep);
}
*/
  

//! Retourne la valeur de l'integrale
double
GLInteg2D::ComputeIntegral(double xmin, double xmax, double ymin, double ymax)
{
  if (!mXPos) ComputeWeights();
  double s = 0;
  // set<double>::iterator i=mBounds.begin();
  // set<double>::iterator j=i; j++;
  // while (j != mBounds.end()) {
  double x1 = xmin;
  double x2 = xmax;
  double y1 = ymin;
  double y2 = ymax;
  for(int i=0; i<mOrder; i++)  {
    double x = x1 + (x2-x1)*mXPos[i];
    double s1 = 0;
    // On calcule l'integrale sur Y a x fixe 
    for(int j=0; j<mOrder; j++)
      s1 += mWeights[j] * (*mFunc_p)(x, y1 + (y2-y1)*mXPos[j]);
    // On cumule sur x
    s += mWeights[i]*s1;
  }
  s *= (x2-x1)*(y2-y1);
  return s;
}

//! Imprime l'ordre et la valeur des poids sur cout
void
GLInteg2D::Print(int lp)
{
  // Integrator::Print(lp);
 cout<<"GLInteg2D order="<<mOrder<<endl;
 if(lp>0 && mOrder>0) {
   for(int i=0;i<mOrder;i++)
     cout<<" ("<<mXPos[i]<<","<<mWeights[i]<<")";
   cout<<endl;
 }
}

} // Fin du namespace
