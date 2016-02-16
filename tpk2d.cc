/* ----
   Project   LSST/BAO/PhotoZ
   Tests de calcul de P(k)-2D et fct d'auto-correlation 
   R.Ansari - Feb 2015 
                                                     -------  */

//----- c/c++ includes 
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>

#include <typeinfo>

//----- sophya includes 
#include "machdefs.h"
#include "sopnamsp.h"

#include "histats.h"
#include "array.h"
#include "fitsioserver.h"
#include "fiosinit.h"


//---- this modules include files
#include "corfunc.h"
#include "myinteg2d.h"
#include "hsplfit.h"
#include "gpkspec.h"

//-------------------------------------------------------------------------
//      ------------------ MAIN PROGRAM ------------------------------
//-------------------------------------------------------------------------
int main(int narg, const char* arg[])
{
  if ((narg<3)||((narg>1)&&(strcmp(arg[1],"-h")==0))) {
    cout << " Usage: tpk2d [options] Nx,Ny,Nz OutPk_PPFFile \n" 
	 << " options: [-d dx,dy,dz] [-k kmin,kmax,nbin] [-t lowval] [-r rfac] [-p lev] \n"
	 << "   -d dx,dy,dz : define input 3D map cell size (in 100-Mpc) \n"
	 << "   -k nbin,kmin,kmax : define number of bins and k-range for P(k) computation \n"
	 << "   -t lowval : set input map cells with val<lowval to lowval before computing P(k) \n"
      	 << "   -s sigmaz : gaussian smearing along z (in 100-Mpc units) \n"
	 << "   -r rfac : P(k) renormalisation factor \n"
	 << "   -p lev : define print level (0,1,2..) \n" << endl;
    return 1;
  }
  Timer tm("tpk2d");
  int rc = 0;
  try { 
    string outppfname;
    // Taille du cube 
    long int Nx,Ny,Nz; 
    bool fgsetcellsize=false;
    double dx,dy,dz;
    dx=dy=dz=0.05; 
    int nkbin=100;
    double kmin=0.001, kmax=100.;
    bool fgtrunclowval=false;
    double lowval=-9.e19;
    double sigmaz=0.;
    bool fgrfac=false;
    double rfac=1.;
    int prtlev=0;
    //----------------------------------------------------------
    // decodage arguments optionnel 
    bool fgoptarg=true;
    while (fgoptarg&&(narg>1)) {
      string fbo = arg[1];
      if (fbo=="-d")  {   // specification taille de cellules (Mpc) drho/rho ou ngal 
	if (narg<3) { cout << " tpk2d/missing/bad argument, tpk2d -h for help " << endl;  return 2; }
	sscanf(arg[2],"%lf,%lf,%lf",&dx,&dy,&dz); 
	fgsetcellsize=true; arg+=2; narg-=2; 
      }
      else if (fbo=="-k")  {   // specification taille de cellules (Mpc) drho/rho ou ngal 
	if (narg<3) { cout << " tpk2d/missing/bad argument, tpk2d -h for help " << endl;  return 2; }
	sscanf(arg[2],"%d,%lf,%lf",&nkbin,&kmin,&kmax);  arg+=2; narg-=2; 
      }
      else if (fbo=="-t")  { 
	if (narg<3) { cout << " tpk2d/missing/bad argument, tpk2d -h for help " << endl;  return 2; }
	lowval=atof(arg[2]);  fgtrunclowval=true; arg+=2; narg-=2; 
      }
      else if (fbo=="-s")  { 
	if (narg<3) { cout << " tpk2d/missing/bad argument, tpk2d -h for help " << endl;  return 2; }
	sigmaz=atof(arg[2]);  arg+=2; narg-=2; 
      }
      else if (fbo=="-r")  { 
	if (narg<3) { cout << " tpk2d/missing/bad argument, tpk2d -h for help " << endl;  return 2; }
	rfac=atof(arg[2]);  fgrfac=true; arg+=2; narg-=2; 
      }
      else if (fbo=="-p")  { 
	if (narg<3) { cout << " tpk2d/missing/bad argument, tpk2d -h for help " << endl;  return 2; }
	prtlev=atoi(arg[2]);  fgrfac=true; arg+=2; narg-=2; 
      }
      else fgoptarg=false;
    }
    //----------------------------------------------------------
    if (narg<3) { cout << " tpk2d/missing/bad argument, tpk2d -h for help " << endl;  return 2; }
    sscanf(arg[1],"%ld,%ld,%ld",&Nx,&Ny,&Nz);
    outppfname=arg[2];

    cout << "tpk2d[1.a] : Opening output PPF file " << outppfname <<endl;
    POutPersist pof(outppfname);
    
    cout << "tpk2d[1.b] : Creating 3D map Nx Ny Nz= " << Nx<<" x "<<Ny<<" x "<<Nz<<" and GFour3DPk"<<endl;
    TArray<r_4> ingrid(Nx,Ny,Nz);
    /*
    {
      FitsInOutFile fin(infitsname, FitsInOutFile::Fits_RO);
      fin >> ingrid;
    }
    double mean, sigma;
    MeanSigma(ingrid, mean, sigma);
    cout << "tpk2d[1.b] Input grid sizes " << ingrid.InfoString() << endl;
    ingrid.Show(); 
    cout << "... Input grid Mean=" << mean << " Sigma=" << sigma << endl;
    tm.Split(" After read ");
    */
    
    GFour3DPk  gpkc(ingrid);
    gpkc.SetPrtLevel(prtlev);
    cout << "tpk2d[1.c] : setting grid cell size: "<< dx<<" x "<<dy<<" x "<<dz<<endl;    
    gpkc.SetGridCellSize(dx,dy,dz);

    /*
    if (fgtrunclowval) {  // truncating grid value below threshold 
      cout << "tpk2d[2] : calling CleanNegatives() ... " << endl;    
      gpkc.CleanNegatives(lowval);
      MeanSigma(ingrid, mean, sigma);
      cout << "... After CleanNegatives grid Mean=" << mean << " Sigma=" << sigma << endl;
      tm.Split(" After CleanNegatives ");
    }
    */
    cout << "tpk2d[2] : Input power spectrum ... "<<endl;    
    Xsi1 xsi;
    Vector vxsi=xsi.FillVec();
    //DBG     cout<<" *DBG* done Vector vxsi=xsi.FillVec();" << endl;
    PkFrXsi pk(xsi);
    Vector vpk=pk.FillVec();
    //DBG    cout<<" *DBG* done PkFrXsi pk(xsi).FillVec();" << endl;
    InterpPk interpkA(pk,0.,200.,2500,false);
    XsiFrPk xsipk(interpkA);
    Vector vxsipk=xsipk.FillVec();
    //DBG     cout<<" *DBG* done XsiFrPk xsipk(interpkA).FillVec();" << endl;
    InterpPk interpk(pk);
    Vector vpkinterp=interpk.FillVec();
    //DBG    cout<<" *DBG* done InterpPk interpk(pk).FillVec();" << endl;
    XsiFrPk xsipki(interpk);
    Vector vxsipki=xsipki.FillVec();
    //DBG    cout<<" *DBG* done XsiFrPk xsipki(interpk).FillVec();" << endl;
    Histo2D hxsi2d=xsipk.FillHisto();
    //DBG    cout<<" *DBG* done XsiFrPk xsipki(interpk).FillHisto();" << endl;
    XsiFrPk2D xsi2d(interpk);
    Histo2D hxsi2dB=xsi2d.FillHisto();

    pof<<PPFNameTag("vxsi")<<vxsi;
    pof<<PPFNameTag("vpk")<<vpk;
    pof<<PPFNameTag("vpkinterp")<<vpkinterp;
    pof<<PPFNameTag("vxsipk")<<vxsipk;
    pof<<PPFNameTag("vxsipki")<<vxsipki;
    pof<<PPFNameTag("hxsi2d")<<hxsi2d;
    pof<<PPFNameTag("hxsi2dB")<<hxsi2dB;

    cout << "tpk2d[2.b] : saved vxsi vpk vpkinterp vxsipk vxsipki hxsi2d hxsi2dB to output PPF  ... "<<endl;    
    tm.Split(" After Input power spectrum ");

    cout << "tpk2d[3] : Generating Fourier coefficients from power spectrum ... " << endl;
    if (sigmaz>0.) cout << " Smearing along the z-axis with sigma="<<sigmaz<<endl;
    gpkc.generateFourierAmp(interpk, sigmaz);
    tm.Split(" After generateFourierAmp ");
    
    cout << "tpk2d[4.a] : Computing 1D power spectrum ... " << endl;
    HProf recPk=gpkc.ComputePk(200, 0., 75., true);

    cout << "tpk2d[4.b] : Computing 2D power spectrum ... " << endl;
    Histo2D recPk2D=gpkc.ComputePk2D(100,50.);
    tm.Split(" After P(k)-reconstruction ");

    pof<<PPFNameTag("recPk")<<recPk;
    pof<<PPFNameTag("recPk2D")<<recPk2D;
    cout << "tpk2d[4.c] : saved  recPk & recPk2D to output PPF  ... "<<endl;    

    cout << "tpk2d[5] : Computing 2D correlation function from 2D-P(k) ... " << endl;
    //    H2DSpline h2spl(recPk2D);
    H2DLinInterpol h2lint(recPk2D, true);
    XsiFrPk2D xsi2drec(h2lint);
    Histo2D recXsi2D=xsi2drec.FillHisto();
    tm.Split(" After Xsi2D-reconstruction ");
    pof<<PPFNameTag("recXsi2D")<<recXsi2D;
    cout << "tpk2d[5.b] : saved  recXsi2D to output PPF  ... "<<endl;    
    
    /*
    cout << "tpk2d[4] : computing power spectrum ... " << endl;
    HProf hpk = gpkc.ComputePk(nkbin,kmin,kmax,true);
    DataTable dtpk; 
    Histo hrap = gpkc.FillPkDataTable(dtpk, rfac);
    tm.Split(" After ComputePk ");
    if (prtlev>0)  { 
      dtpk.SetShowMinMaxFlag(true);
      cout << dtpk; 
    }
    {   // Saving computed P(k) to text (ascii) file 
      cout << "calcpk[5] : wrting P(k) to text file " << outtextname << endl;
      ofstream tof(outtextname.c_str());
      dtpk.WriteASCII(tof);
    }
    if (outppfname.length()>0)  {   // Saving computed P(k) to PPF file 
      cout << "calcpk[6] : writing profile histo P(k) (hpk) and DataTable P(k) dtpk to PPF file  " << outppfname << endl;
      POutPersist po(outppfname);
      po<<PPFNameTag("hpk")<<hpk;
      po<<PPFNameTag("dtpk")<<dtpk;
    }
    */
  }  // End of try bloc 
  catch (PThrowable & exc) {  // catching SOPHYA exceptions
    cerr << " tpk2d.cc: Catched Exception (PThrowable)" << (string)typeid(exc).name() 
         << "\n...exc.Msg= " << exc.Msg() << endl;
    rc = 99;
  }
  catch (std::exception & e) {  // catching standard C++ exceptions
    cerr << " tpk2d.cc: Catched std::exception "  << " - what()= " << e.what() << endl;
    rc = 98;
  }
  catch (...) {  // catching other exceptions
    cerr << " tpk2d.cc: some other exception (...) was caught ! " << endl;
    rc = 97;
  }
  cout << " ==== End of tpk2d.cc program  Rc= " << rc << endl;
  return rc;    
}
