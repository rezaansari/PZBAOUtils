/*  ------------------------ Projet BAO/PhotoZ/LSST -------------------- 
  Programme de calcul du spectre de puissance (3D) a partir d'un 
  cube de donnees (delta rho/rho ou NGal )
    R. Ansari (pour Adline Choyer) - Feb 2015 
  Usage: grid2pkgrid2pk [options] In3DMap_FitsName OutPk_TextFile [OutPk_PPFFile]
         Options: [-d dx,dy,dz] [-k kmin,kmax,nbin] [-t lowval] [-r rfac] [-p lev]
---------------------------------------------------------------  */

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
#include "fitsioserver.h"
#include "fiosinit.h"

#include "ctimer.h"

//---- this modules include files 
#include "gpkspec.h"

//-------------------------------------------------------------------------
//      ------------------ MAIN PROGRAM ------------------------------
//-------------------------------------------------------------------------
int main(int narg, const char* arg[])
{
  if ((narg<3)||((narg>1)&&(strcmp(arg[1],"-h")==0))) {
    cout << " Usage: grid2pk [options] In3DMap_FitsName OutPk_TextFile [OutPk_PPFFile] \n" 
	 << " options: [-d dx,dy,dz] [-k kmin,kmax,nbin] [-t lowval] [-r rfac] [-p lev] \n"
	 << "   -d dx,dy,dz : define input 3D map cell size (in Mpc) \n"
	 << "   -k nbin,kmin,kmax : define number of bins and k-range for P(k) computation \n"
	 << "   -t lowval : set input map cells with val<lowval to lowval before computing P(k) \n"
	 << "   -r rfac : P(k) renormalisation factor \n"
	 << "   -p lev : define print level (0,1,2..) \n" << endl;
    return 1;
  }
  Timer tm("grid2pk");
  int rc = 0;
  try { 
    string infitsname;
    string outtextname;
    string outppfname;

    bool fgsetcellsize=false;
    double dx=1.,dy=1.,dz=1.;
    int nkbin=100;
    double kmin=0.001, kmax=1.001;
    bool fgtrunclowval=false;
    double lowval=-9.e19;
    bool fgrfac=false;
    double rfac=1.;
    int prtlev=0;
    //----------------------------------------------------------
    // decodage arguments optionnel 
    bool fgoptarg=true;
    while (fgoptarg&&(narg>1)) {
      string fbo = arg[1];
      if (fbo=="-d")  {   // specification taille de cellules (Mpc) drho/rho ou ngal 
	if (narg<3) { cout << " grid2pk/missing/bad argument, grid2pk -h for help " << endl;  return 2; }
	sscanf(arg[2],"%lf,%lf,%lf",&dx,&dy,&dz); 
	fgsetcellsize=true; arg+=2; narg-=2; 
      }
      else if (fbo=="-k")  {   // specification taille de cellules (Mpc) drho/rho ou ngal 
	if (narg<3) { cout << " grid2pk/missing/bad argument, grid2pk -h for help " << endl;  return 2; }
	sscanf(arg[2],"%d,%lf,%lf",&nkbin,&kmin,&kmax);  arg+=2; narg-=2; 
      }
      else if (fbo=="-t")  { 
	if (narg<3) { cout << " grid2pk/missing/bad argument, grid2pk -h for help " << endl;  return 2; }
	lowval=atof(arg[2]);  fgtrunclowval=true; arg+=2; narg-=2; 
      }
      else if (fbo=="-r")  { 
	if (narg<3) { cout << " grid2pk/missing/bad argument, grid2pk -h for help " << endl;  return 2; }
	rfac=atof(arg[2]);  fgrfac=true; arg+=2; narg-=2; 
      }
      else if (fbo=="-p")  { 
	if (narg<3) { cout << " grid2pk/missing/bad argument, grid2pk -h for help " << endl;  return 2; }
	prtlev=atoi(arg[2]);  fgrfac=true; arg+=2; narg-=2; 
      }
      else fgoptarg=false;
    }
    //----------------------------------------------------------
    if (narg<3) { cout << " grid2pk/missing/bad argument, grid2pk -h for help " << endl;  return 2; }
    infitsname=arg[1];
    outtextname=arg[2];
    if (narg>3) outppfname=arg[3];

    cout << "grid2pk[1] : reading 3D map from file " << infitsname << endl;
    TArray<r_4> ingrid;
    {
      FitsInOutFile fin(infitsname, FitsInOutFile::Fits_RO);
      fin >> ingrid;
    }
    double mean, sigma;
    MeanSigma(ingrid, mean, sigma);
    cout << "grid2pk[1.b] Input grid sizes " << ingrid.InfoString() << endl;
    ingrid.Show(); 
    cout << "... Input grid Mean=" << mean << " Sigma=" << sigma << endl;
    tm.Split(" After read ");
    
    GFour3DPk  gpkc(ingrid);
    gpkc.SetPrtLevel(prtlev);
    if (fgsetcellsize) {
      gpkc.SetGridCellSize(dx,dy,dz);
    }
    if (fgtrunclowval) {  // truncating grid value below threshold 
      cout << "grid2pk[2] : calling CleanNegatives() ... " << endl;    
      gpkc.CleanNegatives(lowval);
      MeanSigma(ingrid, mean, sigma);
      cout << "... After CleanNegatives grid Mean=" << mean << " Sigma=" << sigma << endl;
      tm.Split(" After CleanNegatives ");
    }

    cout << "grid2pk[3] : Computing Fourier coefficients ... " << endl;    
    gpkc.doFFT();
    tm.Split(" After doFFT ");

    cout << "calcpk[4] : computing power spectrum ... " << endl;
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
  }  // End of try bloc 
  catch (PThrowable & exc) {  // catching SOPHYA exceptions
    cerr << " grid2pk.cc: Catched Exception (PThrowable)" << (string)typeid(exc).name() 
         << "\n...exc.Msg= " << exc.Msg() << endl;
    rc = 99;
  }
  catch (std::exception & e) {  // catching standard C++ exceptions
    cerr << " grid2pk.cc: Catched std::exception "  << " - what()= " << e.what() << endl;
    rc = 98;
  }
  catch (...) {  // catching other exceptions
    cerr << " grid2pk.cc: some other exception (...) was caught ! " << endl;
    rc = 97;
  }
  cout << " ==== End of grid2pk.cc program  Rc= " << rc << endl;
  return rc;    
}
