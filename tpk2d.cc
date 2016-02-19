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
    cout << " Usage: tpk2d [options] InputPkFile Out_PPFFile \n"
	 << "   InputPkFile : text file with pairs of values k  Pk on each line \n"
	 << "       Specify  InputPkFile = '.' or '-' for default P(k) defined in corfunc.h \n"
	 << " options: [-N Nx,Ny,Nz] [-d dx,dy,dz] [-k kmin,kmax,nbin] [-k2d nbin2d,kmax2d] \n"
	 << "          [-s sigmaz] [-t lowval] [-r rfac] [-p lev] \n"
      	 << "   -N Nx,Ny,Nz : define input 3D grid size default=400x400x400 \n"
	 << "   -d dx,dy,dz : define input 3D grid cell size (in Mpc) default=5x5x5 Mpc\n"
	 << "   -k nbin,kmin,kmax : define number of bins and k-range for P(k) computation (def=200,0.,1.) \n"
      	 << "   -k2d nbin2d,kmax2d : define number of bins and k-max for 2D-P(k) computation (def=100,1.) \n"
	 << "   -s sigmaz : gaussian smearing along z (in Mpc units) default=0.\n"
      	 << "   -p lev : define print level (0,1,2..) \n" 
	 << "   -t lowval : set input map cells with val<lowval to lowval before computing P(k) **NOT USED** \n"
	 << "   -r rfac : P(k) renormalisation factor **NOT USED**  \n"
	 << endl;
    return 1;
  }
  Timer tm("tpk2d");
  int rc = 0;
  try { 
    string inpkname=".";
    string outppfname="toto.ppf";
	
    // Taille du cube 
    long int Nx,Ny,Nz;
    Nx=Ny=Nz=400;
    bool fgsetcellsize=false;
    double dx,dy,dz;
    dx=dy=dz=5; 
    int nkbin=200;
    double kmin=0., kmax=1.;
    int nkbin2d=100;
    double kmax2d=1.;
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
      if (fbo=="-N")  {   // specification taille de cellules (Mpc) drho/rho ou ngal 
	if (narg<3) { cout << " tpk2d/missing/bad argument, tpk2d -h for help " << endl;  return 2; }
	sscanf(arg[1],"%ld,%ld,%ld",&Nx,&Ny,&Nz);
	arg+=2; narg-=2; 
      }
      else if (fbo=="-d")  {   // specification taille de cellules (Mpc) drho/rho ou ngal 
	if (narg<3) { cout << " tpk2d/missing/bad argument, tpk2d -h for help " << endl;  return 2; }
	sscanf(arg[2],"%lf,%lf,%lf",&dx,&dy,&dz); 
	fgsetcellsize=true; arg+=2; narg-=2; 
      }
      else if (fbo=="-k")  {   // specification nb-bin kmin,kmax pour calcul P(k)  
	if (narg<3) { cout << " tpk2d/missing/bad argument, tpk2d -h for help " << endl;  return 2; }
	sscanf(arg[2],"%d,%lf,%lf",&nkbin,&kmin,&kmax);  arg+=2; narg-=2; 
      }
      else if (fbo=="-k2d")  {   // specification nb-bin kmax pour calcul 2D-P(k)  
	if (narg<3) { cout << " tpk2d/missing/bad argument, tpk2d -h for help " << endl;  return 2; }
	sscanf(arg[2],"%d,%lf",&nkbin2d,&kmax2d);  arg+=2; narg-=2; 
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

    inpkname=arg[1];
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
    InterpPk* interpk_p=NULL;
    if ((inpkname != ".")&&(inpkname != "-")) {
      cout << "tpk2d[2.a] : reading input power spectrum from file "<<inpkname<<endl;
      interpk_p= new InterpPk(inpkname);
    }
    else {
      cout << "tpk2d[2.a] : creating default input P(k) from Xsi1 (corfunc.h) "<<endl;
      Xsi1 xsi;
      Vector vxsi=xsi.FillVec();
      PkFrXsi pk(xsi);
      Vector vpk=pk.FillVec();
      interpk_p= new InterpPk(pk);
      pof<<PPFNameTag("vxsi")<<vxsi;
      pof<<PPFNameTag("vpk")<<vpk;
      cout << " ... saved vxsi vpk  to output PPF  ... "<<endl;
    }
    cout << "tpk2d[2.b] : computing xsi from interpolated P(k)  "<<endl;
    InterpPk & interpk = (*interpk_p);
    Vector vpkinterp=interpk.FillVec();
    //DBG    cout<<" *DBG* done InterpPk interpk.FillVec()" << endl;
    XsiFrPk xsipki(interpk);
    Vector vxsipki=xsipki.FillVec();
    //DBG    cout<<" *DBG* done XsiFrPk xsipki(interpk).FillVec();" << endl;
    Histo2D hxsi2d=xsipki.FillHisto();
    //DBG    cout<<" *DBG* done XsiFrPk xsipki(interpk).FillHisto();" << endl;
    XsiFrPk2D xsi2d(interpk,5.,150);
    Histo2D hxsi2dB=xsi2d.FillHisto();

    pof<<PPFNameTag("vpkinterp")<<vpkinterp;
    pof<<PPFNameTag("vxsipki")<<vxsipki;
    pof<<PPFNameTag("hxsi2d")<<hxsi2d;
    pof<<PPFNameTag("hxsi2dB")<<hxsi2dB;

    cout << "tpk2d[2.c] : saved vxsi vpk vpkinterp vxsipk vxsipki hxsi2d hxsi2dB to output PPF  ... "<<endl;    
    tm.Split(" After Input power spectrum ");

    cout << "tpk2d[3] : Generating Fourier coefficients from power spectrum ... " << endl;
    if (sigmaz>0.) cout << " Smearing along the z-axis with sigma="<<sigmaz<<endl;
    gpkc.generateFourierAmp(interpk, sigmaz);
    tm.Split(" After generateFourierAmp ");
    cout << "tpk2d[4.a] : Computing 1D power spectrum ... nbin="<<nkbin<<" kmin="<<kmin<<" kmax="<<kmax<<endl;
    HProf recPk=gpkc.ComputePk(nkbin, kmin, kmax, true);

    cout << "tpk2d[4.b] : Computing 2D power spectrum ... nbin="<<nkbin2d<<" kmax="<<kmax2d<<endl; 
    Histo2D recPk2D=gpkc.ComputePk2D(nkbin2d,kmax2d);
    tm.Split(" After P(k)-reconstruction ");

    pof<<PPFNameTag("recPk")<<recPk;
    pof<<PPFNameTag("recPk2D")<<recPk2D;
    cout << "tpk2d[4.c] : saved  recPk & recPk2D to output PPF  ... "<<endl;    

    cout << "tpk2d[5.c] : Computing 1D correlation function from recP(k) ... " << endl;    
    H1DLinInterpol hlint(recPk);
    XsiFrPk xsirec(hlint);
    Vector recXsi=xsirec.FillVec();
    
    cout << "tpk2d[5.b] : Computing 2D correlation function from 2D-recP(k) ... " << endl;
    //    H2DSpline h2spl(recPk2D);
    H2DLinInterpol h2lint(recPk2D, true);
    XsiFrPk2D xsi2drec(h2lint);
    Histo2D recXsi2D=xsi2drec.FillHisto();
    tm.Split(" After Xsi/Xsi2D-reconstruction ");

    pof<<PPFNameTag("recXsi")<<recXsi;
    pof<<PPFNameTag("recXsi2D")<<recXsi2D;
    cout << "tpk2d[5.c] : saved  recXsi & recXsi2D to output PPF  ... "<<endl;    
    
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
    if (interpk_p) delete interpk_p;
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
