/*--- programme de retraitement des On,Off Amas  ---
 ----------------------------------------------------------*/

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <iostream>
#include <cmath>
#include <complex>

#include "array.h"
#include "ntuple.h"

#include "ctimer.h"
// #include "tarrinit.h"
// #include "histinit.h"

#include "swfitsdtable.h"
#include "fiosinit.h"

using namespace std; 
using namespace SOPHYA; 

int Usage()
{
  cout << " ---- galcatext: fits catalog extraction program ---- \n"
       << " Usage: galcatext InFitsFile OutFitsFile Range [HDU=2] [SegSize=8192] \n"
       << "    Range: start,end,step  (starting from zero) " << endl;
  return 0;
}
/* -- MAIN -- */ 
int main (int narg, char* arg[])
{
  if ((narg < 4)||(strcmp(arg[1],"-h")==0)) {
    return Usage();
  }
  cout << " ------ galcatext : fits catalog extraction program ------- " << endl;

  int rc = 0;
  try {
    SophyaInit();
    Timer tm("galcatext.cc");
    string infile = arg[1];
    string outfile = arg[2];
    long int kstart=0,kend=0,kstep=1;
    sscanf(arg[3],"%ld,%ld,%ld",&kstart,&kend,&kstep);
    int numhdu=2;
    size_t segsz=8192;
    if (narg>4)  numhdu=atoi(arg[4]);
    if (narg>5)  segsz=atol(arg[5]);
    cout << "[1] galcatext: opening fits file: " << infile << " , reading SwFitsDataTable ..."<<endl;
    FitsInOutFile si(infile, FitsInOutFile::Fits_RO);
    // Position the fits file on the first extension (BinTable)
    si.MoveAbsToHDU(numhdu);  
    SwFitsDataTable dt(si, segsz, false);
    // Printing table info 
    cout << dt ;
    if (kstart<0) kstart=0;
    if ((kend<=kstart)||(kend>dt.NRows()))  kend=dt.NRows();
    cout << "[2] extracting rows k, "<<kstart<<" <=k< "<<kend<< " step(k+=):"<<kstep
	 << " Out of NRows="<<dt.NRows()<< ((kend-kstart)*100/kstep)/dt.NRows()<<endl;
    FitsInOutFile so(outfile, FitsInOutFile::Fits_Create);
    vector<size_t> rows;
    for(size_t k=kstart; k<kend; k+=kstep)  rows.push_back(k);
    
    SwFitsDataTable dte(so, segsz/8, true);
    dte.SelectFrom(dt,rows,ProgBarM_Time);
    dte.Info()=dt.Info();
    dte.Info()["GEXTFILE"]=infile;
    dte.Info()["KFIRST"]=(int_8)kstart;
    dte.Info()["KEND"]=(int_8)kend;
    dte.Info()["KSTEP"]=(int_8)kstep;
  }
  catch (PThrowable & exc) {
    cerr << " galcatext.cc: Catched Exception (PThrowable)" << (string)typeid(exc).name() 
         << " - Msg= " << exc.Msg() << endl;
    rc = 99;
  }
  catch (std::exception & e) {
    cerr << " galcatext.cc: Catched std::exception "  
         << " - what()= " << e.what() << endl;
    rc = 98;
  }
  catch (...) {
    cerr << " galcatext.cc: some other exception (...) was caught ! " << endl;
    rc = 97;
  }
  cout << " ------------ END OF galcatext  (Rc= " << rc << ") ------------ " << endl;
  return rc;
}
