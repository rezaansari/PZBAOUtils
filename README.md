PZBAOUtils
==========

small utility programs for BAO/PhotoZ catalog analysis -
Creation du module en Fevrier 2015
MAJ en septembre 2015, Février 2016


-------------------------------------------------------------------------
tpk2d.cc : small utility program to compute correlation function from power spectrum,
generate grid according to a power spectrum, compute back P(k), 2D power spectrum P(k)-2D
and 1D and 2D correlation function

 Usage: tpk2d [options] InputPkFile Out_PPFFile 
   InputPkFile : text file with pairs of values k  Pk on each line 
       Specify  InputPkFile = '.' or '-' for default P(k) defined in corfunc.h 
 options: [-N Nx,Ny,Nz] [-d dx,dy,dz] [-k kmin,kmax,nbin] [-k2d nbin2d,kmax2d] 
          [-s sigmaz] [-t lowval] [-r rfac] [-p lev] 
   -N Nx,Ny,Nz : define input 3D grid size default=400x400x400 
   -d dx,dy,dz : define input 3D grid cell size (in Mpc) default=5x5x5 Mpc
   -k nbin,kmin,kmax : define number of bins and k-range for P(k) computation (def=200,0.,1.) 
   -k2d nbin2d,kmax2d : define number of bins and k-max for 2D-P(k) computation (def=100,1.) 
   -s sigmaz : gaussian smearing along z (in Mpc units) default=0.
   -p lev : define print level (0,1,2..) 
   -t lowval : set input map cells with val<lowval to lowval before computing P(k) **NOT USED** 
   -r rfac : P(k) renormalisation factor **NOT USED**  

##  Compute reconstructed P(k) and xsi(r) starting from Xsi1 test auto-correlation function
# with 7.5 Mpc smearing along z
./Objs/tpk2d -s 7.5 - totos.ppf

## Check and display some of the results using piapp
spiapp -term

setaxesatt 'font=helvetica,bold,20 fixedfontsize minorticks'
openppf totos.ppf
newwin 2 1 1000 500
disp vpkinterp "arra_xlim=0.,${vpkinterp.info.kmax}  red"
disp recPk 'same black'
disp vxsipki "arra_xlim=0.,${vxsipki.info.rmax}  red"
disp recXsi  'same black'

newwin 2 1 1000 500
disp recPk2D 'xylimits=-0.4,0.4,-0.4,0.4 h2disp=img colbr128'
disp recXsi2D  'xylimits=-200,200,-200,200 h2disp=img colbr128'



-------------------------------------------------------------------------
galcatext.cc : to extract rows from a fits bin table (Galaxy catalog) and write 
 the subset to a new fits file
csh> Objs/galcatext -h
 ---- galcatext: fits catalog extraction program ---- 
 Usage: galcatext InFitsFile OutFitsFile Range [HDU=2] [SegSize=8192] 
    Range: start,end,step  (starting from zero) 

##  To extract the first 10000 rows (galaxies)
csh> rm smallcat.fits  
csh> Objs/galcatext bigcat.fits smallcat.fits 0,10000,1 
##  To extract one every 50 rows (galaxies) - 2%
csh> rm smallcat.fits  
csh> Objs/galcatext bigcat.fits smallcat.fits 0,0,50

 
-------------------------------------------------------------------------
grid2pk.cc : compute power spectrum (P(k)) from grids of de rho/ro or ngals

 Objs/grid2pk -h
 SophyaInitiator::SophyaInitiator() BaseTools Init
 PIOPersist::Initialize() Starting Sophya Persistence management service 
 Usage: grid2pk [options] In3DMapFitsName OutPkTextFile [OutPk_PPFFile] 
 options: [-d dx,dy,dz] [-k kmin,kmax,nbin] [-t lowval] [-r rfac] [-p lev] 
   -d dx,dy,dz : define input 3D map cell size (in Mpc) 
   -k nbin,kmin,kmax : define number of bins and k-range for P(k) computation 
   -t lowval : set input map cells with val<lowval to lowval before computing P(k) 
   -r rfac : P(k) renormalisation factor 
   -p lev : define print level (0,1,2..) 



-------------------------------------------------------------------------
Files :
gpkspec.h  gpkspec.cc : class **GFour3DPk**   grid to P(k) and P(k) to grid computations
corfunc.h : correlation function to P(k) and reverse - 1D and 2D computation
myinteg2d.h myinteg2d.cc : Integrator for P(k)-2D <> xsi-2D(r)
hsplfit.h : utility class to interpolate histograms to be represented as functions
Pk_z_1_0.txt : set of points (k,P(k) representing cosmological power spectrum
