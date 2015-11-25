#include <cstdio>
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <complex>
#include <vector>
#include "fftw3.h"
#include "dft_pwpp.h"

void parameter_read(parmtr &param);
void pw_find       (numset &ns, parmtr param, vec2d &pw);
void ck_initialize (numset ns, parmtr param, vec2d &pw, fftw_ *ck, fftw_ *nr, 
		    fftw_ *nk);
void G2_initialize (numset &ns, parmtr param, double *G2);
void ck_fft_cr     (numset ns, parmtr param, fftw_ *ck, fftw_ *cr);
void nr_fft_nk     (numset &ns, parmtr &param,fftw_ *cr, fftw_ *nr, fftw_ *nk,
		    double *Fnew, double *Fold, double **G, double *rho_old);
void nr_fft_nk3     (numset &ns, parmtr &param,fftw_ *cr, fftw_ *nr, fftw_ *nk,
		    double *Fnew, double *Fold, double **G, double *rho_old);
void nr_fft_nk2    (numset &ns, parmtr &param, fftw_ *cr, fftw_ *nr, fftw_ *nk);
void V_HT          (numset ns, parmtr param, double *G2, fftw_ *nk, fftw_ *VH, 
		    status &sys, vec2d pw);
void V_loc         (numset ns, parmtr param, double *G2, fftw_ *Vloc);
void E_loc         (numset ns, parmtr param, double *G2, fftw_ *nk, 
		    fftw_ *Vloc, status &sys);
void V_nl          (numset ns, parmtr param, double *G2, fftw_ **Hnl, 
		    vec2d pw, fftw_ **bkf);
void V_nl2          (numset ns, parmtr param, double *G2, fftw_ **Hnl, 
		    vec2d pw, fftw_ **bkf);
void V_HGH         (numset ns, parmtr param, double *G2, fftw_ **Hnl, 
		    vec2d pw, fftw_ **bkf);
void E_nl          (numset ns, parmtr param, vec2d pw, fftw_ **bkf, fftw_ *ck, 
		    status &sys);
void V_xc          (numset ns, parmtr param, fftw_ *nr, fftw_ *Vxc, 
		    status &sys);
void H_kin         (numset ns, vec2d pw, double *G2, double *Hkin);
void E_kin         (numset ns, parmtr param, vec2d pw, fftw_ *ck, 
		    double *Hkin, status &sys);
void Hermitian_lapack(numset &ns, parmtr param, vec2d pw, fftw_ *VH,fftw_ *Vxc,
		      fftw_ *Vloc, double *Hkin, fftw_ **Hnl, fftw_ *ck, 
		      fftw_ **Hold);

void print_rho(numset ns, parmtr param, fftw_complex *nr)
{
  int i, j, k, id, ii, jj;
  FILE *fp = fopen("rho.dat","w");
  k = 0;// ns.ngrid/2;
  for (i=0;i<ns.ngrid+1;i++) {
    ii = i%ns.ngrid;
    for (j=0;j<ns.ngrid+1;j++) {
      jj = j%ns.ngrid;
      id = k + ns.ngrid*(ii + ns.ngrid*jj);
      fprintf (fp, "%d %d %f\n", j, i, nr[id][0]);
    }    
    fprintf(fp,"\n");
  }
  int fclose(FILE *fp);
  //std::cout << nr[129][0];
}

void check(numset ns, parmtr param, fftw_complex *ck)
{
  int i, id, jd;
  double real = 0.0, imag = 0.0, dv = param.omega/(double) ns.ngrid3;
  for (i=0;i<ns.ngrid3;i++) {
    id = i + 4*ns.ngrid3; jd = i + 4*ns.ngrid3;
    real += ck[id][0]*ck[jd][0] + ck[id][1]*ck[jd][1];
    imag += ck[id][0]*ck[jd][1] - ck[id][1]*ck[jd][0];
  }
  std::cout << "normality check = " << real << " " << imag << std::endl;
}

int radix_decompose(int n)
{
  int n2 = 1, n3 = 0, imax, dump; 
  bool tag = true;
  do { 
    n2 ++;
    if (n < pow(2,n2)) tag = false; 
  } while(tag); 
  tag = true;
  do { 
    n3 ++;
    if (n < 2*pow(3,n3)) tag = false; 
  } while(tag);
  imax = pow(2,n2);
  for (int i=1;i<n2;i++) {
    for (int j=0;j<=n3;j++) {
      dump = pow(2,i)*pow(3,j);
      if (dump >= n) imax = std::min(imax, dump);
    }
  }
  int final = imax;
  return final;
}

int main()
{
  //
  // variable declaration
  int     i, j, n;
  numset    ns;
  parmtr param;
  status   sys;  
  vec2d pw(4, std::vector <int> (0,0));
  double E_old, E_total;
  time_t t0 = time(NULL);
  srand((int) t0);
  //
  // Initial data parsing and setup
  parameter_read(param);
  pw_find (ns, param, pw); ns.nband = NBAND;
  double *G2, *Hkin, *rho_old, *Fnew, *Fold, **G;
  G2   = new double [ns.ngrid3];
  Hkin = new double [ns.npw];
  fftw_complex *ck, *cr, *nr, *nk, *VH, *Vloc, *Vxc, **Hnl, **bkf, **Hold;
  ck   = new fftw_complex [ns.ngrid3*ns.nband];
  cr   = new fftw_complex [ns.ngrid3*ns.nband];
  nr   = new fftw_complex [ns.ngrid3];
  nk   = new fftw_complex [ns.ngrid3];
  VH   = new fftw_complex [ns.ngrid3];
  Vloc = new fftw_complex [ns.ngrid3];
  Vxc  = new fftw_complex [ns.ngrid3];
  Fnew = new double [ns.ngrid3];
  Fold = new double [ns.ngrid3];
  rho_old=new double[ns.ngrid3];
  G    = new double * [ns.ngrid3];
  for (i=0;i<ns.ngrid3;i++) G[i] = new double [ns.ngrid3];
  Hnl  = new fftw_ * [ns.npw]; Hold  = new fftw_ * [ns.npw]; 
  for (i = 0;i<ns.npw;i++) {
    Hnl[i] = new fftw_ [ns.npw];
    Hold[i]= new fftw_ [ns.npw];
  }
  bkf  = new fftw_complex * [Nang];
  for (i = 0;i<Nang;i++) bkf[i] = new fftw_ [ns.npw];
  ck_initialize(ns, param, pw, ck, nr, nk); param.alpha = 1.0;
  G2_initialize(ns, param, G2);
  V_loc(ns, param, G2, Vloc);
  H_kin(ns, pw, G2, Hkin);
  V_nl(ns, param, G2, Hnl, pw, bkf);
  n = 0;
  ns.nstage = 0;
  E_total =0.0;
  for (i=0;i<ns.npw;i++) {
    for (j=0;j<ns.npw;j++) {
      Hold[i][j][0] = Hold[i][j][1] = 0.0;
    }
  }
  do {
    E_old = E_total;
    ck_fft_cr(ns, param, ck, cr);
    nr_fft_nk(ns, param, cr, nr, nk, Fnew, Fold, G, rho_old);
    //nr_fft_nk2(ns, param, cr, nr, nk);    
    V_HT     (ns, param, G2, nk, VH, sys, pw);
    V_xc     (ns, param,     nr, Vxc,  sys);
    E_loc    (ns, param, G2, nk, Vloc, sys);
    E_nl     (ns, param, pw, bkf, ck,  sys);
    E_kin    (ns, param, pw, ck, Hkin, sys);
    Hermitian_lapack(ns, param, pw, VH, Vxc, Vloc, Hkin, Hnl, ck, Hold);
    E_total = sys.E_kin + sys.E_HT + sys.E_xc + sys.E_loc + sys.E_nl;
    //std::cout << "Total E = " << E_total << " with dE = " << E_total - E_old 
    //<< " at " << n <<  std::endl;
    //std::cout << sys.E_kin << " " << sys.E_HT << " " << sys.E_xc << " " << sys.E_nl << std::endl;
    //std::cout << "Total E = " <<  sys.E_kin + sys.E_HT +  sys.E_xc + sys.E_nl + sys.E_loc << std::endl;
    n ++;
    param.alpha = 0.3;    //param.alpha = 0.1;
  } while (n<15);
  /*
  for(i=0;i<ns.npw;i++) {
    for (j =0;j<ns.npw;j++) {
      printf("%5.3f ", Hnl[i][j][1]);
    }
    std::cout << std::endl;
    }
  */
  //check(ns, param, ck);
  print_rho(ns, param, nr);
  std::cout << "E_kin = " << sys.E_kin << std::endl;
  std::cout << "E_HT = " << sys.E_HT << std::endl;
  std::cout << "E_xc = " << sys.E_xc << std::endl;
  std::cout << "E_loc = " << sys.E_loc << std::endl;
  std::cout << "E_nl = " << sys.E_nl << std::endl;
  
  delete(G2); delete(Hkin); delete(Fnew); delete(Fold);
  delete(ck); delete(cr);   delete(nr);   delete(nk); 
  delete(VH); delete(Vloc); delete(Vxc); 
  for (i=0;i<ns.npw;i++)  { delete(Hnl[i]); delete(Hold[i]); }
  for (i=0;i<ns.ngrid3;i++) delete(G[i]);
  for (i=0;i<Nang;i++) delete(bkf[i]); delete(bkf);
  delete(Hnl); delete(G); delete(rho_old); delete(Hold);
  return 0;

}

void parameter_read(parmtr &param)
{
  param.a = 8.0; param.omega = pow(param.a, 3.);
  param.lv[1][0] = param.lv[2][0] = 0.; param.lv[0][0] = 1.0;
  param.lv[0][1] = param.lv[2][1] = 0.; param.lv[1][1] = 1.0;
  param.lv[0][2] = param.lv[1][2] = 0.; param.lv[2][2] = 1.0;
  param.rx[0][0] = param.rx[1][0] = param.rx[2][0] = 0.0; 
  
  param.rv[1][0] = param.rv[2][0] = 0.; param.rv[0][0] = 1.;
  param.rv[0][1] = param.rv[2][1] = 0.; param.rv[1][1] = 1.;
  param.rv[0][2] = param.rv[1][2] = 0.; param.rv[2][2] = 1.;
  param.frc[0] = 2.0; param.frc[1] = 1.0; param.frc[2] = 0.0; 
  param.frc[3] = 0.0; param.frc[4] = 0.0; param.frc[5] = 0.0;
  param.ecut = 2.0; param.Zn = 3.0;param.alpha = 1.0;
}

void pw_find (numset &ns, parmtr param, vec2d &pw)
{
  int i, j, k, npw_, imax;
  int n = 10;
  double k2, kmax;
  double Gcut2 = 2.*param.ecut*param.a*param.a/TWOPI/TWOPI;
  npw_ = 0;
  imax = 0;
  kmax = 0.0;
  for (i=-n;i<=n;i++) {    
    for (j=-n;j<=n;j++) {
      for (k=-n;k<=n;k++) {
        k2 = (double) (i*i + j*j + k*k);
        if ( k2 < Gcut2 ) {
          kmax = std::max(k2,kmax);
	  //std::cout << i << j << k << " " << std::endl;
          pw[0].push_back(i);
          pw[1].push_back(j);
          pw[2].push_back(k);
	  pw[3].push_back(0);
          npw_++;
        }
      }
    }

  }
  std::cout << "Number of plane wave is " << npw_ << " with " << imax 
            << " " << kmax << std::endl;
  ns.npw = npw_; 
  imax = int(sqrt(kmax)*2.);
  ns.ngrid = 12 ; //radix_decompose(2*imax+1);
  if (ns.ngrid != 2*(ns.ngrid/2)) {
    std::cout << "Number of Grid is not even " << ns.ngrid << std::endl;
    exit(1);
  }
  ns.ngrid3 = ns.ngrid*ns.ngrid*ns.ngrid;
  std::cout << "Number of grid is " << ns.ngrid << std::endl; 
}
