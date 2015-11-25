#include <cstdio>
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <complex>
#include <vector>
#include "fftw3.h"
#include "dft_pwpp.h"



extern "C" void zheev_(char *job, char *uplo, int *n,  dcomplex *a, int *lda, 
                       double *w, dcomplex *work,  int *lwork, double *rwork, 
                       int *info);
void gram_schmidt(numset ns, fftw_complex *ck);

//////////////////////////////////////////////////////////////////////
void Hermitian_lapack(numset &ns, parmtr param, vec2d pw, fftw_ *VH,fftw_ *Vxc,
		      fftw_ *Vloc, double *Hkin, fftw_ **Hnl, fftw_ *ck, 
		      fftw_ **Hold)
{

  int i, j, n, ix, iy, iz, jx, jy, jz, kx, ky, kz, id, id2;
  double phase, A, B;
  double *w=NULL, *rwork=NULL;
  dcomplex *a=NULL, *work=NULL;
  int info=0, lwork;
  int npw = ns.npw;
  a = new dcomplex [ns.npw*ns.npw];
  w = new double [ns.npw];
  rwork = new double [3*ns.npw-2];
  
  for (i=0;i<ns.npw;i++) {
    ix = pw[0][i]; iy = pw[1][i]; iz = pw[2][i];
    for (j=0;j<ns.npw;j++) {
      jx = pw[0][j]; jy = pw[1][j]; jz = pw[2][j];
      kx = ix - jx;kx = (kx+ns.ngrid)%ns.ngrid;
      ky = iy - jy;ky = (ky+ns.ngrid)%ns.ngrid;
      kz = iz - jz;kz = (kz+ns.ngrid)%ns.ngrid;
      id = kz + ns.ngrid*(ky + ns.ngrid*kx);

      A = VH[id][0] + Vxc[id][0] + Vloc[id][0] + Hnl[i][j][0];
      B = VH[id][1] + Vxc[id][1] + Vloc[id][1] + Hnl[i][j][1];

      a[i+ns.npw*j].real = A; //(1.-param.alpha)*Hold[i][j][0] + param.alpha*A;
      a[i+ns.npw*j].imag = B; //(1.-param.alpha)*Hold[i][j][1] + param.alpha*B;
      Hold[i][j][0] = A;
      Hold[i][j][1] = B;
      
      //a[i+ns.npw*j].real = A;
      //a[i+ns.npw*j].imag = B;
      //std::cout <<  VH[id][0] <<  Vxc[id][0] << Vloc[id][0] << std::endl;
      if (i == j)  a[i+ns.npw*j].real += Hkin[i];
      //printf("%5.3f ", a[i+ns.npw*j].imag);
    }
    //std::cout << std::endl;
  }
  work = new dcomplex [1];  
  lwork = -1;
  // query
  zheev_( "V", "U", &npw, a, &npw, w, work, &lwork, rwork, &info );
  lwork = (int) work[0].real;
  delete(work);
  work = new dcomplex [lwork];

  // actual solving
  zheev_( "V", "U", &npw, a, &npw, w, work, &lwork, rwork, &info );
  // w - eigen value   If INFO = 0, the eigenvalues in ascending order.
  

  std::cout << "Min energy is " << w[0] << " "  << w[1] << " " << w[2] 
            << " " << w[3] <<  " " << w[4] <<std::endl;

  for (n=0;n<ns.nband;n++) {
    ns.Eband[n] = w[n];
    for (i=0;i<ns.npw;i++) {
      id = pw[3][i] + n*ns.ngrid3;
      id2 = i + n*ns.npw;
      ck[id][0] = a[id2].real;
      ck[id][1] = a[id2].imag;
    }
  }
  gram_schmidt(ns, ck);   
  delete(w); delete(rwork); delete(a); delete(work);
}


