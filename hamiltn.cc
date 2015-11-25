#include <cstdio>
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <complex>
#include <vector>
#include "fftw3.h"
#include "dft_pwpp.h"
#define ZERO 1.e-20

int nint(double x) // Nearest integer function
{
  int y;
  x > 0.0 ? y=int(floor(x+0.5)): y=int(ceil(x-0.5));
  // if x > 0 then return floor or ceil.
  return y;
}

double angle_bw(parmtr param, int n1, int n2, int n3, 
                   int m1, int m2, int m3)
{
  double x = (double) n1 * param.rv[0][0] + (double) n2 * param.rv[0][1] + 
    (double) n3 * param.rv[0][2];
  double y = (double) n1 * param.rv[1][0] + (double) n2 * param.rv[1][1] + 
    (double) n3 * param.rv[1][2];
  double z = (double) n1 * param.rv[2][0] + (double) n2 * param.rv[2][1] + 
    (double) n3 * param.rv[2][2];
  double a = (double) m1 * param.rv[0][0] + (double) m2 * param.rv[0][1] + 
    (double) m3 * param.rv[0][2];
  double b = (double) m1 * param.rv[1][0] + (double) m2 * param.rv[1][1] + 
    (double) m3 * param.rv[1][2];
  double c = (double) m1 * param.rv[2][0] + (double) m2 * param.rv[2][1] + 
    (double) m3 * param.rv[2][2];
  double dot = x*a + y*b + z*c;

  
  double xx = x*x + y*y + z*z;
  double aa = a*a + b*b + c*c;
  if (xx < ZERO) dot = sqrt(xx);
  else if ( aa < ZERO) dot = sqrt(aa);
  return dot;
}

void Ypx_return(parmtr param, int n1, int n2, int n3, 
                    dcomplex &Y11, dcomplex &Y10, dcomplex &Y1_1)
{
  /*
  double x = (double) n1 * param.rv[0][0] + (double) n2 * param.rv[0][1] + 
    (double) n3 * param.rv[0][2];
  double y = (double) n1 * param.rv[1][0] + (double) n2 * param.rv[1][1] + 
    (double) n3 * param.rv[1][2];
  double z = (double) n1 * param.rv[2][0] + (double) n2 * param.rv[2][1] + 
    (double) n3 * param.rv[2][2];
  */
  double x,y,z;
  x = (double) n1; y = (double) n2; z = (double) n3;
  double r = sqrt(x*x + y*y + z*z);
  if (r < ZERO) {
   Y11.real = Y1_1.real = Y11.imag = Y10.imag = Y1_1.imag = 0.0; 
   Y10.real = sqrt(3./4./PI); 
  }
  else {
    double A = sqrt(3./8./PI);    
    Y11.real = -A*x/r; Y11.imag = -A*y/r;
    Y1_1.real = A*x/r; Y1_1.imag = -A*y/r;
    A = A*sqrt(2.);
    Y10.real = A*z/r; Y10.imag = 0.0;
    //std::cout << Y11.real << " " << Y11.imag << std::endl;
  }
}


void G2_initialize(numset &ns, parmtr param, double *G2)
{
  int i, j, k, i2, ix, iy, iz;  
  int id = 0;
  int nhalf = ns.ngrid/2;  


  for (i=0;i<ns.ngrid;i++) {
    if (i > nhalf) ix = i - ns.ngrid;
    else ix = i;
    for (j=0;j<ns.ngrid;j++) {
      if (j > nhalf) iy = j - ns.ngrid;
      else iy = j;
      for (k=0;k<ns.ngrid;k++) {
	if (k > nhalf) iz = k - ns.ngrid;
	else iz = k;
	i2 = ix*ix + iy*iy + iz*iz;
	id = k + ns.ngrid*(j + ns.ngrid*i);
	G2[id] = FOURPISQ*(double) i2/param.a/param.a;
	id ++;
      }
    }
  }
}

void V_HT(numset ns, parmtr param, double *G2, fftw_ *nk, fftw_ *VH, 
	  status &sys, vec2d pw)
{
  int i;
  sys.E_HT = 0.0;
  
  VH[0][0] = VH[0][1] = 0.0;
  for (i=1;i<ns.ngrid3;i++) {
    VH[i][0] = 4.*PI*nk[i][0]/G2[i];
    VH[i][1] = 4.*PI*nk[i][1]/G2[i];
    sys.E_HT += 
      2.*PI*param.omega*(nk[i][0]*nk[i][0] + nk[i][1]*nk[i][1])/G2[i];
  }
  
}

//
//
void V_loc(numset ns, parmtr param, double *G2, fftw_ *Vloc)
{
  int i;
  double rs, c1, c2, A, B;
  rs = 0.45;
  c1 = -6.8340578;
  c2 = 0.0;

  Vloc[0][0] = Vloc[0][1] = 0.0;
  for (i=1;i<ns.ngrid3;i++) {
    A = exp(-G2[i]*rs*rs*0.5);
    B = pow(2.*PI,1.5)*pow(rs,3.)*A*(c1 + c2*(3. - G2[i]*rs*rs))/param.omega
      - 4.*PI*param.Zn*A/G2[i]/param.omega;
    Vloc[i][0] = B;
    Vloc[i][1] = 0.0;
  }  
}
    
//
//
void E_loc(numset ns, parmtr param, double *G2, fftw_ *nk, fftw_ *Vloc, 
	   status &sys)
{
  int i;
  sys.E_loc = 0.0;
  
  for (i=1;i<ns.ngrid3;i++) {
    sys.E_loc += param.omega*(Vloc[i][0]*nk[i][0] + Vloc[i][1]*nk[i][1]);
  }
}

//
//
void V_nl(numset ns, parmtr param, double *G2, fftw_ **Hnl, vec2d pw, 
	  fftw_ **bkf)
{  
  int i, j, k, id, ix, iy, iz;
  double omega, rs, rp, h1s, h2s, h1p, Y00, k2, k2rs2, k2rp2, k2exp;
  dcomplex Y11, Y1m, Y10, A;
  double pi25 = pow(PI,2.5);

  omega = param.omega;
  rs = 0.4654363;
  rp = 0.5462433;
  h1s = 2.8140777;
  h2s = 1.9395165;
  h1p = 1.9160118;
  Y00 = 1./sqrt(4.*PI);
  //h1p = 0.0;
  for (i=0;i<ns.npw;i++) {
    ix = pw[0][i]; iy =  pw[1][i]; iz = pw[2][i]; id = pw[3][i];
    k2 = G2[id];
    k2rs2 = k2*rs*rs; 
    k2rp2 = k2*rp*rp;
    Ypx_return(param, ix, iy, iz, Y11, Y10, Y1m);

    bkf[0][i][0] = exp(-0.5*k2rs2);         bkf[0][i][1] = 0.0;
    bkf[1][i][0] = bkf[0][i][0]*(3.-k2rs2); bkf[1][i][1] = 0.0;

    k2exp = exp(-0.5*k2rp2)*sqrt(k2);
    bkf[2][i][0] = k2exp*Y11.real; bkf[2][i][1] = k2exp*Y11.imag;
    bkf[3][i][0] = k2exp*Y10.real; bkf[3][i][1] = k2exp*Y10.imag;
    bkf[4][i][0] = k2exp*Y1m.real; bkf[4][i][1] = k2exp*Y1m.imag;
  }

  h1s = h1s*Y00*Y00* 32.*pow(rs,3.)*pi25/omega;
  h2s = h2s*Y00*Y00*128.*pow(rs,3.)*pi25/omega/15.;
  h1p = h1p        * 64.*pow(rp,5.)*pi25/omega/3.;
  double alpha[Nang] = { h1s, h2s, h1p, h1p, h1p };

  double phase = 0.0;
  for (i=0;i<ns.npw;i++) {
    for (j=0;j<ns.npw;j++) {
      A.real = A.imag = 0.0;
      for (k=0;k<Nang;k++) {
	A.real += alpha[k]*(bkf[k][i][0]*bkf[k][j][0] + 
			    bkf[k][i][1]*bkf[k][j][1]);
	A.imag += alpha[k]*(bkf[k][i][1]*bkf[k][j][0] - 
			    bkf[k][i][0]*bkf[k][j][1]);
      }
      //if (i==j) A.real = A.imag = 0.0;
      Hnl[i][j][0] = A.real;
      Hnl[i][j][1] = A.imag;
    }
  }
}

void E_nl(numset ns, parmtr param, vec2d pw, fftw_ **bkf, fftw_ *ck, 
	  status &sys)
{
  int i, n, id, k;
  dcomplex F[Nang][ns.nband];
  double omega, rs, rp, h1s, h2s, h1p, Y00;
  
  omega = param.omega;
  rs = 0.4654363;
  rp = 0.5462433;
  h1s = 2.8140777;
  h2s = 1.9395165;
  h1p = 1.9160118;
  sys.E_nl = 0.0;
  Y00 = 1./sqrt(4.*PI);

  for (k=0;k<Nang;k++) {
    for (n=0;n<ns.nband;n++) {
      F[k][n].real = F[k][n].imag = 0.0;
    }
  }

  for (i=0;i<ns.npw;i++) {
    for (n=0;n<ns.nband;n++) {
      id = pw[3][i] + n*ns.ngrid3;
      for (k=0;k<Nang;k++) {
	F[k][n].real += bkf[k][i][0]*ck[id][0] + bkf[k][i][1]*ck[id][1];
	F[k][n].imag += bkf[k][i][1]*ck[id][0] - bkf[k][i][0]*ck[id][1];
      }
    }
  }

  h1s = h1s*Y00*Y00*32.*pow(rs,3.)*pow(PI,2.5)/omega;
  h2s = h2s*Y00*Y00*128.*pow(rs,3.)*pow(PI,2.5)/omega/15.;
  h1p = h1p*64.*pow(rp,5.)*pow(PI,2.5)/omega/3.;
  
  double alpha [Nang] = { h1s, h2s, h1p, h1p, h1p };

  for (n=0;n<ns.nband;n++) {    
    for (k=0;k<Nang;k++ ) {
      sys.E_nl += param.frc[n]*alpha[k]*(F[k][n].real*F[k][n].real + 
					 F[k][n].imag*F[k][n].imag);
    }
  }
}


//
//
void V_xc(numset ns, parmtr param, fftw_ *nr, fftw_ *Vxc, status &sys)
{
  int i;
  double asum, aasum, bsum, bbsum, exc, ndedn, rs;

  double a[4] = {0.4581652932831429, 2.217058676663745, 
                 0.7405551735357053, 0.01968227878617998};
  double b[4] = {1.0, 4.504130959426697, 1.110667363742916, 
                 0.02359291751427506};
  double omega = param.omega; 
  fftw_complex *Wxc; Wxc = new fftw_complex [ns.ngrid3];
  double dv = omega/(double) ns.ngrid3;
  sys.E_xc = 0.0;

  for (i=0;i<ns.ngrid3;i++) {
    //std::cout << nr[i][0] << std::endl;
    if (nr[i][0] > 0.0) {
      rs = pow(3./4./PI/nr[i][0], 1./3.);
      asum  = a[0] +    a[1]*rs +    a[2]*rs*rs +    a[3]*rs*rs*rs;
      aasum = a[1] + 2.*a[2]*rs + 3.*a[3]*rs*rs;
      bsum  = b[0]*rs + b[1]*rs*rs + b[2]*rs*rs*rs + b[3]*pow(rs,4.0);
      bbsum = b[0] + 2.*b[1]*rs + 3.*b[2]*rs*rs + 4.*b[3]*rs*rs*rs;
      exc = -asum/bsum;
      ndedn = rs*(aasum*bsum - asum*bbsum)/bsum/bsum/3.;
      Wxc[i][0] = exc + ndedn; Wxc[i][1] = 0.0;
      sys.E_xc += exc*nr[i][0]*dv;        
    }
    //else exit(1);
  }
  fftw_plan plan;
  plan = fftw_plan_dft_3d(ns.ngrid, ns.ngrid, ns.ngrid, Wxc, 
                          Vxc, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(plan); fftw_destroy_plan(plan); 
  
  for (i=0;i<ns.ngrid3;i++) {
    Vxc[i][0] = Vxc[i][0]/(double) ns.ngrid3;
    Vxc[i][1] = Vxc[i][1]/(double) ns.ngrid3;
  }
  delete(Wxc);
}

//
//
void H_kin(numset ns, vec2d pw, double *G2, double *Hkin)
{
  int i, id;

  for (i=0;i<ns.npw;i++) {
    id = pw[3][i];
    Hkin[i] = 0.5*G2[id];
  }
}

//
//
void E_kin(numset ns, parmtr param, vec2d pw, fftw_ *ck, double *Hkin, status &sys)
{
  int i, n, id;
  sys.E_kin = 0.0;
  
  for (i=0;i<ns.npw;i++) {
    for (n=0;n<ns.nband;n++) {
      id = pw[3][i] + ns.ngrid3*n;
      sys.E_kin += 
	Hkin[i]*param.frc[n]*(ck[id][0]*ck[id][0] + ck[id][1]*ck[id][1]);
    }
  }
}

//
//
//
void V_HGH(numset ns, parmtr param, double *G2, fftw_ **Hnl, vec2d pw, 
	  fftw_ **bkf)
{  
  int i, j, k, l, id, ix, iy, iz;
  double rs, rp, rd, rf, hs[3][3], hp[3][3], hd[3][3], hf, Y00, PI5_4;
  double omega, gsq, g2r0sq, r2r0sq;
  dcomplex Y11, Y1m, Y10, A;

  rs = 0.422738; hs[0][0] = 5.906928; hs[1][1] = 3.258196; hs[2][2] = 0.0;
  rp = 0.484278; hp[0][0] = 2.727013; hp[1][1] = 0.0;      hp[2][2] = 0.0;
  rd = 0.0;      hd[0][0] = 0.0;      hd[1][1] = 0.0;      hd[2][2] = 0.0;
  rf = 0.0;      hf       = 0.0;

  hp[0][0] = hp[1][1] = hp[2][2] = 0.0;

  hs[0][1] = hs[1][0] = -0.5*sqrt(3./5.)   *hs[1][1];
  hs[0][2] = hs[2][0] =  0.5*sqrt(5./21.)  *hs[2][2];
  hs[1][2] = hs[2][1] = -0.5*sqrt(100./63.)*hs[2][2];
  
  hp[0][1] = hp[1][0] = -0.5*sqrt(5./7.)*hp[1][1];
  hp[0][2] = hp[2][0] = sqrt(35./11.)*hp[2][2]/6.;
  hp[1][2] = hp[2][1] = -14.*hp[2][2]/6./sqrt(11.);

  hd[0][1] = hd[1][0] = -0.5*sqrt(7./9.)*hd[1][1];
  hd[0][2] = hd[2][0] = 0.5*sqrt(63./143.)*hd[2][2];
  hd[1][2] = hd[2][1] = -0.5*18./sqrt(143.)*hd[2][2];
  omega = param.omega;
  Y00 = 1./sqrt(4.*PI);
  PI5_4 = pow(PI,5./4.)/sqrt(omega);

  for (i=0;i<ns.npw;i++) {
    ix = pw[0][i]; iy = pw[1][i]; iz = pw[2][i]; id = pw[3][i];
    gsq = G2[id];

    g2r0sq = gsq*rs*rs; 
    Ypx_return(param, ix, iy, iz, Y11, Y10, Y1m);
    
    bkf[0][i][0] = 4.*sqrt(2.)*pow(rs,1.5)*PI5_4*exp(-0.5*r2r0sq);
    bkf[1][i][0] = bkf[0][i][0]*2.*(3.-g2r0sq)/sqrt(15.);
    bkf[2][i][0] = bkf[0][i][0]*4.*(15.-10.*g2r0sq+g2r0sq*g2r0sq)/3./sqrt(10.);
    bkf[0][i][1] = bkf[1][i][1] = bkf[2][i][1] = 0.0;

  }

  for (i=0;i<ns.npw;i++) {
    for (j=0;j<ns.npw;j++) {
      A.real = A.imag = 0.0;
      for (k=0;k<3;k++) {
	for (l=0;l<3;l++) {
	  A.real += hs[k][l]*(bkf[k][i][0]*bkf[l][j][0] + 
			      bkf[k][i][1]*bkf[l][j][1]);
	  A.imag += hs[k][l]*(bkf[k][i][1]*bkf[l][j][0] + 
			      bkf[k][i][0]*bkf[l][j][1]);
	}
      }
      Hnl[i][j][0] = A.real*Y00*Y00;
      Hnl[i][j][1] = A.imag*Y00*Y00;
      printf("%5.3f ", Hnl[i][j][0]);
      //std::cout << Hnl[i][j][0] << " ";
    }
    //std::cout << std::endl;
  } 
}

