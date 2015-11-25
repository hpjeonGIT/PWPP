#include <cstdio>
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <complex>
#include <vector>
#include "fftw3.h"
#include "dft_pwpp.h"

void find_E_fermi(numset ns, double &Efermi, double frc[])
{
  int i;
  double Ecrit, Emin, Emax, fsum, old_E, scaling;
  double kbT = 0.004;
  double Nel = 3.0;
  Emin = ns.Eband[0];
  Emax = ns.Eband[NBAND-1];
  Efermi = 0.5*(Emin + Emax);
  do {
    fsum = 0;
    for (i=0;i<NBAND;i++) {
      fsum += 2./(1. + exp((ns.Eband[i] - Efermi)/kbT));
    }
    Ecrit = fsum - Nel;
    if (fsum  > Nel) {
      Emax = Efermi;
      old_E = Efermi;
      Efermi = 0.5*(Efermi + Emin);
    }
    else {
      Emin = Efermi;
      old_E = Efermi;
      Efermi = 0.5*(Efermi + Emax);
    }
  } while (fabs(old_E - Efermi) > 1.E-8);
  //std::cout << "Fermi energy is " << Efermi << " with " << fsum << std::endl;
  scaling = Nel/fsum;
  for (i=0;i<NBAND;i++) {
    frc[i] =  2./(1. + exp((ns.Eband[i] - Efermi)/kbT));
    //std::cout << frc[i] << " ";
  }
  //std::cout << std::endl;

}

dcomplex dot_product(numset ns, fftw_ *ak, fftw_ *bk)
{
  int i;
  double real, imag;
  dcomplex  ans;
  real = imag = 0.0;
  for (i=0;i<ns.ngrid3;i++) {
    real += ak[i][0]*bk[i][0] + ak[i][1]*bk[i][1];
    imag += ak[i][0]*bk[i][1] - ak[i][1]*bk[i][0];
  }
  ans.real = real;
  ans.imag = imag;
  return ans;
}

//
// Random number generator
double rand_double() // random number generator
{
  double x = (double)(rand())/(double)(RAND_MAX) - 0.5;
  return x;
}

//
// Orthogonalization of reciprocal wave vectors
void gram_schmidt(numset ns, fftw_ *ck)
{
  int i, id, jd, n, m;
  dcomplex ans;
  fftw_complex *bk;
  double norm;
  bk = new fftw_complex [ns.nband*ns.ngrid3];

  for (n=0;n<ns.nband;n++) {
    for (i=0;i<ns.ngrid3;i++) {
      id = i + n*ns.ngrid3;
      bk[id][0] = ck[id][0]; bk[id][1] = ck[id][1];
    }
              
    for (m=0;m<n;m++) {
      id = n*ns.ngrid3;
      jd = m*ns.ngrid3;
      ans = dot_product(ns, &bk[jd], &ck[id]);
      for (i=0;i<ns.ngrid3;i++) {
        id = i + n*ns.ngrid3;
        jd = i + m*ns.ngrid3;
        bk[id][0] -= (ans.real*bk[jd][0] - ans.imag*bk[jd][1]); 
        bk[id][1] -= (ans.real*bk[jd][1] + ans.imag*bk[jd][0]);
      }
    }
    id = n*ns.ngrid3;
    ans = dot_product(ns, &bk[id], &bk[id]);
    norm = sqrt(ans.real);
    for (i=0;i<ns.ngrid3;i++) {
      id = i + n*ns.ngrid3;
      bk[id][0] = bk[id][0]/norm; bk[id][1] = bk[id][1]/norm;
    }
  }
  for (i=0;i<ns.nband*ns.ngrid3;i++) {
    ck[i][0] = bk[i][0]; ck[i][1] = bk[i][1];
  }
  delete(bk); 
}

//
// Initialization of reciprocal wave vectors using random numbers
void ck_initialize (numset ns, parmtr param, vec2d &pw, fftw_ *ck, fftw_ *nr, 
		    fftw_ *nk)
{
  int i, n, kx, ky, kz, id;
  double x;
  for (i=0;i<ns.ngrid3*ns.nband;i++) ck[i][0] = ck[i][1] = 0.0;
  for (i=0;i<ns.ngrid3;i++) nr[i][0] = nr[i][1] = nk[i][0] = nk[i][1] = 0.0;
  
  //
  // pw component -> grid element mapping
  for (i=0;i<ns.npw;i++) {
    kx = pw[0][i]; ky = pw[1][i]; kz = pw[2][i];
    kx = (kx + ns.ngrid)%ns.ngrid; 
    ky = (ky + ns.ngrid)%ns.ngrid; 
    kz = (kz + ns.ngrid)%ns.ngrid; 

    id = kz  + ns.ngrid*(ky + ns.ngrid*kx);
    pw[3][i] = id;
  }

  //
  // Initialization and normalization
  for (n=0;n<ns.nband;n++) {
    x = 0.0;
    for (i=0;i<ns.npw;i++) {
      id = pw[3][i] + n*ns.ngrid3;
      ck[id][0] = rand_double();
      ck[id][1] = rand_double();
      x += ck[id][0]*ck[id][0] + ck[id][1]*ck[id][1];
    }
    x = sqrt(x);
    for (i=0;i<ns.npw;i++) {
      id = pw[3][i] + n*ns.ngrid3;
      ck[id][0] = ck[id][0]/x;
      ck[id][1] = ck[id][1]/x;
    }
  }
  gram_schmidt (ns, ck);
}

//
// |phi> => <r|phi>
void ck_fft_cr(numset ns, parmtr param, fftw_ *ck, fftw_ *cr)
{
  int n;
  fftw_plan plan;
  for (n=0;n<ns.nband;n++) {
    plan = fftw_plan_dft_3d(ns.ngrid, ns.ngrid, ns.ngrid, &ck[ns.ngrid3*n], 
                            &cr[ns.ngrid3*n], FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan); fftw_destroy_plan(plan); 
  }
  double omega = param.omega;
	  
  for (n=0;n<ns.nband*ns.ngrid3;n++) {
    cr[n][0] = cr[n][0]/sqrt(omega);
    cr[n][1] = cr[n][1]/sqrt(omega);
  } 
}
				
//
//
void nr_fft_nk(numset &ns, parmtr &param, fftw_ *cr, fftw_ *rin_m, fftw_ *nk,
	       double *Fm, double *Fm_1, double **G, double *rin_m_1)
{
  int i, j, n, id;
  double *rout_m, GF[ns.ngrid3], dF[ns.ngrid3], Fsum, nel, drhoG[ns.ngrid3];
  double alpha = 0.5;
  double Efermi, kbT = 0.04; // Fermi surface smearing
  //Efermi = ns.Eband[1];
  rout_m = new double [ns.ngrid3];
  double dv = param.omega/(double) ns.ngrid3;
  double frc[NBAND];

  if (ns.nstage == 0) {
    for (i=0;i<ns.ngrid3;i++) {
      for (j=0;j<ns.ngrid3;j++) {
	G[i][j] = 0.0;
      }
      rout_m[i] = 0.0;
      Fm[i] = 0.0;
      rin_m_1[i] = 0.0;
      rin_m[i][0] = 0.0;
    }    
    for(n=0;n<ns.nband;n++) {
      for(i=0;i<ns.ngrid3;i++) {      
	id = i+n*ns.ngrid3;
	rin_m[i][0] += param.frc[n]*(cr[id][0]*cr[id][0] + cr[id][1]*cr[id][1]);
      }
    }
    ns.nstage = 1;
  }
  else if (ns.nstage == 1 ) {
    for (i=0;i<ns.ngrid3;i++) {
      rout_m[i] = 0.0;
    }
    find_E_fermi(ns, Efermi, frc);    
    for(n=0;n<ns.nband;n++) {
      for(i=0;i<ns.ngrid3;i++) {      
	id = i+n*ns.ngrid3;
	rout_m[i] += frc[n]*(cr[id][0]*cr[id][0] + cr[id][1]*cr[id][1]);
      }
    }
    for (i=0;i<ns.ngrid3;i++) {
      Fm_1[i] = Fm[i];
      Fm[i] = rout_m[i] - rin_m[i][0];
      G[i][i] = -alpha;
      rin_m_1[i] = rin_m[i][0];
    }
    for (i=0;i<ns.ngrid3;i++) {
      for (j=0;j<ns.ngrid3;j++) {
	rin_m[i][0] -= G[i][j]*Fm[j];
      }
    }
    ns.nstage = 2;
  }
  /*
    else if (ns.nstage >100) {
    for (i=0;i<ns.ngrid3;i++) rout_m[i] = 0.0;
    
    for(n=0;n<ns.nband;n++) {
      for(i=0;i<ns.ngrid3;i++) {      
	id = i+n*ns.ngrid3;
	rout_m[i] += frc[n]*(cr[id][0]*cr[id][0] + cr[id][1]*cr[id][1]);
      }
    }
    for (i=0;i<ns.ngrid3;i++) {      
      rin_m[i][0]  -=  0.1*(rin_m[i][0] - rout_m[i]);
    }
  } 
  */
  else {
    //std::cout << "density treatment3 " << std::endl;
    for (i=0;i<ns.ngrid3;i++) {
      rout_m[i] = 0.0;
    }
    //std::cout << "stage 1" << std::endl;
    find_E_fermi(ns, Efermi, frc);
    nel = 0.0;
    for(i=0;i<ns.ngrid3;i++) {      
      for(n=0;n<ns.nband;n++) {
	id = i+n*ns.ngrid3;
	rout_m[i] += frc[n]*(cr[id][0]*cr[id][0] + cr[id][1]*cr[id][1]);
      }
      //nel += dv* rout_m[i];
    }
    //std::cout << "Total number of electron " << nel << std::endl;
    Fsum = 0.0;
    for (i=0;i<ns.ngrid3;i++) {
      Fm_1[i] = Fm[i];
      Fm[i] = rout_m[i] - rin_m[i][0];
      dF[i] = Fm[i] - Fm_1[i];
      Fsum += dF[i]*dF[i];
    }
    for (i=0;i<ns.ngrid3;i++) {
      GF[i] = 0.0;
      drhoG[i] = 0.0;
      for (j=0;j<ns.ngrid3;j++) {
	GF[i] += G[i][j]*dF[j];
	drhoG[i] += (rin_m[j][0] - rin_m_1[j])*G[j][i];
      }
    }
    for(i=0;i<ns.ngrid3;i++) {
      for(j=0;j<ns.ngrid3;j++) {
	// Broyden's bad method
	G[i][j] += (rin_m[i][0] - rin_m_1[i] - GF[i])*dF[j]/Fsum;
	// Broyden's good method
	//G[i][j] += (rin_m[i][0] - rin_m_1[i] - GF[i])*drhoG[j]/Fsum;
      }
    }
    //std::cout << "stage 3" << std::endl;
    for (i=0;i<ns.ngrid3;i++) {
      rin_m_1[i] = rin_m[i][0];
    }
    for (i=0;i<ns.ngrid3;i++) {
      for (j=0;j<ns.ngrid3;j++) {
	rin_m[i][0] -= G[i][j]*Fm[j];
      }
      //if (rin_m[i][0] < 0.0) rin_m[i][0] = rin_m_1[i] +  alpha*Fm[i];
    }
    ns.nstage ++;
  }
  
  //std::cout << "stage 4" << std::endl;

  fftw_plan plan;
  plan = fftw_plan_dft_3d(ns.ngrid, ns.ngrid, ns.ngrid, rin_m,
			  nk, FFTW_FORWARD, FFTW_ESTIMATE); 
  fftw_execute(plan); fftw_destroy_plan(plan); 
  for (i=0;i<ns.ngrid3;i++) {
    nk[i][0] = nk[i][0]/(double) ns.ngrid3;
    nk[i][1] = nk[i][1]/(double) ns.ngrid3;
    //std::cout << nk[i][0] << " " << nk[i][1] << std::endl;
  }
}

//
//
void nr_fft_nk2(numset &ns, parmtr &param, fftw_ *cr, fftw_ *nr, fftw_ *nk)
{
  int i, n, id;
  double *rho;
  rho = new double [ns.ngrid3];
  double Efermi, frc[NBAND];

  for (i=0;i<ns.ngrid3;i++) rho[i] = 0.0;
  
  find_E_fermi(ns, Efermi, frc);
  if (ns.nstage == 0) {
    for(n=0;n<ns.nband;n++) {
      for(i=0;i<ns.ngrid3;i++) {      
	id = i+n*ns.ngrid3;
	rho[i] += param.frc[n]*(cr[id][0]*cr[id][0] + cr[id][1]*cr[id][1]);
      }
    }
    ns.nstage = 1;
  }
  else {
   for(n=0;n<ns.nband;n++) {
     frc[n] = param.frc[n];
     for(i=0;i<ns.ngrid3;i++) {      
       id = i+n*ns.ngrid3;
       rho[i] += frc[n]*(cr[id][0]*cr[id][0] + cr[id][1]*cr[id][1]);
     }
   }
  }
  for (i=0;i<ns.ngrid3;i++) {
    nr[i][0]  -=  param.alpha*(nr[i][0] - rho[i]);
  }
  fftw_plan plan;
  plan = fftw_plan_dft_3d(ns.ngrid, ns.ngrid, ns.ngrid, nr,
                          nk, FFTW_FORWARD, FFTW_ESTIMATE); 
  fftw_execute(plan); fftw_destroy_plan(plan); 
  for (i=0;i<ns.ngrid3;i++) {
    nk[i][0] = nk[i][0]/(double) ns.ngrid3;
    nk[i][1] = nk[i][1]/(double) ns.ngrid3;
  }
}
