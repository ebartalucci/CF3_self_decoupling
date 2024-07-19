/*

decouple_cwL.cc

CW decoupling 
additional spin diffusion simulated by Liouville space
relaxation matrix
use MAS small step integration
works for CH2 and CH3 groups

arbitrary orientation of the CSA tensor

*/

#include "gamma.h"
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>

#define MAXSPINS 6
#define NPROP 100

using namespace std;

//spin 0 is detected, spins 1 ... n are irradiated

int main(int argc, char *argv[])

{ spin_system ax;
  gen_op temp, Hrf, Ham, H[5], detect, sigma, sigma1;
  super_op U[NPROP];
  super_op Iop, L, L0;
  space_T Adip[MAXSPINS][MAXSPINS], Acsa[MAXSPINS], Acsa_R[MAXSPINS], Adip_R[MAXSPINS][MAXSPINS];
  double D[MAXSPINS][MAXSPINS], J[MAXSPINS][MAXSPINS], delta_CSA[MAXSPINS], eta_CSA[MAXSPINS], iso_CSA[MAXSPINS];
  double gamB1, offset;
  double k, k1;
  int i,j,m,n,Fnp,count,qu,steps;
  string name, names;
  const double thetam=54.73561032;
  double mas_freq, time_gamB1;
  double ltime, time, time_prop, scale;
  int nspins,index, n_prop;
  double alpha_CSA[MAXSPINS],beta_CSA[MAXSPINS],gamma_CSA[MAXSPINS];
  double alpha_D[MAXSPINS][MAXSPINS],beta_D[MAXSPINS][MAXSPINS],gamma_D[MAXSPINS][MAXSPINS];
  double alpha,beta,gamma;
  double theta;
  struct rusage me;
  char hostname[200];


  int value1[] = {2, 50, 100, 144, 200, 300, 538, 1154};
  int value2[] = {1,  7,  27,  11,  29,  37,  55,  107};
  int value3[] = {1, 11,  41,  53,  79,  61, 229,  271};

  count = 1;
  query_parameter(argc,argv,count++,"Spin System Name            ? ", names);
//setup for the spin system
  ax.read(names);
  nspins = ax.spins();
  if(nspins < 2 || nspins > MAXSPINS)
  { cerr << "This program is written for a two to MAXSPINS spin system.\n";
    cerr << "Please change your spin system definition in \n";
    cerr << "the file " << names << ".\n";
    cerr << "Aborting ....\n\n";
    exit(-1);
  }
  for(i=0;i<nspins-1;++i)
  { for(j=i+1;j<nspins;++j)
    { query_parameter(argc,argv,count++,"J       Coupling Constant   ? ", J[i][j]);
      query_parameter(argc,argv,count++,"Dipolar Coupling Constant   ? ", D[i][j]);
      query_parameter(argc,argv,count++,"Euler angle alpha           ? ", alpha_D[i][j]);
      query_parameter(argc,argv,count++,"Euler angle beta            ? ", beta_D[i][j]);
      query_parameter(argc,argv,count++,"Euler angle gamma           ? ", gamma_D[i][j]);
    }
  }
  for(i=0;i<nspins;++i)
  { query_parameter(argc,argv,count++,"CSA tensor (iso   CSA)      ? ", iso_CSA[i]);
    query_parameter(argc,argv,count++,"CSA tensor (delta CSA)      ? ", delta_CSA[i]);
    query_parameter(argc,argv,count++,"CSA tensor (eta CSA)        ? ", eta_CSA[i]);
    query_parameter(argc,argv,count++,"Euler angle alpha           ? ", alpha_CSA[i]);
    query_parameter(argc,argv,count++,"Euler angle beta            ? ", beta_CSA[i]);
    query_parameter(argc,argv,count++,"Euler angle gamma           ? ", gamma_CSA[i]);
  }
  query_parameter(argc,argv,count++,"Field Strength gamB1        ? ", gamB1);
  query_parameter(argc,argv,count++,"Decoupler Offset            ? ", offset);
  query_parameter(argc,argv,count++,"Rate constant k (1/s)       ? ", k);
  query_parameter(argc,argv,count++,"Powder Quality (cheng)      ? ", qu);
  query_parameter(argc,argv,count++,"Number of time steps        ? ", steps);
  query_parameter(argc,argv,count++,"spinning speed              ? ", mas_freq);
  query_parameter(argc,argv,count++,"Number of sampling per rotor? ", n_prop);
  query_parameter(argc,argv,count++,"Number of Data Points       ? ", Fnp);
  query_parameter(argc,argv,count++,"Output Filename             ? ", name);


  time_prop  = fabs((1.0/mas_freq)/n_prop);
  time       = fabs((1.0/mas_freq)/steps);
  time_gamB1 = 1.0/gamB1;
  theta      = atan2(gamB1,offset);
  k1=k*(3.0*cos(theta)*cos(theta)-1)*(3.0*cos(theta)*cos(theta)-1)/4.0;


  block_1D data(Fnp);
  for(i=0;i<Fnp;++i)
    data(i) = 0;

  (void) gethostname(hostname,199);

  cout << "\n\nSimulation of isotropic chemical shift by dipolar coupling\n";
  cout << "==========================================================\n\n";
  cout << "Program version: " << __FILE__ << " compiled at " << __DATE__ ", "
       << __TIME__ << "\n";
  cout << "running on machine " << hostname << "\n\n";
  cout << "Parameters:\n";
  cout << "rotation angle thetam    : " << thetam << " Degree\n";
  cout << "size of spin system            : " << nspins << " spins\n";
  for(i=0;i<nspins-1;++i)
  { for(j=i+1;j<nspins;++j)
    { cout << "J       coupling constant (" << i << "," << j << ") : " << 
              J[i][j] << " Hz\n";
      cout << "dipolar coupling constant (" << i << "," << j << ") : " << 
              D[i][j] << " Hz\n";
      cout << "relativ orientation of D tensor: (" << alpha_D[i][j] << "," <<
               beta_D[i][j] << "," << gamma_D[i][j] << ")\n";
    }
  }
  for(i=0;i<nspins;++i)
  { cout << "iso   of CSA tensor (" << i << ") : " << 
            iso_CSA[i] << " Hz\n";
    cout << "delta of CSA tensor (" << i << ") : " << 
            delta_CSA[i] << " Hz\n";
    cout << "eta of CSA tensor (" << i << ") : " << 
            eta_CSA[i] << "\n";
    cout << "relativ orientation of CSA tensor: (" << alpha_CSA[i] << "," <<
             beta_CSA[i] << "," << gamma_CSA[i] << ")\n";
  }
  cout << "rf-Field Strength:         " << gamB1 << " Hz\n";
  cout << "rf-Field Offset:           " << offset << " Hz\n";
  cout << "Rate constant (k):         " << k << " 1/s\n";
  cout << "Spin-lock angle:           " << theta*180/PI << " Degree\n";
  cout << "Effective rate constant(k):" << k1 << " 1/s\n";
  cout << "# of sampling points:      " << n_prop << "\n";
  cout << "Powder Quality Number:     " << qu << " (" << value1[qu] << " points)\n";
  cout << "Number of data points:     " << Fnp << " points\n";
  cout << "basic rf pulse frequency:  " << gamB1 << " Hz\n";
  cout << "length of one prop. (dw):  " << time_prop << " s\n";
  cout << "MAS frequency:             " << mas_freq << " Hz\n";
  cout << "time increments:           " << time << "s\n";
  cout << "Output filename:           " << name << "\n";
  cout << "\n";



//these are the constant parts of H and L

//these are the constant parts of H and L
  matrix exchange(pow(2.0,2*nspins),pow(2.0,2*nspins),0);
  for(i=0;i<pow(2.0,2*nspins);++i)
    exchange.put_h(1,i,i);
  Iop = super_op(exchange);

  temp = Fz(ax,"19F");
  L0  = k1*d_commutator(temp);
  temp = Fx(ax,"19F");
  L0 += k1*d_commutator(temp);
  temp = Fy(ax,"19F");
  L0 += k1*d_commutator(temp);

  detect = Im(ax,0);
  sigma  = Ix(ax,0);
  Hrf   = gamB1 * Fx(ax,"19F")+offset * Fz(ax,"19F");

//setup for the space tensor
  matrix help(3,3,0);
  for(i=0;i<nspins-1;++i)
  { for(j=i+1;j<nspins;++j)
    { help.put_h(-1.0/2.0,0,0);
      help.put_h(-1.0/2.0,1,1);
      help.put_h( 1.0,2,2);
      help = - (complex) D[i][j] * help;
      Adip[i][j] = A2(help);
      Adip[i][j] = Adip[i][j].rotate(alpha_D[i][j],beta_D[i][j],gamma_D[i][j]);
    }
  }
  for(i=0;i<nspins;++i)
  { help.put_h(-1.0/2.0*(1.0+eta_CSA[i]),0,0);
    help.put_h(-1.0/2.0*(1.0-eta_CSA[i]),1,1);
    help.put_h( 1.0,2,2);
    help = (complex) delta_CSA[i] * help;
    Acsa[i] = A2(help);
//rotate CSA tensor into the PA system of the dipolar tensor
    Acsa[i] = Acsa[i].rotate(alpha_CSA[i],beta_CSA[i],gamma_CSA[i]);
  }
  string name1 = name+".mat";
  string name2 = name;

//here starts the powder loop
//reference JCP 59 (8), 3992 (1973).

  for(count=1; count<value1[qu]; ++count)
  { beta  = 180.0 * count/value1[qu];
    alpha = 360.0 * ((value2[qu]*count) % value1[qu])/value1[qu];
    gamma = 360.0 * ((value3[qu]*count) % value1[qu])/value1[qu];
    if(count % 10 == 1)
    { getrusage(0, & me);
      cout << count << "\tbeta = " << beta << "\talpha = "
           << alpha << "\tgamma = " << gamma 
           << ",\ttime used: " << me.ru_utime.tv_sec << " seconds\n";
      cout.flush();
    }
    scale = sin(beta/180.0*PI);

//now we rotate the space tensor
    for(i=0;i<nspins-1;++i)
    { for(j=i+1;j<nspins;++j)
      { Adip_R[i][j] = Adip[i][j].rotate(alpha,beta,gamma);
      }
    }
    for(i=0;i<nspins;++i)
    { Acsa_R[i] = Acsa[i].rotate(alpha,beta,gamma);
    }

//now we can combine the hamiltonians for the different side diagonals

//zero all components
    for(i=0;i<5;++i)
      H[i] = gen_op();

  for(i=0;i<n_prop;++i)
    U[i] = Iop;

//this is the dipolar part
  for(j=-2;j<=2;++j)
  { for(m=0;m<nspins-1;++m)
    { for(n=m+1;n<nspins;++n)
      { if(ax.isotope(m) != ax.isotope(n))
        { H[j+2] += Adip_R[m][n].component(2,j) * d2(j,0,thetam) * 1/sqrt(6.0)*2*Iz(ax,m)*Iz(ax,n);
          if(j == 0)
            H[j+2] += J[m][n]*Iz(ax,m)*Iz(ax,n);
        }
        else
        { H[j+2] += Adip_R[m][n].component(2,j) * d2(j,0,thetam) * 1/sqrt(6.0)*
                    (2*Iz(ax,m)*Iz(ax,n)-(Ix(ax,m)*Ix(ax,n)+Iy(ax,m)*Iy(ax,n)));
          if(j == 0)
          { H[j+2] += J[m][n]*(Iz(ax,m)*Iz(ax,n)+Ix(ax,m)*Ix(ax,n)+Iy(ax,m)*Iy(ax,n));
	  }
        }
      }
    }
    for(m=0;m<nspins;++m)
    { H[j+2] += Acsa_R[m].component(2,j) * d2(j,0,thetam) * 1/sqrt(6.0)*2*Iz(ax,m);
      if(j == 0)
      { H[j+2] += iso_CSA[m]*Iz(ax,m);
      }
    }
  }
//now we calculate the propagator over one cycle of the MAS
    for(ltime=0.5*time;ltime<=steps*time;ltime += time)
    { Ham = Hrf;
      for(i=-2;i<=2;++i)
        Ham += exp(complex(0,i*2.0*PI*ltime*mas_freq)) * H[i+2];
      index = int(ltime/time_prop);
      Ham.set_DBR();
      L = complex(0,2.0*PI)*commutator(Ham);
      L = L+L0;
      U[index] &= exp(L,-time);
    }
    for(i=0;i<n_prop;++i)
    { U[i].set_HBR();
    }
    detect.set_DBR();
    sigma1=sigma;
    sigma1.set_DBR();
    for(i=0;i<Fnp;++i)
    { data(i) += trace(detect,sigma1)*scale;
      sigma1 = U[i%n_prop]*sigma1;
    }
  } // end of powder loop
  MATLAB(name1,name2,data,1);
  exit(0);
}
