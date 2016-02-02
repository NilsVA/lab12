#include <iostream>
#include <fstream>
#include <complex>
#include <sstream>
//-----------------------------------
using namespace std;
//-----------------------------------
typedef complex<double> cmplx;
//-----------------------------------
void init( cmplx* const psi0, const double eta, const double sigma, const double dx,
          const int Nx);

void writeToFile(const cmplx* const v, const string s, const double dx,
                 const int Nx, const double xmin);

void linstep(cmplx* const psi0, cmplx* const psi1, const double dt,
	  const double dx, const int Nx);

void nlinstep(cmplx* const psi0, cmplx* const psi1, const double dt,
	  const double dx, const int Nx);
//-----------------------------------
int main(){

	const int Nx = 4000;
	const double L = 800;
	const double xmin = 0;
	const double Tend = 50;
	const double dx = L / (Nx - 1);
	const double dt = dx  / 10;
	const int Na = 10;
	int Nk = int(Tend / Na / dt + 0.5);

	const double eta = 0.2;

	stringstream strm;

	cmplx* psi0 = new cmplx[Nx];
	cmplx* psi1 = new cmplx[Nx];
	cmplx* h; 
	init(psi0, eta, dx, dt,Nx);

	writeToFile(psi0,"psi_0", dx,Nx,xmin);


	for (int i = 1; i <= Na; i++) {

		for (int j = 1; j <= Nk-1; j++) 
		{
		  linstep(psi0,psi1,dt,dx,Nx);
		  
		  h = psi1;
		  psi1 = psi0;
		  psi0 = h;
		  
		  nlinstep(psi0,psi1,dt,dx,Nx);
		  
		  h = psi1;
		  psi1 = psi0;
		  psi0 = h;
		}
		strm.str("");
		strm << "psi_" << i;
		writeToFile(psi0,strm.str(), dx,Nx,xmin);
	}

	return 0;
}
//-----------------------------------
void writeToFile(const cmplx* const v, const string s, const double dx,
                 const int Nx, const double xmin)
{
	ofstream out(s.c_str());
	for(int i=0; i<Nx; i++){
		double x = xmin + i * dx;
		out << x << "\t" << norm(v[i]) << "\t" << v[i].real() << "\t" << v[i].imag() << endl;
	}
	out.close();
}
//-----------------------------------

void init( cmplx* const psi0, const double eta,  const double dx, const double dt,
          const int Nx)
{
	const double x0 = dx*Nx * 0.5;
	const double f = sqrt(2) * eta;
	for(int i=0;i<Nx; i++){
		double x = i*dx - x0;
		psi0[i] = 5*f/cosh(eta * x);
	}
}

void linstep(cmplx* const psi0, cmplx* const psi1, const double dt,
	  const double dx, const int Nx)
{
  
const cmplx alpha = cmplx(0.0,1.0)*dt/dx/dx * (-1.0);
const cmplx alpha2 = alpha*alpha;
cmplx* d=new cmplx[Nx];

for(int i=0;i<Nx;i++) d[i] = 1.0 + 2.0*alpha;

for(int i=1;i<Nx;i++){
  d[i] -= alpha2/d[i-1];
  psi0[i] -= -alpha/d[i-1]*psi0[i-1];
  }
  
  psi1[Nx-1] = psi0[Nx-1]/d[Nx-1];

  for(int i=Nx-2;i>0; i--)

  psi1[i] = (psi0[i] + alpha*psi1[i+1])/d[i];

  delete[] d;
}

void nlinstep(cmplx* const psi0, cmplx* const psi1, const double dt,
	  const double dx, const int Nx)
{
  const cmplx alpha = cmplx(0.0,1.0)*dt/dx/dx * (-1.0);
  for(int i = 0; i< Nx ; i++)
  {
    psi1[i] = psi0[i]*exp((-1.)*cmplx(0.0,1.0)*norm(psi0[i])*dt);
  }
}