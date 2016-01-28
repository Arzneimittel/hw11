#include <iostream>
#include <fstream>
#include <complex>
#include <sstream>
//-----------------------------------
using namespace std;
//-----------------------------------
typedef complex<double> cmplx;
//-----------------------------------
void init( cmplx* const psi0, const double alpha, const double lambda,
           const double dx, const double dt,
           const int Nx, const double xmin);


void writeToFile(const cmplx* const v, const string s, const double dx,
         const int Nx, const double xmin, const double alpha,
         const double lambda, const double omega, const double t);

void step(cmplx* const f1, cmplx* const f0,
          const double dt, const double dx,
          const double kv, const int N, const double xmin );
//-----------------------------------
int main(){

	const int Nx = 300;
	const double xmin =-40 ;
  const double xmax = 40;
	const double Tend = 10.0 * M_PI;
	const double dx = (xmax-xmin)/(Nx-1.0) ;
	const double dt = 0.01 * dx  ;
  double t = 0;
	const int Na = 10;
	int Nk = int(Tend / Na / dt + 0.5);

	const double lambda = 10;
  const double omega = 0.2;
  const double kv  = omega * omega;
  const double alpha = pow(kv,0.25);

  stringstream strm;

	cmplx* psi0 = new cmplx[Nx];
	cmplx* psi1 = new cmplx[Nx];
	cmplx* h = new cmplx[Nx];
	

	init(psi0, alpha, lambda, dx, dt, Nx, xmin);

	writeToFile(psi0,"psi_0", dx,Nx,xmin, alpha, lambda, omega,t);



	for (int i = 1; i <= Na; i++) {
		for (int j = 1; j <= Nk-1; j++) {
		  step(psi1,psi0,dt,dx,kv,Nx,xmin);
		  h = psi0;
		  psi0=psi1;
		  psi1=h;
		  t+=dt;
		}
		strm.str("");
		strm << "psi_" << i;
		writeToFile(psi0,strm.str(), dx,Nx,xmin, alpha, lambda, omega,t);
	}
  cout << "t = " << t << endl;
  delete[] psi0 ;
  delete[] psi1 ;
  delete[] h ;

	return 0;
}
//-----------------------------------
void step(cmplx* const f1, cmplx* const f0,
          const double dt, const double dx,
          const double kv, const int N, const double xmin )
{

  cmplx* d=new cmplx[N];  // Variable d
  cmplx* a=new cmplx[N]; // Konstante a
  cmplx* ak=new cmplx[N]; // komplex konj. a
  cmplx* dk=new cmplx[N]; // komplex konj. d
  cmplx* r=new cmplx[N]; // rechte Seite des Gleichungssystems
  double x;
  

  for(int i=0;i<N;i++){
    x = xmin + i*dx;
    d[i] = cmplx(1.0, dt/(2.0*dx*dx) + (dt/2.0) * 0.5 * x * x * kv);
    dk[i] = cmplx(1.0, -dt/(2.0*dx*dx) - (dt/2.0) * 0.5 * x * x * kv);
  } 
  for(int i=0;i<N;i++) a[i] = cmplx(0.0, - dt/(4*dx*dx));
  for(int i=0;i<N;i++) ak[i] = cmplx(0.0, + dt/(4*dx*dx));
  
  // rechte Seite vor der Anpassung
  r[0] = dk[0] * f0[0] + ak[0] * f0[1];
  r[N-1] = ak[0] * f0[N-2] + dk[N-1] * f0[N-1];
   for(int i=1;i<N-2;i++){
   r[i] = ak[0] * f0[i-1] + dk[i] * f0[i] + ak[0] * f0[i+1];    
  }
  
  // Forward substitution
  cmplx tri, zw1, zw2;
  
  for(int i=0; i<N-2; i++){
    tri = a[0] / d[i] ;
    zw1 = tri * a[i];
    d[i+1] -= zw1;
    zw2 = tri * r[i];
    r[i+1] -=  zw2; // rechte Seite nach Anpassung auf Zeilen-Stufen-Form.   
  }
  // Backward substitution
  f1[N-1] = r[N-1]/d[N-1];
  
  for(int i=N-2; i>=0; i--){
    f1[i] = (r[i] - a[0]*f1[i+1]) / d[i];   
  }
  
  delete[] d;
  delete[] a;
  delete[] dk;
  delete[] ak; 
  delete[] r;
}

//-----------------------------------
void writeToFile(const cmplx* const v, const string s, const double dx,
                 const int Nx, const double xmin, const double alpha,
                 const double lambda, const double omega, const double t)
{
	ofstream out(s.c_str());
  double x, xi, xil;
  double h1, h2, h3;
  cmplx ana;
	for(int i=0; i<Nx; i++){
		x = xmin + i * dx;
    xi = alpha * x;
    xil = alpha * lambda;
    h1 = -0.5 * pow(xi - xil*cos(omega*t),2 );
    h2 = omega*t/2 + xi * xil * sin(omega*t);
    h3 =  - 0.25 * xil*xil* sin(2*omega*t);
    ana = cmplx( h1 , h2 + h3  );
    ana = sqrt(alpha / sqrt(M_PI) ) * exp(ana);
		out << x << "\t" << norm(v[i]) << "\t" << v[i].real() << "\t" << v[i].imag()
         << "\t" << norm(ana) << "\t" << ana.real() << "\t" << ana.imag() <<  endl;
	}
	out.close();
}
//-----------------------------------

void init( cmplx* const psi0, const double alpha, const double lambda,
           const double dx, const double dt,
           const int Nx, const double xmin)
{
	const double x0 = dx*Nx * 0.5;
	for(int i=0;i<Nx; i++){
		double x = xmin + i*dx ;
		psi0[i] = sqrt(alpha/sqrt(M_PI)) * exp(- pow(alpha*(x-lambda),2)/2 );
	}
}
