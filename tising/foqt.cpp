#include <cstring>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <complex>
#include <math.h>

#include <armadillo>
// hi A&M
using namespace arma;
typedef std::complex<double> complex;

double g, hi, hf, h, tau;

// L - SIZE	s - HSCALING	g	u - TAUSCALING		c - CYCLE		d - dtime

class chain {
	public:
		int L, states;
		int PBC;
		Mat<int> spin;
		cx_dvec gs, fes;
		/*sp_*/cx_dmat Ham, rho;//, Ham0, Vp;
		double energy;
		double Mval, Aval, Zval;
		double time;
		
		chain(int _L, int _pbc) {
			L = _L;
			states = std::pow(2, _L);
			PBC = _pbc;
			spin.zeros(states, _L);
			gs.zeros(states);
			fes.zeros(states);
			//Ham0.set_size(states, states);
			Ham.set_size(states, states);
			rho.set_size(states, states);
			//Vp.set_size(states, states);
			energy = 0.;
			time = hi*tau;
			Mval = 0.;
			Zval = 0.;
			for(int i = 1; i <= states; ++i) {
				int z = i-1;
				for(int j=0; j < L; ++j) {
					spin(i-1, j) = int(std::pow(-1, z%2+1));
					z = int(z/2);
				}
			
			}
		}
		
		void BuildHam(void);
		inline void GibbsState(double Temperature);
		//void GroundState(void);
		//void twolevels(void);
		void Measure(int point);
		inline void rhoMeasure(void);
		
		void derivate( double t, const cx_dvec& _gs, cx_dvec& k);
		void RungeKuttaPsi(double _dt);
		
		void SuzukiTrotter(double _dt);
		
};


void chain::BuildHam(void) {

		Ham.zeros();
		double temp;
		for(int i = 1; i <= states; ++i) {
			temp = 0.;
			for(int j = 1; j < L; ++j) {
				int est1 = i - spin(i-1,j-1)*std::pow(2, j-1);
				Ham(est1-1,i-1) += - h;
				//Vp(est1-1,i-1) = - 1.;
				int est2 = est1-spin(est1-1,j+1-1)*std::pow(2, j);
				Ham(est2-1,i-1) += - 1.;
				temp += - g * spin(i-1,j-1);
			}
			int est1 = i - spin(i-1,L-1)*std::pow(2, L-1);
			Ham(est1-1,i-1) += - h;
			//Vp(est1-1,i-1) = - 1.;
			temp += - g*spin(i-1,L-1);
			Ham(i-1,i-1) += temp;
		}
		if(PBC==1) {
			for(int i = 1; i <= states; ++i) {
				int est1 = i - spin(i-1,0);
				int est2 = est1-spin(est1-1,L-1)*std::pow(2, L-1);
				Ham(est2-1,i-1) += - 1.;
			}
		}
		//Ham = Ham0 + h*Vp;
}

void chain::GibbsState(double Temperature){

		if(Temperature == 0){
			//rho = expmat_sym( -Ham0/Temperature );
			dvec eigval;
			cx_dmat eigvec;
			//vec = eigs_sym(HH);
			//std::cout << "ererere\n";
			//eigs_gen(eigval, eigvec, Ham0, states-1, "sr");
			eig_sym(eigval, eigvec, Ham);
			for(int i=0; i<states; ++i)
				for(int j=0; j<states; ++j)
				//gs(i) = (eigvec.col(0))(i);// complex( (eigvec.col(0))(i), 0.);
					rho(i,j) = (eigvec.col(0))(i) * conj((eigvec.col(0))(j));
			rho = rho / trace(rho);
			return;
		}
/*
		//rho = expmat_sym( -Ham0/Temperature );
		dvec eigval;
		cx_dmat eigvec;
		//vec = eigs_sym(HH);
		//std::cout << "ererere\n";
		//eigs_gen(eigval, eigvec, Ham0, states-1, "sr");
		eig_sym(eigval, eigvec, Ham);

		cx_dmat temp(states,states);
		for(int i=0; i<states; ++i){
			temp(i,i) = exp(-eigval(i)/Temperature);
			//std::cout << eigval(i) << "a\n";
		}
		rho = eigvec * (temp * conj(trans(eigvec)));
		//std::cout << temp << "		\n\n" << conj(trans(eigvec)) << "\n";
		//exit(8);
		complex trace_rho = trace(rho);
		rho = rho / trace_rho;

		cx_dmat rho1;
		rho1.set_size(states, states);
*/
		rho = expmat_sym(-Ham/Temperature);
		rho = rho / trace(rho);
		//std::cout << rho << "		\n\n" << rho1 << "\n";
		//exit(8);
		//rho = rho1;
}

/*
*
*
void chain::GroundState(void) {

		//sp_dmat HH = real(Ham);
		//eigs_opts opts;
		//opts.maxiter = 1000;
		cx_dvec eigval;
		cx_dmat eigvec;

		//vec = eigs_sym(HH);
		eigs_gen(eigval, eigvec, Ham, 1, "sr");
		for(int i=0; i<states; ++i)
			gs(i) = (eigvec.col(0))(i);// complex( (eigvec.col(0))(i), 0.);
		gs = gs / norm(gs);

		//std::cout  << eigval <<"8gsgsgs\n" ;

}

void chain::twolevels(void) {

		//sp_dmat HH = real(Ham);
		//eigs_opts opts;
		//opts.maxiter = 1000;
		cx_dvec eigval;
		//dvec eigval;
		cx_dmat eigvec;
		int indexgs;
		//dmat rHam

		//vec = eigs_sym(HH);
		eigs_gen(eigval, eigvec, Ham, 2, "sr");
		//eigs_sym(eigval, eigvec, Ham, 2, "sa");
		if (real(eigval(0)) > real(eigval(1)))
			indexgs = 1;
		else
			indexgs = 0;
		
		for(int i=0; i<states; ++i){
			gs(i) = (eigvec.col(indexgs))(i);// complex( (eigvec.col(0))(i), 0.);
			fes(i) = (eigvec.col(1-indexgs))(i);
		}
		gs = gs / norm(gs);
		fes = fes / norm(fes);
		//std::cout << real(eigval) << "\n";
		energy = real(eigval(1-indexgs) - eigval(indexgs));
		//std::cout  << eigvec.size() <<"\n" ;

}
*/


void chain::Measure(int point) {

	complex obs(0.), obs1(0.);
	assert(point > 0 && point <= L);

	for(int j=1; j <= L; ++j) {
		for(int i=1; i <= states; ++i) {
		
//std::cout << j << " " << i <<  " aa\n";
			//obs +=(gs(i-1)*gs.t()(i-1))*complex(spin(i-1, point-1),0.);
			int est = i - spin(i-1, j-1)*std::pow(2, (j-1));
			//std::cout << est << " " << gs <<  " aa\n";
			obs += conj(gs(i-1))*gs(est-1);
			obs1 += complex(spin(i-1, j-1))* conj(gs(i-1)) * gs(i-1);
		}
	}

	Mval = obs.real() / double(L);
	Zval = obs1.real() / double(L);
	/*
	for(int i=0; i < gs.size(); ++i)
		obs += ( gs(i) * gs.t()(i) ) * complex(spin(i, point-1), 0.);
	Mval = obs.real();
	*/
	//delete &obs;
}

void chain::rhoMeasure(void){

	complex temp(0.), temp1(0.);
	for(int point=1; point <= L; ++point ){
		for(int i=1; i<=states; ++i){
			int est = i - spin(i-1,point-1)*std::pow(2, point-1);
			temp = temp + rho(i-1,est-1);
			temp1 += rho(i-1,i-1) * double(spin(i-1, point-1));
		}
	}
	Mval = temp.real() / double(L);
	Zval = temp1.real() / double(L);

	//std::cout << "Mx" << Zval << '\n';
}

void chain::derivate( double t, const cx_dvec& _gs, cx_dvec& k) {
	//int dim = _gs.size();
	//SparseMatrix<complex> Ham(states, states);
	//std::cout << _Ham  << "\n";
	//std::cout << "a01\n";
	//SparseMatrix<complex> HH(Ham.rows(), Ham.cols());
	//for(int i=0; i<Ham.cols(); ++i)
	//	for(int j=0; j<Ham.rows(); ++j)
	//		HH(i, j) = complex(Ham(i,j), 0.);
	k = complex(0., -1.) * ( Ham * _gs );
	//std::cout << "a02\n";
}

void chain::RungeKuttaPsi(double _dt) {

	int dim(gs.size());
	cx_dvec k1(dim), k2(dim), k3(dim), k4(dim), temp;
	double time(0.);
//std::cout << "aa\n";
	derivate( time         , gs            , k1);
//	std::cout << "aa1\n";
	temp = gs + k1*_dt/2.;
	derivate( time + _dt/2., temp, k2);
	temp = gs + k2*_dt/2.;
	derivate( time + _dt/2., temp, k3);
	temp = gs + k3*_dt;
	derivate( time + _dt   , temp, k4);

	gs = gs + _dt/6.*( k1 + 2.*k2 + 2.*k3 + k4 );
	//time = time + _dt;
	//gs.normalise();
	gs = gs / norm(gs);
	//delete &k1, &k2, &k3, &k3, &temp;

}





int main(int argc, char *argv[]) {

	int Lsize(6); int pbc(0);// int n_cycle(1), yy(0);
	double Temperature(0.001), eta(0.);//double dt(0.02);
	g = 1.; hi = -0.25; hf = 0.25;// tau = 100;
	double upsilon = 0.01, sigma = 1.;//double sigmaf = -7.;// , omega = 1.
	char folder[10] = "data";
	std::string lineterm= "";
	for (int i=1;i<argc;i++)
		lineterm.append(std::string(argv[i]).append(" "));

while( argc > 1 ) {

	switch(argv[1][0]) {
		case 'L':
				Lsize = atoi( &argv[1][1] );
			break;
		case 'g':
				g = atof( &argv[1][1] );
			break;
		case 'T':
				Temperature = atof( &argv[1][1] );
			break;
		case 'e':
				eta = atof( &argv[1][1] );
			break;
		case 'h':
			if(argv[1][1] == 'i')		hi = atof( &argv[1][2] );
			else if(argv[1][1] == 'f')	hf = atof( &argv[1][2] );
			else				hi = atof( &argv[1][1] );
			break;
		case 'b':
				pbc = atoi( &argv[1][1] );
			break;
		case 'd':
				//if(argv[1][1] == 'S')	ds = atof( &argv[1][2] );
				///*else*/				dt = atof( &argv[1][1] );
			break;
		case 'c':
				//n_cycle = atoi( &argv[1][1] );
			break;
		case 't':
				tau = atof( &argv[1][1] );
			break;
		case 'u':
				upsilon = atof( &argv[1][1] );
			break;
		//case 'o':
				//omega = atof( &argv[1][1] );
			//break;
		case 's':
				if(argv[1][1] == 'R')	sigma = -std::sqrt(atof( &argv[1][2] ));
				//else if(argv[1][1] == 'f')	sigmaf = atof( &argv[1][2] );
				else			sigma = atof( &argv[1][1] );
			break;
		case '-':
				std::sprintf(folder, "%s", &argv[1][1]);
			break;
		case 'y':
				//yy = atoi(&argv[1][1]);
			break;
		default:
			std::cerr << "Unlucky: Retry input values\n";
			exit (8);
	}
	++argv;
	--argc;
}

/*
char output[80];
if(yy == 0)
	std::sprintf(output, "data/f5otqt%du%dg%dL%d.dat", int(abs(sigma)), int(1000*upsilon), int(10*g), Lsize);
else
	std::sprintf(output, "%s/fotqu%.0fh%.0fl%d.dat", folder, upsilon*10., 10.*abs(hi), Lsize );

std::ofstream out_file(output, std::ios::out | std::ios::trunc);
out_file.precision(16);
*/
std::cout.precision(16);

if (eta != 0.) Temperature = eta / Lsize;
h = 0.;
g = sigma / std::pow(Lsize, 1.) + 1.;
std::cout << "# Lsize		sigma		ups		eta		L^(1/8) Mx		Mz\n";
for(int hindex=0; hindex<50; ++hindex){
	//h = (upsilon + hindex*0.01 )/ std::pow(Lsize, 15./8.);
	g = (sigma - hindex*0.1)  / std::pow(Lsize, 1.) + 1.;
	class chain C(Lsize, pbc);

	C.BuildHam();
	C.GibbsState(Temperature);
	C.rhoMeasure();

	std::cout << Lsize << "	" << g << "	" << upsilon+ hindex*0.01 << "	" << eta << "	" << std::pow(Lsize, 1./8.)*C.Mval << "		" << C.Zval << '\n';
}
/*
class chain C(Lsize, pbc);
class chain E(Lsize, pbc);

C.BuildHam();
//C.GroundState();
h = 0.000001;

E.BuildHam();
E.Ham = E.Ham0 + h * E.Vp;
//E.GroundState();
E.twolevels();
E.Measure(3);
double Szc(E.Zval), DeltaL ;//= E.energy;
DeltaL = 2.*std::pow(g, Lsize)*(1.- g*g);
double M0g = std::pow(1. - g*g, 1./8.);//M0g = E.Mval;

tau = Lsize * upsilon / std::pow(DeltaL, 2.);

out_file << "#		u = " << upsilon / std::pow(DeltaL, 2.) << "		v = " << upsilon << "		sqrt(u) = " << std::pow(upsilon, 0.5)/DeltaL << "		comm_line = " << lineterm <<'\n';
out_file << "# t/ sqrt(u)		M_long/M00		(M_trans - Mzc)/(Mz00-Mzc)		A		B		Mx/M0g		Mx" << "\n";

//E = C;
cx_dvec temp_gs(C.states);

if( yy == 0) {

	sigma = sigma / std::pow(upsilon, 0.5);

	hi = sigma * DeltaL / Lsize;
	hf = - sigma * DeltaL / Lsize ;
}

C.time = hi*tau;
h = hi;
	//h = 0.5;
C.Ham = C.Ham0 + h * C.Vp;
C.GroundState();
//exit(8);
C.Measure(3);
const double M00(C.Mval), Z00(C.Zval);


temp_gs = C.gs;

double dir(1.), sum(hi);//, t00(hi*tau), tmax(4*n_cycle*hi*tau + t00);
int events = 2*n_cycle*int((tau*hf - tau*hi)/dt) + 1;
int frac = 10, prnt = int(2.*abs(sigma)*20);

if( prnt > events ) prnt = events;

hi = std::abs(hi);


std::cout << M0g << " = M0g		M00 = " << M00 <<'\n';

for(int sweep = 0; sweep < events; ++sweep) {
		if ( sweep%(events/frac) == 0) {
			std::cout << int(sweep/float(events)*100.) << "%   " << output << "\n";
			out_file << std::flush;
		}
		
		C.time += dt;
		sum += dir * dt / tau;
		//h = C.time / tau;
		//if ( int( (C.time-t00)/tau/2./hi ) % 2 == 0 )

		//if( (sweep+1) % (events/2/n_cycle) == 0 ) dir = -dir;

		C.Ham = C.Ham + C.Vp * dir * (dt / tau);
		//h = 0.;

		//C.SuzukiTrotter(dt);
		C.RungeKuttaPsi(dt);
		//C.BuildHam(true);

		//E.time += dt;
			
		//if( sum + dir*dt/tau > hf ) {
		if(sweep%(events/prnt) == 0) {
			C.Measure(3);
			
			E.Ham = C.Ham;
			//E.GroundState();
			E.twolevels();

			out_file << C.time*DeltaL/std::pow(upsilon, 0.5) << "		" << C.Mval/M00 << "		" << (C.Zval-Szc)/(Z00-Szc) << "		" << norm((E.gs.t()*C.gs)) << "		" << norm((E.fes.t()*C.gs)) << "		" << C.Mval/M0g << "	" << C.Mval << '\n';
			// M0g instead of M00
			//std::cout << C.time*DeltaL/std::pow(upsilon, 0.5) << "	" << sum << '\n';
		}
}

*/

return(0.);
}
