#ifndef __QSYSTEM_H
#define __QSYSTEM_H

#include <armadillo>

using namespace arma;
typedef std::complex<double> complex;

/*#####################################################*/

const double YMU = 1., ZETA = 1., MUC = -2.;
double delta(1.), mu(MUC), dgamma(0.);

/*#####################################################*/

class qsystem {
public:
	int L;		//! 0=open bc   1=abc   !!!-1=pbc
	//cx_dvec gr_state;
	cx_dmat ham;
	cx_dmat corr; //cx_dmat temp;
	dvec omega, fdist;
	complex energy;
	//dvec obs;
	double PBC, time, Temper;

	qsystem(int _L, double _pbc, double _Temper) {
		L = _L;
		PBC = _pbc;
		Temper = _Temper;
		//gr_state.set_size(2*L);
		ham.set_size(2*L, 2*L);
		corr.set_size(2*L, 2*L);
		//temp.set_size(2*L, 2*L);
		omega.set_size(L);
		fdist.set_size(L);
		energy = 0.;
		//obs.set_size(3);
		time = 0.;
	}

	void hamBuild(void);
	void spectrum(void);
	void corrMatrix(void);
	void Quench_Temper(const double& _Temper, const double& _mu);
	void pureThermTimeEvol(double ttime, cx_dmat& _corr);
	void genCorrMatrix(double ttime);
	void Measure(const int& xx,const int& yy, cx_dmat& _corr, dvec& obs);
	void derstep(const cx_dmat& _corr, cx_dmat& _k);
	void RKmethod(double _dt);
};

/*#####################################################*/

#endif














