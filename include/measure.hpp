#ifndef __MEASURE_H
#define __MEASURE_H

#include "qsystem.hpp"

dvec qsystem::Measure(const int& xx,const int& yy){
	/*
	 * Measure of the correlation functions C and P in the point
	 * with coordinate (_xx, _yy);
	*/
	int _xx = xx, _yy = yy;

	dvec obss(2);
	complex obs(0.);

	for(int m=0; m<L; ++m){
			obs = obs + conj(ham(m,_xx)) * ham(m,_yy);
			for(int q=0; q<L;++q)
			obs = obs + ( conj(ham(m,_xx+L))* ham(q,_yy+L) - conj(ham(m,_xx))*ham(q,_yy) )*corr(m,q) + conj(ham(m,_xx+L))*ham(q,_yy)*corr(m+L,q+L) + conj(ham(m,_xx))*ham(q,_yy+L)*conj(corr(q+L,m+L));
		}
	obss(0) = 2.*real(obs);
	obs = 0.;

	_xx = xx + L; _yy = yy + L;

	for(int m=0; m<L; ++m){
			obs = obs + conj(ham(m,_xx-L)) * conj(ham(m,_yy));
			for(int q=0; q<L;++q)
			obs = obs + ( conj(ham(m,_xx)) * conj(ham(q,_yy-L)) - conj(ham(m,_xx-L)) * conj(ham(q,_yy)) ) * corr(m,q) + conj(ham(m,_xx)) * conj(ham(q,_yy))*corr(m+L,q+L) + conj(ham(m,_xx-L))*conj(ham(q,_yy-L))*conj(corr(q+L,m+L));
		}
	obss(1) = 2. * real(obs);

	return obss;
}




#endif
