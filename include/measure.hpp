#ifndef __MEASURE_H
#define __MEASURE_H

#include "qsystem.hpp"

void qsystem::Measure(const int& xx,const int& yy, cx_dmat& _corr, dvec& obs){
	/*
	 * Measure of the correlation functions C and P in the point
	 * with coordinate (_xx, _yy);
	*/
	int _xx = xx, _yy = yy;

	complex obss(0.);

	for(int m=0; m<L; ++m){
			obss = obss + conj(ham(m,_xx)) * ham(m,_yy);
			for(int q=0; q<L;++q)
			obss = obss + ( conj(ham(m,_xx+L))* ham(q,_yy+L) - conj(ham(m,_xx))*ham(q,_yy) )*_corr(m,q) + conj(ham(m,_xx+L))*ham(q,_yy)*_corr(m+L,q+L) + conj(ham(m,_xx))*ham(q,_yy+L)*conj(_corr(q+L,m+L));
		}
	obs(0) = 2.*real(obss);
	obss = 0.;

	_xx = xx + L; _yy = yy + L;

	for(int m=0; m<L; ++m){
			obss = obss + conj(ham(m,_xx-L)) * conj(ham(m,_yy));
			for(int q=0; q<L;++q)
			obss = obss + ( conj(ham(m,_xx)) * conj(ham(q,_yy-L)) - conj(ham(m,_xx-L)) * conj(ham(q,_yy)) ) * _corr(m,q) + conj(ham(m,_xx)) * conj(ham(q,_yy))*_corr(m+L,q+L) + conj(ham(m,_xx-L))*conj(ham(q,_yy-L))*conj(_corr(q+L,m+L));
		}
	obs(1) = 2. * real(obss);
	obss = 0.;

	for(int _xsite=0; _xsite<L; ++_xsite){
		for(int m=0; m<L; ++m){
			obss = obss + conj(ham(m,_xsite)) * ham(m,_xsite);
			for(int q=0; q<L;++q)
			obss = obss + ( conj(ham(m,_xsite+L))* ham(q,_xsite+L) - conj(ham(m,_xsite))*ham(q,_xsite) )*_corr(m,q) + conj(ham(m,_xsite+L))*ham(q,_xsite)*_corr(m+L,q+L) + conj(ham(m,_xsite))*ham(q,_xsite+L)*conj(_corr(q+L,m+L));
		}
	}

	obs(2) = real(obss);

	//return obss;
}




#endif
