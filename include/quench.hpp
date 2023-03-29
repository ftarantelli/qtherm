#ifndef __QUENCH_H
#define __QUENCH_H

#include "qsystem.hpp"
#include "init_ham.hpp"

cx_dmat corr0;

void qsystem::Quench_Temper(const double& _Temper, const double& _mu){
	/*
	 * ATTENTION: with this function, we are changing the
	 * initial coordinate correlations to the Bogoliubov base;
	 *
	 * corr0 corresponds with the INITIAL values of the
	 * Bogoliubov correlations;
	*/

	Temper = _Temper;
	mu = _mu;

	cx_dmat hham = ham;

	hamBuild();
	spectrum();

	//cx_dmat temp(2*L, 2*L);
	corr0.set_size(2*L,2*L);
	corr0.zeros();
	for(int i=0; i<L; ++i)
		for(int j=i; j<L; ++j){
			for(int m=0; m<L; ++m){
				corr0(j,i) = corr0(j,i) +  ham(j,m) * conj(ham(i,m));
				for(int q=0; q<L; ++q)
					corr0(j,i) = corr0(j,i) + ( conj(ham(j,m+L)) * ham(i,q+L) - ham(j,m) * conj(ham(i,q)) )*corr(m,q) + conj(ham(j,m+L))*conj(ham(i,q))*corr(m+L,q+L) +ham(j,m)*ham(i,q+L)*conj(corr(q+L,m+L));
			}
			if(i != j) corr0(i, j) = conj(corr0(j,i));
		}

	for(int i=L; i<2*L; ++i)
		for(int j=i; j<2*L; ++j){
			for(int m=0; m<L; ++m){
				corr0(j,i) = corr0(j,i) + ham(j-L,m) * conj(ham(i-L,m+L));
				for(int q=0; q<L; ++q)
					corr0(j,i) = corr0(j,i) + ( conj(ham(j-L,m+L)) * ham(i-L,q) - ham(j-L,m) * conj(ham(i-L,q+L)) ) * corr(m,q) + conj(ham(j-L,m+L))*conj(ham(i-L,q+L))*corr(m+L,q+L) + ham(j-L,m)*ham(i-L,q)*conj(corr(q+L,m+L));
			}
			if(i == j) corr0(i, i) = 0.;
			corr0(i,j) = - corr0(j,i);
		}

	//corr = temp;
	//corr0 = ( b'^\dagger b' )( b'^\dagger b'^\dagger );

	if (Temper == 0.)
		for(int k=0; k<L; ++k)
			fdist(k) = 0.;
	else
		for(int k=0; k<L; ++k)
			fdist(k) = 1./(1. + exp(-omega(k)/Temper));
	/*
	 * at the end, we change the ham, the spectrum, the fdist
	 * and the Bogoliubov correlations
	*/
}


void qsystem::pureThermTimeEvol(double ttime){
	/*
	 * compute the time evolution in the Bogoliubov base
	 * for pure dissipation dynamics with
	 * HOMOGENEOUS IDENTICAL THERMAL BATHS
	 * ATTENTION: before use the function Quench_Temper
	*/

	time = ttime;
	//ATTENTION: corr = Bogoliubov correlations;
	for(int kk=0; kk<L; ++kk){
		corr(kk,kk) = fdist(kk)*(1-exp(-2.*dgamma*ttime)) + corr0(kk,kk)*exp(-2.*dgamma*ttime);
		for(int qq=kk+1; qq<L; ++qq){
			corr(kk, qq) = corr0(kk,qq) * exp(complex(0.,1.)*(omega(qq)-omega(kk))*ttime - 2.*dgamma*ttime);
			corr(qq, kk) = conj(corr(kk,qq));
			corr(kk+L,qq+L) = corr0(kk+L,qq+L) * exp(complex(0.,-1.)*(omega(qq)+omega(kk))*ttime - 2.*dgamma*ttime );
			corr(qq+L,kk+L) = - corr(kk+L,qq+L);
		}
	}
}

void qsystem::genCorrMatrix(double ttime){
	/*
	 * this function transform the correlations from the
	 * Bogoliubov base to the coordinate base
	*/

	pureThermTimeEvol(ttime);
	cx_dmat temp(2*L,2*L);
	temp.zeros();
	// correlations c^\dagger_j c_i
	for(int i=0; i<L; ++i)
		for(int j=i; j<L; ++j){
			for(int m=0; m<L; ++m){
				temp(j,i) = temp(j,i) + conj(ham(m,j)) * ham(m,i);
				for(int q=0; q<L;++q)
				temp(j,i) = temp(j,i) + ( conj(ham(m,j+L))* ham(q,i+L) - conj(ham(m,j))*ham(q,i) )*corr(m,q) + conj(ham(m,j+L))*ham(q,i)*corr(m+L,q+L) + conj(ham(m,j))*ham(q,i+L)*conj(corr(q+L,m+L));
			}
			if(i != j) temp(i, j) = temp(j,i);
		}

	// correlations c^\dagger_j c^\dagger_i
	for(int i=L; i<2*L; ++i)
		for(int j=i; j<2*L; ++j){
			for(int m=0; m<L; ++m){
				temp(j,i) = temp(j,i) + conj(ham(m,j-L)) * conj(ham(m,i));
				for(int q=0; q<L;++q)
				temp(j,i) = temp(j,i) + ( conj(ham(m,j)) * conj(ham(q,i-L)) - conj(ham(m,j-L)) * conj(ham(q,i)) ) * corr(m,q) + conj(ham(m,j)) * conj(ham(q,i))*corr(m+L,q+L) + conj(ham(m,j-L))*conj(ham(q,i-L))*conj(corr(q+L,m+L));
			}
			if(i == j) temp(i, i) = 0.;
			temp(i,j) = - temp(j,i);
		}

	corr = temp;

	/*
	 * we erase the corr in the Bogoliubov base in favor of the
	 * corr in the coordinate space
	*/

}

#endif
