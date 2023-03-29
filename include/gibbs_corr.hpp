#ifndef __GIBBS_CORR_H
#define __GIBBS_CORR_H

#include "qsystem.hpp"

void qsystem::corrMatrix(void){
	corr.zeros();
	if (Temper == 0.)
		for(int k=0; k<L; ++k)
			fdist(k) = 0.;
	else
		for(int k=0; k<L; ++k)
			fdist(k) = 1./(1. + exp(-omega(k)/Temper));

	//cx_dmat temp(2*L,2*L);
	////for(int i=0; i<2*L; ++i) for(int j=0; j<2*L; ++j) temp(i, j) = conj(ham(j,i));
	//temp = ham;
	////temp = conj(trans(ham));

	// correlations c^\dagger_j c_i
	for(int i=0; i<L; ++i)
		for(int j=i; j<L; ++j){
			for(int m=0; m<L; ++m)
				corr(j,i) = corr(j,i) + conj(ham(m,j+L)) * fdist(m) * ham(m,i+L) - conj(ham(m,j)) * fdist(m) * ham(m,i) + conj(ham(m,j)) * ham(m,i);
			if(i != j) corr(i, j) = conj(corr(j,i));
		}
	// correlations c^\dagger_j c^\dagger_i
	for(int i=L; i<2*L; ++i)
		for(int j=i; j<2*L; ++j){
			for(int m=0; m<L; ++m)
				corr(j,i) = corr(j,i) + conj(ham(m,j)) * fdist(m) * conj(ham(m,i-L)) - conj(ham(m,j-L)) * fdist(m) * conj(ham(m,i)) + conj(ham(m,j-L)) * conj(ham(m,i));
			if(i == j) corr(i, i) = 0.;
			corr(i,j) = - corr(j,i);
		}
	//ham = temp;



////////////////////////////////////////////////////////////
/*
	ham = temp;		//cx_dmat temp(2*L,2*L);
	temp.zeros();

	for(int i=0; i<L; ++i)
		for(int j=0; j<L; ++j){
			for(int m=0; m<L; ++m){
				temp(j,i) = temp(j,i) +  ham(j,m) * conj(ham(i,m));
				for(int q=0; q<L; ++q)
					temp(j,i) = temp(j,i) + ( conj(ham(j,m+L)) * ham(i,q+L) - ham(j,m) * conj(ham(i,q)) )*corr(m,q) + conj(ham(j,m+L))*conj(ham(i,q))*corr(m+L,q+L) +ham(j,m)*ham(i,q+L)*conj(corr(q+L,m+L));
			}
			//if(i != j) temp(i, j) = conj(temp(j,i));
		}

	for(int i=L; i<2*L; ++i)
		for(int j=L; j<2*L; ++j){
			for(int m=0; m<L; ++m){
				temp(j,i) = temp(j,i) + ham(j-L,m) * conj(ham(i-L,m+L));
				for(int q=0; q<L; ++q)
					temp(j,i) = temp(j,i) + ( conj(ham(j-L,m+L)) * ham(i-L,q) - ham(j-L,m) * conj(ham(i-L,q+L)) ) * corr(m,q) + conj(ham(j-L,m+L))*conj(ham(i-L,q+L))*corr(m+L,q+L) + ham(j-L,m)*ham(i-L,q)*conj(corr(q+L,m+L));
			}
			//if(i == j) temp(i, i) = 0.;
			//temp(i,j) = - temp(j,i);
		}

	//std::cout << real(temp) << "\n\n\n";
	//std::cout << real(ham) << "\n\n\n";
	std::cout << temp(4, 5) << "	" << temp(L-1,L-1) << "	" << temp(7+L, 10+L) << " " << fdist(L-1) << '\n';
	exit(8);
*/
}

/*
	// correlations c^\dagger_j c_i
	for(int i=0; i<L; ++i)
		for(int j=i; j<L; ++j){
			for(int m=0; m<L; ++m)
				corr(j,i) = corr(j,i) + conj(temp(m,j+L)) * fdist(m) * temp(m,i+L) - conj(temp(m,j)) * fdist(m) * temp(m,i) + conj(temp(m,j)) * temp(m,i);
			if(i != j) corr(i, j) = conj(corr(j,i));
		}
	// correlations c^\dagger_j c^\dagger_i
	for(int i=L; i<2*L; ++i)
		for(int j=i; j<2*L; ++j){
			for(int m=0; m<L; ++m)
				corr(j,i) = corr(j,i) + conj(temp(m,j)) * fdist(m) * conj(temp(m,i-L)) - conj(temp(m,j-L)) * fdist(m) * conj(temp(m,i)) + conj(temp(m,j-L)) * conj(temp(m,i));
			if(i == j) corr(i, i) = 0.;
			corr(i,j) = - corr(j,i);
		}
*/

#endif
