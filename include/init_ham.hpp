#ifndef __INIT_HAM_H
#define __INIT_HAM_H

#include "qsystem.hpp"

void qsystem::hamBuild(void){
	ham.zeros();
	for(int i=0;i<L-1;++i){
		ham(i,i) = ham(i,i) -  mu;
		ham(i,i+1) = ham(i,i+1) - 1.;
		ham(i+1,i) = ham(i+1,i) - 1.;

		ham(i,i+L+1) = ham(i,i+L+1) + delta;
		ham(i+1,i+L) = ham(i+1,i+L) - delta;

		ham(i+L,i+1) = ham(i+L,i+1) - delta;
		ham(i+L+1,i) = ham(i+L+1,i) + delta;

		ham(i+L,i+L) = ham(i+L,i+L) +  mu;
		ham(i+L,i+L+1) = ham(i+L,i+L+1) + 1.;
		ham(i+L+1,i+L) = ham(i+1+L,i+L) + 1.;
	}

	ham(L-1,L-1) = ham(L-1,L-1) -  mu;
	ham(2*L-1,2*L-1) = ham(2*L-1,2*L-1) +  mu;

	if (abs(PBC) == 1){
		ham(L-1,0) = ham(L-1,0) + PBC*1.;
		ham(0,L-1) = ham(0,L-1) + PBC*1.;

		ham(L-1,L) = ham(L-1,L) - PBC*delta;
		ham(0,L+L-1) = ham(0,L+L-1) + PBC*delta;

		ham(L+L-1,0) = ham(L+L-1,0) + PBC*delta;
		ham(L,L-1) = ham(L,L-1) - PBC*delta;

		ham(L+L-1,L) = ham(L+L-1,L) - PBC*1.;
		ham(L,L+L-1) = ham(L,L+L-1) - PBC*1.;
	}

}

void qsystem::spectrum(void) {

		//sp_dmat HH = real(Ham);
		//eigs_opts opts;
		//opts.maxiter = 1000;
		dvec eigval;
		cx_dmat eigvec;

		eig_sym(eigval, eigvec, ham);
		//eigs_gen(eigval, eigvec, ham, 1, "sr");
		for(int i=0; i<L; ++i)
			omega(i) = eigval(i);
		//gr_state = gr_state / norm(gr_state);
		//ham = eigvec; //trans(eigvec);
		ham = conj(trans(eigvec));
		//std::cout  << eigvec(1,3) <<"8gsgsgs\n" ;

}

#endif
