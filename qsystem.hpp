#include <cstring>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <complex>
#include <math.h>

#include <armadillo>

using namespace arma;
typedef std::complex<double> complex;

const double YMU = 1., ZETA = 1., MUC = -2.;
double delta(1.), mu(MUC);

/*#####################################################*/

class qsystem {
public:
	int L;		//! 0=open bc   1=abc   !!!-1=pbc
	//cx_dvec gr_state;
	cx_dmat ham;
	cx_dmat corr;
	dvec omega, fdist;
	complex energy;
	double PBC, time, Temper;

	qsystem(int _L, double _pbc, double _Temper) {
		L = _L;
		PBC = _pbc;
		Temper = _Temper;
		//gr_state.set_size(2*L);
		ham.set_size(2*L, 2*L);
		corr.set_size(2*L, 2*L);
		omega.set_size(L);
		fdist.set_size(L);
		time = 0.;
	}

	void hamBuild(void);
	void spectrum(void);
	void corrMatrix(void);
	void derstep(const cx_dmat& _corr, cx_dmat& _k);
	void RKmethod(double _dt);
};

/*#####################################################*/

void qsystem::hamBuild(void){

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
			omega(i) = eigval(i+L);
		//gr_state = gr_state / norm(gr_state);
		ham = eigvec; //trans(eigvec);
		//std::cout  << eigvec(1,3) <<"8gsgsgs\n" ;

}

void qsystem::corrMatrix(void){
	if (Temper == 0.)
		for(int k=0; k<L; ++k)
			fdist(k) = 0.;
	else
		for(int k=0; k<L; ++k)
			fdist(k) = 1./(1. + exp(omega(k)/Temper));

	cx_dmat temp(2*L,2*L);
	for(int i=0; i<2*L; ++i)
		for(int j=0; j<2*L; ++j)
			temp(i, j) = conj(ham(j,i));

	// correlations c^\dagger_j c_i
	for(int i=0; i<L; ++i)
		for(int j=i; j<L; ++j){
			for(int m=0; m<L; ++m)
				corr(j,i) = corr(j,i) + conj(temp(m,j+L)) * fdist(m) * temp(m,i+L) - conj(temp(m,j)) * fdist(m) * temp(m,i) + conj(temp(m,j)) * temp(m,i);
			if(i != j) corr(i, j) = conj(corr(j,i));
		}

	for(int i=L; i<2*L; ++i)
		for(int j=i; j<2*L; ++j){
			for(int m=0; m<L; ++m)
				corr(j,i) = corr(j,i) + conj(temp(m,j)) * fdist(m) * conj(temp(m,i-L)) - conj(temp(m,j-L)) * fdist(m) * conj(temp(m,i)) + conj(temp(m,j-L)) * conj(temp(m,i));
			if(i == j) corr(i, i) = 0.;
			corr(i,j) = - corr(j,i);
		}
	//complex obs(0.);
	//for(int i=0; i<L; ++i) obs += corr(i,i);
	//std::cout << corr(6 + L, 23+L) << "	" << obs << '\n';
}

void qsystem::derstep(const cx_dmat& _corr, cx_dmat& _k){

	for(int j=0; j<L; ++j) {
		for(int m=j; m<L; ++m){
			_k(j,m) = 0.;
			_k(j +L,m +L) = 0.;


			//j=pump   and   m=pump
			//k(j,m) = k(j,m) - u * _corr(j,m)
			//if(j==m) k(j,m) = k(j,m) + u
			//k(j +L,m +L) = k(j +L,m +L) - u * _corr(j +L,m +L)

			//j=loss   and   m=loss
			//k(j,m) = k(j,m) - u * _corr(j,m)
			//k(j +L,m +L) = k(j +L,m +L) - u * _corr(j +L,m +L)

			//DISSIPATION AT THE SITE L/b --- all loss
			/*
			if( (j-1) % bb == 0 ) {
				k(j,m) = k(j,m) - u/2. * _corr(j,m)
				k(j +L,m +L) = k(j +L,m +L) - u/2. * _corr(j +L,m +L)
			}

         if( mod(m-1,bb)==0 ) then
            k(j,m) = k(j,m) - u/2. * _corr(j,m)
            k(j +L,m +L) = k(j +L,m +L) - u/2. * _corr(j +L,m +L)
         endif
			*/

		if (j < L-1) {
			_k(j,m) = _k(j,m) - complex(0.,1.) * _corr(j+1,m) + complex(0.,1.) * delta * conj(_corr(m +L,j+1 +L));
			_k(j +L,m +L) = _k(j +L,m +L) - complex(0.,1.) * _corr(j+1 +L,m +L) - complex(0.,1.) * delta * _corr(m,j+1);
			if ( m == (j+1) )
				_k(j +L,m +L) = _k(j +L,m +L) + complex(0.,1.) * delta * 1.;
		} else {
			if (abs(PBC) == 1){
				_k(j,m) = _k(j,m) + PBC*complex(0.,1.) * _corr(0,m) - PBC*complex(0.,1.) * delta * conj(_corr(m +L,L));
				_k(j +L,m +L) = _k(j +L,m +L) + PBC*complex(0.,1.)*_corr(L,m +L) + PBC*complex(0.,1.)* delta * _corr(m,0);
                if ( m == 0 )
					_k(j +L,m +L) = _k(j +L,m +L) - PBC*complex(0.,1.) * delta * 1.;
			}
		}



		if (m < L-1){

			_k(j,m) = _k(j,m) + complex(0.,1.) * _corr(j,m+1) - complex(0.,1.) * delta * _corr(j +L,m+1 +L);
			_k(j +L,m +L) = _k(j +L,m +L) - complex(0.,1.) * _corr(j +L,m+1 +L) + complex(0.,1.) * delta * _corr(j,m+1);

        } else{

			if (abs(PBC) == 1){
				_k(j,m) = _k(j,m) - PBC*complex(0.,1.) * _corr(j,0) + PBC*complex(0.,1.) * delta * _corr(j +L,L);
				_k(j +L,m +L) = _k(j +L,m +L) + PBC*complex(0.,1.)*_corr(j +L,L) - PBC*complex(0.,1.)* delta * _corr(j,0);
			}

		}



		if (j > 0){

			_k(j,m) = _k(j,m) - complex(0.,1.) * _corr(j-1,m) - complex(0.,1.) * delta * conj(_corr(m +L,j-1 +L));
			_k(j +L,m +L) = _k(j +L,m +L) - complex(0.,1.) * _corr(j-1 +L,m +L) + complex(0.,1.) * delta * _corr(m,j-1);
			if ( m==(j-1) ) _k(j +L,m +L) = _k(j +L,m +L) - complex(0.,1.) * delta * 1.;

		} else {

			if (abs(PBC)==1){
				_k(j,m) = _k(j,m) + PBC*complex(0.,1.) * _corr(L-1,m) + PBC*complex(0.,1.) * delta * conj(_corr(m +L,L +L-1));
				_k(j +L,m +L) = _k(j +L,m +L) + PBC*complex(0.,1.)*_corr(L +L-1,m +L) - PBC*complex(0.,1.)* delta * _corr(m,L-1);
				if ( m==(L-1) ) _k(j +L,m +L) = _k(j +L,m +L) + PBC*complex(0.,1.) * delta * 1.;
			}

		}



		if (m > 0){

            _k(j,m) = _k(j,m) + complex(0.,1.) * _corr(j,m-1) + complex(0.,1.) * delta * _corr(j +L,m-1 +L);
            _k(j +L,m +L) = _k(j +L,m +L) - complex(0.,1.) * _corr(j +L,m-1 +L) - complex(0.,1.) * delta * _corr(j,m-1);

		} else {

             if (abs(PBC)==1){
				_k(j,m) = _k(j,m) - PBC*complex(0.,1.) * _corr(j,L-1) - PBC*complex(0.,1.) * delta * _corr(j +L,L +L-1);
				_k(j +L,m +L) = _k(j +L,m +L) + PBC*complex(0.,1.)*_corr(j +L,L +L-1) + PBC*complex(0.,1.)*delta * _corr(j,L-1);
			}

		}


		_k(j +L,m +L) = _k(j +L,m +L) - 2. * complex(0.,1.) * mu * _corr(j +L,m +L);



		_k(m +L,j +L) = - _k(j +L,m +L);

		if (m == j)
			_k(j +L,m +L) = 0.;
		else
			_k(m,j) = conj(_k(j,m));
//DELETE THE DISSIPATION OR COMMENT IT
		}
	}

}

void qsystem::RKmethod(double _dt){

	cx_dmat k1(2*L,2*L), k2(2*L,2*L), k3(2*L,2*L), k4(2*L,2*L);

	derstep(corr, k1);					// time
	derstep(corr + k1*_dt/2., k2);		// time + _dt/2
	derstep(corr + k2*_dt/2., k3);	// time + _dt/2
	derstep(corr + k3*_dt, k4);			//  time + _dt

	corr = corr + _dt/6.*( k1 + 2.*k2 + 2.*k3 + k4 );
	time = time + _dt;
}



















