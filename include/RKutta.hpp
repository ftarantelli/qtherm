#ifndef __RKUTTA_H
#define __RKUTTA_H

#include "qsystem.hpp"

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


#endif
