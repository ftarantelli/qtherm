#include "mpo.hpp"

using namespace itensor;



int main()
{
    // Numero di spin
    const int N = 10;
	double eta = 100.;

    // Costanti del modello di Ising
    const double J = 4.0*1.;
    const double h = 1.*J/2.;

    // Creazione degli spazi di Hilbert per ogni spin
    auto sites = SpinHalf(N);

    // Costruzione dell'Hamiltoniana di Ising con campo magnetico esterno
    auto ampo = AutoMPO(sites);
    for (int i = 1; i < N; i++)
    {
        ampo += -J, "Sz", i, "Sz", i + 1;
    }
    ampo += -h, "Sx", 1;
    ampo += -h, "Sx", N;
    auto H = toMPO(ampo);

    // Costruzione dello stato di Gibbs a temperatura T = 1
    auto temp = 1.0;
	if (eta != 0.) temp = eta / N;

    auto beta = 1.0 / temp;
    auto state = MPS(sites);
    auto noise = [](Real){ return 1.0; };
    randomize(state);
    auto psi = expHermitian(H, -beta, "Method=Exp", "Noise", noise);
    psi /= norm(psi);

    // Costruisce l'operatore di magnetizzazione lungo z
    auto sigmaz = op(sites,"Sz",1);
    for(auto j = 2; j <= N; ++j) {
        sigmaz *= op(sites,"Sz",j);
    }

    // Calcola il valore di aspettazione dell'operatore su psi
    auto Mz = std::real(inner(psi,sigmaz,psi));

	std::cout << "# Lsize		eta		Mz \n";
	std::cout << N << "	" << eta << "	" << Mz << '\n';
    return 0;
}
