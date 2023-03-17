#include "qsystem.hpp"
#include <chrono>



int main(int argc, char *argv[]){

	int Lsize(30);
	double Temperature(1.*std::pow(Lsize, ZETA));
	double pbc(0.), tmax(5.*Lsize), dtime(0.02);

	int point(Lsize/2), events;
	double kappai(0.), kappa(1.), eta(0.);

	std::string lineterm= "";
	for (int i=1;i<argc;i++)
		lineterm.append(std::string(argv[i]).append(" "));

	while( argc > 1 ) {

		switch(argv[1][0]) {
			case 'L':
					Lsize = atoi( &argv[1][1] );
				break;
			case 'm':
					mu = atof( &argv[1][1] );
				break;
			case 't':
					if(argv[1][1] == 'M')	tmax = atof( &argv[1][2] );
					else					dtime = atof( &argv[1][1] );
				break;
			case 'T':
					if(argv[1][1] == 'z')	eta = atof( &argv[1][2] );
					else			Temperature = atof( &argv[1][1]);
				break;
			case 'e':
					events = atoi( &argv[1][1] );
				break;
			case 'p':
					pbc = atof( &argv[1][1] );
				break;
			case 'k':
					if(argv[1][1] == 'i')	kappai= atof( &argv[1][2] );
					else					kappa = atof( &argv[1][1] );
				break;
			default:
				std::cerr << "Unlucky: Retry input values\n";
				std::cerr << " L - size \n m - mu\n t - delta time\n tM - time max\n T - Temperature\n Tz - eta\n e - num events\n p - PBC\n ki - kappai\n k - kappa\n";
				exit (8);
		}
		++argv;
		--argc;
	}

	//###data###data###data###data###data###data###data###data###data##
	char output[80];
	std::sprintf(output, "data/temp.dat");//qthk%.0fq%.0fe%.0fl%d.dat", kappai, kappa, eta, Lsize);
	std::ofstream out_file(output, std::ios::out | std::ios::trunc);
	out_file.precision(12);
	std::cout.precision(16);
	//###data###data###data###data###data###data###data###data###data##

	//###chrono###chrono###chrono###chrono###chrono###chrono###chrono##
	char runlog[15];
	std::sprintf(runlog, "run.log");
	std::ofstream run_log(runlog, std::ios::out | std::ios::app);
	auto start = std::chrono::system_clock::now();
	std::time_t start_time = std::chrono::system_clock::to_time_t(start);
	run_log << output << " with INPUT: " << lineterm << "	at	" << std::ctime(&start_time) << std::flush;
	//###chrono###chrono###chrono###chrono###chrono###chrono###chrono##

	//###initial##condition###initial##condition###initial##condition##
	tmax = 5.*Lsize;
	events =  int(tmax/dtime);
	point = int(Lsize/2);
	mu = kappai * std::pow(Lsize, -YMU) + MUC;
	if(eta != 0.) Temperature = eta * std::pow(Lsize, -ZETA);
	//###initial##condition###initial##condition###initial##condition##

	class qsystem KC(Lsize, pbc, Temperature);

	KC.hamBuild();
	KC.spectrum();
	KC.corrMatrix();

	//###Parameters###Parameters###Parameters###Parameters##Parameters#
	out_file << "#		T = " << Temperature << "		mu = " << mu << "		ypoint = " << point << "		dtime = " << dtime << "		L = " << Lsize << "		comm_line = " << lineterm <<'\n';
	out_file << "# t		LC(x,y)		P(x,y)		D(y)" << "\n";
	//###Parameters###Parameters###Parameters###Parameters##Parameters#

	//###quench###quench###quench###quench###quench###quench###quench##
	mu = kappa * std::pow(Lsize, -YMU) + MUC;
	//###quench###quench###quench###quench###quench###quench###quench##

	// * EQUILIBRIUM SCALING
	for(int hindex=0; hindex<50; ++hindex){
		mu = (kappa + hindex*0.1 )/ std::pow(Lsize, YMU) + MUC;

		class qsystem KC(Lsize, pbc, Temperature);

		KC.hamBuild();
		KC.spectrum();
		KC.corrMatrix();

		double dens(0.);
		for(int j=0; j<Lsize; ++j) dens += real(KC.corr(j,j));

		std::cout << Lsize << "	" << mu << "	" << eta << "	" << real(double(Lsize)*(KC.corr(0, point) + KC.corr(point, 0))) << "		" <<  2.*real(KC.corr(point+Lsize, 0+Lsize)) << "		" << 1-2*dens/double(Lsize) << '\n';
	}


/*
	//###time###Evolution###time###Evolution###time###Evolution##time##
	for(int i=0; i<events; ++i){

		if ( i % (events/20) == 0)
			out_file << std::flush;

		KC.RKmethod(dtime);

		double dens(0.);
		for(int j=0; j<Lsize; ++j) dens += real(KC.corr(j,j));

		out_file << KC.time << "		" << double(Lsize)*real(KC.corr(0, point) + KC.corr(point, 0)) << "		" <<  2.*real(KC.corr(point+Lsize, 0+Lsize)) << "		" << dens << '\n';

	}

	//###time###Evolution###time###Evolution###time###Evolution##time##
*/



	//###chrono###chrono###chrono###chrono###chrono###chrono###chrono##
	auto end = std::chrono::system_clock::now();
	std::time_t end_time = std::chrono::system_clock::to_time_t(end);
	run_log << output << " END with " << lineterm << "   " << std::ctime(&end_time) << std::flush;
	//###chrono###chrono###chrono###chrono###chrono###chrono###chrono##

	return(0);
}
