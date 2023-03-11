#include "qsystem.hpp"
#include <chrono>



int main(int argc, char *argv[]){

	int Lsize(30);
	double Temperature(1.*std::pow(Lsize, 1.));
	double pbc(0.), tmax(10.), dtime(0.01);

	int point(Lsize/2), events;

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
					Temperature = atof( &argv[1][1]);
				break;
			case 'e':
					events = atoi( &argv[1][1] );
				break;
			case 'p':
					pbc = atof( &argv[1][1] );
				break;
			default:
				std::cerr << "Unlucky: Retry input values\n";
				std::cerr << " L - size \n m - mu\n t - delta time\n tM - time max\n T - Temperature\n e - num events\n p - PBC\n";
				exit (8);
		}
		++argv;
		--argc;
	}

	//###data###data###data###data###data###data###data###data###data##
	char output[80];
	std::sprintf(output, "data/temp.dat");
	//std::sprintf(output, "data/Swireu%.0fk%.0ff%.0fl%d.dat", 100*upsilon, abs(kappa), abs(kappaf), L);
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


	class qsystem KC(Lsize, pbc, Temperature);

	KC.hamBuild();
	KC.spectrum();
	KC.corrMatrix();

	tmax = 5.*Lsize;
	events =  int(tmax/dtime);

	//###Parameters###Parameters###Parameters###Parameters##Parameters#
	out_file << "#		T = " << Temperature << "		mu = " << mu << "		ypoint = " << point << "		dtime = " << dtime << "		L = " << Lsize << "		comm_line = " << lineterm <<'\n';
	out_file << "# t		C(x,y)		P(x,y)		D(y)" << "\n";
	//###Parameters###Parameters###Parameters###Parameters##Parameters#

	mu = -1.;
	//###time###Evolution###time###Evolution###time###Evolution##time##
	for(int i=0; i<events; ++i){

		if ( i % (events/20) == 0)
			out_file << std::flush;

		KC.RKmethod(dtime);

		double dens(0.);
		for(int j=0; j<Lsize; ++j) dens += real(KC.corr(j,j));

		out_file << KC.time << "		" << real(KC.corr(0, point) + KC.corr(point, 0)) << "		" <<  2.*real(KC.corr(point+Lsize, 0+Lsize)) << "		" << dens << '\n';

	}
	//###time###Evolution###time###Evolution###time###Evolution##time##




	//###chrono###chrono###chrono###chrono###chrono###chrono###chrono##
	auto end = std::chrono::system_clock::now();
	std::time_t end_time = std::chrono::system_clock::to_time_t(end);
	run_log << output << " END with " << lineterm << "   " << std::ctime(&end_time) << std::flush;
	//###chrono###chrono###chrono###chrono###chrono###chrono###chrono##

	return(0);
}
