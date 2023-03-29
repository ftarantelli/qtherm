#include "include/all.h"



int main(int argc, char *argv[]){

	int Lsize(30);
	double Temperature(1.*std::pow(Lsize, ZETA));
	double pbc(0.), tmax(5.*Lsize), dtime(0.02);

	int point(2), events;
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
			case 'g':
					dgamma = atof( &argv[1][1]);
				break;
			case 'T':
					if(argv[1][1] == 'z')	eta = atof( &argv[1][2] );
					else			Temperature = atof( &argv[1][1]);
				break;
			case 'e':
					eta = atof( &argv[1][1] );
				break;
			case 'p':
					pbc = atof( &argv[1][1] );
				break;
			case 'y':
					point = atoi( &argv[1][1] );
				break;
			case 'k':
					if(argv[1][1] == 'i')	kappai= atof( &argv[1][2] );
					else					kappa = atof( &argv[1][1] );
				break;
			default:
				std::cerr << "Unlucky: Retry input values\n";
				std::cerr << " L - size \n m - mu\n t - delta time\n tM - time max\n g - gamma\n T - Temperature\n Tz - eta\n e - num events\n p - PBC\n ki - kappai\n k - kappa\n";
				exit (8);
		}
		++argv;
		--argc;
	}

	//###data###data###data###data###data###data###data###data###data##
	char output[80];
	std::sprintf(output, "data/temp.dat");//"data/qthk%.0fq%.0fe%.0fg%.0fl%d.dat", kappai, kappa, eta, dgamma, Lsize);
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
	//point = int(Lsize/point);
	mu = kappai * std::pow(Lsize, -YMU) + MUC;
	if(eta != 0.) Temperature = eta * std::pow(Lsize, -ZETA);
	//###initial##condition###initial##condition###initial##condition##

	class qsystem KC(Lsize, pbc, Temperature);

	KC.hamBuild();
	KC.spectrum();
	KC.corrMatrix();

	//###Parameters###Parameters###Parameters###Parameters##Parameters#
	out_file << "#		T = " << Temperature << "		mu = " << mu << "		ypoint = " << point << "		dtime = " << dtime << "		L = " << Lsize << "		comm_line = " << lineterm <<'\n';
	out_file << "# t/L		LC(x,y)		LP(x,y)		D" << "\n";
	//###Parameters###Parameters###Parameters###Parameters##Parameters#

	//###quench###quench###quench###quench###quench###quench###quench##
	mu = kappa * std::pow(Lsize, -YMU) + MUC;
	//###quench###quench###quench###quench###quench###quench###quench##

/*
	// * EQUILIBRIUM SCALING
	//##equilibrium##equilibrium##equilibrium##equilibrium#equilibrium#
	std::cout << "#		T = " << Temperature << "		mu = " << mu << "		ypoint = " << point << "		dtime = " << dtime << "		L = " << Lsize << "		comm_line = " << lineterm <<'\n';
	std::cout << "# L	k	eta	LC(x,y)		LP(x,y)		1-2D/L" << "\n";
	for(int hindex=0; hindex<100; ++hindex){

		if ( hindex % (100/20) == 0)
			std::cout << std::flush;
		mu = (kappa + hindex*0.2 )/ std::pow(Lsize, YMU) + MUC;

		class qsystem KC(Lsize, pbc, Temperature);

		KC.hamBuild();
		KC.spectrum();
		KC.corrMatrix();

		double dens(0.);
		for(int j=0; j<Lsize; ++j) dens += real(KC.corr(j,j));

		std::cout << Lsize << "	" << kappa + hindex*0.2 << "	" << eta << "	" << real(double(Lsize)*(KC.corr(0, point) + KC.corr(point, 0))) << "		" << double(Lsize)*2.*real(KC.corr(point+Lsize, 0+Lsize)) << "		" << 1-2*dens/double(Lsize) << '\n';
	}
	//##equilibrium##equilibrium##equilibrium##equilibrium#equilibrium#
*/

	//###time###Evolution###time###Evolution###time###Evolution##time##
	int xx = int(Lsize/2 - Lsize/2/point) - 1;
	int yy = int(Lsize/2 + Lsize/2/point) - 1;

	out_file << KC.time/double(Lsize) << "		" << double(Lsize)*real(KC.corr(xx, yy) + KC.corr(yy, xx)) << "		" <<  2.*double(Lsize)*real(KC.corr(yy+Lsize, xx+Lsize)) << '\n';
	/* CHECK
	std::cout << KC.corr(3,4) << "   " << KC.fdist(Lsize-1) <<'\n';
	mu = kappai * std::pow(Lsize, -YMU) + MUC;
	KC.Quench_Temper(Temperature, mu);
	std::cout << corr0(3,4) << "   " << corr0(Lsize-1, Lsize-1) << '\n';
	KC.genCorrMatrix(0);
	//KC.corr = corr0; KC.genCorrMatrix(dtime);
	std::cout << KC.corr(3,4) << '\n';
	exit(8);
	*/

	double Temper2 = 1. * std::pow(Lsize, -ZETA);
	KC.Quench_Temper(0., mu);
	for(int i=0; i<events; ++i){

		if ( i % (events/20) == 0 ){
			out_file << std::flush;
			std::cout << int(100. / events * i) << "%\n";
		}

		//KC.RKmethod(dtime);
		//if ( i % (events/80) == 0 ){
			//KC.genCorrMatrix(double(i+1)*dtime);
			KC.pureThermTimeEvol(double(i+1)*dtime);
			out_file << KC.time/double(Lsize) << "		" << double(Lsize)*KC.Measure(xx,yy)(0) << "		" << double(Lsize)*KC.Measure(xx,yy)(1) << '\n';

			/*
			double dens(0.);
			for(int j=0; j<Lsize; ++j) dens += real(KC.corr(j,j));

			out_file << KC.time/double(Lsize) << "		" << double(Lsize)*real(KC.corr(xx, yy) + KC.corr(yy, xx)) << "		" <<  2.*double(Lsize)*real(KC.corr(yy+Lsize, xx+Lsize)) << "		" << dens << '\n';
			*/
		//}

	}
	//###time###Evolution###time###Evolution###time###Evolution##time##


	//###chrono###chrono###chrono###chrono###chrono###chrono###chrono##
	auto end = std::chrono::system_clock::now();
	std::time_t end_time = std::chrono::system_clock::to_time_t(end);
	run_log << output << " END with " << lineterm << "   " << std::ctime(&end_time) << std::flush;
	//###chrono###chrono###chrono###chrono###chrono###chrono###chrono##

	return(0);
}
