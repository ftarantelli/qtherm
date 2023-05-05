#include "include/all.h"


int main(int argc, char *argv[]){

	int Lsize(30);
	double Temperature(1.*std::pow(Lsize, -ZETA)), Tempb(0.);
	double pbc(0.), tmax(5.), dtime(0.02);

	int point(2), range(1), events;
	double kappai(0.), kappa(1.), eta(0.), etab(1.), sgamma(0.);

	std::string lineterm = "";
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
					else			dtime = atof( &argv[1][1] );
				break;
			case 'g':
					dgamma = atof( &argv[1][1] );
				break;
			case 'T':
					if(argv[1][1] == 'b')	Tempb = atof( &argv[1][2] );
					else 				Temperature = atof( &argv[1][1]);
				break;
			case 'e':
					if(argv[1][1] == 'b')	etab = atof( &argv[1][2] );
					else					eta = atof( &argv[1][1] );
				break;
			case 'r':
					range = atoi( &argv[1][1] );
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
			case 's':
					sgamma = atof( &argv[1][1] );
				break;
			default:
				std::cerr << "Unlucky: Retry input values\n";
				std::cerr << " L - size \n m - mu\n t - delta time\n tM - time max\n g - gamma\n T - Temperature\n Tb - Temp_bath\n e - eta = TL\n eb - eta_bath\n r - power of size time range\n p - PBC\n ki - kappai\n k - kappa\n s - scaling gamma\n";
				exit (8);
		}
		++argv;
		--argc;
	}

	//###data###data###data###data###data###data###data###data###data##
	char output[80];
	//temp.dat"); //

	if(sgamma != 0.)
		std::sprintf(output, "data/eqk%.0fl%d.dat", kappai, Lsize);
	else
		std::sprintf(output, "data/eqk%.0fl%d.dat", kappai, Lsize);

	//std::sprintf(output, "test.dat");
	FILE * test;
	test = fopen(output, "wt");

	//std::ofstream out_file(output, std::ios::out | std::ios::trunc);
	//out_file.precision(12);
	//std::cout.precision(16);

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
	tmax = tmax*double(std::pow(Lsize, range));
	events =  int(tmax/dtime);
	//point = int(Lsize/point);
	if(eta != 0.) Temperature = eta * std::pow(Lsize, -ZETA);
	//##critical#value##of##the#density#operator###########
	mu = MUC;
	class qsystem Eq(Lsize, pbc, 0.);
	Eq.hamBuild();
	Eq.spectrum();
	Eq.corrMatrix();
	double dens(0.), etav(eta);
	for(int j=0; j<Lsize; ++j) dens += real(Eq.corr(j,j));
	std::cout << dens/double(Lsize) << '\n';
	//##critical#value##of##the#density#operator###########
	mu = kappai * std::pow(Lsize, -YMU) + MUC;
	//###initial##condition###initial##condition###initial##condition##

	int xx = int(Lsize/2 - Lsize/2/point) - 1;
	int yy = int(Lsize/2 + Lsize/2/point) - 1;

	//###Parameters###Parameters###Parameters###Parameters##Parameters#
	//out_file << "#		T = " << Temperature << "		mu = " << mu << "		ypoint = " << point << "		dtime = " << dtime << "		L = " << Lsize << "		comm_line = " << lineterm <<'\n';
	//out_file << "# t/L		LC(x,y)		LP(x,y)		D-D(0)" << "\n";

	fprintf(test, "# %s# T = %.9f		mu = %.9f		ypoint = %d		dtime = %.9f	L = %d		comm_line = %s\n# eta		LC(x,y)		LP(x,y)		D-D(0)\n", std::ctime(&start_time), Temperature, mu, point, dtime, Lsize, lineterm.c_str());


	fprintf(test, "%.9f	%.9f	%.9f	%.9f\n", 0., double(Lsize)*real(Eq.corr(xx, yy) + Eq.corr(yy, xx)), 2.*double(Lsize)*real(Eq.corr(xx+Lsize, yy+Lsize)), 0.);

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

	class qsystem KC(Lsize, pbc, Temperature);

	//cx_dmat _corr_(2*Lsize, 2*Lsize);
	//dvec obs(3);

	KC.hamBuild();
	KC.spectrum();

	//#pragma omp parallel
	//#pragma omp parallel for
	for(int i=0; i<events; ++i){

		if ( i % (events/10) == 0 ){
			//out_file << std::flush;
			std::cout << int(100. / events * i) << "%	" << lineterm << "	" << output << "\n";
		}

		etav = eta*double(i+1)/double(events);
		KC.Temper = etav * std::pow(Lsize, -ZETA);

		KC.corrMatrix();
		double dens0(0.);
		for(int j=0; j<Lsize; ++j) dens0 += real(KC.corr(j,j));
		//KC.Measure(xx,yy, _corr_, obs);

		fprintf(test, "%.9f	%.9f	%.9f	%.9f\n", etav, double(Lsize)*real(KC.corr(xx, yy) + KC.corr(yy, xx)), 2.*double(Lsize)*real(KC.corr(xx+Lsize, yy+Lsize)), dens0 -dens);

	}
	//###time###Evolution###time###Evolution###time###Evolution##time##


	//###chrono###chrono###chrono###chrono###chrono###chrono###chrono##
	auto end = std::chrono::system_clock::now();
	std::time_t end_time = std::chrono::system_clock::to_time_t(end);
	run_log << output << " END with " << lineterm << "   " << std::ctime(&end_time) << std::flush;
	//###chrono###chrono###chrono###chrono###chrono###chrono###chrono##

	fclose(test);
	return(0);
}
