#include "modified.h"

using namespace std;

int main(int argc, char* argv[]) {
	int L;
	int ntherm, nmeasure, Corrlength;
	double D, MinBeta, MaxBeta, DeltaBeta;
	double energy, mag;

	
	if(argc != 9) {
		cout << "Wrong usage!" << endl; }
	
	sscanf(argv[1], "%d", &L);
	sscanf(argv[2], "%d", &ntherm);
	sscanf(argv[3], "%d", &nmeasure);
	sscanf(argv[4], "%d", &Corrlength);
	sscanf(argv[5], "%lf", &D);
	sscanf(argv[6], "%lf", &MinBeta);
	sscanf(argv[7], "%lf", &MaxBeta);
	sscanf(argv[8], "%lf", &DeltaBeta);
	

	
	srand(time(NULL));
	
	FILE* output;
	std::stringstream output_name;
	output_name << "XY" << "Delta" << D << "L=" << L << ".dat";	
	const string &temp=output_name.str();
	const char *name=temp.c_str();
	output = fopen(name, "w");
	
	for(double beta = MinBeta; beta < MaxBeta; beta += DeltaBeta) {
		
		energy = 0; mag = 0;
		
		XYSimulation sim(D, beta, L);
		
		for(int iter = 0; iter < L*L*ntherm; iter ++) {
			for(int i = 0; i < Corrlength; i ++)
				sim.Spinflip();  }
		
		sim.thermalize();
		
		for(int iter = 0; iter < L*L*nmeasure; iter ++) {
			for(int i = 0; i < Corrlength; i ++)
				sim.Spinflip();
			
			if(iter%200==0)
				sim.Reset();
			
			energy += sim.Show_energy()/nmeasure/L/L;
			
		}
		
		
		fprintf(output, "%lf\t%lf\n", beta, energy/L/L);
		fflush(output);
		
	}
	fclose(output);
}

