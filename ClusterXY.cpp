#include "ClusterXY.h"

using namespace std;

int main(int argc, char* argv[]) {
	srand(time(NULL));
	int L;
	int ntherm, nmeasure;
	double MinBeta, MaxBeta, DeltaBeta;
	double energy;
	
	if(argc != 7) {
		cout << "Wrong usage!" << endl; }
	
	sscanf(argv[1], "%d", &L);
	sscanf(argv[2], "%d", &ntherm);
	sscanf(argv[3], "%d", &nmeasure);
	sscanf(argv[4], "%lf", &MinBeta);
	sscanf(argv[5], "%lf", &MaxBeta);
	sscanf(argv[6], "%lf", &DeltaBeta);
	
	FILE* output;
	
	std::stringstream output_name;
	output_name << "XY" << "L=" << L << ".dat";	
	const string &temp=output_name.str();
	const char *name=temp.c_str();
	output = fopen(name, "w");
	
	for(double beta = MinBeta; beta < MaxBeta; beta += DeltaBeta) {
		
		energy = 0;
		
		XYSimulation sim(beta, L);
		
		for(int iter = 0; iter < ntherm; iter ++) {
			sim.Clusterflip();
			sim.Reset(); 
	/*		if(iter < 5) {
				for(int i = 0; i < L; i ++) {
					for(int j = 0; j < L; j ++) {
						printf("%lf  ",sim.Show_spin(j,i));
					}
					printf("\n");
				} 
				cout << "T = " << 1/beta << ", Energy = " << sim.Show_energy()/L/L << endl;
				printf("------------\n");} */
			
				
		}
		
		sim.thermalize();
		
		for(int iter = 0; iter < nmeasure; iter ++) {
			sim.Clusterflip();
			
			sim.Reset();
			
			energy += sim.Show_energy()/double(nmeasure);
			
			
		}
		
		cout << "T = " << 1/beta << ", Energy = " << energy/L/L << endl;
		
		fprintf(output, "%lf\t%lf\n", beta, energy/L/L);
		
	}
	fclose(output);
}

