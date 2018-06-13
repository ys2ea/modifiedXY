#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <sstream>
#include <fstream>

const double PI = 3.14159;

class XYSimulation {
private:
	bool Istherm;                                              //whether it's thermalized
	double T_;
	double Beta_;
	double Delta_;
	int L_;
	int N_;
	double* Spinconfig_;                                         //Spin configuration L x L array ,  y*L + x, weird(spiral) boundary condition;
																 //spin 0-2*pi
	double Energy_;
	
public:
	XYSimulation(double Delta, double T, int L);
	void Reset();
	double Show_energy();
	double Show_spin(int, int);
	void thermalize();
	
	void Spinflip();
};

XYSimulation::XYSimulation(double Delta, double Beta, int L)
:	Istherm(0)
,	Delta_(Delta)
,	Beta_(Beta)
,	L_(L)
{
	T_ = 1./Beta;
	N_ = L*L;
	Spinconfig_ = (double*)malloc(sizeof(double)*L*L);
	
	for(int i=0; i<L*L; i++) 
		Spinconfig_[i] = 2*PI*(double)rand()/(double)RAND_MAX;	  
	
	Reset();
	
}
void XYSimulation::Spinflip() {
	int current = L_*(rand()%L_) + rand()%L_;
	double oldspin = Spinconfig_[current];
	double spinchange = 2*PI*(double)rand()/(double)RAND_MAX;
	double Echange = 0.0;
	int neighbor;

	if((neighbor = current + 1) > N_-1) neighbor -= N_;
	Echange += (-cos(oldspin+spinchange-Spinconfig_[neighbor]) + cos(oldspin-Spinconfig_[neighbor]))*(1-Delta_)
		+Delta_*(-cos(2*oldspin+2*spinchange-2*Spinconfig_[neighbor]) + cos(2*oldspin-2*Spinconfig_[neighbor]));

	if((neighbor = current - 1) < 0) neighbor += N_;
	Echange += (-cos(oldspin+spinchange-Spinconfig_[neighbor]) + cos(oldspin-Spinconfig_[neighbor]))*(1-Delta_)
		+Delta_*(-cos(2*oldspin+2*spinchange-2*Spinconfig_[neighbor]) + cos(2*oldspin-2*Spinconfig_[neighbor]));

	if((neighbor = current + L_) > N_-1) neighbor -= N_;
	Echange += (-cos(oldspin+spinchange-Spinconfig_[neighbor]) + cos(oldspin-Spinconfig_[neighbor]))*(1-Delta_)
		+Delta_*(-cos(2*oldspin+2*spinchange-2*Spinconfig_[neighbor]) + cos(2*oldspin-2*Spinconfig_[neighbor]));

	if((neighbor = current - L_) < 0) neighbor += N_;	
	Echange += (-cos(oldspin+spinchange-Spinconfig_[neighbor]) + cos(oldspin-Spinconfig_[neighbor]))*(1-Delta_)
		+Delta_*(-cos(2*oldspin+2*spinchange-2*Spinconfig_[neighbor]) + cos(2*oldspin-2*Spinconfig_[neighbor]));

	double r = (double)rand()/(double)RAND_MAX;

	if(Echange < 0.0 || exp(-Beta_*Echange) > r) {
		Spinconfig_[current] = oldspin + spinchange > 2*PI ? oldspin + spinchange - 2*PI : oldspin + spinchange;
		Energy_ += Echange;  }
}

//calculate the energy
void XYSimulation::Reset() {
	double e = 0;
	double m = 0;
	int ns;
	for( int i = 0; i < N_; i++) {
		ns = i+1;
		if(ns > N_ -1) ns -= N_;
		e += -cos(Spinconfig_[i]-Spinconfig_[ns])*(1-Delta_)+Delta_*cos(2*Spinconfig_[i]-2*Spinconfig_[ns]);
		ns = i+L_;
		if(ns > N_-1) ns -= N_;
		e += -cos(Spinconfig_[i]-Spinconfig_[ns])*(1-Delta_)+Delta_*cos(2*Spinconfig_[i]-2*Spinconfig_[ns]);
			
	}
	Energy_ = e;
}

double XYSimulation::Show_energy() {
	return Energy_; }

void XYSimulation::thermalize() {
	Istherm = true;
}

double XYSimulation::Show_spin(int x, int y) {
	return Spinconfig_[x+y*L_];
}

