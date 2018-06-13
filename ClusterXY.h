/*
 *  ClusterXY.h
 *  XY model with cluster update algorithm.
 *
 *  Created by Yifei Shi on 10/24/13.
 *  Copyright 2013 __MyCompanyName__. All rights reserved.
 *
 */

#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <sstream>
#include <fstream>
#include "MersenneTwister.h"

const double PI = 3.14159;

class XYSimulation {
private:
	bool Istherm;                                              //whether it's thermalized
	double T_;
	double Beta_;
	int L_;
	int N_;
	double* Spinconfig_;                                         //Spin configuration L x L array ,  y*L + x, periodic boundary condition;
																 //spin 0-2*pi
	double Energy_;
	MTRand ran;
	
public:
	XYSimulation(double T, int L);
	void Reset();
	double Show_energy();
	double Show_spin(int, int);
	void thermalize();
	
	void Clusterflip();
};

XYSimulation::XYSimulation(double Beta, int L)
:	Istherm(0)
,	Beta_(Beta)
,	L_(L)
{
	T_ = 1./Beta;
	N_ = L*L;
	Spinconfig_ = (double*)malloc(sizeof(double)*L*L);
	
	for(int i=0; i<L*L; i++) 
		Spinconfig_[i] = 0.0;	  
	
	Reset();
	
}


//Wolf cluster flip.
void XYSimulation::Clusterflip() {
	int stack[N_];
	bool Instack[N_];
	for(int ii = 0; ii < N_; ii++) Instack[ii] = 0;
	
	int i, nstack, neighbor, current; 
	double dir, Padd, oldspin;
	double sign;
	
	i = ran.randDblExc()*L_*L_;
	
	stack[0] = i;
	nstack = 1;
	dir = PI*ran.randDblExc();
	if(fabs(Spinconfig_[i] - dir) < PI/2. || fabs(Spinconfig_[i] - dir) > 3*PI/2.)
		sign = 1.;
	else sign = -1.;
	
	Instack[i] = true;
	oldspin = Spinconfig_[i];
	if((2*dir + PI - oldspin) < 0) Spinconfig_[i] = 2*dir + 3*PI - oldspin;
	else if((2*dir + PI - oldspin) > 2*PI) Spinconfig_[i] = 2*dir - PI - oldspin;
	else Spinconfig_[i] = 2*dir + PI - oldspin;
	double r;
	
	while(nstack) {
		current = stack[--nstack];
		neighbor = (current%L_ == L_-1) ? current - L_ + 1:  current + 1;         //PBC
		if((fabs(PI-fabs(dir-Spinconfig_[neighbor]))-PI/2.)*sign > 0 && (!Instack[neighbor])) {
			r = ran.randDblExc();
			oldspin = Spinconfig_[neighbor];
			Padd = 1. - exp(2.*Beta_*cos(dir-Spinconfig_[current])*cos(dir-oldspin));
			if(r < Padd) {
				stack[nstack++] = neighbor;
				
				//fit the anguler to 0~2pi
				if((2*dir + PI - oldspin) < 0) Spinconfig_[neighbor] = 2*dir + 3*PI - oldspin;
				else if((2*dir + PI - oldspin) > 2*PI) Spinconfig_[neighbor] = 2*dir - PI - oldspin;
				else Spinconfig_[neighbor] = 2*dir + PI - oldspin;
				Instack[neighbor] = true;
			}
		}
		
		neighbor = (current%L_ == 0) ? current + L_ -1:  current - 1;
		if((fabs(PI-fabs(dir-Spinconfig_[neighbor]))-PI/2.)*sign > 0 && (!Instack[neighbor])) {
			r = ran.randDblExc();
			oldspin = Spinconfig_[neighbor];
			Padd = 1. - exp(2.*Beta_*cos(dir-Spinconfig_[current])*cos(dir-oldspin));
			if(r < Padd) {
				stack[nstack++] = neighbor;
				if((2*dir + PI - oldspin) < 0) Spinconfig_[neighbor] = 2*dir + 3*PI - oldspin;
				else if((2*dir + PI - oldspin) > 2*PI) Spinconfig_[neighbor] = 2*dir - PI - oldspin;
				else Spinconfig_[neighbor] = 2*dir + PI - oldspin;
				Instack[neighbor] = true;
			}
		}
		
		neighbor = (current/L_ == L_-1) ? current%L_ :  current + L_;
		if((fabs(PI-fabs(dir-Spinconfig_[neighbor]))-PI/2.)*sign > 0 && (!Instack[neighbor])) {
			r = ran.randDblExc();
			oldspin = Spinconfig_[neighbor];
			Padd = 1. - exp(2.*Beta_*cos(dir-Spinconfig_[current])*cos(dir-oldspin));
			if(r < Padd) {
				stack[nstack++] = neighbor;
				if((2*dir + PI - oldspin) < 0) Spinconfig_[neighbor] = 2*dir + 3*PI - oldspin;
				else if((2*dir + PI - oldspin) > 2*PI) Spinconfig_[neighbor] = 2*dir - PI - oldspin;
				else Spinconfig_[neighbor] = 2*dir + PI - oldspin;
				Instack[neighbor] = true;
			}
		}
		
		neighbor = (current/L_ == 0) ? current%L_ + L_*(L_-1) :  current - L_;
		if((fabs(PI-fabs(dir-Spinconfig_[neighbor]))-PI/2.)*sign > 0 && (!Instack[neighbor])) {
			r = ran.randDblExc();
			oldspin = Spinconfig_[neighbor];
			Padd = 1. - exp(2.*Beta_*cos(dir-Spinconfig_[current])*cos(dir-oldspin));
			if(r < Padd) {
				stack[nstack++] = neighbor;
				if((2*dir + PI - oldspin) < 0) Spinconfig_[neighbor] = 2*dir + 3*PI - oldspin;
				else if((2*dir + PI - oldspin) > 2*PI) Spinconfig_[neighbor] = 2*dir - PI - oldspin;
				else Spinconfig_[neighbor] = 2*dir + PI - oldspin;
				Instack[neighbor] = true;
			}
		}
	}
}

//calculate the energy and magnetization
void XYSimulation::Reset() {
	double e = 0;
	int ns;
	for( int i = 0; i < N_; i++) {
		ns = (i%L_ == L_-1) ? i - i%L_ :  i + 1;
		e += -cos(Spinconfig_[i]-Spinconfig_[ns]);
		ns = (i + L_ > N_ -1) ? i + L_ - N_ :  i + L_;
		e += -cos(Spinconfig_[i]-Spinconfig_[ns]);
				
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

