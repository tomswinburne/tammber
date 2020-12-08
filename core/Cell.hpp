/*
   Copyright (c) 2016, Los Alamos National Security, LLC
   All rights reserved.
   Copyright 2016. Los Alamos National Security, LLC. This software was produced under U.S. Government contract DE-AC52-06NA25396 for Los Alamos National Laboratory (LANL), which is operated by Los Alamos National Security, LLC for the U.S. Department of Energy. The U.S. Government has rights to use, reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is modified to produce derivative works, such modified software should be clearly marked, so as not to confuse it with the version available from LANL.

   Additionally, redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
   1.      Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
   2.      Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
   3.      Neither the name of Los Alamos National Security, LLC, Los Alamos National Laboratory, LANL, the U.S. Government, nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LOS ALAMOS NATIONAL SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef Cell_hpp
#define Cell_hpp

#include "Constants.hpp"

#include <Eigen/Dense>
#include <array>

#include <boost/serialization/vector.hpp>
#include <boost/serialization/array.hpp>
#include <iostream>

using Eigen::MatrixXd;

class Cell {
public:
std::array< std::array<double,NDIM>, NDIM > rsp;
std::array< std::array<double,NDIM>, NDIM > csp;
std::array<double,NDIM> periodic;
std::array<double,NDIM> origin;

void print(){

	std::cout<<"ORIGIN"<<std::endl;
	for(int i=0; i<NDIM; i++) {
		std::cout<<origin[i]<<" ";
	}
	std::cout<<std::endl;

	std::cout<<"PERIODICITY"<<std::endl;
	for(int i=0; i<NDIM; i++) {
		std::cout<<periodic[i]<<" ";
	}
	std::cout<<std::endl;

	std::cout<<"RSP"<<std::endl;
	for(int i=0; i<NDIM; i++) {
		for(int j=0; j<NDIM; j++) {
			std::cout<<rsp[i][j]<<" ";
		}
		std::cout<<std::endl;
	}

	std::cout<<"CSP"<<std::endl;
	for(int i=0; i<NDIM; i++) {
		for(int j=0; j<NDIM; j++) {
			std::cout<<csp[i][j]<<" ";
		}
		std::cout<<std::endl;
	}

};


void update(){
	//compute the projectors here.
	MatrixXd rspe(NDIM,NDIM);
	MatrixXd cspe(NDIM,NDIM);

	for(int i=0; i<NDIM; i++) {
		for(int j=0; j<NDIM; j++) {
			rspe(i,j)=rsp[i][j];
		}
	}
	cspe=rspe.inverse();
	for(int i=0; i<NDIM; i++) {
		for(int j=0; j<NDIM; j++) {
			csp[i][j]=cspe(i,j);
		}
	}
};


bool inside(std::array<double,NDIM> r){

	std::array<double, NDIM> s;
	//std::cout<<"r"<<std::endl;
	for(int i=0; i<NDIM; i++) {
		r[i]-=origin[i];
		//std::cout<<r[i]<<" ";
	}
	//std::cout<<std::endl;

	//project to cell-space
	for(int i=0; i<NDIM; i++) {
		s[i]=0;
		for(int j=0; j<NDIM; j++) {
			s[i]+=csp[i][j]*r[j];
		}
	}

	bool in=true;
	//std::cout<<"s"<<std::endl;
	for(int i=0; i<NDIM; i++) {
		//std::cout<<s[i]<<" ";
		if( (s[i]<0 or s[i]>=1)and not periodic[i] ) {
			in=false;
			break;
		}
	}
	//std::cout<<std::endl;
	return in;
};

bool inside(std::array<double,NDIM> r, Cell globalCell,std::array<double, NDIM> &x){

	x=r;
	std::array<double, NDIM> s;

	for(int i=0; i<NDIM; i++) {
		r[i]-=origin[i];

		if(r[i]<0. and globalCell.periodic[i]) {
			r[i]+=globalCell.rsp[i][i];
			x[i]+=globalCell.rsp[i][i];
		}
		if(r[i]>=rsp[i][i] and globalCell.periodic[i]) {
			r[i]-=globalCell.rsp[i][i];
			x[i]-=globalCell.rsp[i][i];
		}

		//std::cout<<r[i]<<" ";
		//std::cout<<x[i]<<" ";
	}
	//std::cout<<std::endl;

	//project to cell-space
	for(int i=0; i<NDIM; i++) {
		s[i]=0;
		for(int j=0; j<NDIM; j++) {
			s[i]+=csp[i][j]*r[j];
		}
		//std::cout<<s[i]<<" ";
	}
	//std::cout<<std::endl;
	//std::cout<<std::endl;

	bool in=true;
	for(int i=0; i<NDIM; i++) {
		if( (s[i]<0 or s[i]>=1)and not periodic[i] ) {
			in=false;
			break;
		}
	}
	return in;
};

// wrap a vector to lie in cell
void wrap(std::array<double,NDIM> &dr) {
	std::array<double,NDIM> ds;

	//project to cell-space
	for(int i=0; i<NDIM; i++) {
		ds[i] = 0.0;
		for(int j=0; j<NDIM; j++) ds[i] += csp[i][j]*dr[j];
	}

	//wrap if necessary where allowed
	for(int i=0; i<NDIM; i++) ds[i] -= round(ds[i])*periodic[i];

	//project to R-space
	for(int i=0; i<NDIM; i++) {
		dr[i] = 0.0;
		for(int j=0; j<NDIM; j++) dr[i] += rsp[i][j]*ds[j];
	}
};

// wrap a vector to lie in [-L/2,L/2]
void wrapc(std::array<double,NDIM> &dr) {
	std::array<double,NDIM> ds;

	//project to cell-space
	for(int i=0; i<NDIM; i++) {
		ds[i] = 0.0;
		for(int j=0; j<NDIM; j++) ds[i] += csp[i][j]*dr[j];
	}

	//wrap to lie in [-0.5,0.5]
	for(int i=0; i<NDIM; i++) if(periodic[i]) {
		while(ds[i]>=+0.5) ds[i] -= 1.0;
		while(ds[i]<=-0.5) ds[i] += 1.0;
	}

	//project to R-space
	for(int i=0; i<NDIM; i++) {
		dr[i] = 0.0;
		for(int j=0; j<NDIM; j++) dr[i] += rsp[i][j]*ds[j];
	}
};


void wrapin(std::array<double,NDIM> &dr) {
	std::array<double,NDIM> ds;

	//project to cell-space
	for(int i=0; i<NDIM; i++) {
		ds[i] = 0.0;
		for(int j=0; j<NDIM; j++) ds[i] += csp[i][j]*dr[j];
	}

	//wrap if necessary where allowed
	for(int i=0; i<NDIM; i++) {
		while(ds[i]<0.0) ds[i] += 1.0;
		while(ds[i]>1.0) ds[i] -= 1.0;
	}

	//project to R-space
	for(int i=0; i<NDIM; i++) {
		dr[i] = 0.0;
		for(int j=0; j<NDIM; j++) dr[i] += rsp[i][j]*ds[j];
	}
};

// pbc wrapping with some overloading
void minimumImageVector(std::vector<double> &x,std::vector<double> &y,std::vector<double> &dx, double &mag, double &maxmag) {
	int ne = x.size();
	assert(ne%NDIM==0); // required....
	int n = ne/NDIM;
	double amag;
	std::array<double,NDIM> dr;

	dx.resize(ne,0.0);

	mag=0.0;
	maxmag=0.0;
	for(int i=0;i<n;i++) {
		for(int j=0;j<NDIM;j++) dr[j] = y[NDIM*i+j] - x[NDIM*i+j];
		wrap(dr);
		amag = 0.0;
		for(int j=0;j<NDIM;j++) {
			dx[NDIM*i+j] = dr[j];
			amag += dr[j] * dr[j];
		}
		mag += amag;
		maxmag = std::max(maxmag,amag);
	}
};
void minimumImageVector(std::vector<double> &x,std::vector<double> &y,std::vector<double> &dx) {
	double mag,maxmag;
	minimumImageVector(x,y,dx,mag,maxmag);
};

void minimumImageVector(std::vector<double> &x,std::vector<double> &y, double &mag, double &maxmag) {
	std::vector<double> dx;
	minimumImageVector(x,y,dx,mag,maxmag);
};

double msd(std::vector<double> &x,std::vector<double> &y,bool maxd=false) {
	double mag,mmag;
	minimumImageVector(x,y,mag,mmag);
	if(maxd) return sqrt(mmag);
	else return sqrt(mag);
};


void minimumImageVector(std::vector<double> &dx, double &mag, double &maxmag) {
	int ne = dx.size();
	assert(ne%NDIM==0); // required....
	int n = ne/NDIM;
	double amag;
	std::array<double,NDIM> dr;
	mag=0.0;
	maxmag=0.0;
	for(int i=0;i<n;i++) {
		for(int j=0;j<NDIM;j++) dr[j] = dx[NDIM*i+j];
		wrap(dr);
		amag = 0.0;
		for(int j=0;j<NDIM;j++) {
			dx[NDIM*i+j] = dr[j];
			amag += dr[j] * dr[j];
		}
		mag += amag;
		maxmag = std::max(maxmag,amag);
	}
};
void minimumImageVector(std::vector<double> &dx) {
	double mag,maxmag;
	minimumImageVector(dx,mag,maxmag);
};


template<class Archive>
void serialize(Archive & ar, const unsigned int version) {
	ar & rsp;
	ar & csp;
	ar & periodic;
	ar & origin;
};

};

#endif
