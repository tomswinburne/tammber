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



#ifndef BoundaryConditions_hpp
#define BoundaryConditions_hpp

#include <stdio.h>
#include <array>
#include <cmath>

#include "AbstractSystem.hpp"
#include "Constants.hpp"



inline static std::pair<std::array<double, NDIM>, std::array<double, NDIM> > primaryCellPosition(std::array<double, NDIM> &r, Cell &cell){

	std::array<double, NDIM> s;
	for(int i=0; i<NDIM; i++) {
		r[i]-=cell.origin[i];
		//std::cout<<r[i]<<std::endl;
	}

	//project to cell-space
	for(int i=0; i<NDIM; i++) {
		s[i]=0;
		for(int j=0; j<NDIM; j++) {
			s[i]+=cell.csp[i][j]*r[j];
		}
	}

	//wrap if necessary
	for(int i=0; i<NDIM; i++) {
		s[i]-=floor(s[i])*cell.periodic[i];
	}

	//project to R-space
	for(int i=0; i<NDIM; i++) {
		r[i]=0;
		for(int j=0; j<NDIM; j++) {
			r[i]+=cell.rsp[i][j]*s[j];
		}
		r[i]+=cell.origin[i];
	}

	return std::make_pair(r,s);
};



inline static std::pair<std::array<double, NDIM>, std::array<double, NDIM> > primaryCellPosition(int iAtom, AbstractSystem &system){

	Cell &c=system.cell;
	std::array<double, NDIM> s;
	std::array<double, NDIM> r;
	for(int i=0; i<NDIM; i++) {
		r[i]=system.getPosition(iAtom,i);
		//std::cout<<"-"<<r[i]<<std::endl;
	}

	return primaryCellPosition(r,c);
};

inline static double nearestImageDistance(int iAtom, int jAtom, AbstractSystem &system, Cell &cell){
	double r2=0;

	std::array<double, NDIM> ds;
	std::array<double, NDIM> dr;

	//obtain the distance in R-space

	for(int i=0; i<NDIM; i++) {
		dr[i]=system.getPosition(iAtom,i)-system.getPosition(jAtom,i);
	}

	//project to cell-space
	for(int i=0; i<NDIM; i++) {
		ds[i]=0;
		for(int j=0; j<NDIM; j++) {
			ds[i]+=cell.csp[i][j]*dr[j];
		}
	}


	//wrap if necessary (where allowed)
	for(int i=0; i<NDIM; i++) {
		ds[i]-=round(ds[i])*cell.periodic[i];
	}

	//project to R-space
	r2=0;
	for(int i=0; i<NDIM; i++) {
		dr[i]=0;
		for(int j=0; j<NDIM; j++) {
			dr[i]+=cell.rsp[i][j]*ds[j];
		}
		r2+=dr[i]*dr[i];
	}
	return sqrt(r2);
};

inline static double nearestImageDistance(std::array<double, NDIM> positions, int jAtom, AbstractSystem &system, Cell &cell){
	double r2=0;

	std::array<double, NDIM> ds;
	std::array<double, NDIM> dr;

	//obtain the distance in R-space

	for(int i=0; i<NDIM; i++) {
		dr[i]=positions[i]-system.getPosition(jAtom,i);
	}

	//project to cell-space
	for(int i=0; i<NDIM; i++) {
		ds[i]=0;
		for(int j=0; j<NDIM; j++) {
			ds[i]+=cell.csp[i][j]*dr[j];
		}
	}


	//wrap if necessary (where allowed)
	for(int i=0; i<NDIM; i++) {
		ds[i]-=round(ds[i])*cell.periodic[i];
	}

	//project to R-space
	r2=0;
	for(int i=0; i<NDIM; i++) {
		dr[i]=0;
		for(int j=0; j<NDIM; j++) {
			dr[i]+=cell.rsp[i][j]*ds[j];
		}
		r2+=dr[i]*dr[i];
	}
	return sqrt(r2);
};

#endif /* BoundaryConditions_hpp */
