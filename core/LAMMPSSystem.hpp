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



#ifndef LAMMPSSystem_hpp
#define LAMMPSSystem_hpp

#include <stdio.h>
#include "AbstractSystem.hpp"


#include "Constants.hpp"
#include <vector>
#include <array>
#include <map>
#include <set>
#include <algorithm>

#include <boost/serialization/vector.hpp>
#include <boost/serialization/array.hpp>

class LAMMPSSystem : public AbstractSystem {
public:
// data used by LAMMPS
// public so LAMMPSEngine can access vectors directly

int64_t ntimestep;              // timestep counter
int natoms;                     // # of atoms in system
int qflag;                      // 1/0 for yes/no charges defined
int fflag;
int mflag;
//boxabc = simulation box vectors abc, one per row
// boxlo,boxhi,xy,xz,yz = LAMMPS representation of simulation box

double energy;

std::array< std::array<double, NDIM>, NDIM > boxabc;
std::array<int, NDIM> periodicity;

std::array<double, NDIM> boxlo;
std::array<double, NDIM> boxhi;
double xy,xz,yz;

// required per-atom attributes

std::vector<int> id;
std::vector<int> species;
std::vector<double> x;
std::vector<double> v;

//std::vector<int> image;

// optional per-atom attributes
std::vector<double> q;
std::vector<double> f;
//std::vector<int> molecule;

// required methods

LAMMPSSystem();
virtual ~LAMMPSSystem();

virtual int getNAtoms();
virtual void setNAtoms(int n);

virtual void setCell(Cell c);
//virtual Cell getCell();

virtual int getUniqueID(int iAtom);
virtual void setUniqueID(int iAtom, int value);

virtual int getSpecies(int iAtom);
virtual void setSpecies(int iAtom, int value);

virtual double getPosition(int iAtom, int iDim);
virtual void setPosition(int iAtom, int iDim, double value);

virtual double getVelocity(int iAtom, int iDim);
virtual void setVelocity(int iAtom, int iDim, double value);

virtual double getForce(int iAtom, int iDim);
virtual void setForce(int iAtom, int iDim, double value);

virtual double getEnergy();
virtual void setEnergy(double value);

virtual double msd(LAMMPSSystem &s,bool maxd);

virtual void minimumImage(LAMMPSSystem &s);

virtual void remap(std::map<int,int> newUniqueID);

virtual void copyAtom(AbstractSystem &s, int iAtom);


std::vector<double> extractCompressibleData();

void restoreCompressibleData(std::vector<double> &data);


template<class Archive>
void serialize(Archive & ar, const unsigned int version) {
	//ar & label;
	//ar & type;
	//ar & necessity;
	ar & energy;
	ar & epoch;
	ar & ntimestep;
	ar & natoms;
	ar & qflag;
	ar & fflag;
	ar & boxabc;
	ar & periodicity;
	ar & boxlo;
	ar & boxhi;
	ar & xy;
	ar & xz;
	ar & yz;
	ar & id;
	ar & species;
	ar & x;
	ar & v;
	ar & f;
	//ar & image;
	ar & q;
	//ar & molecule;
	ar & cell;
};


// extra methods added by LAMMPS

void resize(int n);
int getNTimestep();
void setNTimestep(int64_t n);
void updateCell();

virtual int addAtom(std::array<double,NDIM> positions, int species, std::map<std::string, std::string> attributes);

private:



void clearAtoms();
};

template< class T >
void reorder(std::vector<T> &v, std::vector<int> const &order )  {
	for ( int s = 1, d; s < order.size(); ++s ) {
		for ( d = order[s]; d < s; d = order[d] ) ;
		if ( d == s ) while ( d = order[d], d != s ) std::swap( v[s], v[d] );
	}
}


#endif /* LAMMPSSystem_hpp */
