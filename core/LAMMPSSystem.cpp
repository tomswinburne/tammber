//
//  LAMMPSSystem.cpp

//
//  Created by Danny Perez on 11/17/16.
//  Copyright © 2016 dp. All rights reserved.
//
//  LAMMPS-specific content added by Steve Plimpton, Dec 2016
//


#include <cassert>
#include <map>
#include <numeric>
#include <algorithm>
#include <iostream>

#include <boost/lexical_cast.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/device/back_inserter.hpp>


#include "LAMMPSSystem.hpp"


// --------------------------------------------------------
// required methods
// --------------------------------------------------------

LAMMPSSystem::LAMMPSSystem() {
	epoch = 0;
	qflag = 0;
	fflag = 0;
	ntimestep = 0;              // timestep counter
	natoms = 0;                     // # of atoms in system
};

// --------------------------------------------------------

LAMMPSSystem::~LAMMPSSystem(){
};


int LAMMPSSystem::getNAtoms()
{
	return natoms;
};




void LAMMPSSystem::setNAtoms(int n)
{
	clearAtoms();
	natoms = n;
	id = std::vector<int>(natoms,0);
	species = std::vector<int>(natoms,0);
	x = std::vector<double>(NDIM*natoms,0);
	v = std::vector<double>(NDIM*natoms,0);
	if(fflag) f = std::vector<double>(NDIM*natoms,0);
	if(qflag) q = std::vector<double>(natoms,0);
};


void LAMMPSSystem::resize(int n)
{
	natoms = n;
	id.resize(n);
	species.resize(n);
	x.resize(NDIM*n);
	v.resize(NDIM*n);
	if(qflag) q.resize(n);
	if(fflag) f.resize(NDIM*n);
};


void LAMMPSSystem::clearAtoms()
{
	natoms = 0;
	id.clear();
	species.clear();
	x.clear();
	v.clear();
	f.clear();
	q.clear();
};

// --------------------------------------------------------




void LAMMPSSystem::setCell(Cell c){
	cell=c;
	cell.update();
	for(int i=0; i<NDIM; i++) {
		boxlo[i]=cell.origin[i];
		periodicity[i]=cell.periodic[i];
	}
	boxhi[0]=cell.rsp[0][0]+boxlo[0];
	boxhi[1]=cell.rsp[1][1]+boxlo[1];
	boxhi[2]=cell.rsp[2][2]+boxlo[2];
	xy=cell.rsp[1][0];
	xz=cell.rsp[2][0];
	yz=cell.rsp[2][1];
};
// --------------------------------------------------------

double LAMMPSSystem::msd(LAMMPSSystem &s,bool maxd) {

	return cell.msd(x,s.x,maxd);
}

void LAMMPSSystem::minimumImage(LAMMPSSystem &s) {
	std::vector<double> dx;
	cell.minimumImageVector(x,s.x,dx);
	for(int i=0;i<NDIM*natoms;i++) s.x[i] = x[i] + dx[i];
}

// --------------------------------------------------------

void LAMMPSSystem::setUniqueID(int iAtom, int value)
{
	assert(iAtom < natoms);
	id[iAtom] = value;
}

int LAMMPSSystem::getUniqueID(int iAtom)
{
	assert(iAtom < natoms);
	return id[iAtom];
}

// --------------------------------------------------------

void LAMMPSSystem::setSpecies(int iAtom, int value)
{
	assert(iAtom < natoms);
	species[iAtom] = value;
}

int LAMMPSSystem::getSpecies(int iAtom)
{
	assert(iAtom < natoms);
	return species[iAtom];
}

// --------------------------------------------------------

void LAMMPSSystem::setPosition(int iAtom, int iDim, double value)
{
	int index = iAtom*NDIM + iDim;
	assert(index < natoms*NDIM);
	x[index] = value;
}

double LAMMPSSystem::getPosition(int iAtom, int iDim)
{
	int index = iAtom*NDIM + iDim;
	assert(index < natoms*NDIM);
	return x[index];
}

// --------------------------------------------------------

void LAMMPSSystem::setVelocity(int iAtom, int iDim, double value)
{
	int index = iAtom*NDIM + iDim;
	assert(index < natoms*NDIM);
	v[index] = value;
}

double LAMMPSSystem::getVelocity(int iAtom, int iDim)
{
	int index = iAtom*NDIM + iDim;
	assert(index < natoms*NDIM);
	return v[index];
}

// --------------------------------------------------------

void LAMMPSSystem::setForce(int iAtom, int iDim, double value)
{
	int index = iAtom*NDIM + iDim;
	assert(index < natoms*NDIM);
	f[index] = value;
}

double LAMMPSSystem::getForce(int iAtom, int iDim)
{
	int index = iAtom*NDIM + iDim;
	assert(index < natoms*NDIM);
	return f[index];
}

// --------------------------------------------------------

void LAMMPSSystem::setEnergy(double value)
{
	energy = value;
}

double LAMMPSSystem::getEnergy()
{
	return energy;
}



// --------------------------------------------------------
// extra methods
// --------------------------------------------------------

int LAMMPSSystem::getNTimestep()
{
	return ntimestep;
}

void LAMMPSSystem::setNTimestep(int64_t n)
{
	ntimestep = n;
}

// --------------------------------------------------------

void LAMMPSSystem::updateCell(){
	for(int i=0; i<NDIM; i++) {
		cell.origin[i]=boxlo[i];
		cell.periodic[i]=periodicity[i];
	}

	cell.rsp[0][0] = boxhi[0]-boxlo[0];
	cell.rsp[0][1] = 0.0;
	cell.rsp[0][2] = 0.0;

	cell.rsp[1][0] = xy;
	cell.rsp[1][1] = boxhi[1]-boxlo[1];
	cell.rsp[1][2] = 0.0;

	cell.rsp[2][0] = xz;
	cell.rsp[2][1] = yz;
	cell.rsp[2][2] = boxhi[2]-boxlo[2];

	cell.update();
}

void LAMMPSSystem::copyAtom(AbstractSystem &s, int iAtom){

	//cast into a LAMMPSSystem
	LAMMPSSystem *other=dynamic_cast<LAMMPSSystem*>(&s);

	natoms++;
	resize(natoms);
	for(int k=0; k<NDIM; k++) {
		setPosition(natoms-1,k,other->getPosition(iAtom,k));
		setVelocity(natoms-1,k,other->getVelocity(iAtom,k));
		if(fflag) {
			if(other->fflag) setForce(natoms-1,k,other->getForce(iAtom,k));
			else setForce(natoms-1,k,0.0);
		}
	}
	setUniqueID(natoms-1,other->getUniqueID(iAtom));
	setSpecies(natoms-1,other->getSpecies(iAtom));

	if (qflag != other->qflag) LOGGERA("copyAtom charge mismatch! "<<iAtom)
	else if (qflag && other->qflag) q[natoms-1]=other->q[iAtom];
};

int LAMMPSSystem::addAtom(std::array<double,NDIM> positions, int species, std::map<std::string, std::string> attributes){

	if(attributes.count("q")) qflag=1;

	natoms++;
	resize(natoms);

	setUniqueID(natoms-1,natoms);
	// this was a bug right?
	for(int k=0; k<NDIM; k++) {
		setPosition(natoms-1,k,positions[k]);
		setVelocity(natoms-1,k,0);
		if(fflag) setForce(natoms-1,k,0.0);
	}
	setSpecies(natoms-1,species);
	setUniqueID(natoms-1,natoms);
	if(qflag) q[natoms-1]=boost::lexical_cast<double>(attributes["q"]);
	return natoms;
};

void LAMMPSSystem::remap(std::map<int,int> newUniqueID){


	int nAtoms=getNAtoms();

	for(int i=0; i<nAtoms; i++) setUniqueID(i, newUniqueID[getUniqueID(i)]);
	std::vector<int> asort(id.size()), ids=id, tid=id, tsp=species;
	std::vector<double> tx=x, tv=v, tf=f, tq=q;

	std::iota(asort.begin(), asort.end(), 0);
	auto comparator = [&ids](int a, int b){ return ids[a] < ids[b]; };
	std::sort(asort.begin(), asort.end(), comparator);


	if(qflag) for(int i = 0; i < ids.size(); ++i) q[i]=tq[asort[i]];

	for(int i = 0; i < ids.size(); ++i)  {
		id[i]=tid[asort[i]];
		species[i]=tsp[asort[i]];
		for(int j=0; j<NDIM; j++) {
			x[NDIM*i+j] = tx[ NDIM*asort[i]+j ];
			v[NDIM*i+j] = tv[ NDIM*asort[i]+j ];
		}
	}
	if(fflag)
		for(int i = 0; i < ids.size(); ++i)
			for(int j=0; j<NDIM; j++)
				f[NDIM*i+j] = tf[ NDIM*asort[i]+j ];
};


std::vector<double> LAMMPSSystem::extractCompressibleData(){
	std::vector<double> data;
	data.reserve(2*x.size());

	data.insert(data.end(),x.begin(),x.end());
	data.insert(data.end(),q.begin(),q.end());

	x.clear();
	q.clear();
	v.clear();
	f.clear();
	return data;
};



void LAMMPSSystem::restoreCompressibleData(std::vector<double> &data){
	int nAtoms=getNAtoms();
	x.reserve(NDIM*nAtoms);
	if(qflag) q.reserve(nAtoms);
	x.insert(x.end(),data.begin(),data.begin()+NDIM*nAtoms);
	if(qflag) q.insert(q.end(),data.begin()+NDIM*nAtoms,data.begin()+(NDIM+1)*nAtoms);
	v.insert(v.end(),NDIM*nAtoms,0);
	if(fflag) f.insert(f.end(),NDIM*nAtoms,0);
};
