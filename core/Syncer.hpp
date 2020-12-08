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


#ifndef Syncer_hpp
#define Syncer_hpp

#include <boost/property_tree/ptree.hpp>
#include <boost/function.hpp>
#include <boost/functional/factory.hpp>
#include <boost/bind.hpp>
#include <boost/lexical_cast.hpp>

#include <list>
#include <array>

#include "AbstractSystem.hpp"
#include "Constants.hpp"
#include "BoundaryConditions.hpp"

class AbstractSyncer {
public:
virtual void initialize(std::unordered_map<std::string,std::string> &parameters)=0;
virtual void sync(AbstractSystem &reference, std::list<AbstractSystem* > &neighbors, AbstractSystem &synchedSystem)=0;
};

class NullSyncer : public AbstractSyncer {
virtual void initialize(std::unordered_map<std::string,std::string> &parameters){
};
virtual void sync(AbstractSystem &reference, std::list<AbstractSystem* > &neighbors, AbstractSystem &synchedSystem){
	synchedSystem=reference;
};
};


class SLSyncer : public AbstractSyncer {
public:
Cell globalCell;
std::array<double,NDIM> activePadding;
std::array<double,NDIM> bufferPadding;

virtual void initialize(std::unordered_map<std::string,std::string> &parameters){

	{
		double xlo,xhi,ylo,yhi,zlo,zhi;
		std::stringstream ss;
		std::string s=parameters["GlobalCell"];
		ss<<s;
		ss>>xlo>>xhi>>ylo>>yhi>>zlo>>zhi;

		globalCell.origin={xlo,ylo,zlo};
		globalCell.rsp[0][0] = xhi-xlo;
		globalCell.rsp[0][1] = 0.0;
		globalCell.rsp[0][2] = 0.0;

		globalCell.rsp[1][0] = 0;
		globalCell.rsp[1][1] = yhi-ylo;
		globalCell.rsp[1][2] = 0.0;

		globalCell.rsp[2][0] = 0;
		globalCell.rsp[2][1] = 0;
		globalCell.rsp[2][2] = zhi-zlo;
	}
	{
		std::string s=parameters["Periodicity"];
		std::stringstream ss;
		ss<<s;
		ss>>globalCell.periodic[0]>>globalCell.periodic[1]>>globalCell.periodic[2];
	}
	{
		std::stringstream ss;
		std::string sp=parameters["ActivePadding"];
		ss<<sp;
		ss>>activePadding[0]>>activePadding[1]>>activePadding[2];
	}
	{
		std::stringstream ss;
		std::string sp=parameters["BufferPadding"];
		ss<<sp;
		ss>>bufferPadding[0]>>bufferPadding[1]>>bufferPadding[2];
	}
	globalCell.update();
};


virtual void sync( AbstractSystem &reference, std::list<AbstractSystem* > &neighbors, AbstractSystem &synchedSystem){


	//std::cout<<"SYNCING "<<neighbors.size()<<std::endl;

	std::set<int> usedId;
	Cell c=reference.getCell();
	/*
	     std::cout<<"REFERENCE CELL"<<std::endl;
	     c.print();
	 */

	Cell activeRegion=reference.getCell();

	for(int i=0; i<NDIM; i++) {
		activeRegion.origin[i]+=activePadding[i];
		activeRegion.rsp[i][i]-=2*activePadding[i];
		if(activePadding[i]<1e-12) {
			activeRegion.periodic[i]=1;
		}
		else{
			activeRegion.periodic[i]=0;
		}

		c.origin[i]+=bufferPadding[i];
		c.rsp[i][i]-=2*bufferPadding[i];
		if(bufferPadding[i]<1e-12) {
			c.periodic[i]=1;
		}
		else{
			c.periodic[i]=0;
		}
	}
	activeRegion.update();
	c.update();

	synchedSystem.clearAtoms();
	//synchedSystem.setCell(c);
	synchedSystem.setCell(reference.getCell());

/*
        std::cout<<"ACTIVE CELL"<<std::endl;
        activeRegion.print();


        std::cout<<"NON-BUFFER CELL"<<std::endl;
        c.print();
 */
	//keep all atoms in the active region of the reference system

	int nAtoms=reference.getNAtoms();
	//std::cout<<nAtoms<<" ATOMS IN REFERENCE"<<std::endl;
	for(int i=0; i<nAtoms; i++) {
		int uniqueId=reference.getUniqueID(i);
		//std::cout<<uniqueId<<std::endl;
		if(usedId.count(uniqueId)==0) {
			std::array<double,NDIM> x;
			for(int k=0; k<NDIM; k++) {
				x[k]=reference.getPosition(i,k);
				//std::cout<<x[k]<<" "<<reference.getPosition(i,k)<<" "<<std::endl;
			}
			//std::cout<<std::endl;
			//check if this atoms is inside the active region, if it is, add it
			if(activeRegion.inside(x)) {
				usedId.insert(uniqueId);
				synchedSystem.copyAtom(reference,i);
			}
		}
	}
	//std::cout<<"KEPT "<<usedId.size()<<" ATOMS FROM REFERENCE "<<std::endl;

	//push the reference to the back of the list

	//neighbors.clear();
	//neighbors.push_back(&reference);


	for(auto it=neighbors.begin(); it!=neighbors.end(); it++) {
		//loop on every atom
		int nAtoms=(*it)->getNAtoms();
		for(int i=0; i<nAtoms; i++) {
			int uniqueId=(*it)->getUniqueID(i);
			if(usedId.count(uniqueId)==0) {
				std::array<double,NDIM> x;
				std::array<double,NDIM> xx;
				for(int k=0; k<NDIM; k++) {
					x[k]=(*it)->getPosition(i,k);
				}
				//check if this atoms is inside the box, if it is, add it
				if(c.inside(x,globalCell,xx)) {
					//Move this atom to xx

					for(int k=0; k<NDIM; k++) {
						(*it)->setPosition(i,k,xx[k]);
					}

					usedId.insert(uniqueId);
					synchedSystem.copyAtom( *(*it),i);
				}
			}
		}
	}
	//std::cout<<"KEPT "<<usedId.size()<<" ATOMS IN TOTAL "<<std::endl;


	//std::this_thread::sleep_for(std::chrono::seconds(12000));

};
};


typedef boost::function< std::shared_ptr<AbstractSyncer>() > SyncerFactory;

const std::map<std::string, SyncerFactory> syncerFactory={
	{"NullSyncer", boost::factory<std::shared_ptr<NullSyncer> >() },
	{"", boost::factory<std::shared_ptr<NullSyncer> >() },
	{"SLSyncer", boost::factory<std::shared_ptr<SLSyncer> >() }
};

#endif
