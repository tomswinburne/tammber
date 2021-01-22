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

#ifndef TransitionFilter_hpp
#define TransitionFilter_hpp

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

class AbstractTransitionFilter {
public:
virtual void initialize(std::unordered_map<std::string,std::string> &parameters){
};
virtual bool isValid(AbstractSystem &before, AbstractSystem &after, std::unordered_map<std::string,std::string> &parameters)=0;
};

class NullTransitionFilter : public AbstractTransitionFilter {
public:
virtual bool isValid(AbstractSystem &before, AbstractSystem &after, std::unordered_map<std::string,std::string> &parameters){
	return true;
}
};


class dXTransitionFilter : public AbstractTransitionFilter {
public:
	virtual bool isValid(AbstractSystem &before, AbstractSystem &after, std::unordered_map<std::string,std::string> &parameters) {

		Cell bc = before.getCell();
		bool maxd=true;
		double thresh = 0.6;
		if (parameters.find("MaxNorm")!=parameters.end()) maxd = boost::lexical_cast<bool>(parameters["MaxNorm"]);
		if (parameters.find("Thresh")!=parameters.end()) thresh = boost::lexical_cast<double>(parameters["Thresh"]);
		int natoms = before.getNAtoms();
		double mag=0., max_mag=0.;

		// need to copy over to respect AbstractSystem format
		std::vector<double> dx(NDIM*natoms,0.0);
		for(int i=0; i<natoms; i++) for(int j=0;j<NDIM;j++)
				dx[NDIM*i+j] = after.getPosition(i,j)-before.getPosition(i,j);

		bc.minimumImageVector(dx,mag,max_mag);
		parameters["dXmax"] = boost::lexical_cast<std::string>(sqrt(max_mag));
		parameters["dX"] = boost::lexical_cast<std::string>(sqrt(mag));
		if(maxd) return (sqrt(max_mag) > thresh);
		else return (sqrt(mag) > thresh);
	}
};

class DomainTransitionFilter : public AbstractTransitionFilter {
public:

virtual void initialize(std::unordered_map<std::string,std::string> &parameters){
	std::stringstream ss;
	std::string sp=parameters["DomainPadding"];
	ss<<sp;
	ss>>padding[0]>>padding[1]>>padding[2];
};

virtual bool isValid(AbstractSystem &before, AbstractSystem &after, std::unordered_map<std::string,std::string> &parameters){
	Cell cell=before.getCell();

	//create this using parameters
	Cell activeDomain=cell;

/*
   //TODO: fix me for general cells
        for(int k=0; k<NDIM; k++) {
                double l=0;
                for(int kk=0; kk<NDIM; kk++) {
                        l+=activeDomain.rsp[kk][k]*activeDomain.rsp[kk][k];
                }
                l=sqrt(l);
                double f=(l-2*pad)/l;
                for(int kk=0; kk<NDIM; kk++) {
                        activeDomain.origin[k]+=padding[k]*activeDomain.rsp[kk][k]/l;
                        activeDomain.rsp[kk][k]*=f;
                }
        }
 */
	for(int k=0; k<NDIM; k++) {
		activeDomain.origin[k]+=padding[k];
		activeDomain.rsp[k][k]-=2*padding[k];
		if(padding[k]<1e-12) {
			activeDomain.periodic[k]=1;
		}
		else{
			activeDomain.periodic[k]=0;
		}
	}

	activeDomain.update();

	std::array<double,NDIM> cm = {{0, 0, 0}};
	std::array<double,NDIM> sd = {{0, 0, 0}};
	std::array<double,NDIM> xb = {{0, 0, 0}};
	std::array<double,NDIM> xa = {{0, 0, 0}};

	int nAtoms = before.getNAtoms();
//find the center of mass of the displacement
	for(int i=0; i<nAtoms; i++) {
		for(int k=0; k<NDIM; k++) {
			xb[k]=before.getPosition(i,k);
			xa[k]=after.getPosition(i,k);
		}
		double dr=nearestImageDistance(xa,i,before,cell);
		/*
		            if(dr>0.1) {
		                    std::cout<<"D*: "<<dr<<" "<<i<<" "<<before.getUniqueID(i)<<" "<<after.getUniqueID(i)<<std::endl;
		                    for(int j=0; j<NDIM; j++) {
		                            std::cout<<xb[j]<<" ";
		                    }
		                    std::cout<<std::endl;
		                    for(int j=0; j<NDIM; j++) {
		                            std::cout<<xa[j]<<" ";
		                    }
		                    std::cout<<std::endl;
		            }
		 */
		for(int k=0; k<NDIM; k++) {
			cm[k]+=dr*xb[k];
			sd[k]+=dr;
			//cm[k]+=dr*xa[k];
		}
	}
	//activeDomain.print();
	for(int k=0; k<NDIM; k++) {
		//cm[k]/=nAtoms;
		cm[k]/=sd[k];
	}
	std::cout<<std::endl;
//is the center of mass in the active domain?
	return activeDomain.inside(cm);
};
std::array<double,NDIM> padding;
};

typedef boost::function< std::shared_ptr<AbstractTransitionFilter>() > TransitionFilterFactory;

const std::map<std::string, TransitionFilterFactory> transitionFilterFactory={
	{"NullTransitionFilter", boost::factory<std::shared_ptr<NullTransitionFilter> >() },
	{"dXTransitionFilter", boost::factory<std::shared_ptr<dXTransitionFilter> >() },
	{"", boost::factory<std::shared_ptr<NullTransitionFilter> >() },
	{"DomainTransitionFilter", boost::factory<std::shared_ptr<DomainTransitionFilter> >() }
};

#endif
