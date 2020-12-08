#ifndef SystemModifier_hpp
#define SystemModifier_hpp


#include <cstdint>
#include <limits>
#include <boost/property_tree/ptree.hpp>
#include <Eigen/Dense>
#include <boost/random/random_device.hpp>
#include <boost/random.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/function.hpp>
#include <boost/functional/factory.hpp>
#include <boost/bind.hpp>

#include "AbstractSystem.hpp"
#include "NeighborList.hpp"
#include "BoundaryConditions.hpp"

using Eigen::MatrixXd;


class AbstractSystemModifier {
public:

virtual void initialize(boost::property_tree::ptree &config){
	timeDelay=config.get<int>("Configuration.SystemModifier.Delay",0);
};


virtual bool modificationNeeded (uint64_t previousTime, uint64_t currentTime){
	if(timeDelay==0) {
		return false;
	}
	return previousTime/timeDelay != currentTime/timeDelay;
};

virtual bool modify ( AbstractSystem &s, uint64_t previousTime, uint64_t currentTime)=0;

protected:
uint64_t timeDelay;
};



class NullSystemModifier : public AbstractSystemModifier {
public:

virtual void initialize(boost::property_tree::ptree &config){
	AbstractSystemModifier::initialize(config);
};

virtual bool modify ( AbstractSystem &s, uint64_t previousTime, uint64_t currentTime){
	return false;
};

};




class DisplacementSystemModifier : public AbstractSystemModifier {
public:
virtual void initialize(boost::property_tree::ptree &config) {
	AbstractSystemModifier::initialize(config);


	dr=std::vector<double>(NDIM,0);

	std::string s=config.get<std::string>("Configuration.SystemModifier.Displacement");
	std::stringstream ss;
	ss<<s;
	ss>>dr[0]>>dr[1]>>dr[2];

	{
		s=config.get<std::string>("Configuration.SystemModifier.MovingAtoms");
		std::vector< std::string > splitVec;
		boost::split( splitVec, s, boost::is_any_of("-"), boost::token_compress_on );
		if(splitVec.size()==2) {
			mLow=boost::lexical_cast<int>(splitVec[0]);
			mHigh=boost::lexical_cast<int>(splitVec[1]);
		}
		//assume there is only one moving atom
		else{
			mLow=mHigh=config.get<int>("Configuration.SystemModifier.MovingAtoms");
		}
	}
};

virtual bool modify ( AbstractSystem &s, uint64_t previousTime, uint64_t currentTime){

	int nMoves= currentTime/this->timeDelay -previousTime/this->timeDelay;

	std::cout<<"DisplacementSystemModifier::modify "<<nMoves<<std::endl;

	if(nMoves==0) {
		return false;
	}

	for(int i=0; i<s.getNAtoms(); i++) {
		if(s.getUniqueID(i)>=mLow and s.getUniqueID(i)<=mHigh) {
			for(int j=0; j<NDIM; j++) {
				s.setPosition(i,j,s.getPosition(i,j)+nMoves*dr[j]);
			}
		}
	}

	return true;
};


std::vector<double> dr;
//int movingIndex;
int mLow;
int mHigh;
};




class AtomInserterSystemModifier : public AbstractSystemModifier  {
public:

void initialize(boost::property_tree::ptree &config){
	AbstractSystemModifier::initialize(config);
	boost::random::random_device rd;
	rng.seed(rd());

	x0=std::vector<double>(NDIM,0);
	rcreate=config.get<double>("Configuration.SystemModifier.RSphere");
	std::string scenter=config.get<std::string>("Configuration.SystemModifier.Center");
	std::stringstream ss;
	ss<<scenter;
	ss>>x0[0]>>x0[1]>>x0[2];


	BOOST_FOREACH(boost::property_tree::ptree::value_type &v, config.get_child("Configuration.SystemModifier.Atoms")) {
		int sp=v.second.get<int>("Label");
		double fraction=v.second.get<double>("Fraction");
		double rclose=v.second.get<double>("RClose");
		double rfar=v.second.get<double>("RFar");

		insertionRanges[sp]=std::make_pair(rclose,rfar);
		fractions[sp]=fraction;

		BOOST_FOREACH(boost::property_tree::ptree::value_type &vv, v.second.get_child("Attributes")) {
			std::string key=vv.first;
			std::string data=vv.second.data();
			boost::trim(key);
			boost::trim(data);
			//std::cout<<key<<" "<<data<<std::endl;
			attributes[sp][key]=data;
		}
	}
};



virtual bool modify ( AbstractSystem &s, uint64_t previousTime, uint64_t currentTime){
	std::cout<<"AtomInserterSystemModifier::modify "<<previousTime<<" "<<currentTime<<std::endl;
	/*
	   if( (not this->modificationNeeded(previousTime,currentTime)) or fractions.size()==0) {
	        return false;
	   }
	 */
	int nInsertions= currentTime/this->timeDelay -previousTime/this->timeDelay;
	//int nInsertions=1;

	std::cout<<"AtomInserterSystemModifier::modify "<<nInsertions<<std::endl;
	std::cout<<"NUMBER OF ATOMS BEFORE INSERT: "<<s.getNAtoms()<<std::endl;
	if(nInsertions==0) {
		return false;
	}

	Cell bc=s.getCell();

	boost::random::uniform_01<> uniform;


	for(int n=0; n<nInsertions; n++) {
		std::cout<<"AtomInserterSystemModifier::modify Inserting"<<std::endl;

		int nAtoms=s.getNAtoms();
		double rs=uniform(rng);
		auto it=fractions.begin();
		double ss=it->second;

		while(ss<rs) {
			it++;
			ss+=it->second;
		}
		int sp=it->first;

		bool inRange;
		double rmin;
		std::array<double,NDIM> xs;
		std::array<double,NDIM> xr;

		//this is slow and dumb. We can do better if needed

		//randomly try to insert until atom is within range
		do {
			inRange=false;
			rmin=1e15;
			//generate a random point in cell space
			for(int k=0; k<NDIM; k++) {
				xs[k]=uniform(rng);
			}
			//project to cell space
			for(int k=0; k<NDIM; k++) {
				xr[k]=0;
				for(int l=0; l<NDIM; l++) {
					xr[k]+=bc.rsp[k][l]*xs[l];
				}
			}
			//std::cout<<xr[0]<<" "<<xr[1]<<" "<<xr[2]<<std::endl;
			//compute min and max distances
			for(int i=0; i<nAtoms; i++) {
				double r=nearestImageDistance(xr, i, s, bc);
				rmin=std::min(rmin,r);
			}
			//std::cout<<"RMIN: "<<rmin<<" "<<insertionRanges[sp].first<<" "<<insertionRanges[sp].second<<std::endl;
			if(rmin>insertionRanges[sp].first and rmin<insertionRanges[sp].second) {
				inRange=true;
				s.addAtom(xr,sp,attributes[sp]);
				std::cout<<"INSERTION ACCEPTED"<<std::endl;
				std::cout<<"NUMBER OF ATOMS AFTER INSERT: "<<s.getNAtoms()<<std::endl;
			}

		}
		while(!inRange);

	}


	std::cout<<"AtomInserterSystemModifier::modify DONE "<<std::endl;

	return true;
};

private:

double rcreate;
std::vector<double> x0;
std::map<int,double> fractions;
std::map<int,std::pair<double,double> > insertionRanges;
std::map<int,std::map<std::string,std::string> > attributes;
boost::random::mt11213b rng;
};



typedef boost::function< std::shared_ptr<AbstractSystemModifier>() > ModifierFactory;



const std::map<std::string, ModifierFactory> modifierFactory={
	{"NullSystemModifier", boost::factory<std::shared_ptr<NullSystemModifier> >() },
	{"", boost::factory<std::shared_ptr<NullSystemModifier> >() },
	{"DisplacementSystemModifier", boost::factory<std::shared_ptr<DisplacementSystemModifier> >() },
	{"AtomInserterSystemModifier", boost::factory<std::shared_ptr<AtomInserterSystemModifier> >() }
};





#endif
