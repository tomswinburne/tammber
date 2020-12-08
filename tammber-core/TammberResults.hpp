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



#ifndef TammberTask_hpp
#define TammberTask_hpp

#include "Task.hpp"
#include "Pack.hpp"
#include "Types.hpp"

//#include "CustomTypedefs.hpp"

#include <map>
#include <string>
#include <boost/property_tree/ptree.hpp>

template<class BundleTaskMapper>
class TammberResultBundle : public TaskResultBundle {
public:
TammberResultBundle() : TaskResultBundle(){
};

TammberResultBundle(int workerId_,Label bundleId_) : TaskResultBundle(workerId_,bundleId_){
};

virtual int size() {
	return int(segments.size()+pathways.size());
};

virtual void clear(){
	segments.clear();
	pathways.clear();
	TaskResultBundle::clear();
};

virtual bool empty(){
	return segments.size()==0 and pathways.size()==0 and TaskResultBundle::empty();
};

void assimilate(TaskResultBundle &t){
	//std::cout<<"TammberResultBundle:: ASSIMILATING BUNDLE "<<t.completedTasks.size()<<std::endl;
	TammberResultBundle &tt=dynamic_cast<TammberResultBundle&>(t);
	// no need to group, just a list is fine.....
	//std::cout<<"TammberResultBundle:: segments.size()=="<<tt.segments.size()<<"\n";
	for(auto it=tt.segments.begin();it!=tt.segments.end();it++) segments.push_back(*it);
	//std::cout<<"TammberResultBundle:: pathways.size()=="<<tt.pathways.size()<<"\n";
	for(auto it=tt.pathways.begin();it!=tt.pathways.end();it++) {
		std::cout<<"ASSIMILATING NEB / PAIRMAP: "<<it->info_str()<<std::endl;
		pathways.push_back(*it);
	}
	TaskResultBundle::assimilate(t);
};


virtual void assimilate(GenericTask &t) {
	//std::cout<<"ASSIMILATING TASK "<<t.type<<std::endl;
	if(t.type == mapper.type("TASK_SEGMENT")) {
		std::list<TADSegment> new_segments;
		::extract("TADSegment",t.returns,new_segments);
		for(auto seg: new_segments)	segments.push_back(seg);
	} else if(t.type == mapper.type("TASK_NEB")) {
		std::list<NEBPathway> new_pathways;
		::extract("NEBPathway",t.returns,new_pathways);
		for(auto path: new_pathways)	pathways.push_back(path);
	} else if(t.type == mapper.type("TASK_PAIRMAP")) {
		std::list<NEBPathway> new_pathways;
		::extract("NEBPathway",t.returns,new_pathways);
		for(auto path: new_pathways)	pathways.push_back(path);
	} else {
		TaskResultBundle::assimilate(t);
	}
};

template<class Archive>
void serialize(Archive & ar, const unsigned int version){
	ar & boost::serialization::base_object<TaskResultBundle>(*this);
	ar & segments;
	ar & pathways;
};

std::list<TADSegment> segments;
std::list<NEBPathway> pathways;
BundleTaskMapper mapper;
};



#endif
