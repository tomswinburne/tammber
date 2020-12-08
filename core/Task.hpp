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



#ifndef Task_hpp
#define Task_hpp

#include <map>
#include <set>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <random>
#include <limits>
#include <stdio.h>
#include <memory>

#include <boost/serialization/vector.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/unordered_set.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/bimap.hpp>
#include <boost/random.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/functional/hash/hash.hpp>

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/identity.hpp>
#include <boost/multi_index/member.hpp>


#include "Types.hpp"
#include "AbstractSystem.hpp"
#include "Pack.hpp"
#include "Data.hpp"
#include "Log.hpp"

#define COMPLETED_TASK_TAG 1
#define TASK_REQUEST_TAG 2


/*
	Abstract task mapping class
*/
#define FIRST_TASK_ID -1
#define NOT_FOUND_ID -2
#define NOT_FOUND_STR "NONE"

class AbstractTaskMapper {
public:
AbstractTaskMapper() {
	insert("TASK_DIE");
	insert("TASK_NOTHING");
};

int type(std::string str) {
	boost::bimap<std::string,int>::left_const_iterator it = type_id_bimap.left.find(str);
	if (it == type_id_bimap.left.end()) return NOT_FOUND_ID;
	else return type_id_bimap.left.at(str);
};

std::string type(int id) {
	boost::bimap<std::string,int>::right_const_iterator it = type_id_bimap.right.find(id);
	if (it == type_id_bimap.right.end()) return NOT_FOUND_STR;
	else return type_id_bimap.right.at(id);
};

void insert(std::string str) {
	int id = FIRST_TASK_ID + type_id_bimap.size();
	type_id_bimap.insert( boost::bimap<std::string,int>::value_type(str,id) );
};

int size() {
	return type_id_bimap.size();
};

private:
	boost::bimap< std::string, int > type_id_bimap;
};



/**
 * TODO: Need mechanism to route results to different processes
 */
class AbstractTask {
public:

/**
 * Identify the producer of the request, not used yet
 */
int producer;

/**
 * Identify the consumer of the results, not used yet
 */
int consumer;

/**
 * Identifies the task type to be carried out
 */
int type;

/**
 * Identifies the task flavor
 */
int flavor;

/**
 * Identifies the task instance
 */
uint64_t id;

/**
 * Number of instances to execute
 */
int nInstances;

/**
 * Should the results of this task be consumed according to a preset order?
 */
bool imposeOrdering;

/**
 * Is this task optional?
 */
bool optional;

/**
 * Did the execution fail?
 */
bool failed;

virtual void clear(){
	type=-1;
	id=-1;
	flavor=-1;
	nInstances=0;
	imposeOrdering=false;
	optional=false;
	failed=false;
};

template<class Archive>
void serialize(Archive & ar, const unsigned int version){
	ar & type;
	ar & id;
	ar & flavor;
	ar & nInstances;
	ar & imposeOrdering;
	ar & optional;
	ar & consumer;
	ar & producer;
	ar & failed;
};


};


/**
 * Encapsulates a task to be run
 *
 * Conventions:
 * -after each execution of a task, the inputs will be empty
 */

template < class StoredDataIn > class TaskT :  public AbstractTask {
public:

TaskT(){
	imposeOrdering=false;
	optional=false;
	failed=false;
};

TaskT(AbstractTask &t) : AbstractTask(t) {
};


/**
 * Arguments. This is provided by the originator of the task.
 */
std::multimap<std::string, DataItem> arguments;

/**
 * Returns. This will be returned back to the originator.
 */
std::multimap<std::string, DataItem> returns;

/**
 * Data to be passed to the task. This will be assembled from the data-stores.
 */
//std::list<Data> inputData;
std::multimap<std::string, StoredDataIn> inputData;
/**
 * Data to be extracted from the task. This will be inserted into the data-stores.
 */
//std::list<StoredDataItem> outputData;
std::multimap<std::string, StoredDataItem> outputData;



/**
 * Assess if dependencies are met. Assume this is true, except if using DataPlaceholder, which gets a specialized implementation
 */
bool isRunnable(std::unordered_set<StoredDataIn> &available){
	return true;
};

void clearInputs(){
	arguments.clear();
	inputData.clear();
};

void clearOutputs(){
	returns.clear();
	outputData.clear();
};


/**
 * Transfer the outputs of t to the inputs of *this.
 */

void chain( TaskT &t, bool consume=true){
	inputData=t.outputData;
	arguments=t.arguments;
	if(consume) {
		t.clearOutputs();
	}
};



/**
 * Lists the Data that are necessary for this task to execute
 */
void requiredData(std::unordered_set<StoredDataIn> &deps){
	//loop over dependencies
	for(auto it=inputData.begin(); it!=inputData.end(); it++) {
		//this system must be available
		if(it->second.necessity==NECESSITY::REQUIRED) {
			//auto t=deps.insert(it->second);
			//t.first->shared=true;
			StoredDataIn t=it->second;
			t.shared=true;
			deps.insert(t);
		}
	}
};


virtual void clear(){
	AbstractTask::clear();
	arguments.clear();
	returns.clear();
	inputData.clear();
	outputData.clear();
};

template<class Archive>
void serialize(Archive & ar, const unsigned int version){
	//ar & boost::serialization::base_object<AbstractTask>(*this);
	ar & type;
	ar & id;
	ar & flavor;
	ar & nInstances;
	ar & imposeOrdering;
	ar & optional;
	ar & consumer;
	ar & producer;
	ar & arguments;
	ar & returns;
	ar & inputData;
	ar & outputData;
	ar & failed;
};

};

/**
 * Assess if dependencies on required Systems are met
 */
template<> inline bool TaskT<DataPlaceholder>::isRunnable(std::unordered_set<DataPlaceholder> &available){

	//loop over dependencies
	for(auto it=inputData.begin(); it!=inputData.end(); it++) {
		//this item must be available
		if(it->second.necessity==NECESSITY::REQUIRED) {
			DataPlaceholder s(it->second);
			//by convention, all systems are marked as OPTIONAL in 'available'.
			s.necessity=NECESSITY::OPTIONAL;
			if(available.count(s)==0) {
				return false;
			}
		}
	}
	return true;
};





typedef TaskT<DataPlaceholder> TaskDescriptor;
typedef TaskT<StoredDataItem> GenericTask;


/**
 * A bundle of task descriptors.
 */
class TaskDescriptorBundle {
public:

TaskDescriptorBundle(){
}

bool empty(){
	return tasks.empty();
}

void transferFrom(TaskDescriptorBundle &t){
	tasks.splice(tasks.end(),t.tasks);
	t.clear();
};

void transferTo(TaskDescriptorBundle &t, int nInstancesMax){
	int nTransfer=0;
	while(tasks.size()>0 and nTransfer<nInstancesMax) {
		t.insert(tasks.back());
		nTransfer+=tasks.back().nInstances;
		tasks.pop_back();
	}
};

int count(){
	int nInstances=0;
	for(auto it=tasks.begin(); it!=tasks.end(); it++) {
		nInstances+=it->nInstances;
	}
	return nInstances;
};

void copy(TaskDescriptorBundle &t){
	tasks.insert(tasks.end(),t.tasks.begin(),t.tasks.end());
};

void insert( TaskDescriptor &t){
	tasks.push_back(t);
};

bool extractRunnable(std::unordered_set<DataPlaceholder> &available, TaskDescriptor &t){
	for(auto it=tasks.begin(); it!=tasks.end(); it++) {
		if(it->isRunnable(available)) {
			t=*it;
			tasks.erase(it);
			return true;
		}
	}
	return false;
};

void clear(){
	tasks.clear();
};

template<class Archive>
void serialize(Archive & ar, const unsigned int version){
	ar & producerId;
	ar & tasks;
};

void requiredData( std::unordered_set<DataPlaceholder> &required){
	for(auto it=tasks.begin(); it!=tasks.end(); it++) {
		it->requiredData(required);
	}
};

std::list<TaskDescriptor>::iterator begin(){
	return tasks.begin();
};

std::list<TaskDescriptor>::iterator end(){
	return tasks.end();
};

int producerId;

std::list<TaskDescriptor> tasks;

};



/**
 * Manage the generation of tasks
 */

class TaskDescriptorRelay {
public:

virtual void initialize(int nProducers, int nConsumers){
	c=nConsumers;
	p=nProducers;
	r=0;
};

virtual void openBatch(){
	//std::cout<<"OPEN BATCH"<<std::endl;
	required.clear();
	optional.clear();
};

virtual void insert(int iProducer, TaskDescriptor &t){
	//std::cout<<"INSERTING TASK "<<t.type<<" "<<t.id<<" "<<t.optional<<" "<<t.imposeOrdering<<std::endl;
	int i=int(iProducer/double(p)*c);
	//if i is negative, assign in round-robin fashion
	if(i<0) {
		i=r%c;
		r++;
	}

	//std::cout<<"INSERTING TASK "<<t.type<<" "<<t.id<<" "<<t.optional<<" "<<t.imposeOrdering<<" "<<i<<" "<<c<<" "<<p<<" "<<r<<" "<<iProducer<<std::endl;

	if(t.optional) {
		optional[i].insert(t);
	}
	else{
		required[i].insert(t);
	}
};

virtual void insert(int iProducer,DataItem &data){
	int i=int(iProducer/double(p)*c);
	//if i is negative, assign in round-robin fashion
	if(i<0) {
		i=r%c;
		r++;
	}
	TaskDescriptorBundle bb;
	unpack(data.data,bb);
	for(auto it=bb.tasks.begin(); it!=bb.tasks.end(); it++) {
		if(it->optional) {
			optional[i].insert(*it);
		}
		else{
			required[i].insert(*it);
		}
	}
};

virtual bool release(int consumerId, DataItem &d){
	TaskDescriptorBundle bb;
	bool released=false;
	if(required.count(consumerId)) {
		bb.copy(required[consumerId]);
		released=true;
	}
	if(optional.count(consumerId)) {
		bb.copy(optional[consumerId]);
		released=true;
	}
	if(released) {
		pack(d.data,bb);
		return true;
	}
	else{
		return false;
	}

};

std::map<int,TaskDescriptorBundle> required;
std::map<int,TaskDescriptorBundle> optional;
//TaskDescriptorBundle b;
//TaskDescriptorBundle bFifo;
int c;
int p;
int r;
};



/**
 * Packs the results of completed tasks
 */
class TaskResultBundle {
public:
virtual ~TaskResultBundle() = default;

TaskResultBundle(){
	nReservations=0;
	nPendingReservations=0;
	queuedStamp=std::chrono::high_resolution_clock::now();
};

TaskResultBundle(int workerId_,Label bundleId_){
	workerId=workerId_;
	bundleId=bundleId_;
	nReservations=0;
	nPendingReservations=0;
	queuedStamp=std::chrono::high_resolution_clock::now();
};

virtual bool empty(){
	return completedTasks.size()==0;
};

virtual void clear(){
	completedTasks.clear();
};

virtual void assimilate(std::list<TaskResultBundle> &t) {
	for(auto it=t.begin(); it!=t.end(); it++) {
		assimilate(*it);
	}
};

virtual void assimilate(TaskResultBundle &t){
	//std::cout<<"TaskResultBundle:: ASSIMILATING BUNDLE "<<t.completedTasks.size()<<std::endl;
	for(auto itm=t.completedTasks.begin(); itm!=t.completedTasks.end(); itm++) {
		//std::cout<<itm->first<<std::endl;
		for(auto it=itm->second.begin(); it!=itm->second.end(); it++) {
			assimilate(*it);
		}
	}
};

virtual void assimilate(GenericTask &t){
	//std::cout<<"TaskResultBundle:: ASSIMILATING TASK "<<t.type<<std::endl;
	completedTasks[t.type].push_back(t);
	completedTasks[t.type].back().outputData.clear();
};

virtual void extract(int type, std::list<GenericTask> &tasks){
	if(completedTasks.count(type)) {
		tasks.splice(tasks.begin(),completedTasks[type]);
		completedTasks.erase(type);
	}
};


template<class Archive>
void serialize(Archive & ar, const unsigned int version){
	ar & bundleId;
	ar & workerId;
	ar & nReservations;
	ar & nPendingReservations;
	ar & completedTasks;
};


int workerId;
int nReservations;
int nPendingReservations;
Label bundleId;

std::map<int,std::list<GenericTask> > completedTasks;
std::chrono::high_resolution_clock::time_point queuedStamp;

};



#define RESULT_RESERVATIONS 1
#define RESULT_BUNDLES 2


template <class Bundle> class ResultManager {
public:

ResultManager(){
	currentBundle=0;
}

void setBundleSize(int size){
	maxBundleSize=size;
};

void initialize(boost::property_tree::ptree &config, int id_){
	timeout=config.get<unsigned>("Configuration.ResultManager.Timeout",1000000);
	maxBundleSize=config.get<int>("Configuration.ResultManager.BundleSize",1);
	id=id_;
	currentBundle=0;
	passthroughBundle.workerId=id;
	passthroughBundle.bundleId=0;
};

/**
 * Reserve a spot in a bundle
 */
void reserve(int workerId, Label promiseId){
	if(inProcessBundles.size()==0 or inProcessBundles.rbegin()->second.nReservations>=maxBundleSize ) {
		currentBundle++;
		Bundle b(id,currentBundle);
		bundleIterators[currentBundle]=inProcessBundles.insert(inProcessBundles.end(),  std::make_pair(currentBundle,b)  );

		//we have to register that bundle with the parent
		//unreportedBundles.push_back(std::pair<int,Label>(id,currentBundle));
		unreportedBundles.insert(std::pair<int,Label>(id,currentBundle));
	}

	inProcessBundles.rbegin()->second.nReservations++;
	inProcessBundles.rbegin()->second.nPendingReservations++;

	reservations[std::pair<int,Label>(workerId,promiseId)]=currentBundle;
};

void release(int workerId, Label promiseId){
	auto p=std::pair<int,Label>(workerId,promiseId);
	auto itr=reservations.find(p);
	if(itr!=reservations.end() ) {
		auto it=bundleIterators[itr->second];
		it->second.nPendingReservations--;
		it->second.nReservations--;
		reservations.erase(itr);
	}
};


/**
 * Deliver the results corresponding to a reservation
 */
template <class Result> void checkin(int workerId, Label promiseId, Result &t){
	//std::cout<<"CHECKIN BUNDLE "<<workerId<<" "<<promiseId<<std::endl;
	//this was an unreserved bundle. Pass straigth through
	auto p=std::pair<int,Label>(workerId,promiseId);
	auto itr=reservations.find(p);

	if(itr==reservations.end() ) {
		//std::cout<<"UNRESERVED BUNDLE "<<std::endl;
		passthroughBundle.assimilate(t);
	}
	else{
		//std::cout<<"RESERVED BUNDLE"<<std::endl;
		auto it=bundleIterators[itr->second];
		it->second.assimilate(t);
		it->second.nPendingReservations--;

		reservations.erase(itr);
	}
};

/**
 * Release content to parent
 */
//virtual void release(std::list<DataItem> &out){
virtual void release(std::multimap<std::string,DataItem> &out){

	out.clear();

	//std::cout<<"RELEASE"<<std::endl;
	if(unreportedBundles.size() > 0) {
		insert("RESERVATIONS",out,unreportedBundles);
		unreportedBundles.clear();
	}


	//std::list<Bundle> completedBundles;
	for(auto it=inProcessBundles.begin(); it!=inProcessBundles.end(); ) {

		if(it->second.nPendingReservations==0) {
			insert("COMPLETED_BUNDLE",out,it->second);
			bundleIterators.erase(it->first);
			it=inProcessBundles.erase(it);
		}
		else{
			//it++;
			break;
		}
	}

	if(not passthroughBundle.empty()) {
		//std::cout<<"RELEASING PASSTHROUGH BUNDLE "<<passthroughBundle.completedTasks.size()<<std::endl;
		insert("COMPLETED_BUNDLE",out,passthroughBundle);
		passthroughBundle.clear();
	}
};

virtual void extract(Bundle &out){
	//out.clear();

	unreportedBundles.clear();
	for(auto it=inProcessBundles.begin(); it!=inProcessBundles.end(); ) {
		if(it->second.nPendingReservations==0) {
			out.assimilate(it->second);
			bundleIterators.erase(it->first);
			it=inProcessBundles.erase(it);
		}
		else{
			it++;
		}
	}

	if(not passthroughBundle.empty()) {
		out.assimilate(passthroughBundle);
		passthroughBundle.clear();
	}

};


/**
 * Assimilate the released content of a child
 */
void assimilate(std::multimap<std::string,DataItem> &in){

	std::set<std::pair<int,Label> > newReservations;
	::extract("RESERVATIONS",in,newReservations);
	for(auto it=newReservations.begin(); it!=newReservations.end(); it++) {
		//std::cout<<"ASSIMILATING RESERVATION: "<<it->first<<" "<<it->second<<std::endl;
		reserve(it->first,it->second);
	}

	std::list<Bundle> bundles;
	::extract("COMPLETED_BUNDLE",in,bundles);
	//std::cout<<"ASSIMILATING BUNDLES "<<bundles.size()<<std::endl;
	for(auto it=bundles.begin(); it!=bundles.end(); it++) {
		std::cout<<"ASSIMILATING BUNDLE: "<<it->completedTasks.size()<<std::endl;
		checkin(it->workerId, it->bundleId, *it);
	}
};


void report(){
	#ifdef USE_BOOST_LOG
	boost::log::sources::severity_logger< boost::log::trivial::severity_level > lg;
	BOOST_LOG_SEV(lg, boost::log::trivial::trace) <<"RESULT MANAGER REPORT: "<<inProcessBundles.size()<<" BUNDLES, "<<reservations.size()<<" RESERVATIONS";
	#else
	std::cout<<"RESULT MANAGER REPORT: "<<inProcessBundles.size()<<" BUNDLES, "<<reservations.size()<<" RESERVATIONS"<<std::endl;
	#endif
};

private:

std::unordered_map< Label, typename std::map< Label, Bundle  >::iterator > bundleIterators;
std::map< Label, Bundle  > inProcessBundles;
//std::map<std::pair<int,Label>, Label > reservations;
std::set<std::pair<int,Label> > unreportedBundles;
//std::list< std::pair<int,Label> > unreportedBundles;
std::unordered_map<std::pair<int,Label>, Label, boost::hash< std::pair<int, Label> > > reservations;

Bundle passthroughBundle;
Label currentBundle;
int timeout;
int maxBundleSize;
int id;

};


/*
   inline std::map<std::string,std::string> extractParameters(int type, int flavor, int defaultFlavor,std::map< std::pair<int,int>, std::map<std::string,std::string> > &taskParameters){

        if(taskParameters.count(std::make_pair(type,flavor))) {
                return taskParameters[std::make_pair(type,flavor)];
        }
        else if(taskParameters.count(std::make_pair(type,defaultFlavor))) {
                return taskParameters[std::make_pair(type,defaultFlavor)];
        }
        else{
                //std::cerr<<"ERROR: extractParameters: NO SUCH PARAMETER SET"<<std::endl;
                return std::map<std::string,std::string>();
        }
   };
 */

inline std::unordered_map<std::string,std::string> extractParameters(int type, int flavor, int defaultFlavor,std::unordered_map< std::pair<int,int>, std::unordered_map<std::string,std::string>,  boost::hash< std::pair<int,int> > > &taskParameters){

	auto search = taskParameters.find(std::make_pair(type,flavor));

	if(search!=taskParameters.end()) {
		return search->second;
	}
	search = taskParameters.find(std::make_pair(type,defaultFlavor));
	if(search!=taskParameters.end()) {
		return search->second;
	}
	return std::unordered_map<std::string,std::string>();


	/*
	   else if{
	   search = taskParameters.find(std::make_pair(type,defaultFlavor));
	   if(search!=taskParameters.end())
	   {
	       return search->second;
	   }
	   }
	   else{
	   //std::cerr<<"ERROR: extractParameters: NO SUCH PARAMETER SET"<<std::endl;
	   return std::unordered_map<std::string,std::string>();
	   }
	 */

/*
        if(taskParameters.count(std::make_pair(type,flavor))) {
                return taskParameters[std::make_pair(type,flavor)];
        }
        else if(taskParameters.count(std::make_pair(type,defaultFlavor))) {
                return taskParameters[std::make_pair(type,defaultFlavor)];
        }
        else{
                //std::cerr<<"ERROR: extractParameters: NO SUCH PARAMETER SET"<<std::endl;
                return std::unordered_map<std::string,std::string>();
        }
 */
};


#endif /* Task_hpp */
