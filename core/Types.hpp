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



#ifndef Types_hpp
#define Types_hpp

#include <stdio.h>
#include <vector>
#include <list>
#include <deque>
#include <memory>
#include <set>
#include <atomic>
#include <algorithm>
#include <future>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/unordered_collections_save_imp.hpp>
#include <boost/serialization/unordered_collections_load_imp.hpp>
#include <boost/serialization/split_free.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/unordered_set.hpp>
#include <boost/serialization/deque.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/deque.hpp>
#include <boost/functional/hash.hpp>
#include <boost/random/random_device.hpp>
#include <boost/random.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/type_index.hpp>
#include <iostream>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/count.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <chrono>

#include "Constants.hpp"

#include "Cell.hpp"

#include "Data.hpp"
#include "Log.hpp"


/*
template <typename T1, typename T2>
T1 cast_wrapper(T2 val) {
  try {
    return boost::lexical_cast<T1>(val);
  } catch (boost::bad_lexical_cast &) {
    std::cout<<"FAILED CASTING "<<boost::typeindex::type_id<T2>().pretty_name()<<" AS "<<boost::typeindex::type_id<T1>().pretty_name()<<" WITH VALUE "<<val<<std::endl;
    std::cout<<"RETURNING INSTEAD 1. CAST AS "<<boost::typeindex::type_id<T1>().pretty_name()<<std::endl;
    return boost::lexical_cast<T1>(1.);
  }

}

template <typename T1, typename T2, typename T3>
T1 cast_wrapper_tag(T2 val,T3 tag) {
  try {
    return boost::lexical_cast<T1>(val);
  } catch (boost::bad_lexical_cast &){
    std::cout<<"FAILED CASTING "<<boost::typeindex::type_id<T2>().pretty_name()<<" AS "<<boost::typeindex::type_id<T1>().pretty_name()<<" WITH VALUE "<<val<<std::endl;
    std::cout<<"CAST tag: "<<tag<<std::endl;
    std::cout<<"RETURNING INSTEAD 1. CAST AS "<<boost::typeindex::type_id<T1>().pretty_name()<<std::endl;
    return boost::lexical_cast<T1>(1.);
  }

}
*/

template <typename T1>
T1 safe_extractor(std::unordered_map<std::string,std::string> &para, std::string key, T1 defaultValue) {
	auto val = para.find(key);
	if(val==para.end()) {
		LOGGER("COULD NOT FIND "<<key<<"; RETURNING "<<defaultValue)
		return defaultValue;
	} else {
		try {
			LOGGER("CASTING "<<key<<" AS "<<boost::typeindex::type_id<T1>().pretty_name()<<" WITH VALUE "<<val->second)
			return boost::lexical_cast<T1>(val->second);
		} catch (boost::bad_lexical_cast &) {
			LOGGER("FAILED CASTING "<<key<<" AS "<<boost::typeindex::type_id<T1>().pretty_name()<<" WITH VALUE "<<val->second<<"; RETURNING "<<defaultValue)
		  return defaultValue;
		}
	}
}



class Timer {
public:
Timer(){
	t0 = std::chrono::high_resolution_clock::now();
};
void start(){
	t0 = std::chrono::high_resolution_clock::now();
};
int stop(){
	auto t1 = std::chrono::high_resolution_clock::now();
	return std::chrono::duration_cast<std::chrono::microseconds>(t1-t0).count();
};

int lapInSeconds(){
	auto t1 = std::chrono::high_resolution_clock::now();
	return std::chrono::duration_cast<std::chrono::seconds>(t1-t0).count();
}
private:
std::chrono::time_point<std::chrono::high_resolution_clock> t0;
};


class ExecutionTimer {
public:
ExecutionTimer(){
	logTimer.start();
};

void start(std::string name){
	timers[name]=Timer();
};

void stop(std::string name){
	int duration=timers[name].stop();
	processes[name](duration);
};

void push(std::string name, int duration){
	processes[name](duration);
};


void report(std::string title="EXECUTION"){
	int logInterval=logTimer.stop();
	#ifdef USE_BOOST_LOG
	boost::log::sources::severity_logger< boost::log::trivial::severity_level > lg;
	BOOST_LOG_SEV(lg, boost::log::trivial::trace) <<"==================== BEGIN "<< title <<" TIMER REPORT =====================";
	for(auto it=processes.begin(); it!=processes.end(); it++) {
		BOOST_LOG_SEV(lg, boost::log::trivial::trace) <<std::setw(15)<< it->first<<" Interval: "<<logInterval<<" count: "<< boost::accumulators::count( it->second ) <<" min: "<<boost::accumulators::min( it->second )<<" mean: "<<boost::accumulators::mean( it->second )<<" max: "<<boost::accumulators::max( it->second );
	}
	BOOST_LOG_SEV(lg, boost::log::trivial::trace) <<"==================== END TIMER REPORT =====================";
	#else
	LOGGER("==================== BEGIN "<< title <<" TIMER REPORT =====================")
	for(auto it=processes.begin(); it!=processes.end(); it++) {
		LOGGER(std::setw(15)<< it->first<<" Interval: "<<logInterval<<" count: "<< boost::accumulators::count( it->second ) <<" min: "<<boost::accumulators::min( it->second )<<" mean: "<<boost::accumulators::mean( it->second )<<" max: "<<boost::accumulators::max( it->second ))
	}
	LOGGER("==================== END TIMER REPORT =====================")
	#endif
	logTimer.start();
	processes.clear();
};

private:
Timer logTimer;
std::map<std::string, Timer> timers;
std::map<std::string, boost::accumulators::accumulator_set< int, boost::accumulators::features< boost::accumulators::tag::min, boost::accumulators::tag::max, boost::accumulators::tag::mean > > > processes;

};



class AbstractSpinLock {

public:
virtual void lock()=0;
virtual void unlock()=0;
};

class SpinLock : public AbstractSpinLock {
std::atomic_flag locked = ATOMIC_FLAG_INIT;
public:
void lock() {
	while (locked.test_and_set(std::memory_order_acquire)) {; }
}
void unlock() {
	locked.clear(std::memory_order_release);
}
};

class NoSpinLock : public AbstractSpinLock {

public:
void lock(){
};
void unlock(){
};
};

struct SyncData {
	Label label;
	uint64_t batchId;
	SyncData(){
		label=0;
		batchId=0;
	};
	template<class Archive>
	void serialize(Archive & ar, const unsigned int version){
		ar & label;
		ar & batchId;
	};
};

// ParSplice specific types

struct StateStatistics {
	StateStatistics();
	StateStatistics(uint64_t label_);

	void update(uint64_t finalState);

	void update(StateStatistics s);

	void sampleEscapeBKL(double unknownPrior, boost::random::mt11213b &rand, boost::random::uniform_01<> &uniform, Label &final,  bool &absorbed,int &nSegments);
	void sampleSegmentBKL(double unknownPrior, boost::random::mt11213b &rand, boost::random::uniform_01<> &uniform, Label &final,  bool &absorbed);


	void clear();

	template<class Archive>
	void serialize(Archive & ar, const unsigned int version){
		ar & label;
		ar & nSegments;
		ar & nTransitions;
		ar & counts;
	};

	uint64_t label;
private:
	int nSegments;
	int nTransitions;
	std::unordered_map< uint64_t, int > counts;
};

struct TransitionStatistics {

	TransitionStatistics();

	void clear();

	void update(Label initialState, Label finalState);

	// overloading for TAD

	void assimilate(TransitionStatistics &s);

	void sampleSegmentBKL(Label &lb,double unknownPrior, Label &final, bool &absorbed);
	void sampleEscapeBKL(Label &lb,double unknownPrior, Label &final, bool &absorbed, int &nConsumed);


	int size();

	template<class Archive>
	void serialize(Archive & ar, const unsigned int version){
		ar & statistics;
	};


	std::unordered_map< uint64_t, StateStatistics > statistics;
private:
	boost::random::mt11213b rng;
	boost::random::uniform_01<> uniform;

};

struct Visit {
	Label label;
	unsigned int duration;
	Visit() {
		duration=0;
	};
	template<class Archive>
	void serialize(Archive & ar, const unsigned int version){
		ar & label;
		ar & duration;
	};
};

struct Trajectory {
	Trajectory();

	void clear();

	//append a visit to a trajectory. Extends the last visit if possible
	void appendVisit(Visit &v, bool extend=true);

	//splice two trajectories
	bool splice(Trajectory &t);

	bool empty();

	void truncate(int maxDuration, Trajectory &leftover);

	uint64_t duration();

	uint64_t& overhead();

	void print();

	void log();

	Visit& back();

	Visit& front();

	void pop_back();

	void pop_front();

	void reduce_length(int n);

	template<class Archive>
	void serialize(Archive & ar, const unsigned int version){
		ar & visits;
		ar & length;
		ar & index;
		ar & nSplice;
		ar & overhead_;
	};

	std::list<Visit> visits;
	uint64_t length;
	unsigned nSplice;
	uint64_t overhead_;
private:
	unsigned int index;

};

class SegmentDatabase {

public:

void clear();

void add(int flavor, Trajectory &t);

void merge(SegmentDatabase &supp);

int size();

int count(Label lb, int flavor);

bool front(Label lb, int flavor, Trajectory &t);

void pop_front(Label lb, int flavor);

void print();

void log();

int duration();

int duration(Label lb, int flavor);

template<class Archive>
void serialize(Archive & ar, const unsigned int version){
	ar & db;
};

std::unordered_map< std::pair<Label, int>, std::deque<Trajectory>, boost::hash <std::pair <Label, int> > > db;

};

#endif /* Types_h */
