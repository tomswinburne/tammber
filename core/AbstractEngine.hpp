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



#ifndef AbstractEngine_h
#define AbstractEngine_h

#include "Task.hpp"
#include "Graph.hpp"
#include "Types.hpp"
#include "Log.hpp"

#include <mpi.h>
#include <functional>

#include <boost/lexical_cast.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/random/random_device.hpp>
#include <boost/random.hpp>
#include <boost/timer/timer.hpp>
#include <boost/filesystem.hpp>
#include <boost/regex.hpp>


class ParameterParser {
public:

ParameterParser() {
	boost::random::random_device rd;
	rng.seed(rd());
	count=1;
}

void seed(int s){
	rng.seed(s);
};

std::string parse(std::string r, std::unordered_map<std::string, std::string> &parameters){
	std::string raw=r;
	//replacements based on parameters
	for(auto it=parameters.begin(); it!=parameters.end(); it++ ) {
		std::string key=it->first;
		key="%"+key+"%";
		boost::trim(key);
		boost::replace_all(raw, key, it->second);
	}

	return parse(raw);
};

std::string parse(std::string r){
	std::string raw=r;

	//standard replacements
	boost::random::uniform_01<> uniform;
	boost::random::uniform_int_distribution<> d(1,1000000);

	std::string key="%RANDU%";
	while(not boost::find_first(raw, key).empty()) {
		int r=d(rng);
		std::string s=boost::str(boost::format("%1%" ) % r );
		boost::replace_first(raw, key, s);
	}
	key="%COUNT%";
	while(not boost::find_first(raw, key).empty()) {
		std::string s=boost::str(boost::format("%1%" ) % count );
		boost::replace_first(raw, key, s);
		count++;
	}
	key="%FILE_INCREMENT%";
	if(!boost::find_first(raw, key).empty())
	{
		std::vector< std::string > all_matching_files;
		int max=0;
		boost::replace_all(r, key, "(\\d+)");
		boost::filesystem::path p( r );
		boost::filesystem::path parent=p.parent_path();
		std::string filep=p.filename().string();
		const boost::regex filter( filep );


		bool match=false;
		try{
			boost::filesystem::directory_iterator end_itr; // Default ctor yields past-the-end
			for( boost::filesystem::directory_iterator i( parent ); i != end_itr; ++i )
			{
				// Skip if not a file
				if( !boost::filesystem::is_regular_file( i->status() ) ) continue;

				boost::smatch what;
				// For V3:
				if( !boost::regex_match( i->path().filename().string(), what, filter ) ) continue;

				// File matches, extract number
				for(int ii=1; ii<what.size(); ii++) {
					match=true;
					int index=boost::lexical_cast<int>(what[ii]);
					max=std::max(max,index);
				}
			}
			if(match) {
				max++;
			}
		}
		catch(...) {

		}
		std::string s=boost::str(boost::format("%010d" ) % max );
		boost::replace_all(raw, key, s);
	}
	return raw;
};

std::vector<std::string> splitLines(std::string r){
	std::vector< std::string > sc;
	boost::split( sc, r, boost::is_any_of("\n"), boost::token_compress_off );
	return sc;
};

private:
boost::random::mt11213b rng;
int count;

};

/**
	Abstract base class that defines an interface to a simulation engine.
	Each engine needs to provide an Engine class derived from this one.
*/
template <class EngineTaskMapper>
class AbstractEngine {
public:
AbstractEngine(boost::property_tree::ptree &config, MPI_Comm localComm_, int seed_){
	local_rank=0;
	impls["TASK_DIE"] = die_impl;
	impls["TASK_NOTHING"] = nothing_impl;
};

std::function<void(GenericTask&)> die_impl = [this](GenericTask &task) {
	/* Derived class can overwrite (but needs repeat of
		impls["TASK_DIE"]=die_impl;	in constructor) */
};

std::function<void(GenericTask&)> nothing_impl = [this](GenericTask &task) {
	/* Derived class can overwrite (but needs repeat of
		impls["TASK_NOTHING"]=nothing_impl;	in constructor) */
};

virtual void process(GenericTask &task) {
	std::string str_id = mapper.type(task.type);
	if(impls.find(str_id)!=impls.end()) {
		LOGGER("EXECUTING "<<str_id<<", "<<task.type)
		impls[str_id](task);

	}
};

virtual bool failed() {
	return false;
};

int local_rank;

protected:
EngineTaskMapper mapper;
std::map< const std::string, std::function<void(GenericTask&)> > impls;
#ifdef USE_BOOST_LOG
boost::log::sources::severity_logger< boost::log::trivial::severity_level > lg;
#endif


};

#endif /* AbstractEngine_h */
