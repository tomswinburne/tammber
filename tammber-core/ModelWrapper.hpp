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


#ifndef __Tammber__Wrapper__
#define __Tammber__Wrapper__

#include <new>
#include <stdio.h>
#include <limits>
#include <memory>
#include <atomic>
#include <future>
#include <deque>
#include <chrono>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <thread>

#include <mpi.h>

#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>
#include <boost/optional.hpp>
#include <boost/bimap/bimap.hpp>
#include <boost/bimap/unordered_set_of.hpp>
#include <boost/bimap/unordered_multiset_of.hpp>
#include <boost/bimap/vector_of.hpp>
#include <boost/bimap/multiset_of.hpp>
#include <boost/bimap/support/lambda.hpp>
#include <boost/timer/timer.hpp>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/identity.hpp>
#include <boost/multi_index/sequenced_index.hpp>
#include <boost/random/random_device.hpp>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/functional/hash/hash.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/unordered_set.hpp>
#include <boost/serialization/list.hpp>
#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/info_parser.hpp>
#include <boost/optional.hpp>
#include <thread>

#include "CustomTypedefs.hpp"
#include "Types.hpp"
#include "TammberTypes.hpp"
#include "Constants.hpp"
#include "TammberModel.hpp"
#include "Log.hpp"


class ModelWrapper {
public:
ModelWrapper(boost::property_tree::ptree &config) {
	LOGGER("TammberModelWrapper")

	start=std::chrono::high_resolution_clock::now();
	carryOverTime=0;
	jobcount = 0;

	bool rs=false;
	if(boost::filesystem::exists("./TammberModel.chk")) {
		rs=true;
		std::ifstream ifs("TammberModel.chk");
		boost::archive::text_iarchive ia(ifs);
		// read class state from archive
		ia >> *this;
	}

	deleteVertex = config.get<uint64_t>("Configuration.MarkovModel.DeleteVertex",0);
	modelStr = config.get<bool>("Configuration.MarkovModel.ModelStr",true);

	// initialize the model
	markovModel.initialize(config,rs);

};


template<class Archive>
void save(Archive & ar, const unsigned int version) const {
	// note, version is always the latest when saving
	long int co=std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now()-start).count()+carryOverTime;
	ar & jobcount;
	ar & co;
	ar & markovModel;
};

template<class Archive>
void load(Archive & ar, const unsigned int version){
	ar & jobcount;
	ar & carryOverTime;
	ar & markovModel;
};
BOOST_SERIALIZATION_SPLIT_MEMBER()


virtual void full_print() {
	LOGGERA("PullMMbuilder::full_print()")
	// Save XML file
	if(deleteVertex!=0) {
		LOGGERA("PullMMbuilder::full_print DELETING "<<deleteVertex)
		markovModel.deleteVertex(deleteVertex);
		{
			std::ofstream ofs("./TammberModelNew.chk");
			{
				boost::archive::text_oarchive oa(ofs);
				oa << *this;
			}
		}
	}
	// Print analysis to screen
	std::cout<<markovModel.info_str(true)<<std::endl;
	std::list<TADjob> jobs;
	markovModel.generateTADs(jobs,100,true);
	if(modelStr) markovModel.write_model("MarkovModel.xml");
};

virtual void add_pathway(NEBPathway path) {
	markovModel.add_pathway(path);
};

virtual void save() {
	LOGGERA("PullMMbuilder::save()")
	{
		std::ofstream ofs("./TammberModelNewNEBS.chk");
		{
			boost::archive::text_oarchive oa(ofs);
			oa << *this;
		}
	}
};

protected:

TammberModel markovModel;
unsigned long carryOverTime;
unsigned long jobcount;
Label deleteVertex;
bool modelStr;
std::chrono::high_resolution_clock::time_point start;
};

#endif
