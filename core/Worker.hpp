


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




#include <chrono>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/info_parser.hpp>






#include <stdio.h>
#include <vector>
#include <random>

#include <mpi.h>
#include <unistd.h>
#include <iostream>
#include <chrono>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/info_parser.hpp>


#include "Task.hpp"
#include "Log.hpp"
#include "Constants.hpp"
#include "Pack.hpp"
#include "CustomTypedefs.hpp"


void worker(MPI_Comm localComm, MPI_Comm interComm, int seed){

	TaskMapperType TaskMapper;


	// Create empty property tree object
	boost::property_tree::ptree tree;
	// Parse the XML into the property tree.
	boost::property_tree::read_xml("./input/ps-config.xml", tree, boost::property_tree::xml_parser::no_comments);

	std::chrono::high_resolution_clock::time_point start=std::chrono::high_resolution_clock::now();
	std::chrono::minutes runTime=std::chrono::minutes( tree.get<unsigned>("Configuration.RunTime",1000000)+1 );


	std::default_random_engine generator(seed);
	std::bernoulli_distribution distribution(0.0001);


	{
		DriverTaskManagerType taskManager(interComm,localComm);
		GenericTask t;
		GenericTask tt=t;
		tt.clearInputs();
		tt.clearOutputs();
		//EngineType engine(tree,localComm,seed);
		std::shared_ptr<EngineType> engine=std::make_shared<EngineType>(tree,localComm,seed);
		bool healthy=true;

		while(true) {
			//receive a task

			if(healthy) {
				if(std::chrono::high_resolution_clock::now() - start > runTime) {
					//break;
				}
				LOGGERA("WAITING FOR TASK")
				auto t0 = std::chrono::high_resolution_clock::now();
				healthy=taskManager.pullTask(t);
				auto t1 = std::chrono::high_resolution_clock::now();
				auto d = std::chrono::duration_cast<std::chrono::milliseconds>(t1-t0).count();
				LOGGERA("PROCESSING TASK "<<t.type<<" ")
			}
			if(t.type == TaskMapper.type("TASK_DIE") || not healthy ) {
				//TODO: add cleanup
				break;
			}


#ifdef INJECT_FAULT

			{
				if (distribution(generator)) {
					std::cout<<"SUICIDE!!!"<<std::endl;
					*(char *)0 = 0;
					//__builtin_trap();
				}
			}

#endif

			//process the task
			try{
				{
					auto t0 = std::chrono::high_resolution_clock::now();
					engine->process(t);
					auto t1 = std::chrono::high_resolution_clock::now();
					auto d = std::chrono::duration_cast<std::chrono::milliseconds>(t1-t0).count();

					if(engine->failed()) {
						tt.failed=true;
						auto t0 = std::chrono::high_resolution_clock::now();
						healthy=taskManager.pushTask(tt);
						auto t1 = std::chrono::high_resolution_clock::now();
						auto d = std::chrono::duration_cast<std::chrono::milliseconds>(t1-t0).count();
						LOGGERA("WORKER:  TASK EXECUTION FAILED")
						healthy=false;
					}

					LOGGER("WAITED "<<d<<" ms to process the task "<<t.type)
				}
				//send the results back
				if(healthy) {
					auto t0 = std::chrono::high_resolution_clock::now();
					healthy=taskManager.pushTask(t);
					auto t1 = std::chrono::high_resolution_clock::now();
					auto d = std::chrono::duration_cast<std::chrono::milliseconds>(t1-t0).count();
					if(!healthy) {
						break;
					}
					LOGGER("WAITED "<<d<<" ms for the push")
				}
			}
			catch(...) {
				LOGGERA("WORKER: ENGINE HAS THROWN EXCEPTION")

				if(healthy) {
					tt.failed=true;
					auto t0 = std::chrono::high_resolution_clock::now();
					healthy=taskManager.pushTask(tt);
					auto t1 = std::chrono::high_resolution_clock::now();
					auto d = std::chrono::duration_cast<std::chrono::milliseconds>(t1-t0).count();
					//boost::log::sources::severity_logger< boost::log::trivial::severity_level > lg;
					//BOOST_LOG_SEV(lg, boost::log::trivial::error) <<
					LOGGERA("WORKER:  TASK EXECUTION FAILED");
					healthy=false;
				}
				else{
					LOGGERA("WORKER:  BREAKING")
					break;
				}
				//recreate the engine
				//engine=std::make_shared<EngineType>(tree,localComm,seed);
			}
		}
	}
};
