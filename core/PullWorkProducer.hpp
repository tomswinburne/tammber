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


#ifndef __WorkProducer__
#define __WorkProducer__



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
#include <boost/optional.hpp>


#include "Types.hpp"
#include "Log.hpp"
#include "Task.hpp"
#include "Pack.hpp"
#include "Constants.hpp"
#include "DDS.hpp"
#include "CustomTypedefs.hpp"
#include "TaskManager.hpp"



class AbstractPullWorkProducer {
public:
AbstractPullWorkProducer(MPI_Comm comm_, AbstractDDS *sharedStore_,std::set<int> children_, boost::property_tree::ptree &config) {
	comm=comm_;
	sharedStore=sharedStore_;
	children=children_;
	int rank;
	MPI_Comm_rank(comm, &rank);

	completedTasks.initialize(config,rank);
	completedTasks.setBundleSize(1);
	start=std::chrono::high_resolution_clock::now();
	die=false;
	killed=false;
	runTime=std::chrono::minutes( config.get<unsigned>("Configuration.RunTime",1000000) );


	relay.openBatch();
	//TODO: FIX ME
	maxTaskMesgSize=10000000;

	sendCompleted=1;
	recvCompleted=0;
};


virtual void processCompletedTasks()=0;
virtual TaskDescriptorBundle generateTasks(int consumerID, int nTasks)=0;
virtual void checkpoint_impl()=0;
virtual void report_impl()=0;



~AbstractPullWorkProducer(){
};

virtual void initialize()=0;
//virtual void processTasks()=0;
//virtual void generateTasks()=0;
//virtual void checkpoint_impl()=0;
//virtual void report_impl()=0;


void server(){
	//BOOST_LOG_SEV(lg,debug) <<"#Splicer::server ";
	lastReport=std::chrono::high_resolution_clock::now();
	lastCheckpoint=std::chrono::high_resolution_clock::now();

	//TODO: do this right
	rbufSize=100000000;
	rbuf=std::vector<char>(rbufSize,' ');

	//post a  non-blocking receive
	MPI_Irecv(&(rbuf[0]),rbufSize,MPI_BYTE,MPI_ANY_SOURCE,MPI_ANY_TAG,comm,&incomingRequest);
	recvCompleted=0;


	initialize();

	while(true) {

		if(std::chrono::high_resolution_clock::now() - start > runTime  ) {
			LOGGERA("WP: TIME TO DIE")
			die=true;
		}

		//update databases
		timer.start("UpdateDB");
		updateStores();
		timer.stop("UpdateDB");

		//process non-blocking receives (completed segments from WM)
		timer.start("Recv");
		processRecv();
		timer.stop("Recv");

		//process non-blocking sends (tasks to WM)
		timer.start("Send");
		processSend();
		timer.stop("Send");

		//write outputs if needed
		timer.start("Report");
		report();
		timer.stop("Report");

		//checkpoint if needed
		checkpoint(false);


		//if(die and recvCompleted and sendCompleted) {
		if(die) {
			//force a checkpoint
			//checkpoint(true);
			break;
			LOGGERA("WP: DONE")
		}
	}
};

void updateStores(){
	int imax=1000;
	int i=0;
	while(sharedStore->singleServe()>0 and i<imax) {
		i++;
	}
};

void processRecv(){
	//did we receive new completed tasks?

	if(not recvCompleted) {
		MPI_Test(&incomingRequest,&recvCompleted,&recvStatus);
	}

	if(not recvCompleted and die) {
		MPI_Cancel(&incomingRequest);
	}

	if(recvCompleted and not die) {
		//timer.start("Recv-Mesg");
		int count=0;
		MPI_Get_count( &recvStatus, MPI_BYTE, &count );
		int from=recvStatus.MPI_SOURCE;

		int requests;
		std::multimap<std::string, DataItem> results;
		unpack(rbuf,results,requests,count);

		//assimilate the results
		timer.start("Recv-Assimilate");
		completedTasks.assimilate(results);
		timer.stop("Recv-Assimilate");

		//release completed results
		timer.start("Recv-Extract");
		completedTasks.extract(ready);
		timer.stop("Recv-Extract");

		timer.start("Recv-Process");
		processCompletedTasks();
		timer.stop("Recv-Process");

		if(requests>0) {
			taskRequests[from]+=requests;
		}

		LOGGER("SPLICER RECEIVED "<<requests<<" TASK REQUESTS FROM "<<from)

		//post a new receive for tasks
		MPI_Irecv(&(rbuf[0]),rbufSize,MPI_BYTE,MPI_ANY_SOURCE,MPI_ANY_TAG,comm,&incomingRequest);
		recvCompleted=0;
		timer.stop("Recv-Mesg");
	}

};



bool processSend(){
	//fulfill requests

	if(not sendCompleted) {
		MPI_Test(&outgoingRequest,&sendCompleted,&sendStatus);
	}

	//according to the MPI standard, TEST after CANCEL should eventually succeed.
	//in practice, it does not seem to be the case (see https://www.mail-archive.com/users@lists.open-mpi.org/msg30905.html)
	if(not sendCompleted and die) {
		MPI_Cancel(&outgoingRequest);
	}

	if(sendCompleted and not die) {
		if(not taskRequests.empty()) {

			// auto it=taskRequests.begin();
			for(auto it=taskRequests.begin(); it!=taskRequests.end(); ) {

				// DataItem d;
				timer.start("Generate");
				TaskDescriptorBundle tasks = generateTasks(it->first,it->second);

				timer.stop("Generate");
				bool deleteRequest=false;
				if(tasks.count()>0) {
					timer.start("Send-Mesg");
					pack(sbuf,tasks);

					//pack(d.data,tasks);
					//TaskDescriptorBundle tt;
					//unpack(sbuf,tt);

					LOGGER("WORK PRODUCER FULFILLING "<<it->second<<" TASK REQUESTS FROM "
					<<it->first<<" "<<tasks.count()<<" "<< sbuf.size()<<" "<<boost::hash_value(sbuf))


					MPI_Issend( &(sbuf[0]),sbuf.size(),MPI_BYTE,it->first,TASK_REQUEST_TAG, comm, &outgoingRequest);
					sendCompleted=0;
					//taskRequests.erase(it);
					deleteRequest=true;
					timer.stop("Send-Mesg");
					break;
				}
				if(deleteRequest){
					it=taskRequests.erase(it);
				}
				else{
					it++;
				}
			}
		}
	}

	return true;

};


void report(){
	//report at given intervals
	if(std::chrono::high_resolution_clock::now() - lastReport> reportDelay  ) {
		Timer t;
		timer.report();

		LOGGERA("PENDING TASK REQUESTS: ")
		for(auto tr : taskRequests) LOGGERA(tr.first<<" : "<<tr.second)

		report_impl();
		lastReport=std::chrono::high_resolution_clock::now();
	}
};

void checkpoint(bool now=false){
	//checkpoint at given intervals
	if(now or (std::chrono::high_resolution_clock::now() - lastCheckpoint > checkpointDelay)  ) {
		//checkpoint
		checkpoint_impl();
		lastCheckpoint=std::chrono::high_resolution_clock::now();
	}
};


protected:
AbstractDDS *sharedStore;

TaskMapperType mapper;

BundleType ready;
TaskDescriptorRelay relay;
ResultManager<BundleType> completedTasks;

std::map<int,int> taskRequests;
std::set<int> children;
MPI_Comm comm;
int recvCompleted;
int sendCompleted;
MPI_Status recvStatus;
MPI_Status sendStatus;
MPI_Request incomingRequest;
MPI_Request outgoingRequest;

std::vector<char> rbuf;
std::vector<char> sbuf;
int rbufSize;

unsigned maxTaskMesgSize;

std::vector<std::vector<char> > sBuf;

std::chrono::high_resolution_clock::time_point start;
std::chrono::milliseconds reportDelay;
std::chrono::milliseconds checkpointDelay;
std::chrono::minutes runTime;

std::chrono::high_resolution_clock::time_point lastReport;
std::chrono::high_resolution_clock::time_point lastCheckpoint;


bool die;
bool killed;

ExecutionTimer timer;


};



#endif
