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



#ifndef NodeManager_h
#define NodeManager_h

#include <unistd.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <set>
#include <unordered_set>
#include <chrono>
#include <string>
#include <exception>
#include <iostream>

#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>
#include <boost/serialization/list.hpp>
#include <boost/functional/hash/hash.hpp>
#include <boost/filesystem.hpp>
#include <thread>
#include <chrono>

#include <mpi.h>


#include "AbstractSystem.hpp"
#include "Task.hpp"
#include "TaskManager.hpp"
#include "LocalStore.hpp"
#include "DDS.hpp"
#include "Types.hpp"
#include "Constants.hpp"
#include "Pack.hpp"
#include "TaskQueue.hpp"



template <class Bundle, class DriverHandle, class WMTaskMapper> class WorkManager {

public:

WorkManager(MPI_Comm comm_, int parent_, AbstractDDS *sharedStore_, AbstractLocalDataStore *localStore_, boost::property_tree::ptree &config, std::set<MPI_Comm> workerComms){
	comm=comm_;
	parent=parent_;
	MPI_Comm_rank(comm,&rank);
	die=false;
	everybodyIsDead=false;

	completedTasks.initialize(config,rank);

	//maximum size of a TASK_TAG message

	defaultFlavor=0;

	//TODO: move these away
	rbufSize=config.get<int>("Configuration.MaximumTaskMesgSize",1000000);

	start=std::chrono::high_resolution_clock::now();
	runTime=std::chrono::minutes( config.get<unsigned>("Configuration.RunTime",1000000) );


	sharedStore=sharedStore_;
	localStore=localStore_;
	uint64_t localBufferSize=config.get<uint64_t>("Configuration.WorkManager.LocalCacheSize",1000000000);
	uint64_t sharedBufferSize=config.get<uint64_t>("Configuration.WorkManager.SharedCacheSize",1000000000);
	localStore->setMaximumSize(localBufferSize);
	sharedStore->setMaximumBufferSize(sharedBufferSize);



	lastReport = std::chrono::high_resolution_clock::now();
	reportDelay = std::chrono::milliseconds( config.get<unsigned>("Configuration.WorkManager.ReportDelay",100000) );
	timeout = config.get<int>("Configuration.WorkManager.Timeout",2000);

	std::set<int> clients;

	for(auto it=workerComms.begin(); it!=workerComms.end(); it++) {
		workers.push_back( DriverHandle( *it ) );
		tasks.push_back(GenericTask());
		taskDescriptors.push_back(TaskDescriptor());
		clients.insert(workers.size()-1);
		reservationIds.push_back(1);
	}

	nWorkers=workers.size();

	//std::cout<<"N WORKERS: "<<nWorkers<<std::endl;

	batchId=1;

	//initialize communications
	recvCompleted=0;
	sendCompleted=1;

	rbuf=std::vector<RawData>(rbufSize,' ');
	//post a non-blocking recv for new tasks
	MPI_Irecv(&(rbuf[0]),rbufSize,MPI_BYTE,MPI_ANY_SOURCE,TASK_REQUEST_TAG,comm,&incomingTaskRequest);

	pendingTaskRequest=false;
};




~WorkManager(){
};


/*
* Manager server
*
* Manager will directly:
* -receive tasks from parent
* -send segments to parent
* -send task requests to workers
* -receive completed tasks from workers
*
* Node manager will indirectly:
* -request states to the database
* -receive states from the database
*/
void server(){


	while(true) {

		if(std::chrono::high_resolution_clock::now() - start > runTime) {
			//std::cout<<"WM: TIME TO DIE"<<std::endl;
			die=true;
		}

		//process non-blocking receives
		timer.start("Recv");
		processRecv();
		timer.stop("Recv");

		//process non-blocking sends
		timer.start("Send");
		processSend();
		timer.stop("Send");

		//promote prefetches to ready state
		timer.start("Promote");
		promoteTasks();
		timer.stop("Promote");

		//release reservations that are no longer needed
		timer.start("Release");
		releaseStaleReservations();
		timer.stop("Release");

		//process completed tasks
		timer.start("Process");
		processCompletedTasks();
		timer.stop("Process");

		//assign new tasks to idle workers
		timer.start("Assign");
		assignTasks();
		timer.stop("Assign");

		//update the db's
		timer.start("UpdateDB");
		updateStores();
		timer.stop("UpdateDB");

		report();


		/*
		   if(die) {
		        std::cout<<sendCompleted<<" "<<recvCompleted<<" "<<everybodyIsDead<<std::endl;
		   }
		 */

		//if(sendCompleted and recvCompleted and everybodyIsDead ) {
		if( everybodyIsDead ) {
			std::cout<<"WM: DONE"<<std::endl;
			break;
		}
	}
};


void processSend(){


	//test for send completion
	if(not sendCompleted) {
		//std::cout<<"WORK MANAGER TESTING FOR SEND"<<std::endl;
		MPI_Test(&outgoingRequest,&sendCompleted,&sendStatus);
		//std::cout<<"WORK MANAGER DONE TESTING FOR SEND "<<sendCompleted<<std::endl;
	}

	if(not sendCompleted and die) {
		MPI_Cancel(&outgoingRequest);
	}

	//TODO use a timer to avoid sending too often
	if(sendCompleted and not die) {
		//std::cout<<"WORK MANAGER: SENDING"<<std::endl;

		std::multimap<std::string,DataItem> done;
		completedTasks.release(done);

		int nTaskRequests=0;


		if( (not pendingTaskRequest)and (workQueue.count()+prefetch.count() < nWorkers) ) {
			nTaskRequests=nWorkers;
			pendingTaskRequest=true;
		}

		if(done.size()>0 || nTaskRequests>0) {
			timer.start("Send-Mesg");
			pack(sbuf,done,nTaskRequests);
			int count=sbuf.size();
			assert(count<=rbufSize);
			//std::cout<<"WORK MANAGER: processSend SENDING RESULTS "<<done.size()<<" "<< nTaskRequests <<" "<<pendingTaskRequest<<" "<<workQueue.count()<<" "<<prefetch.count()<<std::endl;
			MPI_Issend( &(sbuf[0]),count,MPI_BYTE,parent,COMPLETED_TASK_TAG, comm, &outgoingRequest);
			sendCompleted=0;
			timer.stop("Send-Mesg");
		}
	}
};

/**
 *
 */
void processRecv(){


	//did we receive new tasks?
	if(not recvCompleted) {
		//std::cout<<"WORK MANAGER TESTING FOR RECV"<<std::endl;
		MPI_Test(&incomingTaskRequest,&recvCompleted,&recvStatus);
		//std::cout<<"WORK MANAGER DONE TESTING FOR RECV "<<recvCompleted<<std::endl;
	}

	if(not recvCompleted and die) {
		MPI_Cancel(&incomingTaskRequest);
	}


	if(recvCompleted and not die) {
		timer.start("Recv-Mesg");
		Timer t;

		//extract
		int count=0;
		MPI_Get_count (&recvStatus, MPI_BYTE, &count);
		std::vector<char> bb(rbuf.begin(), rbuf.begin() + count);
		//std::cout<<"WORK MANAGER "<<rank<<" RECEIVED TASKS "<<count<<" "<<boost::hash_value(bb)<<std::endl;


		TaskDescriptorBundle incomingTasks;
		//std::cout<<"WORK MANAGER "<<rank<<" UNPACKING "<<count<<std::endl;
		unpack(rbuf,incomingTasks,std::size_t(count));
		//std::cout<<"INCOMING: "<<rank<<" "<<incomingTasks.tasks.size()<<" "<<incomingTasks.count()<<std::endl;


		//std::cout<<"PREFETCH: "<<prefetch.count()<<std::endl;
		prefetch.transferFrom(incomingTasks);
		//std::cout<<"PREFETCH: "<<rank<<" "<<prefetch.count()<<std::endl;


		//issue reservations for the data we need to fulfill these tasks
		std::unordered_set<DataPlaceholder>  requiredData;
		prefetch.requiredData(requiredData);
		for(auto it=requiredData.begin(); it!=requiredData.end(); it++) {
			if(activeReservations.count(std::make_pair(it->location,it->label))==0) {
				activeReservations[std::make_pair(it->location,it->label)]=sharedStore->reserve(it->location,it->label);
				std::cout<<"WORK MANAGER "<<rank<<" ISSUING RESERVATION FOR DATA ITEM "<<it->location<<" "<<it->label<<std::endl;
			}
		}
		batchId++;
		//std::cout<<"WORK MANAGER: DONE RECEIVING "<<t.stop()<<std::endl;

		MPI_Irecv(&(rbuf[0]),rbufSize,MPI_BYTE,MPI_ANY_SOURCE,TASK_REQUEST_TAG,comm,&incomingTaskRequest);
		recvCompleted=0;
		pendingTaskRequest=false;
		timer.stop("Recv-Mesg");
	}

};


/**
 * Release stale reservations
 */
void releaseStaleReservations(){

	Timer t;
	std::unordered_set<DataPlaceholder> required;
	//std::cout<<"REQUIRED DATA"<<std::endl;

	prefetch.requiredData(required);
	workQueue.requiredData(required);
	backfillQueue.requiredData(required);


/*
        std::cout<<"ACTIVE"<<std::endl;
        for(auto it=activeReservations.begin(); it!=activeReservations.end(); it++ ) {
                std::cout<<it->first.second<<" "<<it->first.first<<std::endl;
        }
        std::cout<<"REQUIRED"<<std::endl;
        for(auto it=required.begin(); it!=required.end(); it++) {
                DataPlaceholder s;
                s.necessity=NECESSITY::REQUIRED;
                s.label=it->label;
                s.location=it->location;
                s.shared=true;
                std::cout<<"it: "<<it->label<<" "<<it->location<<" "<<(it->necessity==NECESSITY::REQUIRED)<<" "<<required.count(s)<<std::endl;
                std::cout<<"s: "<<s.label<<" "<<s.location<<" "<<(s.necessity==NECESSITY::REQUIRED)<<" "<<(s==*it)<<std::endl;

        }
 */


	for(auto it=activeReservations.begin(); it!=activeReservations.end(); ) {
		DataPlaceholder s;
		s.necessity=NECESSITY::REQUIRED;
		s.label=it->first.second;
		s.location=it->first.first;
		s.shared=true;
		if(required.count(s)==0) {
			//std::cout<<"WORK MANAGER "<<rank<<" RELEASING RESERVATION ON ITEM "<<s.location<<" "<<s.label<<std::endl;
			sharedStore->release(it->second);
			it=activeReservations.erase(it);
		}
		else{
			it++;
		}
	}
	//std::cout<<"NODE MANAGER: DONE RELEASING: "<<t.stop()<<std::endl;
	localStore->purge();
};






void assignTasks(){
	Timer t;
	everybodyIsDead=true;

	for(int i=0; i<nWorkers; i++) {
		if(workers[i].isIdle() and not workers[i].isDead()) {
			if(die) {
				TaskDescriptor td;
				td.type = mapper.type("TASK_DIE");
				taskDescriptors[i]=td;
				tasks[i]=promoteTaskDescriptor(td);
				workers[i].assign(tasks[i]);
				workers[i].die();
				std::cout<<"KILLED WORKER "<<i<<" "<<workers[i].isDead()<<std::endl;
			}
			else {
				//std::cout<<"WORKER "<<i<<" HERE "<<workQueue.size()<<" "<<backfillQueue.size()<<std::endl;
				TaskDescriptor td;

				// always put non

				if(workQueue.pop(td) || backfillQueue.pop(td) ) {
					taskDescriptors[i]=td;
					tasks[i]=promoteTaskDescriptor(td);
					workers[i].assign(tasks[i]);
					if(tasks[i].imposeOrdering) {
						completedTasks.reserve(i,reservationIds[i]);
					}
					if(td.optional) {
						backfillQueue.push(td);
						backfillQueue.resize(nWorkers);
					}
					taskTimer[i].start();
					std::cout<<"WORK MANAGER "<<rank<<" ASSIGNED TASK TYPE "<< mapper.type(tasks[i].type) <<" TO WORKER "<<i<<" "<<tasks[i].arguments.size()<<" "<<td.arguments.size()<<" "<<tasks[i].imposeOrdering<<" "<<tasks[i].optional<<" "<<tasks[i].id<< std::endl;
				}
			}
		}
		everybodyIsDead = everybodyIsDead and workers[i].isDead();
		//std::cout<<"DEAD? "<<i<<" "<< workers[i].isDead()<<std::endl;
	}
	//std::cout<<"NODE MANAGER: DONE ASSIGNING: "<<t.stop()<<" "<<everybodyIsDead<<std::endl;
}

void extractData(GenericTask &t){
	//std::cout<<"WORK MANAGER "<<rank<<" EXTRACTING DATA "<<std::endl;

	for(auto it=t.outputData.begin(); it!=t.outputData.end(); it++) {
		std::cout<<"EXTRACTING: "<<it->first<<" "<<it->second.location<<" "<<it->second.label<<" "<<it->second.shared<<" "<<it->second.data.size()<<std::endl;
		if(it->second.shared==true) {
			std::cout<<"SHARED"<<std::endl;
			sharedStore->put(it->second.location,it->second.label,it->second.data);
		}
		else{
			//std::cout<<"LOCAL"<<std::endl;
			localStore->put(it->second.location,it->second.label,it->second.data);
		}

	}
};

void populateData(TaskDescriptor &td, GenericTask &t){
	t.inputData.clear();

	//std::cout<<"POPULATING DATA"<<std::endl;
	//loop over systems in the task descriptor

	for(auto it=td.inputData.begin(); it!=td.inputData.end(); it++) {
		StoredDataItem data(it->second);
		bool dataIsAvailable=false;
		//std::cout<<it->second.location<<" "<<it->second.label<<std::endl;
		//std::cout<<data.location<<" "<<data.label<<std::endl;

		if(data.shared==true) {
			if(sharedStore->count(data.location,data.label)>0) {
				sharedStore->get(data.location,data.label,data.data);
				dataIsAvailable=true;
			}
		}
		else{
			if(localStore->count(data.location,data.label)>0) {
				localStore->get(data.location,data.label,data.data);
				dataIsAvailable=true;
			}
		}
		//std::cout<<"POPULATING: "<<it->first<<" "<<data.location<<" "<<data.label<<" "<<data.shared<<" "<<data.data.size()<<" "<<dataIsAvailable<<std::endl;


		//std::cout<<data.location<<" "<<data.label<<" "<<data.tag<<" "<<data.shared<<" "<<it->shared<<" "<<dataIsAvailable<<std::endl;
		if(dataIsAvailable) {
			insert(it->first,t.inputData,data);
		}
	}
}

/**
 * Promote prefetches to ready state
 */
void promoteTasks(){


	//std::cout<<"NODE MANAGER:  PROMOTING: "<<std::endl;
	//Timer t;

	//promote required tasks
	for(auto it=prefetch.tasks.begin(); it!=prefetch.tasks.end(); ) {
		std::unordered_set<DataPlaceholder> deps;
		it->requiredData(deps);
		bool ready=true;
		for(auto itd=deps.begin(); itd!=deps.end(); itd++) {
			ready= ready and  (sharedStore->count(itd->location,itd->label) > 0 );
			/*
			#ifdef VERBOSE
			std::cout<<"PROMOTION: "<<itd->label<<" "<<itd->location<<" "<<sharedStore->count(itd->location,itd->label)<<" "<<ready<<std::endl;
			#endif
			*/
			if(not ready) {
				break;
			}
		}
		if(ready) {
			/*
			#ifdef VERBOSE
			std::cout<<"PROMOTION DONE "<<" "<<it->id<<std::endl;
			#endif
			*/
			workQueue.push(*it);
			it=prefetch.tasks.erase(it);
		}
		else{
			it++;
		}
	}

	//std::cout<<"NODE MANAGER: DONE PROMOTING: "<<t.stop()<<std::endl;
};


GenericTask promoteTaskDescriptor(TaskDescriptor &td){
	GenericTask t(td);
	t.arguments=td.arguments;
	populateData(td,t);
	return t;
};

void processCompletedTasks(){
	Timer t;

	for(int i=0; i<nWorkers; i++) {
		if( (not workers[i].isDead()) and workers[i].probe(tasks[i])) {
			//std::cout<<"WORK MANAGER "<<rank<<" PROCESSING COMPLETED TASK FROM WORKER "<<i<<" "<<tasks[i].type<<std::endl;
			processCompletedTask(i,tasks[i]);
			//std::cout<<"WORK MANAGER "<<rank<<" PROCESSED COMPLETED TASK FROM WORKER "<<i<<std::endl;
		}

		if( (workers[i].isDead() and deadWorkers.count(i)==0)or ( deadWorkers.count(i)==0 and taskTimer[i].lapInSeconds()>timeout  ) ) {
			//this worker just died. We won't get a reply.
			#ifdef USE_BOOST_LOG
			boost::log::sources::severity_logger< boost::log::trivial::severity_level > lg;
			BOOST_LOG_SEV(lg, boost::log::trivial::error) <<"WORKER "<<i<<" DIED";
			#else
			std::cout<<"WORKER "<<i<<" DIED"<<std::endl;
			#endif


			//release the reservation
			completedTasks.release(i,reservationIds[i]);

			//put task back in the queue
			prefetch.insert(taskDescriptors[i]);

			deadWorkers.insert(i);
			workers[i].die();
		}
	}
	//std::cout<<"NODE MANAGER: DONE PROCESSING: "<<t.stop()<<std::endl;
};

void processCompletedTask(int workerID, GenericTask &t){
	//uint64_t taskID=t.id;

	//extract the data
	extractData(t);

	if(t.failed) {
		//release the reservation
		completedTasks.release(workerID,reservationIds[workerID]);

		//put task back in the queue
		prefetch.insert(taskDescriptors[workerID]);

		//pretend that the worker is dead
		deadWorkers.insert(workerID);

		//pretend that the worker is dead
		workers[workerID].die();
		#ifdef USE_BOOST_LOG
		boost::log::sources::severity_logger< boost::log::trivial::severity_level > lg;
		BOOST_LOG_SEV(lg, boost::log::trivial::error) <<"WORKER "<<workerID<<" FAILED";
		#else
		std::cout<<"WORKER "<<workerID<<" FAILED"<<std::endl;
		#endif
	}
	else{
		//process the results
		completedTasks.checkin(workerID,reservationIds[workerID],t);
	}
	timer.push("Task",taskTimer[workerID].stop());
};

void updateStores(){
	Timer t;
	int imax=1000;
	int i=0;
	while(sharedStore->singleServe()>0 and i<imax) {
		i++;
	}
	//std::cout<<"NODE MANAGER: DONE UPDATING STORES: "<<t.stop()<<std::endl;
}


virtual void report(){

	if(std::chrono::high_resolution_clock::now() - lastReport> reportDelay  ) {
		timer.report();
		completedTasks.report();
		#ifdef USE_BOOST_LOG
		boost::log::sources::severity_logger< boost::log::trivial::severity_level > lg;
		BOOST_LOG_SEV(lg, boost::log::trivial::trace) <<workQueue.count()<<" TASKS IN WORK QUEUE";
		BOOST_LOG_SEV(lg, boost::log::trivial::trace) <<backfillQueue.count()<<" TASKS IN BACKFILL QUEUE";
		#else
		std::cout<<workQueue.count()<<" TASKS IN WORK QUEUE"<<std::endl;
		std::cout<<backfillQueue.count()<<" TASKS IN BACKFILL QUEUE"<<std::endl;
		#endif
		lastReport=std::chrono::high_resolution_clock::now();
	}
};



private:
int rank;
int parent;

uint64_t batchId;
std::chrono::high_resolution_clock::time_point lastReport;
std::chrono::milliseconds reportDelay;
std::chrono::high_resolution_clock::time_point start;
std::chrono::minutes runTime;

//Communication variables
MPI_Comm comm;
int recvCompleted;
int sendCompleted;
MPI_Status recvStatus;
MPI_Status sendStatus;
MPI_Request incomingTaskRequest;
MPI_Request outgoingRequest;

unsigned int validatorTimeout;
int timeout;

int nWorkers;
bool everybodyIsDead;
//int nRanksDriver;


std::vector<GenericTask> tasks;
std::vector<TaskDescriptor> taskDescriptors;
std::vector<DriverHandle> workers;
std::vector<Label> reservationIds;
std::set<int> deadWorkers;

//task bundle that contains task requests that are in prefetch state
//TaskDescriptorBundle prefetches;

TaskDescriptorBundle prefetch;

WMTaskMapper mapper;

//store the results of completed tasks
ResultManager<Bundle> completedTasks;

//fifo queue for tasks that are in a ready state
FIFOQueue workQueue;
FIFOQueue backfillQueue;


//active reservations to the data-store
std::map< std::pair<unsigned int, Label>, Label> activeReservations;

AbstractDDS *sharedStore;
AbstractLocalDataStore *localStore;

TransitionStatistics statistics;

std::vector<char> rbuf;
std::vector<char> sbuf;

int rbufSize;
bool die;
int defaultFlavor;

bool pendingTaskRequest;


ExecutionTimer timer;


std::map<int, Timer> taskTimer;

};





#endif 
