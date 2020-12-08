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



#ifndef MPICache3_hpp
#define MPICache3_hpp



#include <iostream>
#include <ostream>
#include <stdio.h>
#include <set>
#include <map>
#include <vector>
#include <list>
#include <tuple>

#include <mpi.h>

#include "HCDS.hpp"
#include "Types.hpp"
#include "Constants.hpp"
#include "Pack.hpp"

#include "DDS.hpp"




class MsgChannel {
public:
MsgChannel( ){
	pending_=false;
	req=MPI_REQUEST_NULL;
};

bool pending(){
	return pending_;
};

void cancel(){
	MPI_Cancel(&req);
	MPI_Wait(&req,&status);
};

MPI_Request req;
MPI_Status status;
std::vector<char> buf;
bool pending_;

};



class RecvChannel : public MsgChannel {
public:


bool test(){

	if(!pending_) {
		return false;
	}


	int complete=0;
	MPI_Test(&req,&complete,&status);
	if(complete) {
		pending_=false;
	}
	//pending_=!bool(complete);
	//std::cout<<"TEST "<<complete<<" "<<pending_<<std::endl;
	return bool(complete);
};

bool postRecv(size_t maxMesgSize, int sourceTag, int msgTag, MPI_Comm &communicator){
	//std::cout<<"RECV PENDING "<<pending_<<std::endl;
	if(pending_) {
		return false;
	}
	buf.resize(maxMesgSize);
	MPI_Irecv(&(buf[0]),maxMesgSize,MPI_BYTE,sourceTag,msgTag,communicator,&req);
	pending_=true;
	//std::cout<<"RECV POSTED "<<pending_<<std::endl;
	return pending_;
};

};

class SendChannel : public MsgChannel {
public:
bool test(){



	int complete=0;
	MPI_Test(&req,&complete,&status);
	if(complete) {
		pending_=false;
	}
	//pending_=!bool(complete);
	//std::cout<<"TEST "<<complete<<" "<<pending_<<std::endl;
	return bool(complete);
};

bool postSend(size_t dataSize, RawDataVector data,int destinationTag, int msgTag, MPI_Comm &communicator){
	if(pending_) {
		std::cout<<"WARNING: MESSAGE WAS ALREADY PENDING"<<std::endl;
		return false;
	}

	buf.swap(data);
	MPI_Issend( &(buf[0]),dataSize,MPI_BYTE,destinationTag, msgTag, communicator, &req);

	pending_=true;
	//std::cout<<"SEND POSTED "<<destinationTag<<std::endl;
	return pending_;
};
};


template <class LocalStoreType> class HDDS3 : public AbstractDDS {

public:
/**
 * Initialize
 */
HDDS3(MPI_Comm comm_, int parent_, std::string homeDir, std::string baseName, uint64_t maxDataSize_, uint64_t maxCacheSize_, std::set<int> dbKeys, std::unordered_map< unsigned int, std::pair<bool,bool> > &dbAttributes, bool purge_){
	comm=comm_;
	//recvPending=false;
	//sendPending=false;
	purge=purge_;
	parent=parent_;
	maxDataSize=maxDataSize_;
	maxCacheSize=maxCacheSize_;

	memoryUsage=0;
	//maxMemoryUsage=maxCacheSize;


	//initialize a receive buffer
	//rbuf=std::vector<char>(maxDataSize);


	MPI_Comm_rank(comm, &self);

	//initialize the databases
	lstore.initialize(homeDir,baseName);
	for(auto it=dbKeys.begin(); it!=dbKeys.end(); it++) {
		bool allowDuplicates=dbAttributes[*it].first;
		bool eraseOnGet=dbAttributes[*it].second;
		lstore.createDatabase(*it,allowDuplicates,eraseOnGet);
	}


	reservationID=1;
	nActiveHolds=0;
};


/**
 * Server process. Does not return. Use in multithreaded scenarios.
 */

void serve(){
	while(true) {
		singleServe();
	}
};


void printStatus(){
	std::cout<<" HDDS on RANK: "<<self<<std::endl;
	std::cout<<self<<" "<<maxDataSize<<" "<<std::endl;
	std::cout<<self<<" MEMORY USAGE: "<<memoryUsage<<" "<<maxCacheSize<<std::endl;;

	std::cout<<self<<" HOLDS: "<<nActiveHolds<<std::endl;
	std::cout<<self<<" RESERVATION ID: "<<reservationID<<" "<<reservations.size()<<std::endl;

};

void setMaximumBufferSize(uint64_t maxSize){
	maxCacheSize=maxSize;
};


/**
 * Single step of the server loop
 */
int singleServe(){




	//post a receive (does nothing if none is pending)
	recvChannel.postRecv(maxDataSize,MPI_ANY_SOURCE,DDS_MESG_TAG,comm);

	//receive incoming messages
	if(recvChannel.test()) {
		int size;
		MPI_Get_count(&recvChannel.status,MPI_BYTE,&size);

		std::list<Transaction> transactions;
		unpack(recvChannel.buf,transactions,std::size_t(size));


		/*
		   std::cout<<"N TRANSACTIONS "<<transactions.size()<<std::endl;
		   for(auto it=transactions.begin(); it!=transactions.end(); it++) {
		        it->print();
		   }
		 */

		incomingTransactions.splice(incomingTransactions.begin(),transactions);



		recvChannel.postRecv(maxDataSize,MPI_ANY_SOURCE,DDS_MESG_TAG,comm);
		//std::cout<<"DATA RECEIVED "<<size<<std::endl;
		//std::cout<<recvChannel.pending_ <<std::endl;
	}


	//purge if using too much memory
	if(purge and memoryUsage>maxCacheSize) {
		//std::cout<<"RANK: "<<self<<" PURGING" <<" "<<memoryUsage<<" "<<maxCacheSize<<std::endl;
		//flush the least recently touched items. Make sure we don't flush elements that are still needed
		for(auto it=mru.begin(); it!=mru.end(); ) {
			unsigned int dbKey=it->first;
			Label key=it->second;

			bool onHold=requestCounter[dbKey][key].size()>0;


			if(not onHold) {
				std::cout<<"RANK: "<<self<<" PURGING ITEM "<<key<<std::endl;
				assert(lstore.count(dbKey,key)>0);
				int size=lstore.erase(dbKey,key);
				memoryUsage-=size;

				if(lstore.count(dbKey,key)==0) {
					it=mru.erase(it);
				}
				else{
					it++;
				}
				//std::cout<<"RANK: "<<self<<" PURGING ITEM "<<key<<std::endl;
			}
			else{
				it++;
			}

			if(memoryUsage<maxCacheSize) {
				break;
			}
		}
		if(memoryUsage>maxCacheSize) {
			std::cout<<"RANK: "<<self<<" WARNING: UNABLE TO PURGE ENOUGH DATA TO MEET MEMORY LIMIT "<<std::endl;
		}
	}





	//process an incoming transactions

	//NOTE: this is a coarse lock. Can we do better?
	//transactionLock.lock();
	if(incomingTransactions.size()>0) {
		auto it=incomingTransactions.begin();
		Transaction &transaction=*it;

		//std::cout<<"RANK: "<<self<<" INCOMING: "<<std::endl;
		//transaction.print();

		if(transaction.type==DDS_GET or transaction.type==DDS_PREFETCH) {
			//GET requests should only come from children
			assert(transaction.source!=parent);

			//bool transactionCompleted=false;
			bool inStore=(lstore.count(transaction.dbKey,transaction.key)>0);

			//std::cout<<"RANK: "<<self<<" GET OR PREFETCH INSTORE: "<<inStore<<" PENDING GETS: "<<pendingGets.count( std::make_pair(transaction.dbKey, transaction.key) )<<" "<<std::endl;


			//we have the item in store
			if(inStore) {
				//the source is a remote child. Push a PUT transaction in the queue
				if(transaction.source!=self) {
					Transaction reply;
					reply.type=DDS_PUT;
					reply.source=self;
					reply.destination=transaction.source;
					reply.dbKey=transaction.dbKey;
					reply.key=transaction.key;

					//push to the outgoing queue
					outgoingTransactions[reply.destination].push_back(reply);

					//Register a hold on this item until we process the corresponding DDS_PUT back to the requester
					if(transaction.type==DDS_GET and not transaction.pending) {
						requestCounter[transaction.dbKey][transaction.key].insert(transaction.source);
						nActiveHolds++;
						//std::cout<<"INSERTED A HOLD ON "<<transaction.dbKey<<" "<<transaction.key<<" "<<transaction.source<<std::endl;
					}
				}
				//update the MRU log
				mru.touch(std::make_pair(transaction.dbKey,transaction.key));

				//we are done with this transaction
				incomingTransactions.erase(it);
			}
			//we don't have the item
			else{
				//if we do not have a pending GET request on this item, we need to forward it to parent
				if( not transaction.pending and pendingGets.count( std::make_pair(transaction.dbKey, transaction.key) ) == 0 and parent>=0) {
					Transaction forward;
					forward.type=DDS_GET;
					forward.dbKey=transaction.dbKey;
					forward.key=transaction.key;
					forward.source=self;
					forward.destination=parent;
					outgoingTransactions[forward.destination].push_back(forward);
					//log that we already requested this item
					//this has to be cleared when we receive the item
					pendingGets.insert(std::make_pair(transaction.dbKey, transaction.key));
				}

				//Register a hold on this item until we process the corresponding DDS_PUT back to the requester
				if(transaction.type==DDS_GET and not transaction.pending and transaction.source!=self) {
					requestCounter[transaction.dbKey][transaction.key].insert(transaction.source);
					nActiveHolds++;

					//std::cout<<"INSERTED A HOLD ON "<<transaction.dbKey<<" "<<transaction.key<<" "<<transaction.source<<std::endl;
					transaction.pending=true;


					//we could not fully fulfill this transaction; mark it as pending and move it to the pending set
					pendingTransactions[ std::make_pair(transaction.dbKey, transaction.key) ].push_back(transaction);
				}

				incomingTransactions.erase(it);
			}
		}


		if(transaction.type==DDS_PUT) {
			bool inStore=(lstore.count(transaction.dbKey, transaction.key)>0);

			//std::cout<<"RANK: "<<self<<" PROCESSING INCOMING PUT "<<transaction.key<<" "<<inStore<<" "<<std::endl;

			//add it to the store if we don't have it already
			if( not inStore ) {
				//std::cout<<"RANK: "<<self<<" ADDING  "<<transaction.dbKey<<" "<<transaction.key<<" TO STORE "<<std::endl;
				lstore.put(transaction.dbKey,transaction.key,transaction.data);
				memoryUsage+=transaction.data.size();

				//reinject pending transactions waiting on this item
				if(pendingTransactions.count( std::make_pair(transaction.dbKey, transaction.key) ) > 0 ) {
					incomingTransactions.splice(incomingTransactions.begin(), pendingTransactions[ std::make_pair(transaction.dbKey, transaction.key) ] );
					pendingTransactions.erase( std::make_pair(transaction.dbKey, transaction.key) );
				}
			}

			//parent fulfilled our request
			if(transaction.source==parent) {
				pendingGets.erase( std::make_pair(transaction.dbKey, transaction.key) );
			}
			//Request from child. Foward item to parent only if we never did before
			else{
				if(forwardedPuts.count(std::make_pair(transaction.dbKey, transaction.key))==0 and parent>=0) {
					Transaction forward;
					forward.type=DDS_PUT;
					forward.dbKey=transaction.dbKey;
					forward.key=transaction.key;
					forward.source=self;
					forward.destination=parent;
					forward.data=transaction.data;

					//Register a hold on this item until we process the corresponding DDS_PUT to parent
					requestCounter[forward.dbKey][forward.key].insert(parent);
					//std::cout<<"INSERTED A HOLD ON "<<forward.dbKey<<" "<<forward.key<<" "<<parent<<std::endl;
					nActiveHolds++;
					outgoingTransactions[forward.destination].push_back(forward);
				}
			}
			//no need to forward this item to parent again
			forwardedPuts.insert(std::make_pair(transaction.dbKey, transaction.key));

			//we are done with this transaction
			incomingTransactions.erase(it);

			//update the MRU log
			mru.touch(std::make_pair(transaction.dbKey,transaction.key));
		}
	}

	//transactionLock.unlock();



	//process outgoing transactions
	for(auto it=outgoingTransactions.begin(); it!=outgoingTransactions.end(); ) {

		bool pending=!sendChannels[it->first].test();

		//initiate send of transactions
		if(not pending) {
			//std::cout<<"OUTGOING TRANSACTIONS TO "<<it->first<<std::endl;
			//gather the data
			std::list<Transaction> outgoing;
			int nt=0;
			int ntMax=10;
			while(nt<=ntMax and it->second.size()>0) {
				Transaction &transaction=it->second.front();

				if(transaction.type==DDS_PUT) {

					std::cout<<DDS_PUT<<" "<<transaction.destination<<" "<<transaction.dbKey<<" "<<transaction.key<<std::endl;
					bool inStore=(lstore.count(transaction.dbKey,transaction.key)>0);
					//this item should be in store
					assert(inStore);

					lstore.get(transaction.dbKey,transaction.key,transaction.data);

					//we have fulfilled the request. Release the hold
					assert( requestCounter[transaction.dbKey][transaction.key].count(transaction.destination) > 0 );
					auto itr=requestCounter[transaction.dbKey][transaction.key].find(transaction.destination);
					requestCounter[transaction.dbKey][transaction.key].erase(itr);
					//std::cout<<"RELEASED A HOLD ON "<<transaction.dbKey<<" "<<transaction.key<<" "<<transaction.destination<<std::endl;
					nActiveHolds--;

					outgoing.push_back(transaction);
				}

				if(transaction.type==DDS_GET) {
					//std::cout<<DDS_GET<<" "<<transaction.destination<<" "<<transaction.dbKey<<" "<<transaction.key<<std::endl;
					outgoing.push_back(transaction);
				}

				nt++;
				//done with this
				it->second.pop_front();
			}

			//initiate send of the data
			sbuf.clear();
			pack(sbuf,outgoing);



			int count=sbuf.size();

			bool ret=sendChannels[it->first].postSend(count,sbuf,it->first,DDS_MESG_TAG,comm);
			//std::cout<<"Sending "<<count<<" TO "<<it->first<<" "<<ret<<std::endl;
		}

		if(it->second.size()==0) {
			it = outgoingTransactions.erase(it);
		}
		else{
			it++;
		}


	}




	return incomingTransactions.size();

};

virtual void prefetch(int dbKey, Label key){
	Transaction transaction;
	transaction.type=DDS_PREFETCH;
	transaction.dbKey=dbKey;
	transaction.key=key;
	transaction.source=self;
	transaction.destination=self;

	transactionLock.lock();
	incomingTransactions.push_back(transaction);
	transactionLock.unlock();
};

//a returned reservationID of 0 means the reservation was not accepted
virtual uint64_t reserve(int dbKey, Label key){

	if(memoryUsage>maxCacheSize) {
		//reject the reservation until we have more room available
		//WARNING: disabled since the rest of the code does not deal with this yet
		//return 0;
	}

	Transaction transaction;
	transaction.type=DDS_GET;
	transaction.dbKey=dbKey;
	transaction.key=key;
	transaction.source=self;
	transaction.destination=self;

	transactionLock.lock();
	incomingTransactions.push_back(transaction);
	transactionLock.unlock();


	//Register a hold on this item until the requester calls release
	requestCounter[dbKey][key].insert(self);
	//Register the reservation
	reservations[reservationID]=std::make_pair(dbKey,key);
	nActiveHolds++;

	return reservationID++;
};

virtual int count(int dbKey, Label key){
	return lstore.count(dbKey,key);
};


virtual bool release(uint64_t id){
	if(id==0) {
		return false;
	}

	//make sure the reservation exists
	if(reservations.count(id)==0) {
		return false;
	}
	unsigned int dbKey=reservations[id].first;
	Label key=reservations[id].second;

	if(requestCounter[dbKey][key].count(self)==0) {
		return false;
	}

	//clear the hold
	assert( requestCounter[dbKey][key].count(self)!=0 );
	auto it=requestCounter[dbKey][key].find(self);
	requestCounter[dbKey][key].erase(it);
	nActiveHolds--;

	//release the reservation
	reservations.erase(id);

	return true;
};

virtual void sync(){
	lstore.sync();
}

virtual bool get(int dbKey, Label key,RawDataVector &data){
#ifdef LOG_DDS
	boost::log::sources::severity_logger< boost::log::trivial::severity_level > lg;
	BOOST_LOG_SEV(lg, boost::log::trivial::trace) <<" HDDS GET "<<dbKey<<" "<<key;
#endif
	if(count(dbKey, key)==0) {
		return false;
	}
	lstore.get(dbKey,key,data);
	return true;
};

virtual void put(int dbKey, Label key, RawDataVector &data){
#ifdef LOG_DDS
	boost::log::sources::severity_logger< boost::log::trivial::severity_level > lg;
	BOOST_LOG_SEV(lg, boost::log::trivial::trace) <<" HDDS PUT "<<dbKey<<" "<<key;
#endif
	Transaction transaction;
	transaction.type=DDS_PUT;
	transaction.dbKey=dbKey;
	transaction.key=key;
	transaction.data=data;
	transaction.source=self;
	transaction.destination=self;

	transactionLock.lock();
	incomingTransactions.push_back(transaction);
	transactionLock.unlock();
};


virtual std::set<uint64_t> availableKeys(unsigned int dbKey){
	return lstore.availableKeys(dbKey);
};


virtual void cancelCommunications(){
	recvChannel.cancel();
	for(auto it=sendChannels.begin(); it!=sendChannels.end(); it++) {
		it->second.cancel();
	}
}

private:

bool purge;
int parent;
int self;
uint64_t maxDataSize;
uint64_t maxCacheSize;
uint64_t memoryUsage;
//uint64_t maxMemoryUsage;
//bool recvPending;
//bool sendPending;
MPI_Comm comm;
//MPI_Comm dbComm;
//MPI_Status recvStatus;
//MPI_Request recvReq;
//MPI_Status sendStatus;
//MPI_Request sendReq;

//std::vector<char> rbuf;
std::vector<char> sbuf;

int nActiveHolds;

LocalStoreType lstore;
MRU< std::pair<unsigned int, uint64_t> > mru;

/**
 * Queue of incoming transactions
 */
std::list<Transaction> incomingTransactions;

/**
 * Queue of outgoing transactions
 */
//std::list<Transaction> outgoingTransactions;


std::map<int,std::list<Transaction> > outgoingTransactions;


RecvChannel recvChannel;
std::map<int,SendChannel> sendChannels;

/**
 * Store pending transitions that are awaiting the arrival of a data item
 */
std::map< std::pair<unsigned int, Label>, std::list<Transaction> > pendingTransactions;

/**
 * PUT requests that were already forwarded to parent
 */
std::set< std::pair<unsigned int, Label> > forwardedPuts;

/**
 * GET requests from children we could not immediately fulfil
 */
std::map< std::pair<unsigned int, Label>, std::set<int> > pendingGetReplies;

/**
 * GET requests that we already issued to parent for which we didn't get an answer yet
 */
std::set< std::pair<unsigned int, Label> > pendingGets;

/**
 * PUT requests we already issued to parent
 */
std::set< std::pair<unsigned int, Label> > completedPutsToParent;

/**
 * Store holds on items. requestCounter[dbKey][key]={requestors}
 */
std::map< unsigned int, std::map<Label,std::multiset<int> > > requestCounter;

/**
 * Next available reservation ID
 */
uint64_t reservationID;

/**
 * Store open reservations. reservations[id]=(dbKey,key)
 */
std::map<uint64_t, std::pair<unsigned int, Label > >  reservations;

SpinLock transactionLock;
};


#endif /* MPICache_hpp */
