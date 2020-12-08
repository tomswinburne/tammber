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



#ifndef TaskManager_hpp
#define TaskManager_hpp

#include <stdio.h>
#include <vector>

#include <mpi.h>
#include <unistd.h>
#include <iostream>
#include <chrono>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/info_parser.hpp>


#include "Task.hpp"
#include "Constants.hpp"
#include "Pack.hpp"

class AbstractDriverTaskManager {
virtual bool pullTask( GenericTask &t )=0;
virtual bool pushTask( GenericTask &t )=0;
};


class MPIDriverTaskManager : public AbstractDriverTaskManager {
public:
MPIDriverTaskManager(MPI_Comm parentComm_,MPI_Comm localComm_){
	parentComm=parentComm_;
	//MPI_Comm_get_parent(&parentComm);
	MPI_Comm_rank(localComm_,&rank);
};

bool pullTask( GenericTask &t){
	GenericTask tt;
	//std::cout<<"PULLING TASK"<<std::endl;
	//receive size
	int count;
	MPI_Bcast(&(count), 1, MPI_INT, 0, parentComm);
	buffer.resize(count);
	//std::cout<<"pull - size: "<<count<<std::endl;
	//receive data
	MPI_Bcast(&(buffer[0]), count, MPI_BYTE, 0, parentComm);

	//unpack
	unpack(buffer,tt,std::size_t(count));
	t=tt;

	return true;
};

bool pushTask( GenericTask &t){
	//std::cout<<"PUSHING TASK"<<std::endl;
	//Send the task back to parent
	if(rank==0) {
		std::vector<char> b;
		pack(b,t);
		int count=b.size();
		//send the size
		MPI_Send(&count, 1, MPI_INT, 0, TASK_MANAGER_SIZE_TAG, parentComm);
		//send the count
		MPI_Send(&(b[0]), b.size(), MPI_BYTE, 0, TASK_MANAGER_DATA_TAG, parentComm);
	}

	return true;
};


private:
MPI_Comm parentComm;
int rank;
std::vector<char> buffer;
};

/**
 * This needs a better constructor interface that is generic. Maybe pass a map of properties?
 */

class AbstractDriverHandle {
public:
virtual bool isIdle()=0;
virtual bool assign(GenericTask &t)=0;
virtual bool probe(GenericTask &t)=0;
virtual void die()=0;
virtual bool isDead()=0;
};




class MPIDriverHandle :   AbstractDriverHandle {
public:
MPIDriverHandle(MPI_Comm parentComm_){
	parentComm=parentComm_;
	idle=true;
	dead=false;
};

virtual void die(){
	dead=true;
};

virtual bool isDead(){
	return dead;
};

virtual bool isIdle(){
	return idle;
};


virtual bool assign(GenericTask &t){

	if(!idle or dead) {
		return false;
	}

	std::vector<char> b;
	pack(b,t);

	int count=int(b.size());
	//std::cout<<"ASSIGN "<<count<<std::endl;
	//send the size
	MPI_Bcast(&count, 1, MPI_INT, MPI_ROOT, parentComm);
	//send the data
	MPI_Bcast(&(b[0]), count, MPI_BYTE, MPI_ROOT, parentComm);

	idle=false;

	//post a receive for the size of the result
	MPI_Irecv(&(resultCount), 1, MPI_INT, MPI_ANY_SOURCE, TASK_MANAGER_SIZE_TAG, parentComm, &resultSizeRequest );

	//std::cout<<"DONE ASSIGN "<<std::endl;

	return true;
};


virtual bool probe(GenericTask &t){

	GenericTask tt;

	//std::cout<<"PROBE "<<std::endl;

	//there is no task pending on this driver
	if(idle or dead) {
		return false;
	}

	int flag;
	MPI_Status status;
	MPI_Test(&resultSizeRequest, &flag, &status);
	if(flag==0) {
		return false;
	}
	//if we make it here, the size of the result was received, so the task is complete
	//std::cout<<"probe - size: "<<resultCount<<std::endl;
	//receive the data
	buffer.resize(resultCount);
	MPI_Recv(&(buffer[0]), resultCount, MPI_BYTE, MPI_ANY_SOURCE, TASK_MANAGER_DATA_TAG, parentComm, &status);
	unpack(buffer,tt,std::size_t(resultCount));

	//we are ready to process another task
	idle=true;

	t=tt;

	//std::cout<<"PROBE SUCCESSFUL "<<resultCount<<std::endl;

	return true;
};


private:
std::string executable;
std::string args;
MPI_Comm parentComm;
MPI_Request resultSizeRequest;


bool idle;
bool dead;
int resultCount;
std::vector<char> buffer;

};





#endif /* TaskManager_hpp */
