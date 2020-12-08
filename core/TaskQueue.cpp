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

#include "TaskQueue.hpp"


FIFOQueue::FIFOQueue(){
	nInstances=0;
};


void FIFOQueue::push(TaskDescriptor &t){
	tasks.push_front(t);
	nInstances+=t.nInstances;
};


bool FIFOQueue::pop(TaskDescriptor &t){

	if(not tasks.empty()) {

		// return first non-optional job.... order doesn't matter here
		for(auto tit=tasks.begin();tit!=tasks.end();++tit) {
			if(!tit->optional) {
				t = *tit;
				t.nInstances=1;
				if(--tit->nInstances==0) tasks.erase(tit);
				nInstances-=1;
				return true;
			}
		}

		t=tasks.back();
		t.nInstances=1;
		if(--tasks.back().nInstances==0) {
			tasks.pop_back();
		}
		nInstances-=1;
		return true;
	}
	return false;
};

void FIFOQueue::requiredData( std::unordered_set<DataPlaceholder> &required){
	for(auto it=tasks.begin(); it!=tasks.end(); it++) {
		it->requiredData(required);
	}
};

bool FIFOQueue::empty(){
	return tasks.empty();
};

void FIFOQueue::resize(int size){
	TaskDescriptor t;
	while(this->size()>=size) {
		this->pop(t);
	}
};

int FIFOQueue::size(){
	return nInstances;
};

int FIFOQueue::count(){
	return nInstances;
};

BackfillQueue::BackfillQueue(){
	maxBatchId=0;
	boost::random::random_device rd;
	rng.seed(rd());
};

void BackfillQueue::push(TaskDescriptor &t, uint64_t batchId){
	tasks.push_back(std::make_pair(batchId,t));
	weigths.push_back( t.nInstances);

	if(batchId>maxBatchId) {
		maxBatchId=batchId;
		purge();
	}
};

bool BackfillQueue::pop(TaskDescriptor &t){
	if(tasks.size()>0) {
		boost::random::discrete_distribution<> d(weigths.begin(),weigths.end());
		auto it=tasks.begin();
		std::advance(it, d(rng));
		t=it->second;
		return true;
	}
	return false;
};


void BackfillQueue::purge(){
	weigths.clear();
	auto it=tasks.begin();
	while(it!=tasks.end()) {
		if(it->first<maxBatchId) {
			it=tasks.erase(it);
		}
		else{
			weigths.push_back(it->second.nInstances);
			it++;
		}
	}
};


void BackfillQueue::requiredData( std::unordered_set<DataPlaceholder> &required){
	for(auto it=tasks.begin(); it!=tasks.end(); it++) {
		it->second.requiredData(required);
	}
};


bool BackfillQueue::empty(){
	return tasks.empty();
};
