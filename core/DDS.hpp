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



#ifndef MPICache_hpp
#define MPICache_hpp



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
#include "Data.hpp"


class AbstractDDS {
public:

virtual ~AbstractDDS(){
};

/**
 * Hint that you might need a given item in the future. This does not protect the item from
 * being purged before it is read.
 */
virtual void prefetch(int dbKey, Label key)=0;

/**
 * Get an item. This does not release any reservation.
 */
virtual bool get(int dbKey, Label key, RawDataVector &data)=0;

/**
 * Request an item and places a hold so it remains in store until it is released. Return the reservation ID.
 */
virtual uint64_t reserve(int dbKey, Label key)=0;

/**
 * Count the number of items with a given key currently in store.
 */
virtual int count(int dbKey, Label key)=0;

/**
 * Release a reservation.
 */
virtual bool release(uint64_t id)=0;

/**
 * Put an element in the store
 */
virtual void put(int dbKey, Label key, RawDataVector &data)=0;

/**
 * Return the keys of the items that are currently available.
 */
virtual std::set<Label> availableKeys(unsigned int dbKey)=0;

virtual void setMaximumBufferSize(uint64_t maxSize)=0;

virtual void serve()=0;

/**
 * Perform an iteration of the communication/management layer. This should be constantly called.
 */
virtual int singleServe()=0;


virtual void printStatus()=0;
virtual void sync()=0;
};



#define DDS_GET 1
#define DDS_PREFETCH 2
#define DDS_PUT 3




#endif /* MPICache_hpp */
