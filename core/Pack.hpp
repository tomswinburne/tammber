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



#ifndef Pack_hpp
#define Pack_hpp

#include <boost/serialization/vector.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/list.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/serialization/deque.hpp>


template <typename T1> void pack(std::vector<char> &v, T1 &t1){
	v.clear();
	//v=std::vector<char>();
	boost::iostreams::back_insert_device<std::vector<char> > sink {v};
	boost::iostreams::stream<boost::iostreams::back_insert_device<std::vector<char> > > os {sink};
	boost::archive::binary_oarchive oa(os, boost::archive::no_header | boost::archive::no_codecvt);
	oa << t1;
	os.flush();
	//os.close();
};

template <typename T1, typename T2> void pack(std::vector<char> &v, T1 &t1, T2 &t2){
	v.clear();
	//v=std::vector<char>();
	boost::iostreams::back_insert_device<std::vector<char> > sink {v};
	boost::iostreams::stream<boost::iostreams::back_insert_device<std::vector<char> > > os {sink};
	boost::archive::binary_oarchive oa(os, boost::archive::no_header | boost::archive::no_codecvt);
	oa << t1;
	oa << t2;
	os.flush();
	//os.close();
};

template <typename T1, typename T2, typename T3> void pack(std::vector<char> &v, T1 &t1, T2 &t2, T3 &t3){
	v.clear();
	//v=std::vector<char>();
	boost::iostreams::back_insert_device<std::vector<char> > sink {v};
	boost::iostreams::stream<boost::iostreams::back_insert_device<std::vector<char> > > os {sink};
	boost::archive::binary_oarchive oa(os, boost::archive::no_header | boost::archive::no_codecvt);
	oa << t1;
	oa << t2;
	oa << t3;
	os.flush();
	//os.close();
};

template <typename T1> void unpack(std::vector<char> &v, T1 &t1, std::size_t size=0){
	t1=T1();
	if(size<=0) {
		size=v.size();
	}
	boost::iostreams::array_source source {v.data(), size};
	boost::iostreams::stream<boost::iostreams::array_source> is {source};
	boost::archive::binary_iarchive ia(is, boost::archive::no_header | boost::archive::no_codecvt);
	ia >> t1;
};

template <typename T1, typename T2> void unpack(std::vector<char> &v, T1 &t1, T2 &t2, std::size_t size=0){
	t1=T1();
	t2=T2();
	if(size<=0) {
		size=v.size();
	}
	boost::iostreams::array_source source {v.data(), size};
	boost::iostreams::stream<boost::iostreams::array_source> is {source};
	boost::archive::binary_iarchive ia(is, boost::archive::no_header | boost::archive::no_codecvt);
	ia >> t1;
	ia >> t2;
};

template <typename T1, typename T2, typename T3> void unpack(std::vector<char> &v, T1 &t1, T2 &t2, T3 &t3,std::size_t size=0){
	t1=T1();
	t2=T2();
	t3=T3();
	if(size<=0) {
		size=v.size();
	}
	boost::iostreams::array_source source {v.data(), size};
	boost::iostreams::stream<boost::iostreams::array_source> is {source};
	boost::archive::binary_iarchive ia(is,  boost::archive::no_header | boost::archive::no_codecvt);
	ia >> t1;
	ia >> t2;
	ia >> t3;
};



 #endif
