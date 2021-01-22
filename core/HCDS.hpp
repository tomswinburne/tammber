//
//  Hcds.hpp

//
//  Created by Danny Perez on 1/9/17.
//  Copyright © 2017 dp. All rights reserved.
//

#ifndef Hcds_h
#define Hcds_h


#include <boost/multi_index_container.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/identity.hpp>
#include <boost/multi_index/sequenced_index.hpp>
#include <boost/multi_index/member.hpp>

#include <boost/serialization/vector.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/map.hpp>


#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/device/back_inserter.hpp>


#include "Types.hpp"
#include "Data.hpp"
#include "Log.hpp"
#include <iostream>


/**
 * Wraps a single database transaction.
 */
struct Transaction {
	Transaction(){
		pending=false;
	};
	Transaction(unsigned int type_, unsigned int dbKey_, Label &key_, RawDataVector &data_,int source_, int destination_){
		type=type_;
		dbKey=dbKey_;
		key=key_;
		data=data_;
		source=source_;
		destination=destination_;
		pending=false;
	};


	void print(){
		LOGGER(" type: "<<type<<" dbKey: "<<dbKey<<" key: "<<key<<" source: "<<source<<" destination: "<<destination<<" pending: "<<pending)
	};

	unsigned int dbKey;
	uint64_t key;
	RawDataVector data;
	int source;
	int destination;
	unsigned int type;
	bool pending;

	template<class Archive>
	void serialize(Archive & ar, const unsigned int version){
		ar & type;
		ar & dbKey;
		ar & key;
		ar & data;
		ar & source;
		ar & destination;
		ar & pending;
	}

};



/**
 * Maintains an ordered list of the most recently touched items. Oldest element comes first.
 */
template <class LabelType > class MRU
{
typedef boost::multi_index::multi_index_container<
		std::pair<unsigned int, Label>,
		boost::multi_index::indexed_by<
			boost::multi_index::sequenced<>,
			boost::multi_index::hashed_unique<boost::multi_index::identity<std::pair<unsigned int, Label> > >
			>
		> item_list;

public:
//typedef Item                         item_type;
typedef typename item_list::iterator iterator;
typedef typename item_list::reverse_iterator reverse_iterator;

void touch( LabelType item){
	std::pair<iterator,bool> p=il.push_back(item);

	if(!p.second) {                     /* duplicate item */
		il.relocate(il.end(),p.first); /* put in back */
	}
}

void erase(LabelType item){
	il.get<1>().erase(item);
}

iterator erase(iterator it){
	return il.erase(it);
}


std::pair<bool, LabelType > oldest(){
	LabelType b;
	bool ok=il.size()>0;
	if(ok) {
		b=il.front();
	}
	return std::make_pair(ok,b);
}

bool pop_oldest(){
	LabelType b;
	bool ok=il.size()>0;
	if(ok) {
		il.pop_front();
	}
	return ok;
}

iterator begin(){
	return il.begin();
};

iterator end(){
	return il.end();
};


private:
item_list il;

};






#endif /* Hcds_h */
