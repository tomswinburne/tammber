

#ifndef Data_hpp
#define Data_hpp

#include <vector>
#include <list>
#include <iostream>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/map.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/optional.hpp>

#include "Pack.hpp"

enum class NECESSITY : int {REQUIRED, OPTIONAL};

typedef uint64_t Label;

typedef char RawData;
typedef std::vector<RawData> RawDataVector;


/**
 * Abstract bit of data
 */
struct AbstractData {

	AbstractData(){
	};

	AbstractData(Label label_){

	};
	/*
	   template<class Archive> void serialize(Archive & ar, const unsigned int version){
	   };
	 */

};

/**
 * Abstract bit of data that lives in a data-store
 */
struct AbstractStoredData : public AbstractData {

	AbstractStoredData(){
	};

	AbstractStoredData(Label label_, unsigned int location_, bool shared_ ){
		label=label_;
		location=location_;
		shared=shared_;
		canonical_label=0;
	};

	// overload with canonical label
	AbstractStoredData(Label canonical_label_, Label label_, unsigned int location_, bool shared_ ){
		label=label_;
		canonical_label=canonical_label_;
		location=location_;
		shared=shared_;
	};


	template<class Archive> void serialize(Archive & ar, const unsigned int version){
		//ar & boost::serialization::base_object<AbstractData>(*this);
		ar & location;
		ar & shared;
		ar & label;
		ar & canonical_label;
	};

	Label label;
	Label canonical_label; //
	unsigned int location;
	bool shared;
};

/**
 * Holds a reference to a data item stored in a data-store
 */
struct DataPlaceholder : public AbstractStoredData {
	DataPlaceholder(){
	};

	DataPlaceholder(Label label_, unsigned int location_, bool shared_, NECESSITY necessity_ ){
		label=label_;
		canonical_label=0;
		location=location_;
		shared=shared_;
		necessity=necessity_;
	};

	// overload with canonical label
	DataPlaceholder(Label canonical_label_, Label label_, unsigned int location_, bool shared_, NECESSITY necessity_ ){
		label=label_;
		canonical_label=canonical_label_;
		location=location_;
		shared=shared_;
		necessity=necessity_;
	};

	bool operator ==(const DataPlaceholder &b) const {
		return (label==b.label and location==b.location and shared==b.shared and necessity==b.necessity);
	};

	bool operator !=(const DataPlaceholder &b) const {
		return !(*this==b);
	};

	template<class Archive> void serialize(Archive & ar, const unsigned int version){
		//ar & boost::serialization::base_object<AbstractStoredData>(*this);
		ar & location;
		ar & shared;
		ar & label;
		ar & canonical_label;
		ar & necessity;
	};

	NECESSITY necessity;
};

namespace std
{
template<>
struct hash<DataPlaceholder> {
	size_t operator()(const DataPlaceholder &s) const {
		return std::hash<uint64_t>()(s.label) ^ std::hash<int>()(static_cast<int>(s.necessity)) ^ std::hash<unsigned int>()(static_cast<unsigned int>(s.location));
	}
};
}


struct DataItem : public AbstractData {
	DataItem(){
	};

	template <class T> DataItem(T &t){
		pack(data,t);
	};

	template <class T> bool unpack(T &t){
		try{
			::unpack(data,t);
		}
		catch(...) {
			return false;
		}
		return true;
	};

	void clear(){
		data.clear();
	};

	template<class Archive> void serialize(Archive & ar, const unsigned int version){
		//ar & boost::serialization::base_object<AbstractData>(*this);
		ar & data;
	};

	RawDataVector data;
};



struct StoredDataItem : public AbstractStoredData {
	StoredDataItem(){
	};

	StoredDataItem(DataPlaceholder &d){
		label=d.label;
		canonical_label=d.canonical_label;
		shared=d.shared;
		location=d.location;
	};

	template <class T> StoredDataItem( Label label_, bool shared_, unsigned int location_, T &t){
		label=label_;
		canonical_label = 0;
		shared=shared_;
		location=location_;
		pack(data,t);
	};

	// overload with canonical label
	template <class T> StoredDataItem(Label canonical_label_, Label label_, bool shared_, unsigned int location_, T &t){
		label=label_;
		canonical_label=canonical_label_;
		shared=shared_;
		location=location_;
		pack(data,t);
	};


	template <class T> bool unpack(T &t){
		try{
			unpack(data,t);
		}
		catch(...) {
			return false;
		}
		return true;
	};

	void clear(){
		data.clear();
	};

	template<class Archive> void serialize(Archive & ar, const unsigned int version){
		//ar & boost::serialization::base_object<AbstractStoredData>(*this);
		ar & location;
		ar & shared;
		ar & label;
		ar & canonical_label;
		ar & data;
	};

	RawDataVector data;
};




/**
 * Helper functions to deal with Data maps
 */

/*
   template<class Data, class T>  bool extract(std::string name, std::map<std::string, Data> &map, T &out){
        if(map.count(name)) {
                try{
                        auto it=map.find(name);
                        unpack(it->second.data,out);
                }
                catch(...) {
                        std::cerr<<"EXTRACTION ERROR"<<std::endl;
                        return false;
                }
                return true;
        }
        return false;
   };
 */

template<class Data, class T>  bool extract(std::string name, std::multimap<std::string, Data> &map, T &out){
	if(map.count(name)) {
		try{
			auto it=map.find(name);
			unpack(it->second.data,out);
		}
		catch(...) {
			std::cerr<<"EXTRACTION ERROR"<<std::endl;
			return false;
		}
		return true;
	}
	return false;
};

template<class Data, class T>  void extract(std::string name, std::multimap<std::string, Data> &map, std::list<T> &out){
	out.clear();
	auto itr=map.equal_range(name);
	for(auto it=itr.first; it!=itr.second; it++) {
		try{
			T item;
			unpack(it->second.data,item);
			out.push_back(item);
		}
		catch(...) {
			std::cerr<<"LIST EXTRACTION ERROR"<<std::endl;
		}
	}
};



template<class T> void insert(std::string name, std::multimap<std::string, DataItem> &map, T &in ){
	map.insert(std::make_pair(name,DataItem(in)));
};

inline void insert(std::string name, std::multimap<std::string, StoredDataItem> &map, StoredDataItem &in ){
	map.insert(std::make_pair(name,in));
};

inline void insert(std::string name, Label label, unsigned int location, bool shared, NECESSITY necessity, std::multimap<std::string, DataPlaceholder> &map){
	map.insert(std::make_pair(name,DataPlaceholder(label, location, shared, necessity)));
};
// overload with canonical label
inline void insert(std::string name, Label canonical_label, Label label, unsigned int location, bool shared, NECESSITY necessity, std::multimap<std::string, DataPlaceholder> &map){
	map.insert(std::make_pair(name,DataPlaceholder(canonical_label, label, location, shared, necessity)));
};


template<class T> void insert(std::string name, Label label, unsigned int location, bool shared, std::multimap<std::string, StoredDataItem> &map, T &in){
	map.insert(std::make_pair(name,StoredDataItem(label, location, shared,in)));
};
// overload with canonical label
template<class T> void insert(std::string name, Label canonical_label, Label label, unsigned int location, bool shared, std::multimap<std::string, StoredDataItem> &map, T &in){
	map.insert(std::make_pair(name,StoredDataItem(canonical_label, label, location, shared,in)));
};


template<class T> void insert(std::string name, std::multimap<std::string, StoredDataItem> &map, T &in){
	map.insert(std::make_pair(name,StoredDataItem(0, 0, 0, false,in)));
};

// not used as far as I can see....
class CompressibleContainer {

public:

void compress(){
	compressed=true;
	rawInt.clear();
	compressedInt.clear();
};
void uncompress(Label index){
	treeIndex=index;
	compressed=false;
	compressedInt.clear();
	compressedDouble.clear();
};

template<class Archive>
void serialize(Archive & ar, const unsigned int version){
	ar & compressed;
	ar & treeIndex;
	ar & rawInt;
	ar & rawDouble;
	ar & compressedInt;
	ar & compressedDouble;
};

bool compressed;
Label treeIndex;
std::vector<int> rawInt;
std::vector<double> rawDouble;

std::vector<char> compressedInt;
std::vector<char> compressedDouble;
};




#endif
