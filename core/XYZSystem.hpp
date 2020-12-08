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



#ifndef XYZSystem_hpp
#define XYZSystem_hpp

#include <stdio.h>
#include "AbstractSystem.hpp"

#include "Constants.hpp"
#include <vector>
#include <array>
#include <map>

#include <boost/serialization/vector.hpp>
#include <boost/serialization/array.hpp>


/**
 * Test implementation of MinimalSystem. Contains only members required by the interface specification.
 */
class XYZSystem : public AbstractSystem {
public:
XYZSystem();
virtual ~XYZSystem(){
};
virtual int getNAtoms();
virtual void setNAtoms(int n);

virtual double getPosition(int iAtom, int iDim);
virtual void setPosition(int iAtom, int iDim, double p);

virtual int getSpecies(int iAtom);
virtual void setSpecies(int iAtom, int species);

virtual int getUniqueID(int iAtom);
virtual void setUniqueID(int iAtom, int i);

virtual void setBox(int iDim, int jDim, double box);
virtual double getBox(int iDim, int jDim);

virtual void setPeriodic(int iDim, bool periodic);
virtual int getPeriodic(int iDim);

virtual void pack(std::vector<char> &v);
virtual void unpack(std::vector<char> &v);

//static const bool placeholder=false;

template<class Archive>
void serialize(Archive & ar, const unsigned int version){
	ar & nAtoms;
	ar & b;
	ar & p;
	ar & spe;
	ar & pos;
	ar & id;
};

private:
void clearAtoms();
std::vector<int> spe;
std::vector<double> pos;
std::vector<int> id;


std::array< std::array<double, NDIM>, NDIM >  b;
std::array<int, NDIM> p;
int nAtoms;
};


#endif /* XYZSystem_hpp */
