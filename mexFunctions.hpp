#pragma once
#include "mex.hpp"
#include "mexAdapter.hpp"

using namespace matlab::data;
using matlab::mex::ArgumentList;


bool isfield(StructArray& s, std::string fieldname){
    bool retVal = true;
    try {
        s[0][fieldname].getType();
    }
    catch(const matlab::data::InvalidFieldNameException& e)
    {
        retVal = false;
    }
    return retVal;
}

inline bool isscalar(const Array& x){
    return (1 == x.getNumberOfElements());
}

inline bool iswholenumber(double x){
    return (0. == x-std::floor(x));
}

bool isinteger(const Array& x){
    TypedArray<double> y(x);
    bool retVal = true;
    for (auto& elem : y)
        retVal &= iswholenumber(elem);
    return retVal;
}

bool isscalarinteger(const Array& x){
    TypedArray<double> y(x);
    return (isinteger(y) && isscalar(y));
}

inline bool isstring(const Array& x){
    return ( (ArrayType::CHAR == x.getType()) || (ArrayType::MATLAB_STRING == x.getType()) );
}

inline bool isnumeric(const Array& x){
    bool retVal = false;
    for (int i = 3; i < 23; i++)
        retVal |= (ArrayType(i) == x.getType());
    return retVal;
}

std::string getStringValue(const Array& x) {
    if ( ArrayType::CHAR == x.getType() ) {
        CharArray retVal(std::move(x));
        return retVal.toAscii();
    } else if ( ArrayType::MATLAB_STRING == x.getType() ) {
        StringArray retVal(std::move(x));
        return retVal[0];
    }
    return std::string("");
}