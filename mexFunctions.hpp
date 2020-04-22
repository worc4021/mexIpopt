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

bool ispositive(const Array& x) {
    TypedArray<double> y(x);
    return 0 < y[0];
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

inline bool ishandle(const Array& x) {
    return ( ArrayType::HANDLE_OBJECT_REF == x.getType() );
}

inline bool issparse(const Array& x) {
    return ( ArrayType::SPARSE_DOUBLE == x.getType() );
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

void displayError(  std::shared_ptr<matlab::engine::MATLABEngine>& matlabPtr,
                    ArrayFactory& factory,
                    std::string& message )
                    {
        matlabPtr->feval(
                        matlab::engine::convertUTF8StringToUTF16String("error"),
                        0, 
                        std::vector<Array>({
                            factory.createScalar(message)
                            })
                        );
}