#pragma once
#include "mex.hpp"
#include "mexAdapter.hpp"

/* Buffer class to send string stream into Matlab output */
class Buffer : public std::stringbuf
{
private:
    std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr;
    matlab::data::ArrayFactory factory;
public:
    Buffer(void) : matlabPtr(nullptr) {};
    void setMatlabPtr(std::shared_ptr<matlab::engine::MATLABEngine>& mptr){
        matlabPtr = mptr;
    }
    int sync() {
        matlabPtr->feval(   matlab::engine::convertUTF8StringToUTF16String("fprintf"),
                            0,
                            std::vector<matlab::data::Array>({
                                    factory.createScalar(this->str())
                                
                            })
                        );
        str("");
        return 0;
    }
};

