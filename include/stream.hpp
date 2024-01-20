#pragma once
#include "utilities.hpp"

/* Buffer class to send string stream into Matlab output */
class Buffer 
    : public std::stringbuf
{
private:
    matlab::data::ArrayFactory factory;
public:
    Buffer() = default;
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

    bool operator==(const Buffer& other) const
    {
        return false;
    }

    Buffer(const matlab::data::ArrayFactory& factory) = delete;
};