#pragma once
#include "utilities.hpp"
#include "IpJournalist.hpp"

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

    Buffer(const Buffer&) = delete;
    Buffer(Buffer&&) = delete;
    Buffer operator=(const Buffer&) = delete;
};

class MatlabJournal : 
    public Ipopt::Journal
{
private:
    Buffer buffer;
    std::ostream cout;
public:

    MatlabJournal(const std::string& name,
                        Ipopt::EJournalLevel level) : 
        Ipopt::Journal(name, level),
        cout(&buffer) {}

    virtual ~MatlabJournal() {}

protected:

    virtual std::string Name() {
        return std::string("MatlabJournal");
    }

    // These functions override the functions in the Journal class.
    virtual
    void
    PrintImpl(  [[maybe_unused]]Ipopt::EJournalCategory category,
                [[maybe_unused]]Ipopt::EJournalLevel    level,
                char const *     str) {
        cout << std::string(str);
    }

    virtual
    void
    PrintfImpl( [[maybe_unused]]Ipopt::EJournalCategory category,
                [[maybe_unused]]Ipopt::EJournalLevel    level,
                char const *     pformat,
                va_list ap ) {
                    char _buffer[1000];
                    vsnprintf(_buffer, 1000, pformat, ap);
                    cout << std::string(_buffer);
                }

    virtual
    void FlushBufferImpl()
    {
        cout.flush();
    }

public:

    MatlabJournal() = delete;
    MatlabJournal(const MatlabJournal&) = delete;
    MatlabJournal(const Buffer& buffer, const std::ostream& cout) = delete;

    MatlabJournal operator=(const MatlabJournal& other) = delete;
};