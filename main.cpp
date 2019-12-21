// #include <eigen3/Eigen/Core>

// using namespace Eigen;
// // autodiff include
// #include <autodiff/forward.hpp>
// #include <autodiff/forward/eigen.hpp>
// using namespace autodiff;


#include "mex.hpp"
#include "mexAdapter.hpp"
#include "stream.hpp"
#include "nlp.hpp"


using namespace matlab::data;
using matlab::mex::ArgumentList;

typedef matlab::data::Array mexArray;


class MexFunction : public matlab::mex::Function {
private:
    Buffer buffer;
    std::ostream cout;
    std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr;
    ArrayFactory factory;
    void displayError(std::string errorMessage) {
        matlabPtr->feval(
                        matlab::engine::convertUTF8StringToUTF16String("error"),
                        0, 
                        std::vector<mexArray>({
                            factory.createScalar(errorMessage) 
                            })
                        );
    }

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
        return (ArrayType::CHAR == x.getType());
    }

    inline bool isstringtype(const Array& x){
        return (ArrayType::MATLAB_STRING == x.getType());
    }

    inline bool isnumeric(const Array& x){
        bool retVal = false;
        for (int i = 3; i < 23; i++)
            retVal |= (ArrayType(i) == x.getType());
        return retVal;
    }

    void setSingleOption(SmartPtr<IpoptApplication> app, std::string name, const Array& option){
        if (isnumeric(option)){
            TypedArray<double> opt(std::move(option));
            if (isscalarinteger(opt)){
                app->Options()->SetIntegerValue(name, static_cast<int>(opt[0]));
            } else {
                app->Options()->SetNumericValue(name, opt[0]);
            }
        } else if (isstring(option)) {
            CharArray charopt(std::move(option));
            app->Options()->SetStringValue(name, charopt.toAscii());
        }

    }

    void setOptions(SmartPtr<IpoptApplication> app, StructArray& options) {
        if ( isfield(options, "ipopt") ) {
            StructArray opts = options[0]["ipopt"];
            auto fields = opts.getFieldNames();
            std::string fieldName("");
            Array curField;
            for (auto& field : fields){
                fieldName = field;
                curField = opts[0][field];
                setSingleOption(app, fieldName, curField);
            }
        }
    }

    

public:
    MexFunction() : 
        matlabPtr(getEngine()),
        cout(&buffer)
    {
        buffer.setMatlabPtr(matlabPtr);
    }
    ~MexFunction() {}
    void operator()(ArgumentList outputs, ArgumentList inputs) {
        
        if (inputs.size() != 3)
            displayError("x0, funcs and options are needed.");
        
        TypedArray<double> x0(inputs[0]);

        StructArray funcs(inputs[1]);
        StructArray opts(inputs[2]);

        SmartPtr<myNLP> mynlp = new myNLP();

        mynlp->setup(matlabPtr, x0, funcs, opts);

        SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
        app->RethrowNonIpoptException(true);

        // Set options that were passed
        setOptions(app, opts);

        if (isfield(opts, "ipopt") && (ArrayType::STRUCT == opts[0]["ipopt"].getType()) ) {
            StructArray ipoptopts = opts[0]["ipopt"];
            if (isfield(ipoptopts, "print_level") && isscalarinteger(ipoptopts[0]["print_level"]))
            {
                TypedArray<double> plTemp = ipoptopts[0]["print_level"];
                int pl(static_cast<int>(plTemp[0]));
                EJournalLevel printLevel(static_cast<EJournalLevel>(pl));
                SmartPtr<Journal> console = new MatlabJournal("MatlabJournal",printLevel, matlabPtr);
                app->Jnlst()->AddJournal(console);
            } 
        }


        ApplicationReturnStatus status;
        status = app->Initialize();
        if (status != Solve_Succeeded) {
            cout << std::endl << std::endl << "*** Error during initialization!" << std::endl;
        }

        // Ask Ipopt to solve the problem
        status = app->OptimizeTNLP(mynlp);

        if (status == Solve_Succeeded) {
            cout << std::endl << std::endl << "*** The problem solved!" << std::endl;
        }
        else {
            cout << std::endl << std::endl << "*** The problem FAILED!" << std::endl;
        }

        if (outputs.size() > 0)
            outputs[0] = mynlp->getX();

        if (outputs.size() > 1)
            outputs[1] = mynlp->getInfo();

    }
};