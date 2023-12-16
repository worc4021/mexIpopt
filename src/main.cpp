#include "mex.hpp"
#include "mexAdapter.hpp"
#include "utilities.hpp"
#include "stream.hpp"
#include "nlp.hpp"

struct evals {
    int objective;
    int constraints;
    int gradient;
    int jacobian;
    int hessian;
};

struct MatlabInfo {
    std::unique_ptr<Ipopt::Number[]> x;
    int status;
    std::unique_ptr<Ipopt::Number[]> zl;
    std::unique_ptr<Ipopt::Number[]> zu;
    std::unique_ptr<Ipopt::Number[]> lambda;
    int iter;
    Ipopt::Number cpu;
    Ipopt::Number objective;
    evals eval;
};

class MexFunction 
    : public matlab::mex::Function {
private:
    Buffer buffer;
    std::ostream cout;
    matlab::data::ArrayFactory factory;
    
    void setSingleOption(Ipopt::SmartPtr<Ipopt::IpoptApplication> app, std::string name, const matlab::data::Array& option){
        if (utilities::isnumeric(option)){
            matlab::data::TypedArray<double> opt(std::move(option));
            if (utilities::isscalarinteger(opt)){
                app->Options()->SetIntegerValue(name, static_cast<int>(opt[0]));
            } else {
                app->Options()->SetNumericValue(name, opt[0]);
            }
        } else if ( utilities::isstring(option) ) {
            app->Options()->SetStringValue(name, utilities::getstringvalue(option));
        }

    }

    void setOptions(Ipopt::SmartPtr<Ipopt::IpoptApplication> app, matlab::data::StructArray& options) {
        if ( utilities::isfield(options, "ipopt") ) {
            matlab::data::StructArray opts = options[0]["ipopt"];
            auto fields = opts.getFieldNames();
            std::string fieldName("");
            matlab::data::Array curField;
            for (auto& field : fields){
                fieldName = field;
                curField = opts[0][field];
                setSingleOption(app, fieldName, curField);
            }
        }
    }

    

public:
    MexFunction() : 
        cout(&buffer)
    {
        matlabPtr = getEngine();
    }
    ~MexFunction() {}
    void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) {
        
        if (inputs.size() != 3)
            utilities::error("x0, funcs and options are needed.");
        
        matlab::data::TypedArray<double> x0(inputs[0]);

        matlab::data::StructArray funcs(inputs[1]);
        matlab::data::StructArray opts(inputs[2]);

        Ipopt::SmartPtr<myNLP> mynlp = new myNLP();

        mynlp->setup(x0, funcs, opts);

        Ipopt::SmartPtr<Ipopt::IpoptApplication> app = IpoptApplicationFactory();
        app->RethrowNonIpoptException(true);

        // Set options that were passed
        setOptions(app, opts);

        if (utilities::isfield(opts, "ipopt") && utilities::isstruct(utilities::getfield(opts, "ipopt"))) {
            matlab::data::StructArray ipoptopts = opts[0]["ipopt"];
            if (    utilities::isfield(ipoptopts, "print_level") 
                &&  utilities::isscalarinteger(utilities::getfield(ipoptopts,"print_level"))
                &&  utilities::ispositive(utilities::getfield(ipoptopts,"print_level")))
            {
                matlab::data::TypedArray<double> plTemp = utilities::getfield(ipoptopts,"print_level");
                int pl(static_cast<int>(plTemp[0]));
                Ipopt::EJournalLevel printLevel(static_cast<Ipopt::EJournalLevel>(pl));
                Ipopt::SmartPtr<Ipopt::Journal> console = new MatlabJournal("MatlabJournal",printLevel, matlabPtr);
                app->Jnlst()->AddJournal(console);
            } 
        }

        MatlabInfo infoStruct;

        Ipopt::ApplicationReturnStatus status;
        status = app->Initialize();
        if (status != Ipopt::Solve_Succeeded) {
            cout << std::endl << std::endl << "*** Error during initialization!" << std::endl;
        }

        // Ask Ipopt to solve the problem
        status = app->OptimizeTNLP(mynlp);

        if (outputs.size() > 0)
            outputs[0] = mynlp->getX();

        if (outputs.size() > 1) {
            matlab::data::StructArray info = mynlp->getInfo();
            info[0]["status"] = factory.createScalar<int>(status);
            matlab::data::StructArrayRef evals = info[0]["eval"];

            if (Ipopt::IsValid(app->Statistics())) {
                matlab::data::TypedArrayRef<double> objective = info[0]["objective"];
                objective[0] = app->Statistics()->FinalObjective();
                matlab::data::TypedArrayRef<double> cpu = info[0]["cpu"];
                cpu[0] = app->Statistics()->TotalCpuTime();
                matlab::data::TypedArrayRef<int> iter = info[0]["iter"];
                iter[0] = app->Statistics()->IterationCount();
                matlab::data::TypedArrayRef<int> objectiveCalls = evals[0]["objective"];
                matlab::data::TypedArrayRef<int> constraintCalls = evals[0]["constraints"];
                matlab::data::TypedArrayRef<int> gradientCalls = evals[0]["gradient"];
                matlab::data::TypedArrayRef<int> jacobianCalls = evals[0]["jacobian"];
                matlab::data::TypedArrayRef<int> hessianCalls = evals[0]["hessian"];
                app->Statistics()->NumberOfEvaluations( objectiveCalls[0],
                                                        constraintCalls[0],
                                                        gradientCalls[0],
                                                        jacobianCalls[0],
                                                        hessianCalls[0]
                                                );
            }


            outputs[1] = std::move(info);
        }

    }
};