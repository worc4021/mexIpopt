#include "mex.hpp"
#include "mexAdapter.hpp"
#include "utilities.hpp"
#include "stream.hpp"
#include "nlp.hpp"


std::string getApplicationStatus(Ipopt::ApplicationReturnStatus status) {
    using Ipopt::ApplicationReturnStatus;
    std::string retval{ "unknown" };
    switch (status) {
        case ApplicationReturnStatus::Solve_Succeeded:
            retval = "Solve_Succeeded";
            break;
        case ApplicationReturnStatus::Solved_To_Acceptable_Level:
            retval = "Solved_To_Acceptable_Level";
            break;
        case ApplicationReturnStatus::Infeasible_Problem_Detected:
            retval = "Infeasible_Problem_Detected";
            break;
        case ApplicationReturnStatus::Search_Direction_Becomes_Too_Small:
            retval = "Search_Direction_Becomes_Too_Small";
            break;
        case ApplicationReturnStatus::Diverging_Iterates:
            retval = "Diverging_Iterates";
            break;
        case ApplicationReturnStatus::User_Requested_Stop:
            retval = "User_Requested_Stop";
            break;
        case ApplicationReturnStatus::Feasible_Point_Found:
            retval = "Feasible_Point_Found";
            break;
        case ApplicationReturnStatus::Maximum_Iterations_Exceeded:
            retval = "Maximum_Iterations_Exceeded";
            break;
        case ApplicationReturnStatus::Restoration_Failed:
            retval = "Restoration_Failed";
            break;
        case ApplicationReturnStatus::Error_In_Step_Computation:
            retval = "Error_In_Step_Computation";
            break;
        case ApplicationReturnStatus::Maximum_CpuTime_Exceeded:
            retval = "Maximum_CpuTime_Exceeded";
            break;
        case ApplicationReturnStatus::Maximum_WallTime_Exceeded:
            retval = "Maximum_WallTime_Exceeded";
            break;
        case ApplicationReturnStatus::Not_Enough_Degrees_Of_Freedom:
            retval = "Not_Enough_Degrees_Of_Freedom";
            break;
        case ApplicationReturnStatus::Invalid_Problem_Definition:
            retval = "Invalid_Problem_Definition";
            break;
        case ApplicationReturnStatus::Invalid_Option:
            retval = "Invalid_Option";
            break;
        case ApplicationReturnStatus::Invalid_Number_Detected:
            retval = "Invalid_Number_Detected";
            break;
        case ApplicationReturnStatus::Unrecoverable_Exception:
            retval = "Unrecoverable_Exception";
            break;
        case ApplicationReturnStatus::NonIpopt_Exception_Thrown:
            retval = "NonIpopt_Exception_Thrown";
            break;
        case ApplicationReturnStatus::Insufficient_Memory:
            retval = "Insufficient_Memory";
            break;
        case ApplicationReturnStatus::Internal_Error:
            retval = "Internal_Error";
            break;
    }
    return retval;
        
}

class MexFunction 
    : public matlab::mex::Function {
private:
    Buffer buffer;
    std::ostream cout;
    matlab::data::ArrayFactory factory;
    
    void setSingleOption(Ipopt::SmartPtr<Ipopt::IpoptApplication> app, std::string name, matlab::data::Array& option){
        bool status{false};
        if (utilities::isnumeric(option)){
            matlab::data::TypedArray<double> opt(std::move(option));
            if (utilities::isscalarinteger(opt)){
                status = app->Options()->SetIntegerValue(name, static_cast<int>(opt[0]));
            } else {
                status = app->Options()->SetNumericValue(name, opt[0]);
            }
        } else if ( utilities::isstring(option) ) {
            status = app->Options()->SetStringValue(name, utilities::getstringvalue(option));
        }
        if (!status)
            utilities::printf("Failed to set field {}\n", name);
    }

    void setOptions(Ipopt::SmartPtr<Ipopt::IpoptApplication> app, const matlab::data::StructArray& opts) {
        auto fields = opts.getFieldNames();
        std::string fieldName("");
        matlab::data::Array curField;
        for (auto& field : fields){
            fieldName = field;
            curField = opts[0][field];
            setSingleOption(app, fieldName, curField);
        }

        if (utilities::isfield(opts, "print_level")
            && utilities::isscalarinteger(utilities::getfield(opts, "print_level"))
            && utilities::ispositive(utilities::getfield(opts, "print_level")))
        {
            matlab::data::TypedArray<double> plTemp = utilities::getfield(opts, "print_level");
            int pl(static_cast<int>(plTemp[0]));
            Ipopt::EJournalLevel printLevel(static_cast<Ipopt::EJournalLevel>(pl));
            Ipopt::SmartPtr<Ipopt::Journal> console = new MatlabJournal("MatlabJournal", printLevel);
            app->Jnlst()->AddJournal(console);
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
        
        if (inputs.size() != 1)
            utilities::error("Pass a struct with the fields 'variableInfo', 'funcs' and 'ipopt'");
        
        matlab::data::StructArray problem = std::move(inputs[0]);
        
        if (!utilities::isfield(problem, "variableInfo"))
            utilities::error("Field 'variableInfo' not supplied.");
        if (!utilities::isfield(problem, "funcs"))
            utilities::error("Field 'funcs' not supplied.");
        if (!utilities::isfield(problem, "ipopt"))
            utilities::warning("Field 'ipopt' not supplied.");

        Ipopt::SmartPtr<myNLP> mynlp = new myNLP();
        mynlp->setup(problem);

        
        Ipopt::SmartPtr<Ipopt::IpoptApplication> app = IpoptApplicationFactory();
        app->RethrowNonIpoptException(true);

        // Set options that were passed
        matlab::data::StructArray opts = utilities::getfield(problem, "ipopt");
        setOptions(app, opts);

        
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
            info[0]["appstatus"] = factory.createScalar(getApplicationStatus(status));
            matlab::data::StructArrayRef evals = info[0]["eval"];

            if (Ipopt::IsValid(app->Statistics())) {
                matlab::data::TypedArrayRef<double> objective = info[0]["objective"];
                matlab::data::TypedArrayRef<double> cpu = info[0]["cpu"];
                matlab::data::TypedArrayRef<int> iter = info[0]["iter"];
                matlab::data::TypedArrayRef<int> objectiveCalls = evals[0]["objective"];
                matlab::data::TypedArrayRef<int> constraintCalls = evals[0]["constraints"];
                matlab::data::TypedArrayRef<int> gradientCalls = evals[0]["gradient"];
                matlab::data::TypedArrayRef<int> jacobianCalls = evals[0]["jacobian"];
                matlab::data::TypedArrayRef<int> hessianCalls = evals[0]["hessian"];
                objective[0] = app->Statistics()->FinalObjective();
                cpu[0] = app->Statistics()->TotalCpuTime();
                iter[0] = app->Statistics()->IterationCount();
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
    MexFunction(const MexFunction&) = delete;
    MexFunction(MexFunction&&) = delete;
    MexFunction operator=(const MexFunction&) = delete;
};