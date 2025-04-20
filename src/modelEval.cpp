#include "mex.hpp"
#include "mexAdapter.hpp"

#include "utilities.hpp"
#include "sparse.hpp"
#include "stream.hpp"

#include "IpTNLP.hpp"
#include "IpIpoptApplication.hpp"
#include "IpIpoptCalculatedQuantities.hpp"
#include "IpIpoptData.hpp"
#include "IpTNLPAdapter.hpp"
#include "IpOrigIpoptNLP.hpp"
#include "IpJournalist.hpp"
#include "IpSolveStatistics.hpp"

#include <algorithm>
#include <memory>
#include <numbers>


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

std::string getStatus(Ipopt::SolverReturn status) {
    std::string retVal;
    switch (status)
    {
    case Ipopt::SUCCESS:
        retVal = "SUCCESS";
        break;
    case Ipopt::MAXITER_EXCEEDED:
        retVal = "MAXITER_EXCEEDED";
        break;
    case Ipopt::CPUTIME_EXCEEDED:
        retVal = "CPUTIME_EXCEEDED";
        break;
    case Ipopt::STOP_AT_TINY_STEP:
        retVal = "STOP_AT_TINY_STEP";
        break;
    case Ipopt::STOP_AT_ACCEPTABLE_POINT:
        retVal = "STOP_AT_ACCEPTABLE_POINT";
        break;
    case Ipopt::LOCAL_INFEASIBILITY:
        retVal = "LOCAL_INFEASIBILITY";
        break;
    case Ipopt::USER_REQUESTED_STOP:
        retVal = "USER_REQUESTED_STOP";
        break;
    case Ipopt::FEASIBLE_POINT_FOUND:
        retVal = "FEASIBLE_POINT_FOUND";
        break;
    case Ipopt::DIVERGING_ITERATES:
        retVal = "DIVERGING_ITERATES";
        break;
    case Ipopt::RESTORATION_FAILURE:
        retVal = "RESTORATION_FAILURE";
        break;
    case Ipopt::ERROR_IN_STEP_COMPUTATION:
        retVal = "ERROR_IN_STEP_COMPUTATION";
        break;
    case Ipopt::INVALID_NUMBER_DETECTED:
        retVal = "INVALID_NUMBER_DETECTED";
        break;
    case Ipopt::TOO_FEW_DEGREES_OF_FREEDOM:
        retVal = "TOO_FEW_DEGREES_OF_FREEDOM";
        break;
    case Ipopt::INVALID_OPTION:
        retVal = "INVALID_OPTION";
        break;
    case Ipopt::OUT_OF_MEMORY:
        retVal = "OUT_OF_MEMORY";
        break;
    case Ipopt::INTERNAL_ERROR:
        retVal = "INTERNAL_ERROR";
        break;
    case Ipopt::UNASSIGNED:
        retVal = "UNASSIGNED";
        break;
    case Ipopt::WALLTIME_EXCEEDED:
        retVal = "WALLTIME_EXCEEDED";
        break;
    default:
        retVal = "UNKOWN";
        break;
    }
    return retVal;
}

class modelCache
    : public Ipopt::TNLP
{
public:
    bool useIntermediateCallback{false};
    matlab::data::Array returnStruct;
private:
    
    utilities::Sparse<Ipopt::Number> jac;
    utilities::Sparse<Ipopt::Number> hess;

    std::vector<Ipopt::Number> _x;

    matlab::data::Array problemHandle;
    

    template <typename T>
    matlab::data::TypedArray<T> pullData(std::string propname) const
    {
        matlab::data::TypedArray<T> retval =
            matlabPtr->getProperty(
                problemHandle,
                matlab::engine::convertUTF8StringToUTF16String(propname));
        return retval;
    }

    bool need_modelcall(const Ipopt::Number* x) const
    {
        return !std::equal(_x.cbegin(), _x.cend(), x);
    }

    void update_model(const Ipopt::Number* x) 
    {

        _x.assign(x, x + _x.size());

        matlab::data::ArrayFactory factory;

        matlab::data::TypedArray<double> xAr = factory.createArray<double>({_x.size(), 1}, _x.data(), _x.data() + _x.size());

        matlabPtr->feval(
            matlab::engine::convertUTF8StringToUTF16String("update"),
            0,
            {
                problemHandle,
                xAr
            }
        );

        matlab::data::SparseArray<double> jacMat =
            matlabPtr->getProperty(
                problemHandle,
                matlab::engine::convertUTF8StringToUTF16String("cJacCur"));

        jac.updateValues(jacMat);
    }

    void update_hessian(Ipopt::Number sigma, const Ipopt::Number *lambda) {
        matlab::data::ArrayFactory factory;

        matlab::data::TypedArray<double> lambdaAr = factory.createArray<double>({jac.getNumberOfRows(), 1}, lambda, lambda + jac.getNumberOfRows());

        matlabPtr->feval(
            matlab::engine::convertUTF8StringToUTF16String("constructHessian"),
            0,
            {
                problemHandle,
                factory.createScalar(sigma),
                lambdaAr
            }
        );

        matlab::data::SparseArray<double> hessMat =
            matlabPtr->getProperty(
                problemHandle,
                matlab::engine::convertUTF8StringToUTF16String("HVal"));
        
        hess.updateValues(hessMat);
    
    }

public:
    modelCache(const matlab::data::Array &problem)
        : jac(), hess(), _x(), problemHandle(problem)
    {
        matlab::data::ArrayFactory factory;
        matlab::data::StructArray _returnStruct = factory.createStructArray({ 1,1 }, { "status" });
        _returnStruct[0]["status"] = factory.createScalar("Terminated Early");
        returnStruct = std::move(_returnStruct);

        auto isObjAr = matlabPtr->feval(
            matlab::engine::convertUTF8StringToUTF16String("isa"),
            1,
            {problemHandle,
             factory.createScalar("ipoptInterface")});

        matlab::data::TypedArray<bool> isObj = std::move(isObjAr[0]);

        if (isObj[0] == false)
            utilities::error("Input must be a ipoptInterface object.");

        matlab::data::SparseArray<double> jacMat =
            matlabPtr->getProperty(
                problemHandle,
                matlab::engine::convertUTF8StringToUTF16String("cJacCur"));

        jac.set(jacMat);

        matlab::data::SparseArray<double> hessMat =
            matlabPtr->getProperty(
                problemHandle,
                matlab::engine::convertUTF8StringToUTF16String("HVal"));
        hess.set(hessMat);

        // Check consistency.
        std::size_t nx = jac.getNumberOfColumns();
        std::size_t nc = jac.getNumberOfRows();

        if (nx != hess.getNumberOfColumns() || nx != hess.getNumberOfRows())
            utilities::error("Hessian is {} x {}, but expected to be {} x {}.", hess.getNumberOfRows(), hess.getNumberOfColumns(), nx, nx);

        matlab::data::TypedArray<double> xBnd = pullData<double>("xBnd");
        if (xBnd.getDimensions()[0] != nx)
            utilities::error("xBnd is {} x 2, but expected to be {} x 2.", xBnd.getDimensions()[0], nx);

        matlab::data::TypedArray<double> cBnd = pullData<double>("cBnd");
        if (cBnd.getDimensions()[0] != nc)
            utilities::error("cBnd is {} x 2, but expected to be {} x 2.", cBnd.getDimensions()[0], nc);

        matlab::data::TypedArray<double> x0 = pullData<double>("x0");
        if (x0.getDimensions()[0] != nx)
            utilities::error("xCur is {} x 1, but expected to be {} x 1.", x0.getDimensions()[0], nx);
        
        _x.resize(nx);
        std::fill(_x.begin(), _x.end(), 0.);
    }

    bool get_nlp_info(
        Ipopt::Index &n,
        Ipopt::Index &m,
        Ipopt::Index &nnz_jac_g,
        Ipopt::Index &nnz_h_lag,
        Ipopt::TNLP::IndexStyleEnum &index_style) override
    {
        index_style = Ipopt::TNLP::IndexStyleEnum::C_STYLE;

        n = jac.getNumberOfColumns();
        m = jac.getNumberOfRows();
        nnz_jac_g = jac.getNumberOfNonZeroElements();
        nnz_h_lag = hess.getNumberOfNonZeroElements();

        return true;
    }

    bool get_bounds_info(
        Ipopt::Index n,
        Ipopt::Number *x_l,
        Ipopt::Number *x_u,
        Ipopt::Index m,
        Ipopt::Number *g_l,
        Ipopt::Number *g_u) override
    {

        matlab::data::TypedArray<double> xBnd = pullData<double>("xBnd");
        matlab::data::TypedArray<double> cBnd = pullData<double>("cBnd");

        for (Ipopt::Index i = 0; i < n; ++i)
        {
            x_l[i] = xBnd[i][0];
            x_u[i] = xBnd[i][1];
        }

        for (Ipopt::Index i = 0; i < m; ++i)
        {
            g_l[i] = cBnd[i][0];
            g_u[i] = cBnd[i][1];
        }

        return true;
    }

    bool get_starting_point(
        Ipopt::Index n,
        bool init_x,
        Ipopt::Number *x,
        bool init_z,
        Ipopt::Number *z_L,
        Ipopt::Number *z_U,
        [[maybe_unused]]Ipopt::Index m,
        bool init_lambda,
        Ipopt::Number *lambda) override
    {
        if (init_x)
        {
            matlab::data::TypedArray<double> x0 = pullData<double>("x0");
            std::copy(x0.cbegin(), x0.cend(), x);
        }
        
        if (init_z)
        {
            matlab::data::TypedArray<double> z0 = pullData<double>("z0");
            for (Ipopt::Index i = 0; i < n; ++i)
            {
                z_L[i] = z0[i][0];
                z_U[i] = z0[i][1];
            }
        }

        if (init_lambda)
        {
            matlab::data::TypedArray<double> lambda0 = pullData<double>("lambda0");
            std::copy(lambda0.cbegin(), lambda0.cend(), lambda);
        }

        return true;
    }

    bool eval_f(
        [[maybe_unused]]Ipopt::Index n,
        const Ipopt::Number *x,
        bool new_x,
        Ipopt::Number &obj_value) override
    {
        if (new_x && need_modelcall(x))
            update_model(x);

        auto f = pullData<Ipopt::Number>("fVal");
        obj_value = f[0];

        return true;
    }

    bool eval_grad_f(
        [[maybe_unused]]Ipopt::Index n,
        const Ipopt::Number *x,
        bool new_x,
        Ipopt::Number *grad_f) override
    {
        if (new_x && need_modelcall(x))
            update_model(x);

        auto grad = pullData<Ipopt::Number>("gCur");

        std::copy(grad.cbegin(), grad.cend(), grad_f);

        return true;
    }

    bool eval_g(
        [[maybe_unused]]Ipopt::Index n,
        const Ipopt::Number *x,
        bool new_x,
        [[maybe_unused]]Ipopt::Index m,
        Ipopt::Number *g) override
    {
        if (new_x && need_modelcall(x))
            update_model(x);

        auto cons = pullData<Ipopt::Number>("cVal");

        std::copy(cons.cbegin(), cons.cend(), g);

        return true;
    }

    bool eval_jac_g(
        [[maybe_unused]]Ipopt::Index n,
        const Ipopt::Number *x,
        bool new_x,
        [[maybe_unused]]Ipopt::Index m,
        [[maybe_unused]]Ipopt::Index nele_jac,
        Ipopt::Index *iRow,
        Ipopt::Index *jCol,
        Ipopt::Number *values) override
    {
        if (nullptr == values) {
            jac.iRow(iRow);
            jac.jCol(jCol);
        } else {

            if (new_x && need_modelcall(x))
                update_model(x);
            
            jac.val(values);

        }

        return true;
    }

    bool eval_h(
        [[maybe_unused]]Ipopt::Index n,
        [[maybe_unused]]const Ipopt::Number *x,
        bool new_x,
        Ipopt::Number obj_factor,
        [[maybe_unused]]Ipopt::Index m,
        const Ipopt::Number *lambda,
        bool new_lambda,
        [[maybe_unused]]Ipopt::Index nele_hess,
        Ipopt::Index *iRow,
        Ipopt::Index *jCol,
        Ipopt::Number *values) override
    {
        if (nullptr == values) {
            hess.iRow(iRow);
            hess.jCol(jCol);
        } else {

            if (new_x && need_modelcall(x))
                update_model(x);

            if (new_lambda)
                update_hessian(obj_factor, lambda);

            hess.val(values);

        }
        return true;
    }

    void finalize_solution(
        Ipopt::SolverReturn status,
        [[maybe_unused]]Ipopt::Index n,
        const Ipopt::Number *x,
        const Ipopt::Number *z_L,
        const Ipopt::Number *z_U,
        [[maybe_unused]]Ipopt::Index m,
        [[maybe_unused]]const Ipopt::Number *g,
        const Ipopt::Number *lambda,
        Ipopt::Number obj_value,
        [[maybe_unused]]const Ipopt::IpoptData *ip_data,
        [[maybe_unused]]Ipopt::IpoptCalculatedQuantities *ip_cq) override
    {
        matlab::data::ArrayFactory factory;
        matlab::data::StructArray _retStruct = factory.createStructArray({ 1, 1 }, 
            {   "x",
                "fVal",
                "z_L",
                "z_U",
                "lambda",
                "status"
                });
        
        std::size_t _n = jac.getNumberOfColumns();
        std::size_t _m = jac.getNumberOfRows();

        _retStruct[0]["x"] = factory.createArray<double>({ _n,1 }, x, x + _n);
        _retStruct[0]["fVal"] = factory.createScalar(obj_value);
        _retStruct[0]["z_L"] = factory.createArray<double>({ _n,1 }, z_L, z_L + _n);
        _retStruct[0]["z_U"] = factory.createArray<double>({ _n,1 }, z_U, z_U + _n);
        _retStruct[0]["lambda"] = factory.createArray<double>({ _m,1 }, lambda, lambda + _m);
        _retStruct[0]["status"] = factory.createScalar(getStatus(status));

        returnStruct = std::move(_retStruct);
    }

    bool intermediate_callback(
        [[maybe_unused]]Ipopt::AlgorithmMode mode,
        Ipopt::Index iter,
        Ipopt::Number obj_value,
        Ipopt::Number inf_pr,
        Ipopt::Number inf_du,
        Ipopt::Number mu,
        Ipopt::Number d_norm,
        Ipopt::Number regularization_size,
        Ipopt::Number alpha_du,
        Ipopt::Number alpha_pr,
        Ipopt::Index  ls_trials,
        const Ipopt::IpoptData *ip_data,
        Ipopt::IpoptCalculatedQuantities *ip_cq) override
    {
        if (!useIntermediateCallback)
            return true;

        matlab::data::ArrayFactory factory;
        Ipopt::Index n = jac.getNumberOfColumns();
        Ipopt::Index m = jac.getNumberOfRows();
        
        matlab::data::buffer_ptr_t<Ipopt::Number> x                         = factory.createBuffer<Ipopt::Number>(n);
        matlab::data::buffer_ptr_t<Ipopt::Number> z_L                       = factory.createBuffer<Ipopt::Number>(n);
        matlab::data::buffer_ptr_t<Ipopt::Number> z_U                       = factory.createBuffer<Ipopt::Number>(n);
        matlab::data::buffer_ptr_t<Ipopt::Number> g                         = factory.createBuffer<Ipopt::Number>(m);
        matlab::data::buffer_ptr_t<Ipopt::Number> lambda                    = factory.createBuffer<Ipopt::Number>(m);

        matlab::data::buffer_ptr_t<Ipopt::Number> x_L_violation             = factory.createBuffer<Ipopt::Number>(n);
        matlab::data::buffer_ptr_t<Ipopt::Number> x_U_violation             = factory.createBuffer<Ipopt::Number>(n);
        matlab::data::buffer_ptr_t<Ipopt::Number> compl_x_L                 = factory.createBuffer<Ipopt::Number>(n);
        matlab::data::buffer_ptr_t<Ipopt::Number> compl_x_U                 = factory.createBuffer<Ipopt::Number>(n);
        matlab::data::buffer_ptr_t<Ipopt::Number> grad_lag_x                = factory.createBuffer<Ipopt::Number>(n);
        matlab::data::buffer_ptr_t<Ipopt::Number> nlp_constraint_violation  = factory.createBuffer<Ipopt::Number>(m);
        matlab::data::buffer_ptr_t<Ipopt::Number> compl_g                   = factory.createBuffer<Ipopt::Number>(m);

        bool getScaledValues = true;
        bool status = get_curr_iterate( ip_data, 
                                        ip_cq, 
                                        getScaledValues, 
                                        n, 
                                        x.get(), 
                                        z_L.get(), 
                                        z_U.get(), 
                                        m, 
                                        g.get(), 
                                        lambda.get());

        if (!status)
            utilities::warning("Intermediate callback failed to obtain internal values for current iterate.");

        status = get_curr_violations(   ip_data, 
                                        ip_cq, 
                                        getScaledValues, 
                                        n, 
                                        x_L_violation.get(), 
                                        x_U_violation.get(), 
                                        compl_x_L.get(), 
                                        compl_x_U.get(), 
                                        grad_lag_x.get(), 
                                        m, 
                                        nlp_constraint_violation.get(), 
                                        compl_g.get());

        if (!status)
            utilities::warning("Intermediate callback failed to obtain internal values for current violations.");

        matlab::data::StructArray passVal = factory.createStructArray(
                            {1, 1},
                            {   "x", 
                                "z_L", 
                                "z_U",
                                "g",
                                "lambda",
                                "x_L_violation",
                                "x_U_violation",
                                "compl_x_L",
                                "compl_x_U",
                                "grad_lag_x",
                                "nlp_constraint_violation",
                                "compl_g",
                                "iter",
                                "obj_value",
                                "inf_pr",
                                "inf_du",
                                "mu",
                                "d_norm",
                                "regularization_size",
                                "alpha_du",
                                "alpha_pr",
                                "ls_trials"
                            }
                            );

        std::size_t _n = static_cast<std::size_t>(n);
        std::size_t _m = static_cast<std::size_t>(m);

        passVal[0]["x"]                         = factory.createArrayFromBuffer({_n, 1}, std::move(x));
        passVal[0]["z_L"]                       = factory.createArrayFromBuffer({_n, 1}, std::move(z_L));
        passVal[0]["z_U"]                       = factory.createArrayFromBuffer({_n, 1}, std::move(z_U));
        passVal[0]["g"]                         = factory.createArrayFromBuffer({_m, 1}, std::move(g));
        passVal[0]["lambda"]                    = factory.createArrayFromBuffer({_m, 1}, std::move(lambda));
        passVal[0]["x_L_violation"]             = factory.createArrayFromBuffer({_n, 1}, std::move(x_L_violation));
        passVal[0]["x_U_violation"]             = factory.createArrayFromBuffer({_n, 1}, std::move(x_U_violation));
        passVal[0]["compl_x_L"]                 = factory.createArrayFromBuffer({_n, 1}, std::move(compl_x_L));
        passVal[0]["compl_x_U"]                 = factory.createArrayFromBuffer({_n, 1}, std::move(compl_x_U));
        passVal[0]["grad_lag_x"]                = factory.createArrayFromBuffer({_n, 1}, std::move(grad_lag_x));
        passVal[0]["nlp_constraint_violation"]  = factory.createArrayFromBuffer({_m, 1}, std::move(nlp_constraint_violation));
        passVal[0]["compl_g"]                   = factory.createArrayFromBuffer({_m, 1}, std::move(compl_g));

        passVal[0]["iter"]                      = factory.createScalar(iter);
        passVal[0]["obj_value"]                 = factory.createScalar(obj_value);
        passVal[0]["inf_pr"]                    = factory.createScalar(inf_pr);
        passVal[0]["inf_du"]                    = factory.createScalar(inf_du);
        passVal[0]["mu"]                        = factory.createScalar(mu);
        passVal[0]["d_norm"]                    = factory.createScalar(d_norm);
        passVal[0]["regularization_size"]       = factory.createScalar(regularization_size);
        passVal[0]["alpha_du"]                  = factory.createScalar(alpha_du);
        passVal[0]["alpha_pr"]                  = factory.createScalar(alpha_pr);
        passVal[0]["ls_trials"]                 = factory.createScalar(ls_trials);

        auto retVec = matlabPtr->feval(
            matlab::engine::convertUTF8StringToUTF16String("intermediateCallback"),
            1,
            {
                problemHandle,
                passVal
            }
        );

        matlab::data::TypedArray<bool> ret = std::move(retVec[0]);
        return ret[0];
    }
};


class MexFunction 
    : public matlab::mex::Function {
private:
    Buffer buffer;
    std::ostream cout;
    matlab::data::ArrayFactory factory;
    
    void setSingleOption(Ipopt::SmartPtr<Ipopt::IpoptApplication> app, std::string name, matlab::data::Array& option){
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
        if (inputs.size() < 1)
            utilities::error("Model wrapper of class 'BaseWrapper' must be passed as first arguments.");
        
        matlab::data::ArrayFactory factory;

        matlab::data::StructArray opts = factory.createStructArray({0, 1},{});

        if (inputs.size() > 1)
            opts = std::move(inputs[1]);

        Ipopt::SmartPtr<modelCache> nlp = new modelCache(inputs[0]);
        
        if (utilities::isfield(opts, "useIntermediateCallback"))
        {
            matlab::data::TypedArray<bool> useIntermediateCallback = utilities::getfield(opts, "useIntermediateCallback");
            nlp->useIntermediateCallback = useIntermediateCallback[0];
        }
        
        Ipopt::SmartPtr<Ipopt::IpoptApplication> app = IpoptApplicationFactory();
        app->RethrowNonIpoptException(true);

        if (utilities::isfield(opts, "ipopt"))
        {
            matlab::data::StructArray ipoptOpts = utilities::getfield(opts, "ipopt");
            setOptions(app, ipoptOpts);
        }
        
        Ipopt::ApplicationReturnStatus status;
        status = app->Initialize();
        if (status != Ipopt::Solve_Succeeded) {
            utilities::error("Error during initialization: {}", getApplicationStatus(status));
        }

        // Ask Ipopt to solve the problem
        status = app->OptimizeTNLP(nlp);

        if (outputs.size() > 0)
            outputs[0] = std::move(nlp->returnStruct);

        if (outputs.size() > 1) {
            matlab::data::StructArray statistics = factory.createStructArray({ 1,1 }, {
                "appstatus",
                "num_obj_evals",
                "num_constr_evals",
                "num_obj_grad_evals",
                "num_constr_jac_evals",
                "num_hess_evals",
                "dual_inf",        ///< dual infeasibility (Gradient of Lagrangian not zero)
                "constr_viol",     ///< violation of constraints
                "varbounds_viol",  ///< violation of variable bounds @since 3.14.0
                "complementarity", ///< violation of complementarity
                "kkt_error",
                "totalCpuTime",
                "totalWallclockTime",
                "totalSysTime",
                "iterationCount" });
            
            Ipopt::Index num_obj_evals{ -1 };
            Ipopt::Index num_constr_evals{ -1 };
            Ipopt::Index num_obj_grad_evals{ -1 };
            Ipopt::Index num_constr_jac_evals{ -1 };
            Ipopt::Index num_hess_evals{ -1 };
            Ipopt::Number dual_inf{std::numeric_limits<Ipopt::Number>::quiet_NaN()};
            Ipopt::Number constr_viol{ std::numeric_limits<Ipopt::Number>::quiet_NaN() };
            Ipopt::Number varbounds_viol{ std::numeric_limits<Ipopt::Number>::quiet_NaN() };
            Ipopt::Number complementarity{ std::numeric_limits<Ipopt::Number>::quiet_NaN() };
            Ipopt::Number kkt_error{ std::numeric_limits<Ipopt::Number>::quiet_NaN() };
            Ipopt::Number totalCpuTime{ std::numeric_limits<Ipopt::Number>::quiet_NaN() };
            Ipopt::Number totalWallclockTime{ std::numeric_limits<Ipopt::Number>::quiet_NaN() };
            Ipopt::Number totalSysTime{ std::numeric_limits<Ipopt::Number>::quiet_NaN() };
            Ipopt::Index iterationCount{ -1 };


            if (Ipopt::IsValid(app->Statistics())) {
                app->Statistics()->NumberOfEvaluations(num_obj_evals,num_constr_evals,num_obj_grad_evals,num_constr_jac_evals,num_hess_evals);
                app->Statistics()->Infeasibilities(dual_inf, constr_viol, varbounds_viol, complementarity, kkt_error);
                totalCpuTime = app->Statistics()->TotalCpuTime();
                totalWallclockTime = app->Statistics()->TotalWallclockTime();
                totalSysTime = app->Statistics()->TotalSysTime();
                iterationCount = app->Statistics()->IterationCount();
            }

            statistics[0]["appstatus"] = factory.createScalar(getApplicationStatus(status));
            statistics[0]["num_obj_evals"] = factory.createScalar(num_obj_evals);
            statistics[0]["num_constr_evals"] = factory.createScalar(num_constr_evals);
            statistics[0]["num_obj_grad_evals"] = factory.createScalar(num_obj_grad_evals);
            statistics[0]["num_constr_jac_evals"] = factory.createScalar(num_constr_jac_evals);
            statistics[0]["num_hess_evals"] = factory.createScalar(num_hess_evals);
            statistics[0]["dual_inf"] = factory.createScalar(dual_inf);
            statistics[0]["constr_viol"] = factory.createScalar(constr_viol);
            statistics[0]["varbounds_viol"] = factory.createScalar(varbounds_viol);
            statistics[0]["complementarity"] = factory.createScalar(complementarity);
            statistics[0]["kkt_error"] = factory.createScalar(kkt_error);
            statistics[0]["totalCpuTime"] = factory.createScalar(totalCpuTime);
            statistics[0]["totalWallclockTime"] = factory.createScalar(totalWallclockTime);
            statistics[0]["totalSysTime"] = factory.createScalar(totalSysTime);
            statistics[0]["iterationCount"] = factory.createScalar(iterationCount);

            outputs[1] = std::move(statistics);
        }
        

    }
    MexFunction(const MexFunction&) = delete;
    MexFunction(MexFunction&&) = delete;
    MexFunction operator=(const MexFunction&) = delete;
};