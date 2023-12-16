#pragma once
#include "utilities.hpp"
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

class MatlabJournal : 
    public Ipopt::Journal
{
private:
    Buffer buffer;
    std::ostream cout;
public:

    MatlabJournal(const std::string& name,
                        Ipopt::EJournalLevel level,
                        std::shared_ptr<matlab::engine::MATLABEngine> mptr) : 
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
    PrintImpl( Ipopt::EJournalCategory category,
            Ipopt::EJournalLevel    level,
            char const *     str) {
        cout << std::string(str);
    }

    virtual
    void
    PrintfImpl( Ipopt::EJournalCategory category,
                Ipopt::EJournalLevel    level,
                char const *     pformat,
                va_list ap ) {
                    char buffer[1000];
                    vsnprintf(buffer, 1000, pformat, ap);
                    cout << std::string(buffer);
                }

    virtual
    void FlushBufferImpl()
    {
        cout.flush();
    }

};


class myNLP 
    : public Ipopt::TNLP
{
private:
    
    matlab::data::ArrayFactory factory;
    matlab::data::StructArray funcs;
    matlab::data::StructArray options;
    matlab::data::StructArray retStr;
    std::size_t _n;
    std::size_t _m;
    std::vector<matlab::data::Array> args;
    bool isinitialised;
    bool returnHessian;
    bool intermediateCallback;
    bool diagnosticPrintout;
    Buffer buffer;
    std::ostream stream;

    std::vector<matlab::data::Array> feval(
            const matlab::data::Array &handle, 
            const std::size_t numReturned, 
            const std::vector<matlab::data::Array> &arguments
            ) {
        std::vector<matlab::data::Array> theArgs({handle});
        theArgs.insert(theArgs.end(), arguments.begin(), arguments.end());
        return matlabPtr->feval(
            matlab::engine::convertUTF8StringToUTF16String("feval"),
            numReturned,
            theArgs
        );
    }

    void updateX(const Ipopt::Number *x){
        for (auto i = 0 ; i < args[0].getNumberOfElements(); i++)
            args[0][i] = x[i];
    }

    std::string getStatus(Ipopt::SolverReturn status){
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
        default:
            break;
        }
        return retVal;
    }

public:
    /** default constructor */
    myNLP() : isinitialised(false), returnHessian(true), intermediateCallback(false),
        funcs(factory.createStructArray({0,0},{})),
        options(factory.createStructArray({0,0},{})),
        retStr(factory.createStructArray({1,1}, {"z_L", "z_U", "lambda", "status", "iter","cpu","objective","eval"})), // Make sure the structure is always available
        stream(&buffer), 
        args({factory.createEmptyArray()}),
        diagnosticPrintout(false) 
    {
        retStr[0]["eval"] = factory.createStructArray({1,1},{"objective","constraints","gradient","jacobian","hessian"});
        matlab::data::StructArrayRef evals = retStr[0]["eval"];
        evals[0]["objective"] = factory.createScalar<int>(0);
        evals[0]["constraints"] = factory.createScalar<int>(0);
        evals[0]["gradient"] = factory.createScalar<int>(0);
        evals[0]["jacobian"] = factory.createScalar<int>(0);
        evals[0]["hessian"] = factory.createScalar<int>(0);
        retStr[0]["iter"] = factory.createScalar<int>(0);
        retStr[0]["cpu"] = factory.createScalar<double>(0.);
        retStr[0]["objective"] = factory.createScalar<double>(NAN);
    };

    matlab::data::TypedArray<double> getX(void) {
        return args[0];
    };

    matlab::data::StructArray getInfo(void) {
        return retStr;
    }

    void setup( 
            matlab::data::TypedArray<double>& _x0, 
            matlab::data::StructArray& _funcs, 
            matlab::data::StructArray& _options
        ) {

        args[0] = _x0;

        funcs = _funcs;
        options = _options;

        if (utilities::isfield(options, "ipopt")){
            matlab::data::StructArray ipoptStr = options[0]["ipopt"];
            if (utilities::isfield(ipoptStr, "hessian_approximation")){
                if (utilities::isstring(ipoptStr[0]["hessian_approximation"])  ) {
                    returnHessian = (
                        0 != utilities::getstringvalue(ipoptStr[0]["hessian_approximation"])
                                .compare("limited-memory")
                                );
                } else {
                    utilities::error("limited-memory must be a string type");
                }
            }
        }

        intermediateCallback = utilities::isfield(funcs, "intermediate") && utilities::ishandle(funcs[0]["intermediate"]);
        diagnosticPrintout = utilities::isfield(options, "debug");
    }

    matlab::data::Array fevalWithX(const matlab::data::Array& handle) {
        if (diagnosticPrintout)
            stream << "Enter fevalWithX" << std::endl;

        std::vector<matlab::data::Array> retVal = feval(handle, 1, args);
        if (diagnosticPrintout)
            stream << "Exit fevalWithX" << std::endl;
        return retVal[0];
    }


  /** default destructor */
  ~myNLP() {};

  /**@name Overloaded from TNLP */
  //@{
  /** Method to return some info about the nlp */
  bool get_nlp_info(Ipopt::Index& n, Ipopt::Index& m, Ipopt::Index& nnz_jac_g,
                    Ipopt::Index& nnz_h_lag, IndexStyleEnum& index_style) 
    {
        if (diagnosticPrintout)
            stream << "Enter get_nlp_info" << std::endl;
        index_style = C_STYLE;
        n = args[0].getNumberOfElements();
        _n = n;

        matlab::data::TypedArray<double> cu = options[0]["cu"];
        m = cu.getNumberOfElements();
        _m = m;
        
        // Initialise return arrays in case things fail.
        retStr[0]["z_L"] = factory.createArray<double>({static_cast<size_t>(n),1});
        retStr[0]["z_U"] = factory.createArray<double>({static_cast<size_t>(n),1});
        retStr[0]["lambda"] = factory.createArray<double>({static_cast<size_t>(m),1});
        

        std::vector<matlab::data::Array> nullArgs = {};

        std::vector<matlab::data::Array> jacStrOut = feval(funcs[0]["jacobianstructure"], 1, nullArgs);
        
        if ( !utilities::issparse(jacStrOut[0]) ) {
            utilities::error("Jacobian has to be sparse");
        }

        matlab::data::SparseArray<double> Jac = std::move(jacStrOut[0]);

        nnz_jac_g = Jac.getNumberOfNonZeroElements();

        if (returnHessian){
            std::vector<matlab::data::Array> hesStrOut = feval(funcs[0]["hessianstructure"],1,nullArgs);
            if ( !utilities::issparse(hesStrOut[0]) ) {
                utilities::error("Hessian has to be sparse");
            }
            matlab::data::SparseArray<double> Hes = std::move(hesStrOut[0]);
            nnz_h_lag = Hes.getNumberOfNonZeroElements();
        }
        if (diagnosticPrintout)
            stream << "Exit get_nlp_info" << std::endl;

        return true;
    };

  /** Method to return the bounds for my problem */
  bool get_bounds_info(Ipopt::Index n, Ipopt::Number* x_l, Ipopt::Number* x_u,
                Ipopt::Index m, Ipopt::Number* g_l, Ipopt::Number* g_u) 
    {
        if (diagnosticPrintout)
            stream << "Enter get_bounds_info" << std::endl;
        matlab::data::TypedArray<double> xl = std::move(options[0]["lb"]);
        matlab::data::TypedArray<double> xu = std::move(options[0]["ub"]);
        matlab::data::TypedArray<double> gl = std::move(options[0]["cl"]);
        matlab::data::TypedArray<double> gu = std::move(options[0]["cu"]);

        for (auto i=0; i<n; i++)
            x_l[i] = xl[i];

        for (auto i=0; i<n; i++)
            x_u[i] = xu[i];

        for (auto i=0; i<m; i++)
            g_l[i] = gl[i];

        for (auto i=0; i<m; i++)
            g_u[i] = gu[i];
            
        if (diagnosticPrintout)
            stream << "Exit get_bounds_info" << std::endl;    
        return true;
    };

  /** Method to return the starting point for the algorithm */
  bool get_starting_point(  Ipopt::Index n, bool init_x, Ipopt::Number* x,
                            bool init_z, Ipopt::Number* z_L, Ipopt::Number* z_U,
                            Ipopt::Index m, bool init_lambda,
                            Ipopt::Number* lambda) 
{   
    if (diagnosticPrintout)
        stream << "Enter get_starting_point" << std::endl;
    if (init_x){
        for (auto i=0; i<n; i++)
            x[i] = args[0][i];
    }
    if (init_z){
        for (auto i=0; i<n; i++){
            z_L[i] = 0.;
            z_U[i] = 0.;
        }
        if (init_lambda) {
            for (auto i=0; i<m; i++)
                lambda[i] = 0.;
        }
    }

    if (diagnosticPrintout)
        stream << "Exit get_starting_point" << std::endl;
    
    if (diagnosticPrintout)
        stream << "n: " << n << " init_x: " << init_x << " init_z " << init_z << " m: " << m << " init lambda " << init_lambda << std::endl;

    return true;
};

  /** Method to return the objective value */
  bool eval_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number& obj_value){
    if (diagnosticPrintout)
        stream << "Enter eval_f" << std::endl;
    
    if (new_x) 
        updateX(x);

    matlab::data::TypedArray<double> objOut = fevalWithX(funcs[0]["objective"]);

    obj_value = objOut[0];

    if (diagnosticPrintout)
        stream << "Exit eval_f" << std::endl;

    return true;
  };

  /** Method to return the gradient of the objective */
  bool eval_grad_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number* grad_f) {
    if (diagnosticPrintout)
        stream << "Enter eval_grad_f" << std::endl;

    if (new_x) 
        updateX(x);

    matlab::data::TypedArray<double> gradOut = fevalWithX(funcs[0]["gradient"]);

    for (auto i = 0; i < n; i++)
        grad_f[i] = gradOut[i];

    if (diagnosticPrintout)
        stream << "Exit eval_grad_f" << std::endl;

    return true;
  };

  /** Method to return the constraint residuals */
  bool eval_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Index m, Ipopt::Number* g) {

    if (diagnosticPrintout)
        stream << "Enter eval_g" << std::endl;

    if (new_x) 
        updateX(x);

    matlab::data::TypedArray<double> constOut = fevalWithX(funcs[0]["constraints"]);
    
    for (auto i = 0; i < m; i++)
        g[i] = constOut[i];
    
    if (diagnosticPrintout)
        stream << "Exit eval_g" << std::endl;

    return true;
  };

  /** Method to return:
   *   1) The structure of the jacobian (if "values" is NULL)
   *   2) The values of the jacobian (if "values" is not NULL)
   */
bool eval_jac_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
            Ipopt::Index m, Ipopt::Index nele_jac, Ipopt::Index* iRow, Ipopt::Index *jCol,
            Ipopt::Number* values) 
{
if (diagnosticPrintout)
    stream << "Enter eval_jac_g" << std::endl;
    
    if (nullptr == values){
        matlab::data::Array jacobianStructure = funcs[0]["jacobianstructure"];
        std::vector<matlab::data::Array> args(0);
        std::vector<matlab::data::Array> jacStrOut = feval(jacobianStructure, 1, args);

        matlab::data::SparseArray<double> JacobianStr(std::move(jacStrOut[0]));

        auto i = 0;
        for (matlab::data::TypedIterator<double> it = JacobianStr.begin(); it != JacobianStr.end(); it++){
            iRow[i] = JacobianStr.getIndex(it).first;
            jCol[i] = JacobianStr.getIndex(it).second;
            i++;
        }

    } else {

        if (new_x) 
            updateX(x);

        matlab::data::SparseArray<double> Jacobian = fevalWithX(funcs[0]["jacobian"]);
        auto i = 0;
        for (auto& elem : Jacobian){
            values[i] = elem;
            i++;
        }
    }
    if (diagnosticPrintout)
        stream << "Exit eval_jac_g" << std::endl;
    return true;
};

  /** Method to return:
   *   1) The structure of the hessian of the lagrangian (if "values" is NULL)
   *   2) The values of the hessian of the lagrangian (if "values" is not NULL)
   */
  bool eval_h(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                      Ipopt::Number obj_factor, Ipopt::Index m, const Ipopt::Number* lambda,
                      bool new_lambda, Ipopt::Index nele_hess, Ipopt::Index* iRow,
                      Ipopt::Index* jCol, Ipopt::Number* values) 
        {
        if (diagnosticPrintout)
            stream << "Enter eval_h" << std::endl;

            if (nullptr == values){
                matlab::data::Array hessianStructure = funcs[0]["hessianstructure"];
                std::vector<matlab::data::Array> args(0);
                std::vector<matlab::data::Array> hesStrOut = feval(hessianStructure, 1, args);
                matlab::data::SparseArray<double> HessianStr(std::move(hesStrOut[0]));
                auto i = 0;
                for (matlab::data::TypedIterator<double> it = HessianStr.begin(); it != HessianStr.end(); it++){
                    iRow[i] = HessianStr.getIndex(it).first;
                    jCol[i] = HessianStr.getIndex(it).second;
                    i++;
                }
            } else {
                if (new_x) 
                    updateX(x);

                std::vector<matlab::data::Array> hessianArgs = {
                    args[0], 
                    factory.createScalar(obj_factor), 
                    factory.createArray<double>({static_cast<size_t>(m),1}, lambda, lambda + m)
                };


                std::vector<matlab::data::Array> retVals = feval(funcs[0]["hessian"], 1, hessianArgs);
                matlab::data::SparseArray<double> Hessian = std::move(retVals[0]);
                auto i = 0;
                for (auto& elem : Hessian){
                    values[i] = elem;
                    i++;
                }
            }
            if (diagnosticPrintout)
                stream << "Exit eval_h" << std::endl;
            return true;
        };

  //@}

  bool intermediate_callback(   Ipopt::AlgorithmMode mode, 
                                Ipopt::Index iter, 
                                Ipopt::Number obj_value, 
                                Ipopt::Number inf_pr, 
                                Ipopt::Number inf_du, 
                                Ipopt::Number mu, 
                                Ipopt::Number d_norm, 
                                Ipopt::Number regularization_size, 
                                Ipopt::Number alpha_du, 
                                Ipopt::Number alpha_pr, 
                                Ipopt::Index ls_trials, 
                                const Ipopt::IpoptData* ip_data, 
                                Ipopt::IpoptCalculatedQuantities* ip_cq) 
    {
    if (diagnosticPrintout)
        stream << "Enter intermediate_callback" << std::endl;
      if (intermediateCallback){
            matlab::data::StructArray passVal = factory.createStructArray({0, 0},{});
            Ipopt::TNLPAdapter* tnlp_adapter(nullptr);
            if ( nullptr != ip_cq )
            {
                Ipopt::OrigIpoptNLP* orignlp;
                orignlp = dynamic_cast<Ipopt::OrigIpoptNLP*>(GetRawPtr(ip_cq->GetIpoptNLP()));
                if ( nullptr != orignlp )
                    tnlp_adapter = dynamic_cast<Ipopt::TNLPAdapter*>(GetRawPtr(orignlp->nlp()));
                if ( nullptr != tnlp_adapter ) {
                    passVal = factory.createStructArray(
                            {1, 1},
                            {   "primals", 
                                "duals", 
                                "duallbs", 
                                "dualubs",
                                "iter",
                                "obj_value",
                                "inf_pr",
                                "inf_du",
                                "mu",
                                "d_norm",
                                "regularization_size",
                                "alpha_du",
                                "alpha_pr"
                            }
                            );
                    
                    std::unique_ptr<double[]> buffer_(new double[std::max(2*_n,_m)]);
                    
                    matlab::data::TypedArray<double> primals(factory.createArray<double>({_n,1}));
                    // Weird construct required as x.release.get() somehow kills the mex
                    tnlp_adapter->ResortX(*ip_data->curr()->x(), buffer_.get() );
                    for (auto i = 0; i<_n ; i ++)
                        primals[i] = buffer_[i];
                    passVal[0]["primals"] = std::move(primals);

                    matlab::data::TypedArray<double> duals(factory.createArray<double>({_m,1}));
                    tnlp_adapter->ResortG(*ip_data->curr()->y_c(), *ip_data->curr()->y_d(), buffer_.get() );
                    for (auto i = 0; i<_m; i++)
                        duals[i] = buffer_[i];
                    passVal[0]["duals"]  = std::move(duals);


                    matlab::data::TypedArray<double> duallbs(factory.createArray<double>({_n,1}));
                    matlab::data::TypedArray<double> dualubs(factory.createArray<double>({_n,1}));
                    tnlp_adapter->ResortBnds(   *ip_data->curr()->z_L(), buffer_.get(),
                                                *ip_data->curr()->z_U(), buffer_.get() + _n);
                    for (auto i = 0; i<_n; i++) {
                        duallbs[i] = buffer_[i];
                        dualubs[i] = buffer_[i+_n];
                    }
                    passVal[0]["duallbs"] = std::move(duallbs);
                    passVal[0]["dualubs"] = std::move(dualubs);
                    passVal[0]["iter"] = factory.createScalar(static_cast<double>(iter));
                    passVal[0]["obj_value"] = factory.createScalar(obj_value);
                    passVal[0]["inf_pr"] = factory.createScalar(inf_pr);
                    passVal[0]["inf_du"] = factory.createScalar(inf_du);
                    passVal[0]["mu"] = factory.createScalar(mu);
                    passVal[0]["d_norm"] = factory.createScalar(d_norm);
                    passVal[0]["regularization_size"] = factory.createScalar(regularization_size);
                    passVal[0]["alpha_du"] = factory.createScalar(alpha_du);
                    passVal[0]["alpha_pr"] = factory.createScalar(alpha_pr);

                    std::vector<matlab::data::Array> args({std::move(passVal)});
                    std::vector<matlab::data::Array> retVal = feval(funcs[0]["intermediate"],1,args);

                    matlab::data::TypedArray<bool> retValue = std::move(retVal[0]);
                    
                    if (diagnosticPrintout)
                        stream << "Exit intermediate_callback" << std::endl;
                    
                    return retValue[0];
                }
            }
            if (diagnosticPrintout)
                stream << "Exit intermediate_callback" << std::endl;

            return true;
        } else {
            if (diagnosticPrintout)
                stream << "Exit intermediate_callback" << std::endl;

            return true;
        }

    }

  /** @name Solution Methods */
  //@{
  /** This method is called when the algorithm is complete so the TNLP can store/write the solution */
  void finalize_solution(   Ipopt::SolverReturn status,
                            Ipopt::Index n, const Ipopt::Number* x, const Ipopt::Number* z_L, const Ipopt::Number* z_U,
                            Ipopt::Index m, const Ipopt::Number* g, const Ipopt::Number* lambda,
                            Ipopt::Number obj_value, const Ipopt::IpoptData* ip_data,
				            Ipopt::IpoptCalculatedQuantities* ip_cq) 
    {
        if (diagnosticPrintout)
            stream << "Enter finalize_solution" << std::endl;

        updateX(x);

        matlab::data::TypedArrayRef<double> zL = retStr[0]["z_L"];
        for (auto i = 0; i < n; i++)
            zL[i] = z_L[i];
        matlab::data::TypedArrayRef<double> zU = retStr[0]["z_U"];
        for (auto i = 0; i < n; i++)
            zU[i] = z_U[i];
        matlab::data::TypedArrayRef<double> lam = retStr[0]["lambda"];
        for (auto i = 0; i < m; i++)
            lam[i] = lambda[i];

        if (diagnosticPrintout)
            stream << "Exit finalize_solution" << std::endl;

        };
  //@}
};