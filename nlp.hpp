#pragma once
#include "mexFunctions.hpp"
#include "stream.hpp"
#include "IpTNLP.hpp"
#include "IpIpoptApplication.hpp"
#include "IpIpoptCalculatedQuantities.hpp"
#include "IpIpoptData.hpp"
#include "IpTNLPAdapter.hpp"
#include "IpOrigIpoptNLP.hpp"
#include "IpJournalist.hpp"
#include "mex.hpp"
#include "mexAdapter.hpp"
#include <algorithm>
#include <memory>

using namespace Ipopt;


class MatlabJournal : 
    public Journal 
{
private:
    Buffer buffer;
    std::ostream cout;
public:

    MatlabJournal(const std::string& name,
                        EJournalLevel level,
                        std::shared_ptr<matlab::engine::MATLABEngine> mptr) : 
        Journal(name, level),
        cout(&buffer)
        {
            buffer.setMatlabPtr(mptr);
        }

    virtual ~MatlabJournal(void) {}


protected:

    virtual std::string Name() {
        return std::string("MatlabJournal");
    }

    // These functions override the functions in the Journal class.
    virtual
    void
    PrintImpl( EJournalCategory category,
            EJournalLevel    level,
            char const *     str) {
        cout << std::string(str);
    }

    virtual
    void
    PrintfImpl( EJournalCategory category,
                EJournalLevel    level,
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


class myNLP : public TNLP
{
private:
    std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr;
    ArrayFactory factory;
    StructArray funcs;
    StructArray options;
    StructArray retStr;
    size_t _n;
    size_t _m;
    std::vector<Array> args;
    bool isinitialised;
    bool returnHessian;
    bool intermediateCallback;
    Buffer buffer;
    std::ostream stream;

    void Error(std::string message) {
        displayError(matlabPtr, factory, message);
    }

    std::vector<Array> feval(
            const Array &handle, 
            const size_t numReturned, 
            const std::vector<Array> &arguments
            ) {
        std::vector<Array> theArgs({handle});
        theArgs.insert(theArgs.end(), arguments.begin(), arguments.end());
        return matlabPtr->feval(
            matlab::engine::convertUTF8StringToUTF16String("feval"),
            numReturned,
            theArgs
        );
    }

    void updateX(const Number *x){
        for (auto i = 0 ; i < args[0].getNumberOfElements(); i++)
            args[0][i] = x[i];
    }

    std::string getStatus(SolverReturn status){
        std::string retVal;
        switch (status)
        {
        case SUCCESS:
            retVal = "SUCCESS";
            break;
        case MAXITER_EXCEEDED:
            retVal = "MAXITER_EXCEEDED";
            break;
        case CPUTIME_EXCEEDED:
            retVal = "CPUTIME_EXCEEDED";
            break;
        case STOP_AT_TINY_STEP:
            retVal = "STOP_AT_TINY_STEP";
            break;
        case STOP_AT_ACCEPTABLE_POINT:
            retVal = "STOP_AT_ACCEPTABLE_POINT";
            break;
        case LOCAL_INFEASIBILITY:
            retVal = "LOCAL_INFEASIBILITY";
            break;
        case USER_REQUESTED_STOP:
            retVal = "USER_REQUESTED_STOP";
            break;
        case FEASIBLE_POINT_FOUND:
            retVal = "FEASIBLE_POINT_FOUND";
            break;
        case DIVERGING_ITERATES:
            retVal = "DIVERGING_ITERATES";
            break;
        case RESTORATION_FAILURE:
            retVal = "RESTORATION_FAILURE";
            break;
        case ERROR_IN_STEP_COMPUTATION:
            retVal = "ERROR_IN_STEP_COMPUTATION";
            break;
        case INVALID_NUMBER_DETECTED:
            retVal = "INVALID_NUMBER_DETECTED";
            break;
        case TOO_FEW_DEGREES_OF_FREEDOM:
            retVal = "TOO_FEW_DEGREES_OF_FREEDOM";
            break;
        case INVALID_OPTION:
            retVal = "INVALID_OPTION";
            break;
        case OUT_OF_MEMORY:
            retVal = "OUT_OF_MEMORY";
            break;
        case INTERNAL_ERROR:
            retVal = "INTERNAL_ERROR";
            break;
        case UNASSIGNED:
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
        retStr(factory.createStructArray({0,0},{})),
        stream(&buffer), 
        args({factory.createEmptyArray()}) 
    {};

    TypedArray<double> getX(void) {
        return args[0];
    };

    StructArray getInfo(void) {
        return retStr;
    }

    void setup(
            std::shared_ptr<matlab::engine::MATLABEngine>& mptr, 
            TypedArray<double>& _x0, 
            StructArray& _funcs, 
            StructArray& _options
        ) {
        matlabPtr = mptr;
        buffer.setMatlabPtr(matlabPtr);

        args[0] = _x0;

        funcs = _funcs;
        options = _options;

        if (isfield(options, "ipopt")){
            StructArray ipoptStr = options[0]["ipopt"];
            if (isfield(ipoptStr, "hessian_approximation")){
                if (isstring(ipoptStr[0]["hessian_approximation"])  ) {
                    returnHessian = (
                        0 != getStringValue(ipoptStr[0]["hessian_approximation"])
                                .compare("limited-memory")
                                );
                } else {
                    Error("limited-memory must be a string type");
                }
            }
        }

        isinitialised = isfield(funcs, "intermediate") && ishandle(funcs[0]["intermediate"]);

    }

    Array fevalWithX(const Array& handle) {
        std::vector<Array> retVal = feval(handle, 1, args);
        return retVal[0];
    }


  /** default destructor */
  ~myNLP() {};

  /**@name Overloaded from TNLP */
  //@{
  /** Method to return some info about the nlp */
  bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                    Index& nnz_h_lag, IndexStyleEnum& index_style) 
    {
    
        index_style = C_STYLE;
        n = args[0].getNumberOfElements();
        _n = n;

        TypedArray<double> cu = options[0]["cu"];
        m = cu.getNumberOfElements();
        _m = m;

        std::vector<Array> nullArgs = {};

        std::vector<Array> jacStrOut = feval(funcs[0]["jacobianstructure"], 1, nullArgs);
        
        if ( !issparse(jacStrOut[0]) ) {
            Error("Jacobian has to be sparse");
        }

        SparseArray<double> Jac = std::move(jacStrOut[0]);

        nnz_jac_g = Jac.getNumberOfNonZeroElements();

        if (returnHessian){
            std::vector<Array> hesStrOut = feval(funcs[0]["hessianstructure"],1,nullArgs);
            if ( !issparse(hesStrOut[0]) ) {
                Error("Hessian has to be sparse");
            }
            SparseArray<double> Hes = std::move(hesStrOut[0]);
            nnz_h_lag = Hes.getNumberOfNonZeroElements();
        }
        

        return true;
    };

  /** Method to return the bounds for my problem */
  bool get_bounds_info(Index n, Number* x_l, Number* x_u,
                Index m, Number* g_l, Number* g_u) 
    {
                                   
        TypedArray<double> xl = std::move(options[0]["lb"]);
        TypedArray<double> xu = std::move(options[0]["ub"]);
        TypedArray<double> gl = std::move(options[0]["cl"]);
        TypedArray<double> gu = std::move(options[0]["cu"]);

        int i;

        for (i=0; i<n; i++)
            x_l[i] = xl[i];

        for (i=0; i<n; i++)
            x_u[i] = xu[i];

        for (i=0; i<m; i++)
            g_l[i] = gl[i];

        for (i=0; i<m; i++)
            g_u[i] = gu[i];
            
        return true;
    };

  /** Method to return the starting point for the algorithm */
  bool get_starting_point(  Index n, bool init_x, Number* x,
                            bool init_z, Number* z_L, Number* z_U,
                            Index m, bool init_lambda,
                            Number* lambda) 
{   
    int i;
    if (init_x){
        for (i=0; i<n; i++)
            x[i] = args[0][i];
    }
    if (init_z){
        for (i=0;i<n;i++){
            z_L[i] = 0.;
            z_U[i] = 0.;
        }
    if (init_lambda) {
        for (i=0;i<m;i++)
            lambda[i] = 0.;
    }

    }

    return true;
};

  /** Method to return the objective value */
  bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value){

    if (new_x) 
        updateX(x);

    TypedArray<double> objOut = fevalWithX(funcs[0]["objective"]);

    obj_value = objOut[0];

    return true;
  };

  /** Method to return the gradient of the objective */
  bool eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f) {

    if (new_x) 
        updateX(x);

    TypedArray<double> gradOut = fevalWithX(funcs[0]["gradient"]);

    for (auto i = 0; i < n; i++)
        grad_f[i] = gradOut[i];

    return true;
  };

  /** Method to return the constraint residuals */
  bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g) {

    if (new_x) 
        updateX(x);

    TypedArray<double> constOut = fevalWithX(funcs[0]["constraints"]);
    
    for (auto i = 0; i < m; i++)
        g[i] = constOut[i];
    
    return true;
  };

  /** Method to return:
   *   1) The structure of the jacobian (if "values" is NULL)
   *   2) The values of the jacobian (if "values" is not NULL)
   */
  bool eval_jac_g(Index n, const Number* x, bool new_x,
                Index m, Index nele_jac, Index* iRow, Index *jCol,
                Number* values) 
        {
            int i;
            if (nullptr == values){
                Array jacobianStructure = funcs[0]["jacobianstructure"];
                std::vector<Array> args(0);
                std::vector<Array> jacStrOut = feval(jacobianStructure, 1, args);

                SparseArray<double> JacobianStr(std::move(jacStrOut[0]));

                i = 0;
                for (TypedIterator<double> it = JacobianStr.begin(); it != JacobianStr.end(); it++){
                    iRow[i] = JacobianStr.getIndex(it).first;
                    jCol[i] = JacobianStr.getIndex(it).second;
                    i++;
                }

            } else {

                if (new_x) 
                    updateX(x);

                SparseArray<double> Jacobian = fevalWithX(funcs[0]["jacobian"]);;
                i = 0;
                for (auto elem : Jacobian){
                    values[i] = elem;
                    i++;
                }
            }

            return true;
        };

  /** Method to return:
   *   1) The structure of the hessian of the lagrangian (if "values" is NULL)
   *   2) The values of the hessian of the lagrangian (if "values" is not NULL)
   */
  bool eval_h(Index n, const Number* x, bool new_x,
                      Number obj_factor, Index m, const Number* lambda,
                      bool new_lambda, Index nele_hess, Index* iRow,
                      Index* jCol, Number* values) 
        {
            int i;
            assert(!returnHessian);
            if (nullptr == values){
                Array hessianStructure = funcs[0]["hessianstructure"];
                std::vector<Array> args(0);
                std::vector<Array> hesStrOut = feval(hessianStructure, 1, args);
                SparseArray<double> HessianStr(std::move(hesStrOut[0]));
                i = 0;
                for (TypedIterator<double> it = HessianStr.begin(); it != HessianStr.end(); it++){
                    iRow[i] = HessianStr.getIndex(it).first;
                    jCol[i] = HessianStr.getIndex(it).second;
                    i++;
                }
            } else {
            if (new_x) 
                updateX(x);

            std::vector<Array> hessianArgs = {
                args[0], 
                factory.createScalar(obj_factor), 
                factory.createArray<double>({static_cast<size_t>(m),1}, lambda, lambda + m)
            };


            std::vector<Array> retVals = feval(funcs[0]["hessian"], 1, hessianArgs);
            SparseArray<double> Hessian = std::move(retVals[0]);
            i = 0;
            for (auto elem : Hessian){
                values[i] = elem;
                i++;
            }
            }
            return true;
        };

  //@}

  bool intermediate_callback(   AlgorithmMode mode, 
                                Index iter, 
                                Number obj_value, 
                                Number inf_pr, 
                                Number inf_du, 
                                Number mu, 
                                Number d_norm, 
                                Number regularization_size, 
                                Number alpha_du, 
                                Number alpha_pr, 
                                Index ls_trials, 
                                const IpoptData* ip_data, 
                                IpoptCalculatedQuantities* ip_cq) 
    {

      if (intermediateCallback){
            StructArray passVal = factory.createStructArray({0, 0},{});
            TNLPAdapter* tnlp_adapter(nullptr);
            if ( nullptr != ip_cq )
            {
                OrigIpoptNLP* orignlp;
                orignlp = dynamic_cast<OrigIpoptNLP*>(GetRawPtr(ip_cq->GetIpoptNLP()));
                if ( nullptr != orignlp )
                    tnlp_adapter = dynamic_cast<TNLPAdapter*>(GetRawPtr(orignlp->nlp()));
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
                    
                    TypedArray<double> primals(factory.createArray<double>({_n,1}));
                    // Weird construct required as x.release.get() somehow kills the mex
                    tnlp_adapter->ResortX(*ip_data->curr()->x(), buffer_.get() );
                    for (auto i = 0; i<_n ; i ++)
                        primals[i] = buffer_[i];
                    passVal[0]["primals"] = std::move(primals);

                    TypedArray<double> duals(factory.createArray<double>({_m,1}));
                    tnlp_adapter->ResortG(*ip_data->curr()->y_c(), *ip_data->curr()->y_d(), buffer_.get() );
                    for (auto i = 0; i<_m; i++)
                        duals[i] = buffer_[i];
                    passVal[0]["duals"]  = std::move(duals);


                    TypedArray<double> duallbs(factory.createArray<double>({_n,1}));
                    TypedArray<double> dualubs(factory.createArray<double>({_n,1}));
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

                    std::vector<Array> args({std::move(passVal)});
                    std::vector<Array> retVal = feval(funcs[0]["intermediate"],1,args);

                    TypedArray<bool> retValue = std::move(retVal[0]);

                    return retValue[0];
                }
            }
            return true;
        } else {
            return true;
        }

    }

  /** @name Solution Methods */
  //@{
  /** This method is called when the algorithm is complete so the TNLP can store/write the solution */
  void finalize_solution(   SolverReturn status,
                            Index n, const Number* x, const Number* z_L, const Number* z_U,
                            Index m, const Number* g, const Number* lambda,
                            Number obj_value, const IpoptData* ip_data,
				            IpoptCalculatedQuantities* ip_cq) 
    {

        updateX(x);

        retStr = factory.createStructArray({1,1}, {"z_L", "z_U", "lambda", "solverStatus", "infoStatus"});
        retStr[0]["z_L"] = factory.createArray<double>({static_cast<size_t>(n),1});
        retStr[0]["z_U"] = factory.createArray<double>({static_cast<size_t>(n),1});
        retStr[0]["lambda"] = factory.createArray<double>({static_cast<size_t>(m),1});
        
        retStr[0]["solverStatus"] = factory.createCharArray(getStatus(status));
        TypedArrayRef<double> zL = retStr[0]["z_L"];
        for (auto i = 0; i < n; i++)
            zL[i] = z_L[i];
        TypedArrayRef<double> zU = retStr[0]["z_U"];
        for (auto i = 0; i < n; i++)
            zU[i] = z_U[i];
        TypedArrayRef<double> lam = retStr[0]["lambda"];
        for (auto i = 0; i < m; i++)
            lam[i] = lambda[i];

        };
  //@}
};