#pragma once
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


class myNLP 
    : public Ipopt::TNLP
{
private:
    
    std::size_t _n{};
    std::size_t _m{};
    bool isinitialised{ false };
    bool returnHessian{ false };
    bool intermediateCallback{ false };
    bool diagnosticPrintout{ false };
    matlab::data::ArrayFactory factory;
    matlab::data::StructArray funcs;
    matlab::data::StructArray options;
    matlab::data::StructArray retStr;
    std::vector<matlab::data::TypedArray<double>> args;
    Buffer buffer;
    std::ostream stream;
    utilities::Sparse<double> _jac{};
    utilities::Sparse<double> _hes{};

    matlab::data::TypedArray<double> _x;
    matlab::data::TypedArray<double> _sigma;
    matlab::data::TypedArray<double> _lambda;
    

    void updateX(const Ipopt::Number *x){
        std::size_t i = 0;
        for (auto &elem : _x)
            elem = x[i++];
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
        case Ipopt::WALLTIME_EXCEEDED:
            retVal = "WALLTIME_EXCEEDED";
            break;
        default:
            retVal = "UNKOWN";
            break;
        }
        return retVal;
    }

public:
    
    myNLP(const myNLP&) = delete;
    myNLP operator=(const myNLP&) = delete;

    myNLP() : _n(0), _m(0), isinitialised(false), returnHessian(true), intermediateCallback(false), diagnosticPrintout(false),
        funcs(factory.createStructArray({0,0},{})),
        options(factory.createStructArray({0,0},{})),
        retStr(factory.createStructArray({1,1}, {"z_L", "z_U", "lambda", "status", "iter","cpu","objective","eval"})), // Make sure the structure is always available
        stream(&buffer), 
        _jac(),_hes(),_x(factory.createArray<double>({0,1})), _sigma(factory.createScalar(1.)),_lambda(factory.createArray<double>({0,1}))
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
    }

    matlab::data::TypedArray<double> getX(void) {
        return _x;
    };

    matlab::data::StructArray getInfo(void) {
        return retStr;
    }

    void setup( 
            matlab::data::TypedArray<double>& _x0, 
            matlab::data::StructArray& _funcs, 
            matlab::data::StructArray& _options
        ) {
        
        _x = std::move(_x0);

        funcs = std::move(_funcs);
        options = std::move(_options);

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

        std::vector<matlab::data::Array> retVal = utilities::feval(handle, 1, {_x});
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
      // Robustify this process by a huge margin! Hold an internal representation to verify that the dimensions of all passed arguments work for all function handles.
        if (diagnosticPrintout)
            stream << "Enter get_nlp_info" << std::endl;
        index_style = Ipopt::TNLP::C_STYLE;
        n = static_cast<Ipopt::Index>(_x.getNumberOfElements());
        _n = _x.getNumberOfElements();

        matlab::data::TypedArray<double> cu = utilities::getfield(options,"cu");
        m = static_cast<Ipopt::Index>(cu.getNumberOfElements());
        _m = cu.getNumberOfElements();
        
        // Initialise return arrays in case things fail.
        retStr[0]["z_L"] = factory.createArray<double>({static_cast<size_t>(n),1});
        retStr[0]["z_U"] = factory.createArray<double>({static_cast<size_t>(n),1});
        retStr[0]["lambda"] = factory.createArray<double>({static_cast<size_t>(m),1});
        

        std::vector<matlab::data::Array> nullArgs{};

        std::vector<matlab::data::Array> jacStrOut = utilities::feval(funcs[0]["jacobianstructure"], 1, nullArgs);
        
        if ( !utilities::issparse(jacStrOut[0])) {
            utilities::error("Jacobian pattern must be sparse matrix");
        }

        matlab::data::SparseArray<double> jacobian = std::move(jacStrOut[0]);
        _jac.set(jacobian);
        
        if (_jac.getNumberOfColumns() != _n || _jac.getNumberOfRows() != _m)
            utilities::error("Jacobian structure must be {} x {}, passed matrix is {} x {}", _m, _n, _jac.getNumberOfRows(), _jac.getNumberOfColumns());

        nnz_jac_g = static_cast<Ipopt::Index>(_jac.getNumberOfNonZeroElements());

        if (returnHessian){
            std::vector<matlab::data::Array> hesStrOut = utilities::feval(funcs[0]["hessianstructure"],1,nullArgs);
            if ( !utilities::issparse(hesStrOut[0])) {
                utilities::error("Hessian pattern must be sparse matrix");
            }
            matlab::data::SparseArray<double> hessian = std::move(hesStrOut[0]);
            _hes.set(hessian);
            if (_hes.getNumberOfColumns() != _n || _hes.getNumberOfRows() != _hes.getNumberOfColumns())
                utilities::error("Hessian structure must be {} x {}, passed matrix is {} x {}", _n, _n, _hes.getNumberOfRows(), _hes.getNumberOfColumns());

            nnz_h_lag = static_cast<Ipopt::Index>(_hes.getNumberOfNonZeroElements());
        }
        if (diagnosticPrintout)
            stream << "Exit get_nlp_info" << std::endl;

        return true;
    }

  /** Method to return the bounds for my problem */
  bool get_bounds_info([[maybe_unused]]Ipopt::Index n, Ipopt::Number* x_l, Ipopt::Number* x_u,
                [[maybe_unused]]Ipopt::Index m, Ipopt::Number* g_l, Ipopt::Number* g_u) 
    {
        if (diagnosticPrintout)
            stream << "Enter get_bounds_info" << std::endl;
        matlab::data::TypedArray<double> xl = utilities::getfield(options,"lb");
        matlab::data::TypedArray<double> xu = utilities::getfield(options,"ub");
        matlab::data::TypedArray<double> gl = utilities::getfield(options,"cl");
        matlab::data::TypedArray<double> gu = utilities::getfield(options,"cu");

        if (xl.getNumberOfElements() != _n)
            utilities::error("Size mismatch: Lower bound on x has {} elements whereas x0 has {}.", xl.getNumberOfElements(), _n);
        if (xu.getNumberOfElements() != _n)
            utilities::error("Size mismatch: Upper bound on x has {} elements whereas x0 has {}.", xu.getNumberOfElements(), _n);
        if (gl.getNumberOfElements() != _m)
            utilities::error("Size mismatch: Lower bound on constraints has {} elements whereas {} are expected.", gl.getNumberOfElements(), _m);
        if (gu.getNumberOfElements() != _m)
            utilities::error("Size mismatch: Upper bound on constraints has {} elements whereas {} are expected.", gu.getNumberOfElements(), _m);

        std::copy(xl.cbegin(), xl.cend(), x_l);
        std::copy(xu.cbegin(), xu.cend(), x_u);
        std::copy(gl.cbegin(), gl.cend(), g_l);
        std::copy(gu.cbegin(), gu.cend(), g_u);
            
        if (diagnosticPrintout)
            stream << "Exit get_bounds_info" << std::endl;    
        return true;
    }

  /** Method to return the starting point for the algorithm */
  bool get_starting_point(  Ipopt::Index n, bool init_x, Ipopt::Number* x,
                            bool init_z, Ipopt::Number* z_L, Ipopt::Number* z_U,
                            Ipopt::Index m, bool init_lambda,
                            Ipopt::Number* lambda) 
{   
    if (diagnosticPrintout)
        stream << "Enter get_starting_point" << std::endl;
    if (init_x){
        std::copy(_x.cbegin(), _x.cend(), x);
    }
    if (init_z){
        // We should be able to warmstart duals on the decision variable
        std::fill(z_L, z_L + n, 0.);
        std::fill(z_U, z_U + n, 0.);
        // as well as the constrain multipliers. Figure out a way of passing them in!
        if (init_lambda) {
            std::fill(lambda, lambda + m, 0.);
        }
    }

    if (diagnosticPrintout)
        stream << "Exit get_starting_point" << std::endl;
    
    if (diagnosticPrintout)
        stream << "n: " << n << " init_x: " << init_x << " init_z " << init_z << " m: " << m << " init lambda " << init_lambda << std::endl;

    return true;
}

  /** Method to return the objective value */
  bool eval_f([[maybe_unused]]Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number& obj_value){
    if (diagnosticPrintout)
        stream << "Enter eval_f" << std::endl;
    
    if (new_x) 
        updateX(x);

    matlab::data::TypedArray<double> objOut = fevalWithX(funcs[0]["objective"]);

    obj_value = objOut[0];

    if (diagnosticPrintout)
        stream << "Exit eval_f" << std::endl;

    return true;
  }

  /** Method to return the gradient of the objective */
  bool eval_grad_f([[maybe_unused]]Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number* grad_f) {
    if (diagnosticPrintout)
        stream << "Enter eval_grad_f" << std::endl;

    if (new_x) 
        updateX(x);

    matlab::data::TypedArray<double> gradOut = fevalWithX(funcs[0]["gradient"]);

    if (gradOut.getNumberOfElements() != _n)
        utilities::error("Size mismatch: gradient callback returned {} elements, whereas {} are expected.", gradOut.getNumberOfElements(), _n);

    std::copy(gradOut.cbegin(), gradOut.cend(), grad_f);

    if (diagnosticPrintout)
        stream << "Exit eval_grad_f" << std::endl;

    return true;
  }

  /** Method to return the constraint residuals */
  bool eval_g([[maybe_unused]]Ipopt::Index n, const Ipopt::Number* x, bool new_x, [[maybe_unused]]Ipopt::Index m, Ipopt::Number* g) {

    if (diagnosticPrintout)
        stream << "Enter eval_g" << std::endl;

    if (new_x) 
        updateX(x);

    matlab::data::TypedArray<double> constOut = fevalWithX(funcs[0]["constraints"]);
    
    if (constOut.getNumberOfElements() != _m)
        utilities::error("Size mismatch: constraint callback returned {} elements, whereas {} are expected.", constOut.getNumberOfElements(), _m);

    std::copy(constOut.cbegin(), constOut.cend(), g);
    
    if (diagnosticPrintout)
        stream << "Exit eval_g" << std::endl;

    return true;
  }

  /** Method to return:
   *   1) The structure of the jacobian (if "values" is NULL)
   *   2) The values of the jacobian (if "values" is not NULL)
   */
bool eval_jac_g([[maybe_unused]]Ipopt::Index n, const Ipopt::Number* x, bool new_x,
            [[maybe_unused]]Ipopt::Index m, [[maybe_unused]]Ipopt::Index nele_jac, Ipopt::Index* iRow, Ipopt::Index *jCol,
            Ipopt::Number* values) 
{
if (diagnosticPrintout)
    stream << "Enter eval_jac_g" << std::endl;
    
    if (nullptr == values){
        _jac.iRow(iRow);
        _jac.jCol(jCol);
    } else {

        if (new_x) 
            updateX(x);

        matlab::data::SparseArray<double> jVals = fevalWithX(funcs[0]["jacobian"]);
        _jac.updateValues(jVals);
        _jac.val(values);
    }
    if (diagnosticPrintout)
        stream << "Exit eval_jac_g" << std::endl;
    return true;
}

  /** Method to return:
   *   1) The structure of the hessian of the lagrangian (if "values" is NULL)
   *   2) The values of the hessian of the lagrangian (if "values" is not NULL)
   */
  bool eval_h([[maybe_unused]]Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                      Ipopt::Number obj_factor, [[maybe_unused]]Ipopt::Index m, const Ipopt::Number* lambda,
                      bool new_lambda, [[maybe_unused]]Ipopt::Index nele_hess, Ipopt::Index* iRow,
                      Ipopt::Index* jCol, Ipopt::Number* values) 
        {
        if (diagnosticPrintout)
            stream << "Enter eval_h" << std::endl;

            if (nullptr == values){
                _hes.iRow(iRow);
                _hes.jCol(jCol);
                _lambda = factory.createArray<double>({ _hes.getNumberOfColumns(),1});
            } else {
                if (new_x) 
                    updateX(x);
                // No point providing a more formal mechanism to update things that are only used in here.
                _sigma = factory.createScalar(obj_factor);
                if (new_lambda)
                    std::copy(lambda, lambda + _m, _lambda.begin());

                std::vector<matlab::data::Array> hessianArgs{
                    _x, 
                    _sigma,
                    _lambda
                };

                std::vector<matlab::data::Array> retVals = utilities::feval(funcs[0]["hessian"], 1, hessianArgs);
                matlab::data::SparseArray<double> hessian = std::move(retVals[0]);
                _hes.updateValues(hessian);
                _hes.val(values);
            }
            if (diagnosticPrintout)
                stream << "Exit eval_h" << std::endl;
            return true;
        }

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
                                [[maybe_unused]]Ipopt::Index ls_trials, 
                                [[maybe_unused]]const Ipopt::IpoptData* ip_data, 
                                [[maybe_unused]]Ipopt::IpoptCalculatedQuantities* ip_cq) 
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
                                "alpha_pr",
                                "mode"
                            }
                            );
                    
                    matlab::data::buffer_ptr_t<double> primals_p = factory.createBuffer<double>(_n);
                    matlab::data::buffer_ptr_t<double> duals_p = factory.createBuffer<double>(_m);
                    matlab::data::buffer_ptr_t<double> duallbs_p = factory.createBuffer<double>(_n);
                    matlab::data::buffer_ptr_t<double> dualubs_p = factory.createBuffer<double>(_n);
                    
                    tnlp_adapter->ResortX(*ip_data->curr()->x(), primals_p.get() );
                    tnlp_adapter->ResortG(*ip_data->curr()->y_c(), *ip_data->curr()->y_d(), duals_p.get() );
                    tnlp_adapter->ResortBounds(   *ip_data->curr()->z_L(), duallbs_p.get(),
                                                *ip_data->curr()->z_U(), dualubs_p.get());

                    passVal[0]["primals"] = factory.createArrayFromBuffer<double>({ _n,1 }, std::move(primals_p));
                    passVal[0]["duals"] = factory.createArrayFromBuffer<double>({ _m,1 }, std::move(duals_p));
                    passVal[0]["duallbs"] = factory.createArrayFromBuffer<double>({ _n,1 }, std::move(duallbs_p));
                    passVal[0]["dualubs"] = factory.createArrayFromBuffer<double>({ _n,1 }, std::move(dualubs_p));
                    passVal[0]["iter"] = factory.createScalar(static_cast<double>(iter));
                    passVal[0]["obj_value"] = factory.createScalar(obj_value);
                    passVal[0]["inf_pr"] = factory.createScalar(inf_pr);
                    passVal[0]["inf_du"] = factory.createScalar(inf_du);
                    passVal[0]["mu"] = factory.createScalar(mu);
                    passVal[0]["d_norm"] = factory.createScalar(d_norm);
                    passVal[0]["regularization_size"] = factory.createScalar(regularization_size);
                    passVal[0]["alpha_du"] = factory.createScalar(alpha_du);
                    passVal[0]["alpha_pr"] = factory.createScalar(alpha_pr);
                    passVal[0]["mode"] = factory.createScalar(std::string(mode == Ipopt::AlgorithmMode::RegularMode ? "regular" : "restoration"));

                    std::vector<matlab::data::Array> intermediateArgs{passVal};
                    std::vector<matlab::data::Array> retVal = utilities::feval(funcs[0]["intermediate"],1, intermediateArgs);

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
  void finalize_solution(   [[maybe_unused]]Ipopt::SolverReturn status,
                            [[maybe_unused]]Ipopt::Index n, const Ipopt::Number* x, const Ipopt::Number* z_L, const Ipopt::Number* z_U,
                            [[maybe_unused]]Ipopt::Index m, [[maybe_unused]]const Ipopt::Number* g, const Ipopt::Number* lambda,
                            [[maybe_unused]]Ipopt::Number obj_value, [[maybe_unused]]const Ipopt::IpoptData* ip_data,
				            [[maybe_unused]]Ipopt::IpoptCalculatedQuantities* ip_cq) 
    {
        if (diagnosticPrintout)
            stream << "Enter finalize_solution" << std::endl;

        updateX(x);

        matlab::data::TypedArrayRef<double> zL = retStr[0]["z_L"];
        matlab::data::TypedArrayRef<double> zU = retStr[0]["z_U"];
        matlab::data::TypedArrayRef<double> lam = retStr[0]["lambda"];

        std::copy(z_L, z_L + _n, zL.begin());
        std::copy(z_U, z_U + _n, zU.begin());
        std::copy(lambda, lambda + _m, lam.begin());
        
        if (diagnosticPrintout)
            stream << "Exit finalize_solution" << std::endl;

    }

  //@}

  bool operator==(const myNLP& other) const = delete;

};