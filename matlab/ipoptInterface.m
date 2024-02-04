classdef ipoptInterface < handle
       
    properties (SetAccess = protected)
        fVal (1,1) double       % Objective value
        gCur (:,1) double       % Objective gradient

        cVal (:,1) double       % Constraint values
        cJacCur (:,:) double    % Constraint Jacobian

        xCur (:,1) double       % Current decision variable (Initially x0)

        xBnd (:,2) double       % Bounds on decision variable [lb,ub]
        cBnd (:,2) double       % Bounds on constraints [lb,ub]
        
        % Optional
        zInit (:,2) double      % Initial multipliers for decision variable
        lambdaInit (:,1) double % Initial multipliers for constraints
        HVal (:,:) double       % Hessian matrix (lower triangular)
    end

    methods (Abstract)
        % Single callback to update internal model and set function values
        % and derivatives etc.
        update(obj,x)
        % Callback for updating internal Hessian HVal with appropriate
        % lower triangular Hessian, given objective weight sigma and
        % multipliers lambda.
        constructHessian(obj,sigma,lambda)
    end

    methods 
        % Prototype intermediate callback. Return value can terminate
        % solve.
        function bContinue = intermediateCallback(~,s)
            bContinue = true;
        end
    end

    methods (Static)
        function opt = getDefaultOptions()
            opt.ipopt.tol                               = 1e-8;
            opt.ipopt.max_iter                          = 1500;
            opt.ipopt.mu_strategy                       = 'adaptive';
            opt.ipopt.linear_solver                     = 'ma27';
            opt.ipopt.print_level                       = 5;
            opt.ipopt.print_user_options                = 'yes';
            opt.ipopt.acceptable_tol                    = 10*opt.ipopt.tol;
            opt.ipopt.acceptable_iter                   = 3;
            opt.ipopt.acceptable_dual_inf_tol           = opt.ipopt.acceptable_tol;
            opt.ipopt.hessian_approximation             = "exact";
        end
    end
end