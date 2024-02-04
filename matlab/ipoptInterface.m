classdef ipoptInterface < handle
       
    properties (SetAccess = protected)
        fVal (1,1) double
        gCur (:,1) double

        cVal (:,1) double
        cJacCur (:,:) double

        xCur (:,1) double

        xBnd (:,2) double
        cBnd (:,2) double

        zInit (:,2) double
        lambdaInit (:,1) double
        HVal (:,:) double
    end

    methods (Abstract)
        runmodel(obj, x)
        constructHessian(obj,sigma,lambda)
    end

    methods 
                
        function update(obj,x)
            runmodel(obj,x);
        end

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