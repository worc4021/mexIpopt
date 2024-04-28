classdef BaseProblem < matlab.unittest.TestCase & ProblemData
    
    properties (TestParameter)
        linear_solver = {"ma27", "mumps", "pardisomkl"};
        hessian_approximation = {"exact","limited-memory"};
        warm_start_init_point = {"no","yes"};
    end
    
    methods (Test)
        function testSolve(testCase, linear_solver, hessian_approximation, warm_start_init_point)
            
            problem = struct();
            problem.variableInfo = struct();
            problem.variableInfo.lb = testCase.xBnd(:,1);  % Lower bound on the variables.
            problem.variableInfo.ub = testCase.xBnd(:,2);  % Upper bound on the variables.
            problem.variableInfo.cl = testCase.cBnd(:,1);   % Lower bounds on the constraint functions.
            problem.variableInfo.cu = testCase.cBnd(:,2);   % Upper bounds on the constraint functions.
            problem.variableInfo.x0 = testCase.x0;  % The starting point.
            % Initialize the dual point.
            problem.variableInfo.zl     = testCase.z0(:,1);
            problem.variableInfo.zu     = testCase.z0(:,2);
            problem.variableInfo.lambda = testCase.lambda0;
            
            % Set the IPOPT options.
            problem.ipopt.tol                               = 1e-7;
            problem.ipopt.max_iter                          = 1000;
            problem.ipopt.mu_strategy                       = "adaptive";
            problem.ipopt.linear_solver                     = linear_solver;
            problem.ipopt.print_level                       = 5;
            problem.ipopt.hessian_approximation             = hessian_approximation;
            % problem.ipopt.print_options_documentation       = "yes";
            problem.ipopt.warm_start_init_point             = warm_start_init_point;
            % problem.ipopt.print_user_options                = "yes";
            % problem.ipopt.debug                             = 0;
            problem.ipopt.tol                               = 1e-5;
            problem.ipopt.max_iter                          = 1500;
            problem.ipopt.acceptable_tol                    = 1e-4;
            problem.ipopt.acceptable_iter                   = 3;
            problem.ipopt.acceptable_dual_inf_tol           = problem.ipopt.acceptable_tol;
            
            % The callback functions.
            problem.funcs.objective         = @testCase.objective;
            problem.funcs.constraints       = @testCase.constraints;
            problem.funcs.gradient          = @testCase.gradient;
            problem.funcs.jacobian          = @testCase.jacobian;
            problem.funcs.jacobianstructure = @testCase.jacobianstructure;
            problem.funcs.hessian           = @testCase.hessian;
            problem.funcs.hessianstructure  = @testCase.hessianstructure;
            problem.funcs.intermediate      = @testCase.intermediateCallback;

            % Run IPOPT.
            [~, info] = ipopt(problem);

            testCase.verifyThat(info.status, ...
                                IsMemberOf(["SUCCESS","STOP_AT_ACCEPTABLE_POINT"]), ...
                                "Exit status");

            testCase.verifyThat(info.appstatus, ...
                                 IsMemberOf(["Solve_Succeeded","Solved_To_Acceptable_Level"]), ...
                                 "Application Status");
        end
    end

    methods (Abstract)
        fVal = objective(obj, x)
        fJac = gradient(obj, x)
        cVal = constraints(obj, x)
        cJac = jacobian(obj, x)
        cStr = jacobianstructure(obj)
        % Hessian in lower triangular form
        hVal = hessian(obj, x,sigma,lambda)
        hStr = hessianstructure(obj)
    end

    methods (Hidden)
        function bContinue = intermediateCallback(~, ~)
            % "primals", "duals", "duallbs", "dualubs","iter","obj_value","inf_pr","inf_du","mu","d_norm","regularization_size","alpha_du","alpha_pr","mode"

            bContinue = true;
        end
    end
end

