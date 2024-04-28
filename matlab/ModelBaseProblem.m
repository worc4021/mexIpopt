classdef ModelBaseProblem < matlab.unittest.TestCase & ipoptInterface
    
    properties (TestParameter)
        linear_solver = {"ma27", "mumps", "pardisomkl"};
        hessian_approximation = {"exact","limited-memory"};
        warm_start_init_point = {"no","yes"};
    end
    
    methods(Test)
        function testSolve(testCase, linear_solver, hessian_approximation, warm_start_init_point)
            opt = testCase.getDefaultOptions();
            opt.ipopt.linear_solver         = linear_solver;
            opt.ipopt.hessian_approximation = hessian_approximation;
            opt.ipopt.warm_start_init_point = warm_start_init_point;

            testCase.cJacCur = testCase.jacobianstructure();
            testCase.HVal = testCase.hessianstructure();

            [retStr,stat] = model_eval(testCase, opt);

            testCase.verifyThat(retStr.status, ...
                                IsMemberOf(["SUCCESS","STOP_AT_ACCEPTABLE_POINT"]), ...
                                "Exit status");

            testCase.verifyThat(stat.appstatus, ...
                                 IsMemberOf(["Solve_Succeeded","Solved_To_Acceptable_Level"]), ...
                                 "Application Status");
        end
    end
    methods (Hidden)
        function update(obj,x)
            obj.x0 = x;
            obj.fVal = obj.objective(x);
            obj.gCur = obj.gradient(x);
            obj.cVal = obj.constraints(x);
            obj.cJacCur = obj.jacobian(x);

        end
        
        function constructHessian(obj,sigma,lambda)
            obj.HVal = obj.hessian(obj.x0,sigma,lambda);
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
end