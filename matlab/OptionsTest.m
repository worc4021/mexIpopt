classdef OptionsTest < BaseBinaryFinder
    
    properties 
        outputfile (1,:) char
    end
    
    properties (TestParameter)
        linear_solver = {"mumps", "ma27", "ma57", "ma77", "ma97"};
        hessian_approximation = {"exact","limited-memory"};
    end

    methods (TestMethodSetup)
        function getUniqueFilename(testCase)
            testCase.outputfile = [tempname,'.ipoptout'];
            testCase.addTeardown(@testCase.trydelete,testCase.outputfile);
        end
    end


    methods(Test)
        function baseline(testCase)
            problem = testCase.getProblem();
            
            % Run IPOPT.
            [~, info] = ipopt(problem);
            
            testCase.verifyTrue(ismember(info.status, ...
                                ["SUCCESS","STOP_AT_ACCEPTABLE_POINT"]), ...
                                "Exit status");

            testCase.verifyTrue(ismember(info.appstatus, ...
                                 ["Solve_Succeeded","Solved_To_Acceptable_Level"]), ...
                                 "Application Status");

            testCase.verifyTrue(isfield(info,'build') && isfield(info.build,'number') && isfield(info.build,'hash'),...
                'Build info not returned');
        end

        function intermediateCallbackCanAbort(testCase)
            interCb = @(s)testCase.intermediateCallback(s,"Abort",true,"CorrectFields",true);
            problem = testCase.getProblem("IntermediateCallbackFun",interCb,"UseIntermediateCallback",true);
            
            % Run IPOPT.
            [~, info] = ipopt(problem);
            s = parseOutput(testCase);
            testCase.verifyEqual(info.status,"USER_REQUESTED_STOP",'Intermediate Callback could not abort');
        end

        function warnsOnUnknownOptions(testCase)
            problem = testCase.getProblem();
            problem.ipopt.hello_world = "Wassup";
            testCase.verifyWarning(@()ipopt(problem),'IPOPT:UnknownOption','Unset field does not produce warning');
            s = parseOutput(testCase);
        end

        function canScaleObjectiveByInteger(testCase)
            problem = testCase.getProblem();
            problem.ipopt.obj_scaling_factor = 100;
            [~, info] = ipopt(problem);
            s = parseOutput(testCase);
            testCase.verifyTrue(ismember(info.status, ...
                                ["SUCCESS","STOP_AT_ACCEPTABLE_POINT"]), ...
                                "Exit status");

            testCase.verifyTrue(isfield(s.used_options,'obj_scaling_factor'),'Objective scaling did not make it through');
        end

        function canApplyTolerances(testCase)
            problem = testCase.getProblem();
            problem.ipopt.tol                               = 1e-6; % Relative NLP error
            problem.ipopt.constr_viol_tol                   = 1e-5; % Absolute constraint violation
            problem.ipopt.dual_inf_tol                      = 10*problem.ipopt.constr_viol_tol; % Absolute dual infeasibility
            problem.ipopt.compl_inf_tol                     = 1e-6; % Absolute complementarity
            problem.ipopt.acceptable_tol                    = 10*problem.ipopt.tol; % Relative NLP error
            problem.ipopt.acceptable_constr_viol_tol        = 10*problem.ipopt.constr_viol_tol; % Absolute constraint violation
            problem.ipopt.acceptable_dual_inf_tol           = 10*problem.ipopt.dual_inf_tol; % Absolute dual infeasibility
            problem.ipopt.acceptable_compl_inf_tol          = 10*problem.ipopt.compl_inf_tol; % Absolute complementarity condition
            problem.ipopt.acceptable_obj_change_tol         = 1e-5; 
            problem.ipopt.acceptable_iter                   = 3;
            problem.ipopt.mu_target                         = 1e-6;
            [~, info] = ipopt(problem);
            s = parseOutput(testCase);
            testCase.verifyTrue(ismember(info.status, ...
                                ["SUCCESS","STOP_AT_ACCEPTABLE_POINT"]), ...
                                "Exit status");
            for field = fieldnames(problem.ipopt)'
                testCase.verifyTrue(isfield(s.used_options,field{1}),sprintf('Failed to apply %s',field{1}));
            end
        end

        function matrixOptions(testCase, linear_solver, hessian_approximation)
            problem = testCase.getProblem();
            problem.ipopt.linear_solver = linear_solver;
            problem.ipopt.hessian_approximation = hessian_approximation;
            [~, info] = ipopt(problem);
            s = parseOutput(testCase);
            testCase.verifyTrue(ismember(info.status, ...
                                ["SUCCESS","STOP_AT_ACCEPTABLE_POINT"]), ...
                                "Exit status");
            testCase.verifyTrue(isfield(s.used_options,'linear_solver') ...
                && strcmp(s.used_options.linear_solver,linear_solver),...
                'Failed to set requested linear solver');
            testCase.verifyTrue(isfield(s.used_options,'hessian_approximation') ...
                && strcmp(s.used_options.hessian_approximation,hessian_approximation),...
                'Failed to set requested Hessian approximation methdod');
        end
        
    end

    methods (Hidden)
        function bContinue = intermediateCallback(testCase, sIntermediate, opt)
            arguments
                testCase (1,1) OptionsTest
                sIntermediate (1,1) struct
                opt.Abort (1,1) logical = false;
                opt.CorrectFields (1,1) logical = false;
            end
            if opt.CorrectFields
                for name = ["primals", "duals", "duallbs", "dualubs","iter","obj_value","inf_pr","inf_du","mu","d_norm","regularization_size","alpha_du","alpha_pr","mode"]
                    testCase.verifyTrue(isfield(sIntermediate,name),sprintf('Field "%s" not on intermediate callback structure',name));
                end
            end

            if opt.Abort
                bContinue = false;
            else
                bContinue = true;
            end
        end

        function problem = getProblem(testCase,opt)
            arguments
                testCase (1,1) OptionsTest
                opt.InitialiseDuals (1,1) logical = false;
                opt.UseIntermediateCallback (1,1) logical = false;
                opt.IntermediateCallbackFun (:,1) function_handle = @(~)true;
            end

            problem = struct();
            problem.ipopt = struct();
            problem.ipopt.print_user_options = "yes";
            problem.ipopt.file_print_level = 5;
            problem.ipopt.output_file = testCase.outputfile;
        
            problem.variableInfo = struct();
            problem.variableInfo.lb = ones(4,1);
            problem.variableInfo.ub = ones(4,1)*5;
            problem.variableInfo.cl = [25;40];
            problem.variableInfo.cu = [inf;40];
            problem.variableInfo.x0 = [1; 5; 3; 1];
            
            if opt.InitialiseDuals
                % Initialize the dual point.
                problem.variableInfo.zl     = [0;0;0;1];
                problem.variableInfo.zu     = [0;1;0;0];
                problem.variableInfo.lambda = [-1;-1];
            end

            % The callback functions.
            problem.funcs.objective         = @OptionsTest.objective;
            problem.funcs.constraints       = @OptionsTest.constraints;
            problem.funcs.gradient          = @OptionsTest.gradient;
            problem.funcs.jacobian          = @OptionsTest.jacobian;
            problem.funcs.jacobianstructure = @OptionsTest.jacobianstructure;
            problem.funcs.hessian           = @OptionsTest.hessian;
            problem.funcs.hessianstructure  = @OptionsTest.hessianstructure;
            if ~isempty(opt.UseIntermediateCallback)
                problem.funcs.intermediate      = opt.IntermediateCallbackFun;
            end
        end

        function s = parseOutput(testCase)
            text = fileread(testCase.outputfile);
            s = regexp(text,'.*This\s+is\s+Ipopt\s+version\s(?<ipoptversion>[0-9\.]+),\s+running\s+with\s+linear\s+solver\s+(?<linearsolver>\S+)*','names','dotall');
            lines = strsplit(strtrim(text),newline);
            for i = 2:numel(lines)
                line = strtrim(lines{i});
                if strlength(strtrim(line))
                    m = regexp(line,'^(?<name>\S+)\s+=\s+(?<value>\S+)\s+(?<used>\S+)$','names');
                    if ~isempty(m)
                        if strcmp(m.used,'yes')
                            s.used_options.(m.name) = m.value;
                        else
                            s.unused_options.(m.name) = m.value;
                        end
                    end
                end
            end
        end
    end

    methods (Static, Hidden)
        

        function fVal = objective(x)
            fVal = x(1)*x(4)*sum(x(1:3)) + x(3);
        end
        % ----------------------------------------------------------------------
        function g = gradient(x)
            g = [ x(1)*x(4) + x(4)*sum(x(1:3))
                  x(1)*x(4)
                  x(1)*x(4) + 1
                  x(1)*sum(x(1:3)) ]; 
        end
        
        function c = constraints(x)
            c = [ prod(x); sum(x.^2) ];
        end
      
        function j = jacobian(x)
            j = sparse([ prod(x)./x; 2*x ]);
        end
      
        function cStr = jacobianstructure()
            cStr = sparse(ones(2,4));
        end

        % ----------------------------------------------------------------------
        function H = hessian (x, sigma, lambda)
        
            H = sigma*[ 2*x(4)             0      0   0;
                        x(4)               0      0   0;
                        x(4)               0      0   0;
                        2*x(1)+x(2)+x(3)  x(1)  x(1)  0 ];
        
            H = H + lambda(1)*[    0          0         0         0;
                                x(3)*x(4)     0         0         0;
                                x(2)*x(4) x(1)*x(4)     0         0;
                                x(2)*x(3) x(1)*x(3) x(1)*x(2)     0  ];
            H = H + lambda(2)*diag([2 2 2 2]);
            H = sparse(H);
        end

        function hStr = hessianstructure()
            hStr = sparse(tril(ones(4)));
        end

        function trydelete(filepath)
            if isfile(filepath)
                delete(filepath);
            end
        end
    end
    
end