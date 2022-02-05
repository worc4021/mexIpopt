% Test the "ipopt" Matlab interface on the Hock & Schittkowski test problem
% #71. See: Willi Hock and Klaus Schittkowski. (1981) Test Examples for
% Nonlinear Programming Codes. Lecture Notes in Economics and Mathematical
% Systems Vol. 187, Springer-Verlag.
%
% Copyright (C) 2008 Peter Carbonetto. All Rights Reserved.
% This code is published under the Eclipse Public License.
%
% Author: Peter Carbonetto
%         Dept. of Computer Science
%         University of British Columbia
%         September 18, 2008
clear;

    x0         = [3 3 3 3];  % The starting point.
    options.lb = [1 1 1 1];  % Lower bound on the variables.
    options.ub = [5 5 5 5];  % Upper bound on the variables.
    options.cl = [25  40];   % Lower bounds on the constraint functions.
    options.cu = [inf 40];   % Upper bounds on the constraint functions.
    
    % Initialize the dual point.
%     options.zl     = [1 1 1 1];
%     options.zu     = [1 1 1 1];
%     options.lambda = [1 1];

    
    % Set the IPOPT options.
    options.ipopt.tol                               = 1e-7;
    options.ipopt.max_iter                          = 25;
    options.ipopt.mu_strategy                       = "adaptive";
    options.ipopt.linear_solver                     = "ma27";
    options.ipopt.output_file                       = "ipopt.out";
    options.ipopt.file_print_level                  = 5;
    options.ipopt.print_level                       = 5;
    options.ipopt.hessian_approximation             = "exact";%"limited-memory";
    options.ipopt.limited_memory_update_type        = "damped-bfgs";
    options.ipopt.limited_memory_damping_threshold  = 0.2;
%     options.ipopt.print_options_documentation       = "yes";
    options.ipopt.dependency_detector               = "mumps";
    options.ipopt.print_user_options                = "yes";
%     options.debug = 1;
    
    % The callback functions.
    funcs.objective         = @(x) x(1)*x(4)*sum(x(1:3)) + x(3);
    funcs.constraints       = @constraints;
    funcs.gradient          = @gradient;
    funcs.jacobian          = @jacobian;
    funcs.jacobianstructure = @() sparse(ones(2,4));
    funcs.hessian           = @hessian;
    funcs.hessianstructure  = @() sparse(tril(ones(4)));
    funcs.intermediate      = @intermediate;
    % Run IPOPT.
    [x, info] = ipopt(x0, funcs, options);
  
  % ----------------------------------------------------------------------
  function g = gradient (x)
    g = [ x(1)*x(4) + x(4)*sum(x(1:3))
          x(1)*x(4)
          x(1)*x(4) + 1
          x(1)*sum(x(1:3)) ]; 
  end
    
  function c = constraints(x)
  c = [ prod(x); sum(x.^2) ];
  end

  function b = intermediate(intermediateCallbackStructure)
%     fprintf('iteration %d\n',x.iter);
    b = true;
  end
  
  function j = jacobian(x)
  j = sparse([ prod(x)./x; 2*x ]);
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