classdef hs71 < BaseProblem
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
    methods 
        function obj = hs71()
            obj.x0 = [3; 3; 3; 3];
            obj.xBnd = [1, 5;
                        1, 5;
                        1, 5;
                        1, 5];
            obj.cBnd = [25, inf;
                        40, 40];
            obj.z0 = [1, 1;
                      1, 1;
                      1, 1;
                      1, 1];
            obj.lambda0 = [1;1];
        end
    end

  methods
      function fVal = objective(~, x, ~)
        fVal = x(1)*x(4)*sum(x(1:3)) + x(3);
      end
      % ----------------------------------------------------------------------
      function g = gradient (~,x,~)
        g = [ x(1)*x(4) + x(4)*sum(x(1:3))
              x(1)*x(4)
              x(1)*x(4) + 1
              x(1)*sum(x(1:3)) ]; 
      end
        
      function c = constraints(~,x,~)
        c = [ prod(x); sum(x.^2) ];
      end
      
      function j = jacobian(~,x,~)
        j = sparse([ prod(x)./x; 2*x ]);
      end
      
      function cStr = jacobianstructure(~)
        cStr = sparse(ones(2,4));
      end

      % ----------------------------------------------------------------------
      function H = hessian (~, x, sigma, lambda,~,~)
        
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

      function hStr = hessianstructure(~)
        hStr = sparse(tril(ones(4)));
      end
    end

end