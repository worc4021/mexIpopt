classdef Himmelblau < ProblemData
    methods
        function obj = Himmelblau()
            obj.x0 = [3; 2];
            obj.z0 = [0, 0;
                      0, 0];
            obj.lambda0 = 5;
            obj.cBnd = [-inf,-2];
            obj.xBnd = [-5,5;
                        -5,5];
        end

      function fVal = objective(~, var)
        x = var(1,:);
        y = var(2,:);
        fVal = (x.^2+y-11).^2+(x+y.^2-7).^2;
      end
      % ----------------------------------------------------------------------
      function fJac = gradient (~,var,~)
        x = var(1);
        y = var(2);

        fJac = [2*x + 4*x*(x^2 + y - 11) + 2*y^2 - 14;
                2*y + 4*y*(y^2 + x - 7) + 2*x^2 - 22];
      end
        
      function c = constraints(~,var)
        x = var(1,:);
        y = var(2,:);
        
        c = sin(x)-y;
      end
      
      function j = jacobian(~,var)
        x = var(1);
        j = sparse([cos(x),-1]);
      end
      
      function cStr = jacobianstructure(obj)
        cStr = obj.jacobian(rand(2,1));
      end

      % ----------------------------------------------------------------------
      function hVal = hessian (~, var, sigma, lambda)
        
        x = var(1);
        y = var(2);
            
        H = [12*x^2 + 4*y - 42,4*x + 4*y;
            4*x + 4*y,12*y^2 + 4*x - 26];
        C = [sin(x), 0;
             0, 0];
        hVal = sparse(tril(sigma*H + lambda(1)*C));
      end

      function hStr = hessianstructure(obj)
        fakelambda = rand(1);
        hStr = obj.hessian(rand(2,1),1, fakelambda);
      end
    end

end

