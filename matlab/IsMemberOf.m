classdef IsMemberOf < matlab.unittest.constraints.Constraint
    properties (SetAccess=immutable)
        set (:,1)
    end
    
    methods
        function constraint = IsMemberOf(set)
            constraint.set = set;
        end
        
        function tf = satisfiedBy(constraint,actual)
            tf = constraint.instanceIsMemberOfSet(actual);
        end

        function diagnostic = getDiagnosticFor(constraint,actual)
            import matlab.automation.diagnostics.StringDiagnostic
            if constraint.instanceIsMemberOfSet(actual)
                diagnostic = StringDiagnostic("IsMemberOf passed.");
            else
                diagnostic = StringDiagnostic( ...
                    "IsMemberOf failed." + newline + "Actual value: " ...
                    + string(actual) + newline ...
                    + "Expected Options: [" ...
                    + strjoin(constraint.set,', ') ...
                    + "]");
            end
        end
    end

    methods (Access=private)
        function tf = instanceIsMemberOfSet(constraint,actual)
            tf = ismember(actual,constraint.set);
        end
    end
end

