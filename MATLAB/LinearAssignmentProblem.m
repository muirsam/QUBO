function [isPerm, P, result] = LinearAssignmentProblem(costs, varargin)
    %{ This function takes in a matrix 'costs' of costs and optional arguments
    % 'ts' and 'weight'.
    %
    % It solves the linear assignment problem for the matrix 'costs' using
    % MATLABs built in QUBO solver.
    %    
    % It returns the Boolean 'isPerm' which tells you if the solution is
    % valid, the solution matrix 'mat' and the result of the QUBO solver
    % 'result'
    %
    % The optional argument 'ts' is a tabuSearch object which controls 
    % properties of the algorithm used to find a solution.
    % The optional argument 'weights' is a double which 
    % control the penalty weight.
    %}

    % optional argument stuff
    ip = inputParser;
    addRequired(ip,'A');
    addParameter(ip, 'ts', tabuSearch());
    addParameter(ip, 'weight', 10);
    parse(ip, costs, varargin{:})

    % -------------------------------------

    n = size(costs, 1);

    % Construct C
    C = zeros(2*n, n^2);
    for i=1:n
        for j=1:n
            C(i, n*(i-1) + j) = 1;
        end
    end
    for i=1:n
        for j=i:n:n^2
            C(n + i, j) = 1;
        end
    end

    % construct d
    d = transpose(ones(1, 2*n));
    
    % construct QUBO matrix Q
    Q = diag(reshape(costs, 1, [])) + ip.Results.weight*transpose(C)*C - ip.Results.weight*2*diag(transpose(C)*d);

    qprob = qubo(Q);
    result = solve(qprob, Algorithm=ip.Results.ts);
    
    P = transpose(reshape(result.BestX, n, []));  % transform optimal solution to matrix
    
    % check if solution is valid
    isPerm = true;

    for i=1:n
        isPerm = isPerm & sum(P(:,i)) == true & sum(P(i,:)) == true;
    end

end