% Maximal Independent Set problem

function [isValid, mat, result] = MaximalIndependentSetProblem(A, varargin)
    %{ This function takes in adjacency matrix 'A', and optional arguments 
    % 'ts' and 'weights'.
    % 
    % It then uses MATLABs built-in QUBO solver to solve the maximimal 
    % independent set problem. 
    %
    % It returns Boolean 'isIso' which tells you if the solution satisfies
    % the constraints, the matrix 'mat' which corresponds to the solution and
    % the quboResult 'result' which contains additional information about
    % the solution.
    %
    % The optional argument 'ts' is a tabuSearch object which controls 
    % properties of the algorithm used to find a solution.
    % The optional argument 'weight' controls the penalty weight.
    %}

    % optional argument stuff
    p = inputParser;
    addRequired(p,'A');
    addParameter(p, 'ts', tabuSearch());
    addParameter(p, 'weight', 10);
    parse(p, A, varargin{:})

    n = size(A, 1); % This assumes both adjacency matrices are the same size
    
    Q = zeros(n);

    for i=1:n
        for j=1:n
            if A(i,j) == 1
                Q(i,i) = Q(i,i) - p.Results.weight;
                Q(i,j) = Q(i,j) + p.Results.weight;
            end
        end
    end

    qprob = qubo(Q);
    result = solve(qprob, Algorithm=p.Results.ts);
    
end
