% Graph Matching problem

function [isIso, mat, result] = GraphMatchingProblem(A, B, varargin)
    %{ This function takes in adjacency matrices 'A', 'B' and optional
    % arguments 'ts' and 'weights'.
    % 
    % This function solves the graph matching problem using MATLABs built
    % in QUBO solver. 
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
    addRequired(p,'B');
    addParameter(p, 'ts', tabuSearch());
    addParameter(p, 'weight', 10);
    parse(p, A, B, varargin{:})

    n = size(A, 1); % This assumes both adjacency matrices are the same size

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
    
    % Construct H
    H = (kron(A, eye(n)) - kron(eye(n), B))^2;

    % construct d
    d = transpose(ones(1, 2*n));
    
    % Create QUBO matrix Q
    Q = H + p.Results.weight*(transpose(C)*C - 2*diag(transpose(C)*d));
    qprob = qubo(Q);

    result = solve(qprob, Algorithm=p.Results.ts);
    
    % convert x vector to solution matrix
    mat = reshape(result.BestX, [n,n]);
    
    % Test to see if this is a isomorphic permutation matrix
    isIso = isequal(C*result.BestX, d) & all(H*result.BestX == 0); 

end
