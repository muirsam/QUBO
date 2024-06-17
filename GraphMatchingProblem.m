% Graph Matching problem

function [isIso, mat, result] = GraphMatchingProblem(A, B, ts)
    %{ This function takes in two Adjacency matrices A and B and an
    % (optional) tabuSearch object ts.
    % 
    % This function implements the graph matching problem as a QUBO
    % problem and solves it using MATLABs built in QUBO solver.
    % 
    % It then resurns a Boolean variable isIso which is true when the
    % solution satisfies the constraints, the matrix mat which corresponds
    % to the solution and result which contains additional information
    % about the problem solution.
    %
    % It can take it a tabuSearch object ts in order to control aspects of
    % the solver. More information about ts here,
    % https://uk.mathworks.com/help/matlab/ref/tabusearch.html
    %}
    
    % If ts not included, set it to default value
    if ~exist('ts','var')
      ts = tabuSearch();
    end

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

    % create Q matrix and solve QUBO
    Q = -eye(n*n) + 100*transpose(C)*C + 100)*transpose(H)*H;
    qprob = qubo(Q, 100*(-2)*transpose(d)*C, transpose(d)*d);
    
    result = solve(qprob, Algorithm=ts);
    
    % convert x vector to solution matrix
    mat = reshape(result.BestX, [n,n]);
    
    % Test to see if this is a isomorphic permutation matrix
    isIso = isequal(mat*A, B*mat);

end
