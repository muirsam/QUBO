function [sol, result, varargout] = sudoku(clues, varargin)
    %{ This function takes in a n x 3 array called 'clues'. Each 1 x 3 
    % array in 'clues' repressents a Sudoku position that has already been
    % filled. Each clue has form [row, column, value].
    %
    % This function REQUIRES the file 'sudokuQ.mat' to create the QUBO matrix
    %
    % It will return the boolean the matrix 'sol' which corresponds to the
    % solution and the quboResult 'result' which contains additional
    % information about the solution process.
    %
    % The optional argument 'ts' is a tabuSearch object which controls 
    % properties of the algorithm used to find a solution.
    % It can also check the validity of the result by setting 
    % the paramater 'test' = true
    %}

    % optional argument stuff
    p = inputParser;
    addRequired(p,'clues');
    addParameter(p, 'ts', tabuSearch());
    addParameter(p, 'test', false)
    parse(p, clues, varargin{:})

    n = 9;
    algo = p.Results.ts;

    % converts to index
    conI = @(i, j, k) (i-1)*n^2 + (j-1)*n + k;
    
    % loads constraint matrix
    s = load("sudokuQ.mat"); % This takes ~50ms to run
    Q = s.Q;

    % entries given
    for i=transpose(clues)
        B = transpose(i);
        for j=1:n
            Q(conI(B(1),B(2),j),conI(B(1),B(2),j)) = Q(conI(B(1),B(2),j),conI(B(1),B(2),j)) + 1;
        end
        Q(conI(B(1),B(2),B(3)),conI(B(1),B(2),B(3))) = Q(conI(B(1),B(2),B(3)),conI(B(1),B(2),B(3))) - 2;
    end
    
    % solve QUBO
    qprob = qubo(Q, [], size(clues, 1)+324);
    result = solve(qprob, Algorithm=algo);
    
    % convert solution to vector 
    vec = zeros(n^2,1);
    for l=1:n^2
        for m=1:n
            if result.BestX((l-1)*n + m) == 1
                vec(l) = m;
            end
        end
    end
    
    % obtain matrix form 
    solution = reshape(vec, n, []);
    sol = transpose(solution);
    
    % check solution
    if p.Results.test == true
        isValid = true;

        % Check the solution has each clue in the same place
        for i=transpose(clues)
            B = transpose(i);
            isValid = isValid & sol(B(1),B(2)) == B(3);
        end
    
        % Check each row contains all values
        for i=1:n
            for k=1:n
                isValid = isValid & ismember(k, sol(i,:));
            end
        end
        
        % Check each column contains all values
        for j=1:n
            for k=1:n
                isValid = isValid & ismember(k, sol(:,j));
            end
        end
    
        % Check each box contains all values
        for p=1:sqrt(n)
            for q=1:sqrt(n)
                box = sol((p-1)*3 + 1:(p-1)*3+3, (q-1)*3 + 1:(q-1)*3+3);
                x = reshape(box, 1, []);
                for k = 1:n
                    isValid = isValid & ismember(k, x);
                end
            end
        end
        varargout{1} = isValid;
    end
end