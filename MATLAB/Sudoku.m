function [isValid, sol, result] = Sudoku(clues, varargin)
    %{ This function takes in a n x 3 array called 'clues'. Each 1 x 3 
    % array in 'clues' repressents a Sudoku position that has already been
    % filled. Each clue has form [row, column, value].
    %
    % It solves each 9 by 9 Sudoku using MATLABs built in QUBO solver.
    %
    % It will return the boolean 'isValid' which tells you if the solution
    % satisfies the constraints, the matrix 'sol' which corresponds to the
    % solution and the quboResult 'result' which contains additional
    % information about the solution process.
    %
    % The optional argument 'ts' is a tabuSearch object which controls 
    % properties of the algorithm used to find a solution.
    %
    %}

    % optional argument stuff
    p = inputParser;
    addRequired(p,'clues');
    addParameter(p, 'ts', tabuSearch());
    parse(p, clues, varargin{:})

    n = 9;
    Q = zeros(n^3, n^3);
    algo = p.Results.ts;

    % converts to index
    conI = @(i, j, k) (i-1)*n^2 + (j-1)*n + k;
    
    % entries given
    for i=transpose(clues)
        B = transpose(i);
        Q(conI(B(1),B(2),B(3)),conI(B(1),B(2),B(3))) = Q(conI(B(1),B(2),B(3)),conI(B(1),B(2),B(3))) - 1;
    end
    
    % row penalties
    for i=1:n
        for k=1:n
            x = conI(i,1:n,k);
            for l=1:n
                Q(x(l),x(l)) = Q(x(l),x(l)) - 1;
                for m=l+1:n
                    Q(x(l),x(m)) = Q(x(l),x(m)) + 2;
                end
            end
        end
    end

    % column penalties
    for j=1:n
        for k=1:n
            x = conI(1:n,j,k);
            for l=1:n
                Q(x(l),x(l)) = Q(x(l),x(l)) - 1;
                for m=l+1:n
                    Q(x(l),x(m)) = Q(x(l),x(m)) + 2;
                end
            end
        end
    end

    % each cell must contain a number 
    for i=1:n
        for j=1:n
            x = conI(i,j,1:n);
            for l=1:n
                Q(x(l),x(l)) = Q(x(l),x(l)) - 1;
                for m=l+1:n
                    Q(x(l),x(m)) = Q(x(l),x(m)) + 2;
                end
            end
        end
    end
    
    % box penalties
    for p=1:sqrt(n)
        for q=1:sqrt(n)
            for k = 1:n
                x = [];
                for i=sqrt(n)*p-2:sqrt(n)*p
                    for j=sqrt(n)*q-2:sqrt(n)*q
                        x = [x conI(i,j,k)];
                    end
                end
                for l=1:n
                    Q(x(l),x(l)) = Q(x(l),x(l)) - 1;
                    for m=l+1:n
                        Q(x(l),x(m)) = Q(x(l),x(m)) + 2;
                    end
                end
            end
        end
    end 
    
    % solve QUBO
    qprob = qubo(Q);
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

    % Test if solution is valid
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
end