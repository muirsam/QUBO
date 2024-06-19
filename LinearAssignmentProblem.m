costs = [6 1 2; 2 8 1; 3 4 3];

A = diag(reshape(costs, 1, []));

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

% constraint weight
W = 10;

Q = A + W*transpose(C)*C - W*2*diag(transpose(C)*d);

qprob = qubo(Q);

result = solve(qprob);

P = transpose(reshape(result.BestX, n, []));