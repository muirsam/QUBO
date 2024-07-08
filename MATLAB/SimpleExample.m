% simple example

a = [a_1, a_2, a_3, a_4];
weight = sum(abs(a));

Q = diag(-a - 4*weight) + weight*ones(4);
solve(qubo(Q))