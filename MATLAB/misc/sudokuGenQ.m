% Generate constraint matrix

n = 9;
Q = zeros(n^3, n^3);

% converts to index
conI = @(i, j, k) (i-1)*n^2 + (j-1)*n + k;

% row penalties
for i=1:n
    for k=1:n
        x = conI(i,1:n,k);
        for l=1:n
            Q(x(l),x(l)) = Q(x(l),x(l)) - 1;
            for m=l+1:n
                Q(x(l),x(m)) = Q(x(l),x(m)) + 1;
                Q(x(m),x(l)) = Q(x(m),x(l)) + 1;
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
                Q(x(l),x(m)) = Q(x(l),x(m)) + 1;
                Q(x(m),x(l)) = Q(x(m),x(l)) + 1;
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
                Q(x(l),x(m)) = Q(x(l),x(m)) + 1;
                Q(x(m),x(l)) = Q(x(m),x(l)) + 1;
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
                Q(x(l),x(m)) = Q(x(l),x(m)) + 1;
                Q(x(m),x(l)) = Q(x(m),x(l)) + 1;
                end
            end
        end
    end
end 