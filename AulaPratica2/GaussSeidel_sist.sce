function x = sist(L,b)
n = size(L,1);
x = zeros(n,1);
x(1) = b(1)/L(1,1);
for i = 2:n
    x(i) = (b(i) - sum(L(i,1:i-1)*x(1:i-1)))/L(i,i);
end
endfunction

function [x, d, k, r] =GaussSeidel_sist(A, b, x0, E, M, n) 

L = tril(A);

U = triu(A, 1);

x = sist(L, -U*x0 + b);

k = 0;

while((norm(x-x0, n)>=E) && (k<M))
    x0 = x;
    x = sist(L, -U*x0 + b);
    k = k + 1;
end

d = norm(x-x0, n);

r = norm(b-A*x, n);

endfunction
