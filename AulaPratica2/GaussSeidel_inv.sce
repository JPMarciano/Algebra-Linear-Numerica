function [x, d, k, r] =GaussSeidel_inv(A, b, x0, E, M, n) 

Linv = inv(tril(A));

U = triu(A, 1);

Mg = -Linv*U;

cg = Linv*b;

x = Mg*x0 + cg;

k = 0;

while((norm(x-x0, n)>=E) && (k<M))
    x0 = x;
    x = Mg*x0 + cg;
    k = k + 1;
end

d = norm(x-x0, n);

r = norm(b-A*x, n);

endfunction
