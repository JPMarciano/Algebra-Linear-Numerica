function [x, d, k, r] = Jacobi(A, b, x0, E, M, n) 

invD = diag(1./diag(A));
invD

R = -A+diag(diag(A));
R

db = invD*b;

Mj = invD*R;

x = Mj*x0 + db;

k = 0;

while((norm(x-x0, n)>=E) && (k<M))
    x0 = x;
    x = Mj*x0 + db;
    k = k + 1;
end

d = norm(x-x0, n);

r = norm(b-A*x, n);

endfunction
