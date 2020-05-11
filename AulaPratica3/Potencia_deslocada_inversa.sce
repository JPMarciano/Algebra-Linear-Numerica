function [lambda,x1,k,n_erro] = Potencia_deslocada_inversa(A,x0,epsilon,alfa,M) 

k=1;
x0=x0/norm(x0, 2);
[n, m] = size(A);

while(k<=M)
    [x1, C] = Gaussian_Elimination_4(A - alfa*eye((n,n)), x0);
    x1 = x1/norm(x1, 2);
    lambda = (x1')*A*x1; // Quociente de Rayleigh; x1 é unitário
    n_erro = norm(abs(x1) - abs(x0), 2);
    if (n_erro < epsilon)
        break;
    end 
    x0 = x1;
    k = k + 1;
end

endfunction
