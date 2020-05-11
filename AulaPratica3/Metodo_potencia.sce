function [lambda,x1,k,n_erro] = Metodo_potencia(A,x0,epsilon,M)

k=1;
x0=x0/max(x0);
x1 = A*x0;    // aproximação do autovetor dominante

while(k <= M)
    lambda = max(x1); // aproximação autovalor dominante
    x1 = x1/lambda;
    n_erro = norm(abs(x1) - abs(x0), %inf);
    if n_erro < epsilon
        break; 
    end
    x0 = x1;
    x1 = A*x0;
    k = k + 1;
end

endfunction
