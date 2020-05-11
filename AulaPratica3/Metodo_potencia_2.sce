function [lambda,x1,k,n_erro] = Metodo_potencia_2(A,x0,epsilon,M)

k=1;
x0=x0/norm(x0,2);
x1=A*x0;

while (k<=M)
    lambda = (x0')*x1;  // Quociente de Rayleigh x0 é unitário
    x1 = x1/ norm(x1, 2);
    n_erro = norm(abs(x1) - abs(x0), 2);
    if (n_erro < epsilon) 
        break;
    end
    x0 = x1
    x1 = A*x0
    k = k+1;
end

endfunction
