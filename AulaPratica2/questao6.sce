numelem = [10 100 1000 2000];
t = zeros(2,4);

for i = 1:4
    A = 2*(rand((numelem(1, i),numelem(1, i)), "uniform")-0.5) + diag(rand((numelem(1 ,i),1), "uniform") + numelem(1 ,i) - 1);
    b = 2*(rand((numelem(1, i),1), "uniform") - 0.5);
    tic();
    [x,y,z,w] = GaussSeidel_inv(A, b, zeros(numelem(1, i),1), 10^(-4), 1000, 2);
    t(1,i) = toc();
    tic();
    [x,y,z,w] = GaussSeidel_sist(A, b, zeros(numelem(1, i),1), 10^(-4), 1000, 2);
    t(2,i) = toc();
end

disp(t);
