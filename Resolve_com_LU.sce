function X =Resolve_com_LU(C, B, P)
// C e B são como especificados no enunciado e P é a matriz de permutações

n = size(C,1);
m = size(B,2);
//permuta os elementos dos vetores b_i
B = P*B;

//
Y=zeros(n,m);
// resolve LY=B, ou seja, resolve L y_i = b_i com 1<=i<=m
for j=1:m
    Y(1, j)=B(1,j);
    for i=2:n
        Y(i,j)=(B(i,j)-C(i,1:i)*Y(1:i,j));
    end
end

X=zeros(n,m);

// resolve UX=Y, ou seja, resolve U x_i = y_i com 1<=i<=m
for j=1:m
    X(n, j)=Y(n,j)/C(n,n);
    for i=n-1:-1:1
        X(i,j)=(Y(i,j)-C(i,i:n)*X(i:n,j))/C(i,i);
    end
end



endfunction
