clc;
clear;

% Taking the data inputs
t=input("Enter the x-coordinates of the data points as row vector: ");
y=input("Enter the y-coordinates of the data points as row vector: ");
N1=length(t);
N2=length(y);
p=0;
while (p==0)
    if N1~=N2
        fprintf("Please check the number of data points and re-enter.\n");
        t=input("Enter the x-coordinates of the data points as row vector: ");
        y=input("Enter the y-coordinates of the data points as row vector: ");
        N1=length(t);
        N2=length(y);
    else 
        m=N1;
        p=1;
        fprintf("\n The data is given in a table as: \n\n ");
        disp(cell2mat(compose('%11.5f', [t ; y])));
    end
end
d=input("Enter the degree of the polynomial to fit the above data: ");
n=d+1;

% Finding the coefficient matrix
A=zeros(m,n); 
for i=1:n
    A(:,i)=t'.^(i-1);
end
fprintf("Matrix A is given as:\n ");
disp(A);

%
Q=eye(m);
b=y';
for k=1: min(n,m-1)
    for j=k+1:m
        T=eye(m);
        % a=[A(k,k); A(j,k)];
        a=A([k,j],k);
        if a(2)~=0 
            c=a(1)/norm(a);
            s=-a(2)/norm(a);
            T(k,k)=c;
            T(j,j)=c;
            T(k,j)=-s;
            T(j,k)=s;
            fprintf('The Givens matrix T_(%d,%d) is:\n ', k , j);
            disp(T);
            fprintf('At the current iteration the reduced matrix is:\n ');
            A=T*A;
            b=T*b;
            disp([A, b]);
            Q=T*Q;
        end
    end
end
fprintf("The matrix Q in QR factorization is given as:\n ");
disp(Q');
fprintf("The matrix R in QR factorization is given as:\n ");
disp(A);

% Solution using Householder
r=0;
for i=1:min(m,n)
    if A(i,i)~= 0
        r=r+1;
    end
end
R=A(1:r, :);
disp(R);
c=b(1:r);

% Coefficients of the best fitted polynomial
x=R\c;
fprintf("Coefficients of the best fitted polynomial are :\n ");
disp(x);

% Plotting the best fit polynomial with given data set
grid on
plot(t, y, "k*");
hold on;
h=min(t):.01: max(t);
w=0;
for i=1:n
    w=w+x(i)*h.^(i-1);  
end
plot(h, w, "r-");
title('Best fit plynomial of degree d for the given data set', 'Fontsize', 11)
xlabel('t --->')
ylabel('y --->')
%xlim([(min(t)-1) (max(t)+1)])
%ylim([(min(y)-1) (max(y)+1)])
grid on