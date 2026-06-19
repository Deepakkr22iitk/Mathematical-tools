clc;
clear;
% Taking the data inputs
t=input("Enter the x-coordinates of the data points as row vector: ");
y=input("Enter the y-coordinates of the data points as row vector: ");
m=length(t);
fprintf("\n The data is given in a table as: \n\n ");
disp(cell2mat(compose('%11.5f', [t ; y])));
d=input("Enter the degree of the polynomial to fit the above data: ");
n=d+1;
A=zeros(m,n); % Finding the coefficient matrix
for i=1:n
    A(:,i)=t'.^(i-1);
end
fprintf("Matrix A is given as:\n ");
disp(A);
b=A'*y'; % RHS of the normal system
fprintf("The augmented Normal system of equations associated to the given data is: \n");
disp([A'*A b]);
% Finding Cholesky decomposition of A^T*A
fprintf("Lower triangular matrix of the Cholesky decomposition of the normal system is given as: \n ");
L=chol(A'*A,'lower');
disp(L);
% Solving the system L*L^Tx=A^T*y using backward substitution
z=zeros(1,n);
x=zeros(1,n);
z(1)=b(1)/L(1,1);
for i=2:n
    s=0;
    for j=1:i-1
        s=s+L(i,j)*z(j);
    end
    z(i)=(b(i)-s)/L(i,i);
end
x(n)=z(n)/L(n,n);
for i=n-1:-1:1
    s=0;
    for j=i+1:n
        s=s+L(j,i)*x(j);
    end
    x(i)=(z(i)-s)/L(i,i);
end
fprintf("Solution of the Normal equation is given as: ");
% disp(x);
disp(cell2mat(compose('%14.7f', x)));
% disp(y'-A*x');
e=norm(y'-A*x',2);
E=e^2;
fprintf("Minimized Error for this approximation is: %5.10f \n\n", E);

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
xlim([(min(t)-1) (max(t)+1)])
ylim([(min(y)-1) (max(y)+1)])
grid on