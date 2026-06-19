clc;
clear;
% Taking the data inputs
x=input("Enter the x-coordinates of the data points as row vector: ");
y=input("Enter the y-coordinates of the data points as row vector: ");
N1=length(x);
N2=length(y);

% Checking the consistency of the inputs
p=0;
while (p==0)
    if N1~=N2
        fprintf("Please check the number of data points and re-enter.\n");
        t=input("Enter the x-coordinates of the data points as row vector: ");
        y=input("Enter the y-coordinates of the data points as row vector: ");
        N1=length(x);
        N2=length(y);
    else 
        m=N1;
        p=1;
        fprintf("\n The data is given in a table as: \n\n ");
        disp(cell2mat(compose('%11.5f', [x ; y])));
    end
end
n=5; % number of parameters involved in the equation of the orbit.
A=[y'.^2 x'.*y' x' y' ones(m,1)]; %Finding the coefficient matrix
disp(A);
b=A'*(x').^2; % RHS of the normal system
fprintf("The augmented Normal system of equations associated to the given data is: \n");
disp([A'*A b]);

% Finding Cholesky decomposition of A^T*A
fprintf("Lower triangular matrix of the Cholesky decomposition of the normal system is given as: \n ");
a=A'*A;
L=zeros(n);
L(1,1)=sqrt(a(1,1));
for j=2:n
    L(j,1)=a(j,1)/L(1,1);
end
for i=2:n-1
    s=0;
    for k=1:i-1
        s=s+L(i,k)^2;
    end
    L(i,i)=sqrt((a(i,i)-s));
    for j=i+1:n
        s=0;
        for k=1:i-1
            s=s+L(j,k)*L(i,k);
        end
        L(j,i)=(a(j,i)-s)/L(i,i);
    end
end
s=0;
for k=1:n-1
    s=s+L(n,k)^2;
end
L(n,n)=sqrt((a(n,n)-s));

% L=chol(A'*A,'lower'); % Matlab inbuit command for Cholesky decompostion
disp(L);

% Solving the system A*A^Tu=b using Cholesky decomposition
z=zeros(1,n);
u=zeros(1,n);
z(1)=b(1)/L(1,1);
for i=2:n
    s=0;
    for j=1:i-1
        s=s+L(i,j)*z(j);
    end
    z(i)=(b(i)-s)/L(i,i);
end
u(n)=z(n)/L(n,n);
for i=n-1:-1:1
    s=0;
    for j=i+1:n
        s=s+L(j,i)*u(j);
    end
    u(i)=(z(i)-s)/L(i,i);
end
fprintf("Solution of the Normal equation is given as: ");
% disp(u);
disp(cell2mat(compose('%14.7f', u)));
% disp(x'.^2-A*u');
e=norm(x'.^2-A*u',2);
E=e^2;
fprintf("Minimized Error for this approximation is: %5.10f \n\n", E);

% Plotting the best fit polynomial with given data set
grid on
plot(x, y, "b*");
hold on;

f = @(x,y) u(1)*y.^2+u(2)*x.*y+ u(3)*x+u(4)*y+ u(5)*1-x.^2;
fimplicit(f,[(min(x)-1) (max(x)+1) (min(y)-1) (max(y)+1)])
%fimplicit(f)
title('Orbit of the planet', 'Fontsize', 11)
xlabel('x --->')
ylabel('y --->')
grid on