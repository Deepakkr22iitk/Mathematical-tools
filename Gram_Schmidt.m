clc;
clear;
% Taking the data inputs
n=input("Enter the number of the vectors : ");
m=input("Enter the dimension of the vector space where the vectors lie : ");
v=input("Enter the vectors as row vectors of a matrix : ");
% Displying the vectors
fprintf("Given vectros are displayed as columns in a (%d by %d) matrix given below.\n", m, n);
v=v';
disp(v);
% Gram-Schmidt process for orthogonalization
a=zeros(m,n);
a(:,1)=v(:,1);
for i=2:n
    s=0;
    for k=1:i-1
        s=s+(v(:,i)'*a(:,k))/(a(:,k)'*a(:,k))*a(:,k);
    end
    a(:,i)=v(:,i)-s;
end
% Displying the Orthogonalized vectors
fprintf("Orthogonalized (by Gram-Schmidt process) vectors are given as columns in a matrix given below :\n");
disp(a);
% orthonormalization
b=zeros(m,n);
for i=1:n
    b(:,i)=a(:,i)/norm(a(:,i),2);
end
% Displying the Orthonormalized vectors
fprintf("Orthonormalized (by Gram-Schmidt process) vectors are given as columns in a matrix given below :\n");
disp(b);