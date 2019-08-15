% hello
u=zeros(4,4);
u(2,2)=1;
u(2,3)=1;
u(3,2)=1;
u(3,3)=1;

A=zeros(4,4);
B=zeros(4,4);

phi=1;

%%Jacobi method  %%
for i=1:100

for i=2:3
    for j=2:3
    A(i,j)=1/4*(u(i-1,j)+u(i,j-1)+u(i+1,j)+u(i,j+1)-phi);
    end
end
delta=norm(A-u);
u=A;
if delta<0.001
    break
end
u
end

%Gauss-Seidel method%
u1=zeros(4,4);
u1(2,2)=1;
u1(2,3)=1;
u1(3,2)=1;
u1(3,3)=1;
 
for i=1:100
for i=2:3
    for j=2:3
    u1(i,j)=1/4*(u1(i-1,j)+u1(i,j-1)+u1(i+1,j)+u1(i,j+1)-phi);
    end
end
delta1=norm(B-u1);
B=u1;
if delta1<0.001
    break
end
u1
end 

e=eig(u)
e1=eig(u1)

