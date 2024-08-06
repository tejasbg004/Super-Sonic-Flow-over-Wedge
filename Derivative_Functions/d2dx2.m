function [d2fdx2]=d2dx2(f,dx)
m=size(f,1);
n=size(f,2);
d2fdx2=zeros(m,n);
for j=1:n
    for i=2:m-1
        d2fdx2(i,j)=(f(i+1,j)-2*f(i,j)+f(i-1,j))/(dx^2);
    end
end

k=1;
for j=1:n
    d2fdx2(i,k)=(f(k+1,j)-2*f(k,j)+f(end,j))/dx^2;
end
l=n;
for j=1:n
    d2fdx2(i,l)=(f(1,j)-2*f(l,j)+f(l-1,j))/dx^2;
end
end