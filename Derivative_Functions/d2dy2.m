function [d2fdy2]=d2dy2(f,dy)
m=size(f,1);
n=size(f,2);
d2fdy2=zeros(m,n);
for i=1:m
    for j=2:n-1
        d2fdy2(i,j)=(f(i,j+1)-2*f(i,j)+f(i,j-1))/(dy^2);
    end
end

k=1;
for i=1:m
    d2fdy2(i,k)=(f(i,k+1)-2*f(i,k)+f(i,end))/dy^2;
end
l=n;
for i=1:m
    d2fdy2(i,l)=(f(i,1)-2*f(i,l)+f(i,l-1))/dy^2;
end
end