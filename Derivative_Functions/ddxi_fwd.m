function [derivative_of_f]= ddx_fwd(f,dx)
derivative_of_f=zeros(size(f,1),size(f,2));
for j=1:size(f,2)
    for i=1:size(f,1)-1
        derivative_of_f(i,j)=(f(i+1,j)-f(i,j))/dx;
    end
end
k=size(f,1);
for j=1:size(f,2)
    derivative_of_f(k,j)=(f(k,j)-f(k-1,j))/dx;

end