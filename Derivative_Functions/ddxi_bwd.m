function [derivative_of_f]= ddx_bwd(f,dx)
derivative_of_f=zeros(size(f,1),size(f,2));
for j=1:size(f,2)
    for i=2:size(f,1)
        derivative_of_f(i,j)=(f(i,j)-f(i-1,j))/dx;
    end
end
k=1;
for j=1:size(f,2)
    derivative_of_f(k,j)=(f(k+1,j)-f(k,j))/dx;

end
end