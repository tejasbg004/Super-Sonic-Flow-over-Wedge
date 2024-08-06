function [derivative_of_f]= ddy_bwd(f,dy)
derivative_of_f=zeros(size(f,1),size(f,2));
for i=1:size(f,1)
    for j=2:size(f,2)
        derivative_of_f(i,j)=(f(i,j)-f(i,j-1))/dy;
    end
end
k=1;
for i=1:size(f,1)
    derivative_of_f(i,k)=(f(i,k+1)-f(i,k))/dy;

end