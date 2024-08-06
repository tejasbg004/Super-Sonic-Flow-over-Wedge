function [derivative_of_f]= ddxi_central(f,dx)
derivative_of_f=zeros(size(f,1),size(f,2));
for j=1:size(f,2)
    for i=2:size(f,1)-1
        derivative_of_f(i,j)=(f(i+1,j)-f(i-1,j))/(2*dx);
    end
end
k=1;
for j=1:size(f,2)
    derivative_of_f(k,j)=(4*f(k+1,j)-3*f(k,j)-f(k+2,j))/(2*dx);

end
l=size(f,1);
for j=1:size(f,2)
     derivative_of_f(l,j)=(-4*f(l-1,j)+3*f(l,j)+f(l-2,j))/(2*dx);
    
end
end
