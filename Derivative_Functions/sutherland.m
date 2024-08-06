function [mu] = sutherland(T,u0)
    m=size(T,1);
    n=size(T,2);
    mu=zeros(m,n);
    S1=110.4;
    %u0=1.735*10^-5;
    T0=288.15;
    for i=1:m
        for j=1:n
            mu(i,j)=u0*(((T(i,j)/T0)^(3/2))*((T0+S1)/(T(i,j)+S1)));
        end
    end

            