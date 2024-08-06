function [rho,u,v,T,p,e,Et] = cons2prim(U,R,cv)
    m=size(U,2);
    n=size(U,3);
    rho=zeros(m,n);
    u=zeros(m,n);
    v=zeros(m,n);
    T=zeros(m,n);
    p=zeros(m,n);
    e=zeros(m,n);
    Et=zeros(m,n);
    for k=1:4
        for j=1:n
            for i=1:m
                if k==1
                    rho(i,j)=U(k,i,j);
                elseif k==2
                    u(i,j)=U(k,i,j)/rho(i,j);
                elseif k==3
                    v(i,j)=U(k,i,j)/rho(i,j);
                else
                    T(i,j)=(U(k,i,j)-0.5*rho(i,j)*(u(i,j)^2+v(i,j)^2))/(rho(i,j)*cv);
                end
            end
        end
    end
    for j=1:n
        for i=1:m
            p(i,j)=rho(i,j)*R*T(i,j);
        end
    end

    for j=1:n
        for i=1:m
            e(i,j)=cv*T(i,j);
        end
    end

    for j=1:n
        for i=1:m
            Et(i,j)=U(4,i,j);
        end
    end
end