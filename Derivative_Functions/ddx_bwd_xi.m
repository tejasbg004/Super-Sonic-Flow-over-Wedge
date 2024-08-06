function [dfdx]= ddx_bwd_xi(f,dxi,deta,yy1,xx)
    [nx,ny]     = size(f);
    dxdxi=ddxi_central(xx,dxi);
    dxdeta=ddeta_central(xx,deta);
    dydxi=ddxi_central(yy1,dxi);
    dydeta=ddeta_central(yy1,deta);
    J=dxdxi.*dydeta-dxdeta.*dydxi;

    xix=dydeta./J;
    etax=-dydxi./J;
    
    % determine field size
    

    % allocate return field
    dfdx        = zeros(nx,ny);
    
    % backward difference
    for i=2:nx
        for j=2:ny
            dfdx(i,j) = ((f(i,j)-f(i-1,j))/dxi)*xix(i,j)+((f(i,j)-f(i,j-1))/deta)*etax(i,j);
        end
    end
    
    % forward difference for first collumn
    j = 1;
    for i=2:nx
        dfdx(i,j) = ((f(i,j)-f(i-1,j))/dxi)*xix(i,j)+((f(i,j+1)-f(i,j))/deta)*etax(i,j);
    end

 % forward difference for first row
    i=1;
    for j=2:ny
        dfdx(i,j) = ((f(i+1,j)-f(i,j))/dxi)*xix(i,j)+((f(i,j)-f(i,j-1))/deta)*etax(i,j);
    end

    dfdx(1,1) = ((f(1+1,1)-f(1,1))/dxi)*xix(1,1)+((f(1,1+1)-f(1,1))/deta)*etax(1,1);
end