function [dfdx]= ddx_fwd_xi(f,dxi,deta,yy1,xx)
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
    
    % forward difference
    for i=1:nx-1
        for j=1:ny-1
            dfdx(i,j) = ((f(i+1,j)-f(i,j))/dxi)*xix(i,j)+((f(i,j+1)-f(i,j))/deta)*etax(i,j);
        end
    end
    
    % backward difference for last collumn
    j = ny;
    for i=1:nx-1
        dfdx(i,j) = ((f(i+1,j)-f(i,j))/dxi)*xix(i,j)+((f(i,j)-f(i,j-1))/deta)*etax(i,j);
    end

 % backward difference for last row
    i=nx;
    for j=1:ny-1
        dfdx(i,j) = ((f(i,j)-f(i-1,j))/dxi)*xix(i,j)+((f(i,j+1)-f(i,j))/deta)*etax(i,j);
    end

    dfdx(nx,ny) = ((f(nx,ny)-f(nx-1,ny))/dxi)*xix(nx,ny)+((f(nx,ny)-f(nx,ny-1))/deta)*etax(nx,ny);

end