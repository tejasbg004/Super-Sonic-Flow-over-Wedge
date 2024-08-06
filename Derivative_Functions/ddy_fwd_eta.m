function dfdy = ddy_fwd_eta(f,dxi,deta,yy1,xx)
    [nx,ny]     = size(f);
    dxdxi=ddxi_central(xx,dxi);
    dxdeta=ddeta_central(xx,deta);
    dydxi=ddxi_central(yy1,dxi);
    dydeta=ddeta_central(yy1,deta);
    J=dxdxi.*dydeta-dxdeta.*dydxi;

    xiy=-dxdeta./J;
    etay=dxdxi./J;
    
    % determine field size
    

    % allocate return field
    dfdy        = zeros(nx,ny);
    
    % forward difference
    for i=1:nx-1
        for j=1:ny-1
            dfdy(i,j) = ((f(i+1,j)-f(i,j))/dxi)*xiy(i,j)+((f(i,j+1)-f(i,j))/deta)*etay(i,j);
        end
    end
    
    % backward difference for last collumn
    j = ny;
    for i=1:nx-1
        dfdy(i,j) = ((f(i+1,j)-f(i,j))/dxi)*xiy(i,j)+((f(i,j)-f(i,j-1))/deta)*etay(i,j);
    end

 % backward difference for last row
    i=nx;
    for j=1:ny-1
        dfdy(i,j) = ((f(i,j)-f(i-1,j))/dxi)*xiy(i,j)+((f(i,j+1)-f(i,j))/deta)*etay(i,j);
    end

    dfdy(nx,ny) = ((f(nx,ny)-f(nx-1,ny))/dxi)*xiy(nx,ny)+((f(nx,ny)-f(nx,ny-1))/deta)*etay(nx,ny);

end
