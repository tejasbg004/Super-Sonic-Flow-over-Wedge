function dfdy = ddy_central_eta(f,dxi,deta,yy1,xx)
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
    
    % central difference
    for i=2:nx-1
        for j=2:ny-1
            dfdy(i,j) = ((f(i+1,j)-f(i-1,j))/2/dxi)*xiy(i,j)+((f(i,j+1)-f(i,j-1))/2/deta)*etay(i,j);
        end
    end
    
    % forward difference for first point
    j = 1;
    for i=2:nx-1
        dfdy(i,j) = ((f(i+1,j)-f(i-1,j))/2/dxi)*xiy(i,j)+((-3*f(i,j)+4*f(i,j+1)-f(i,j+2))/2/deta)*etay(i,j);
    end
    
    % backward difference for last point
    j = ny;
    for i=2:nx-1
        dfdy(i,j) = ((f(i+1,j)-f(i-1,j))/2/dxi)*xiy(i,j)+((3*f(i,j)-4*f(i,j-1)+f(i,j-2))/2/deta)*etay(i,j);
    end

    %forward differencing first row
    i = 1;
    for j=2:ny-1
        dfdy(i,j) = ((-3*f(i,j)+4*f(i+1,j)-f(i+2,j))/2/dxi)*xiy(i,j)+((f(i,j+1)-f(i,j-1))/2/deta)*etay(i,j);
    end
    
    % backward difference for last row
    i = nx;
    for j=2:ny-1
        dfdy(i,j) = ((3*f(i,j)-4*f(i-1,j)+f(i-2,j))/2/dxi)*xiy(i,j)+((f(i,j+1)-f(i,j-1))/2/deta)*etay(i,j);
    end





dfdy(1,1) = ((-3*f(1,1)+4*f(1+1,1)-f(1+2,1))/2/dxi)*xiy(1,1)+((-3*f(1,1)+4*f(1,1+1)-f(1,1+2))/2/deta)*etay(1,1);
dfdy(1,ny) = ((-3*f(1,ny)+4*f(1+1,ny)-f(1+2,ny))/2/dxi)*xiy(1,ny)+((3*f(1,ny)-4*f(1,ny-1)+f(1,ny-2))/2/deta)*etay(1,ny);
dfdy(nx,1) = ((3*f(nx,1)-4*f(nx-1,1)+f(nx-2,1))/2/dxi)*xiy(nx,1)+((-3*f(nx,1)+4*f(nx,1+1)-f(nx,1+2))/2/deta)*etay(nx,1);
dfdy(nx,ny) = ((3*f(nx,ny)-4*f(nx-1,ny)+f(nx-2,ny))/2/dxi)*xiy(nx,ny)+((3*f(nx,ny)-4*f(nx,ny-1)+f(nx,ny-2))/2/deta)*etay(nx,ny);





end
