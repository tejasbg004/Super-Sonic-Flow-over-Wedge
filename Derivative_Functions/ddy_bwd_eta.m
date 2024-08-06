function [dfdy]= ddy_bwd_eta(f,dxi,deta,yy1,xx)
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
    
    % backward differencing
    for i=2:nx
        for j=2:ny
            dfdy(i,j) = ((f(i,j)-f(i-1,j))/dxi)*xiy(i,j)+((f(i,j)-f(i,j-1))/deta)*etay(i,j);
        end
    end
    
    % forward difference for first collumn
    j = 1;
    for i=2:nx
        dfdy(i,j) = ((f(i,j)-f(i-1,j))/dxi)*xiy(i,j)+((f(i,j+1)-f(i,j))/deta)*etay(i,j);
    end

 % forward difference for first row
    i=1;
    for j=2:ny
        dfdy(i,j) = ((f(i+1,j)-f(i,j))/dxi)*xiy(i,j)+((f(i,j)-f(i,j-1))/deta)*etay(i,j);
    end

    dfdy(1,1) = ((f(1+1,1)-f(1,1))/dxi)*xiy(1,1)+((f(1,1+1)-f(1,1))/deta)*etay(1,1);


end