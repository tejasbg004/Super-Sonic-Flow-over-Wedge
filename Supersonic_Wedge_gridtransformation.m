function [u,v,T,P,rho,dxi,deta,U,e,etaeta,dydeta,J]=Supersonic_double_cone_gridtransformation(yy1,xx)
    nx=100;
    ny=100;
    %%%%%%sta=50;
    L=5*10^-1;
    H=5*10^-1;
    M=25;
    Tinf=288.15;
    Pinf=101300;
    %rhoinf=1.225;
    mu_o=1.735*10^-5;
    %mu_o=0;
    Prandtl=0.71;
    cp=1005;
    cv=718;
    gamma=cp/cv;
    R=287;
    Uinf=M*(sqrt(gamma*R*Tinf));
    xi=linspace(0,L,nx);
    eta=linspace(0,H,ny);
    [xixi,etaeta]=ndgrid(xi,eta);
    dxi=xixi(2,1)-xixi(1,1);
    deta=etaeta(1,2)-etaeta(1,1);
    addpath('Derivative_Functions')
    dxdxi=ddxi_central(xx,dxi);
    dxdeta=ddeta_central(xx,deta);
    dydxi=ddxi_central(yy1,dxi);
    dydeta=ddeta_central(yy1,deta);
    J=dxdxi.*dydeta-dxdeta.*dydxi;
    
    xiy=-dxdeta./J;
    etay=dxdxi./J;

    
    xix=dydeta./J;
    etax=-dydxi./J;
    dt=2.35*10^-11;
    u=zeros(nx,ny,1);
    v=zeros(nx,ny,1);
    T=zeros(nx,ny,1);
    P=zeros(nx,ny,1);
    %Initial Condition
    u(2:nx,2:ny-1,1)=Uinf;
    %%%%u(1:sta,1,1)=Uinf;%cha
    v(2:nx,2:ny-1,1)=0;
    %%%%v(1:sta,1,1)=0;
    P(2:nx,2:ny-1,1)=Pinf;
    T(2:nx,2:ny-1,1)=Tinf;


   %u(:,1,:)=(4*u(:,2)-u(:,3))/3;%wall cha
    u(:,1,:)=0;%wall cha

   v(:,1,:)=0;%wall cha
    
    T(:,1,:)=Tinf;%wall cha
    %far field u
    u(1,:,:)=Uinf;
    u(:,ny,:)=Uinf;
    %farfield v
    v(1,:,:)=0;
    v(:,ny,:)=0;
    %pressure field(inlet and farfiel)
    P(1,:,:)=Pinf;
    P(:,ny,:)=Pinf;
    %farfield T
    T(1,:,:)=Tinf;
    T(:,ny,:)=Tinf;
    P(:,1,1)=2*P(:,2,1)-P(:,3,1);%wall
    %leading edge u
    u(1,1,:)=0;%cha
    v(1,1,:)=0;
    P(1,1,:)=Pinf;
    T(1,1,:)=Tinf;

    %outflow
    u(nx,:,:)=2*u(nx-1,:,:)-u(nx-2,:,:);
    v(nx,:,:)=2*v(nx-1,:,:)-v(nx-2,:,:);
    T(nx,:,:)=2*T(nx-1,:,:)-T(nx-2,:,:);
    P(nx,:,:)=2*P(nx-1,:,:)-P(nx-2,:,:);
%         u(nx,:,:)=Uinf;
%         v(nx,:,:)=0;
%         T(nx,:,:)=Tinf;
%         P(nx,:,:)=Pinf;


   
    rho=zeros(nx,ny,1);
    for i=1:nx
        for j=1:ny
            rho(i,j)=P(i,j)/(R*T(i,j));
        end
    end

    

   Enew=zeros(4,nx,ny);

   Fnew=zeros(4,nx,ny);

   residual=zeros(1501,1);

   
   n=0;
   y=0;
   figure(1)
   iter=3500;
%     vid=VideoWriter('Wedge','MPEG-4');
%     vid.FrameRate=30;
    %open(vid)
   while n<iter
       U=prim2cons(rho,u,v,T,cv);
       Uinit=U;
       U1=J.*squeeze(U(1,:,:));
       U2=J.*squeeze(U(2,:,:));
       U3=J.*squeeze(U(3,:,:));
       U4=J.*squeeze(U(4,:,:));
       U(1,:,:)=U1;
       U(2,:,:)=U2;
       U(3,:,:)=U3;
       U(4,:,:)=U4;
   
       Uinital=U;
      
       
       mu=sutherland(T,mu_o);
       K=(cp/Prandtl).*mu;
       %updating E
       dudxi=ddx_bwd_xi(u,dxi,deta,yy1,xx);
       dvdeta=ddy_central_eta(v,dxi,deta,yy1,xx);
       dTdxi=ddx_bwd_xi(T,dxi,deta,yy1,xx);
       dudeta=ddy_central_eta(u,dxi,deta,yy1,xx);
       dvdxi=ddx_bwd_xi(v,dxi,deta,yy1,xx);
       tau_xx_new=2*mu.*(dudxi-1/3.*(dudxi+dvdeta));
       tau_xy_new=mu.*(dudeta+dvdxi);
       for k=1:4
        for j=1:ny
            for i=1:nx
                if k==1
                    Enew(k,i,j)=rho(i,j)*u(i,j);
                elseif k==2
                    Enew(k,i,j)=rho(i,j)*u(i,j)*u(i,j)+P(i,j)-tau_xx_new(i,j);
                elseif k==3
                    Enew(k,i,j)=rho(i,j)*u(i,j)*v(i,j)-tau_xy_new(i,j);
                else
                    Enew(k,i,j)=(rho(i,j)*(cv*T(i,j)+0.5*(u(i,j)^2+v(i,j)^2))+P(i,j))*u(i,j)-u(i,j)*tau_xx_new(i,j)-v(i,j)*tau_xy_new(i,j)-K(i,j)*dTdxi(i,j);
                end
            end
        end
       end

       %updating F
       dudeta=ddy_bwd_eta(u,dxi,deta,yy1,xx);
       dvdxi=ddx_central_xi(v,dxi,deta,yy1,xx);
       dvdeta=ddy_bwd_eta(v,dxi,deta,yy1,xx);
       dudxi=ddx_central_xi(u,dxi,deta,yy1,xx);
       dTdeta=ddy_bwd_eta(T,dxi,deta,yy1,xx);
       tau_xy_new=mu.*(dudeta+dvdxi);
       tau_yy_new=2*mu.*(dvdeta-1/3*(dudxi+dvdeta));

       for k=1:4
        for j=1:ny
            for i=1:nx
                if k==1
                    Fnew(k,i,j)=rho(i,j)*v(i,j);
                elseif k==2
                    Fnew(k,i,j)=rho(i,j)*u(i,j)*v(i,j)-tau_xy_new(i,j);
                elseif k==3
                    Fnew(k,i,j)=rho(i,j)*v(i,j)*v(i,j)+P(i,j)-tau_yy_new(i,j);  
                else
                    Fnew(k,i,j)=(rho(i,j)*(cv*T(i,j)+0.5*(u(i,j)^2+v(i,j)^2))+P(i,j))*v(i,j)-u(i,j)*tau_xy_new(i,j)-v(i,j)*tau_yy_new(i,j)-K(i,j)*dTdeta(i,j);
                end
            end
        end
       end

       E_eps(1,:,:)=J.*squeeze(Enew(1,:,:)).*xix+J.*squeeze(Fnew(1,:,:)).*xiy;
       E_eps(2,:,:)=J.*squeeze(Enew(2,:,:)).*xix+J.*squeeze(Fnew(2,:,:)).*xiy;
       E_eps(3,:,:)=J.*squeeze(Enew(3,:,:)).*xix+J.*squeeze(Fnew(3,:,:)).*xiy;
       E_eps(4,:,:)=J.*squeeze(Enew(4,:,:)).*xix+J.*squeeze(Fnew(4,:,:)).*xiy;

       F_eps(1,:,:)=J.*squeeze(Enew(1,:,:)).*etax+J.*squeeze(Fnew(1,:,:)).*etay;
       F_eps(2,:,:)=J.*squeeze(Enew(2,:,:)).*etax+J.*squeeze(Fnew(2,:,:)).*etay;
       F_eps(3,:,:)=J.*squeeze(Enew(3,:,:)).*etax+J.*squeeze(Fnew(3,:,:)).*etay;
       F_eps(4,:,:)=J.*squeeze(Enew(4,:,:)).*etax+J.*squeeze(Fnew(4,:,:)).*etay;
       
       
       
       dEdx1=ddxi_fwd(squeeze(E_eps(1,:,:)),dxi);
       dEdx2=ddxi_fwd(squeeze(E_eps(2,:,:)),dxi);
       dEdx3=ddxi_fwd(squeeze(E_eps(3,:,:)),dxi);
       dEdx4=ddxi_fwd(squeeze(E_eps(4,:,:)),dxi);
       dEdx(1,:,:)=dEdx1;
       dEdx(2,:,:)=dEdx2;
       dEdx(3,:,:)=dEdx3;
       dEdx(4,:,:)=dEdx4;

       dFdy1=deta_fwd(squeeze(F_eps(1,:,:)),deta);
       dFdy2=deta_fwd(squeeze(F_eps(2,:,:)),deta);
       dFdy3=deta_fwd(squeeze(F_eps(3,:,:)),deta);
       dFdy4=deta_fwd(squeeze(F_eps(4,:,:)),deta);
       dFdy(1,:,:)=dFdy1;
       dFdy(2,:,:)=dFdy2;
       dFdy(3,:,:)=dFdy3;
       dFdy(4,:,:)=dFdy4;
       
       Ubar=U-dt.*(dEdx+dFdy);
      
       U1=squeeze(Ubar(1,:,:))./J;
       U2=squeeze(Ubar(2,:,:))./J;
       U3=squeeze(Ubar(3,:,:))./J;
       U4=squeeze(Ubar(4,:,:))./J;
       Ubar(1,:,:)=U1;
       Ubar(2,:,:)=U2;
       Ubar(3,:,:)=U3;
       Ubar(4,:,:)=U4;
   
       [rho,u,v,T,P] = cons2prim(Ubar,R,cv);
       
%        U1=squeeze(Ubar(1,:,:)).*J;
%        U2=squeeze(Ubar(2,:,:)).*J;
%        U3=squeeze(Ubar(3,:,:)).*J;
%        U4=squeeze(Ubar(4,:,:)).*J;
%        Ubar(1,:,:)=U1;
%        Ubar(2,:,:)=U2;
%        Ubar(3,:,:)=U3;
%        Ubar(4,:,:)=U4;
   
    %outflow
    u(nx,:,:)=2*u(nx-1,:,:)-u(nx-2,:,:);
    v(nx,:,:)=2*v(nx-1,:,:)-v(nx-2,:,:);
    T(nx,:,:)=2*T(nx-1,:,:)-T(nx-2,:,:);
    P(nx,:,:)=2*P(nx-1,:,:)-P(nx-2,:,:);
%         u(nx,:,:)=Uinf;
%         v(nx,:,:)=0;
%         T(nx,:,:)=Tinf;
%         P(nx,:,:)=Pinf;

    %u(:,1,:)=(4*u(:,2)-u(:,3))/3;%wall cha
    u(:,1,:)=0;%wall cha

    v(:,1,:)=0;%wall cha
    T(:,1,:)=Tinf;%wall
    P(:,1,:)=2*P(:,2,:)-P(:,3,:);%wall pressure extrapolation
    %far field u
    u(1,:,:)=Uinf;
    u(:,ny,:)=Uinf;
    %farfield v
    v(1,:,:)=0;
    v(:,ny,:)=0;
     %pressure field(inlet and farfiel)
    P(1,:,:)=Pinf;
    P(:,ny,:)=Pinf;
    %farfield T
    T(1,:,:)=Tinf;
    T(:,ny,:)=Tinf;
    %leading edge u
    u(1,1,:)=0;%cha
    v(1,1,:)=0;
    P(1,1,:)=Pinf;
    T(1,1,:)=Tinf;

   Utmp=prim2cons(rho,u,v,T,cv);
   U1=squeeze(Utmp(1,:,:)).*J;
   U2=squeeze(Utmp(2,:,:)).*J;
   U3=squeeze(Utmp(3,:,:)).*J;
   U4=squeeze(Utmp(4,:,:)).*J;
   Ubar(1,:,:)=U1;
   Ubar(2,:,:)=U2;
   Ubar(3,:,:)=U3;
   Ubar(4,:,:)=U4;

     

    for i=1:nx
        for j=1:ny
            rho(i,j)=P(i,j)/(R*T(i,j));
        end
    end
       mu=sutherland(T,mu_o);
       K=(cp/Prandtl).*mu;
       dudxi=ddx_fwd_xi(u,dxi,deta,yy1,xx);
       dvdeta=ddy_central_eta(v,dxi,deta,yy1,xx);
       dTdxi=ddx_fwd_xi(T,dxi,deta,yy1,xx);
       dudeta=ddy_central_eta(u,dxi,deta,yy1,xx);
       dvdxi=ddx_fwd_xi(v,dxi,deta,yy1,xx);
       tau_xx_new=2*mu.*(dudxi-1/3.*(dudxi+dvdeta));
       tau_xy_new=mu.*(dudeta+dvdxi);

       for k=1:4
        for j=1:ny
            for i=1:nx
                if k==1
                    Enew(k,i,j)=rho(i,j)*u(i,j);
                elseif k==2
                    Enew(k,i,j)=rho(i,j)*u(i,j)*u(i,j)+P(i,j)-tau_xx_new(i,j);
                elseif k==3
                    Enew(k,i,j)=rho(i,j)*u(i,j)*v(i,j)-tau_xy_new(i,j);
                else
                    Enew(k,i,j)=(rho(i,j)*(cv*T(i,j)+0.5*(u(i,j)^2+v(i,j)^2))+P(i,j))*u(i,j)-u(i,j)*tau_xx_new(i,j)-v(i,j)*tau_xy_new(i,j)-K(i,j)*dTdxi(i,j);
                end
            end
        end
       end

      %updating F
       dudeta=ddy_fwd_eta(u,dxi,deta,yy1,xx);
       dvdxi=ddx_central_xi(v,dxi,deta,yy1,xx);
       dvdeta=ddy_fwd_eta(v,dxi,deta,yy1,xx);
       dudxi=ddx_central_xi(u,dxi,deta,yy1,xx);
       dTdeta=ddy_fwd_eta(T,dxi,deta,yy1,xx);
       tau_xy_new=mu.*(dudeta+dvdxi);
       tau_yy_new=2*mu.*(dvdeta-1/3*(dudxi+dvdeta)); 

       for k=1:4
        for j=1:ny
            for i=1:nx
                if k==1
                    Fnew(k,i,j)=rho(i,j)*v(i,j);
                elseif k==2
                    Fnew(k,i,j)=rho(i,j)*u(i,j)*v(i,j)-tau_xy_new(i,j);
                elseif k==3
                    Fnew(k,i,j)=rho(i,j)*v(i,j)*v(i,j)+P(i,j)-tau_yy_new(i,j);  
                else
                    Fnew(k,i,j)=(rho(i,j)*(cv*T(i,j)+0.5*(u(i,j)^2+v(i,j)^2))+P(i,j))*v(i,j)-u(i,j)*tau_xy_new(i,j)-v(i,j)*tau_yy_new(i,j)-K(i,j)*dTdeta(i,j);
                end
            end
        end
       end

       E_eps(1,:,:)=J.*squeeze(Enew(1,:,:)).*xix+J.*squeeze(Fnew(1,:,:)).*xiy;
       E_eps(2,:,:)=J.*squeeze(Enew(2,:,:)).*xix+J.*squeeze(Fnew(2,:,:)).*xiy;
       E_eps(3,:,:)=J.*squeeze(Enew(3,:,:)).*xix+J.*squeeze(Fnew(3,:,:)).*xiy;
       E_eps(4,:,:)=J.*squeeze(Enew(4,:,:)).*xix+J.*squeeze(Fnew(4,:,:)).*xiy;

       F_eps(1,:,:)=J.*squeeze(Enew(1,:,:)).*etax+J.*squeeze(Fnew(1,:,:)).*etay;
       F_eps(2,:,:)=J.*squeeze(Enew(2,:,:)).*etax+J.*squeeze(Fnew(2,:,:)).*etay;
       F_eps(3,:,:)=J.*squeeze(Enew(3,:,:)).*etax+J.*squeeze(Fnew(3,:,:)).*etay;
       F_eps(4,:,:)=J.*squeeze(Enew(4,:,:)).*etax+J.*squeeze(Fnew(4,:,:)).*etay;
       
        
       dEdx1=ddxi_bwd(squeeze(E_eps(1,:,:)),dxi);
       dEdx2=ddxi_bwd(squeeze(E_eps(2,:,:)),dxi);
       dEdx3=ddxi_bwd(squeeze(E_eps(3,:,:)),dxi);
       dEdx4=ddxi_bwd(squeeze(E_eps(4,:,:)),dxi);
       dEdx(1,:,:)=dEdx1;
       dEdx(2,:,:)=dEdx2;
       dEdx(3,:,:)=dEdx3;
       dEdx(4,:,:)=dEdx4;

       dFdy1=deta_bwd(squeeze(F_eps(1,:,:)),deta);
       dFdy2=deta_bwd(squeeze(F_eps(2,:,:)),deta);
       dFdy3=deta_bwd(squeeze(F_eps(3,:,:)),deta);
       dFdy4=deta_bwd(squeeze(F_eps(4,:,:)),deta);
       dFdy(1,:,:)=dFdy1;
       dFdy(2,:,:)=dFdy2;
       dFdy(3,:,:)=dFdy3;
       dFdy(4,:,:)=dFdy4;

       U=0.5*((Ubar+Uinital)-dt*(dEdx)-dt*(dFdy));

       U1=squeeze(U(1,:,:))./J;
       U2=squeeze(U(2,:,:))./J;
       U3=squeeze(U(3,:,:))./J;
       U4=squeeze(U(4,:,:))./J;
       U(1,:,:)=U1;
       U(2,:,:)=U2;
       U(3,:,:)=U3;
       U(4,:,:)=U4;

       [rho,u,v,T,P] = cons2prim(U,R,cv);

     %outflow
    u(nx,:,:)=2*u(nx-1,:,:)-u(nx-2,:,:);
    v(nx,:,:)=2*v(nx-1,:,:)-v(nx-2,:,:);
    T(nx,:,:)=2*T(nx-1,:,:)-T(nx-2,:,:);
    P(nx,:,:)=2*P(nx-1,:,:)-P(nx-2,:,:);
%         u(nx,:,:)=Uinf;
%         v(nx,:,:)=0;
%         T(nx,:,:)=Tinf;
%         P(nx,:,:)=Pinf;

    %u(:,1,:)=(4*u(:,2)-u(:,3))/3;%wall cha
    u(:,1,:)=0;%wall cha

    v(:,1,:)=0;%wall cha
    T(:,1,:)=Tinf;%wall
    P(:,1,:)=2*P(:,2,:)-P(:,3,:);%wall pressure extrapolation
    %far field u
    u(1,:,:)=Uinf;
    u(:,ny,:)=Uinf;
    %farfield v
    v(1,:,:)=0;
    v(:,ny,:)=0;
    %pressure field(inlet and farfiel)
    P(1,:,:)=Pinf;
    P(:,ny,:)=Pinf;
    %farfield T
    T(1,:,:)=Tinf;
    T(:,ny,:)=Tinf;
    %leading edge u
    u(1,1,:)=0;%cha
    v(1,1,:)=0;
    P(1,1,:)=Pinf;
    T(1,1,:)=Tinf;
    e=cv*T;
    for i=1:nx
        for j=1:ny
            rho(i,j)=P(i,j)/(R*T(i,j));
        end
    end

    U=prim2cons(rho,u,v,T,cv);
    
    residual(n+1)=norm(squeeze(Uinit(2,:,:))-squeeze(U(2,:,:)));
    %Uinital=U;
    drawnow;


    
    if n==y+100
        subplot(3,2,1)
        pcolor(xx,yy1,(u));
        colormap("Jet")
        colorbar
        shading interp
        axis equal tight
        xlabel('x')
        ylabel('y')
        title('u-velocity','FontSize',15)
        subplot(3,2,2)
        pcolor(xx,yy1,(v));
        colormap("Jet")
        colorbar
        shading interp
        axis equal tight
        xlabel('x')
        ylabel('y')
        title('v-velocity','FontSize',15)
        subplot(3,2,3)
        pcolor(xx,yy1,(rho));
        colormap("Jet")
        colorbar
        shading interp
        axis equal tight
        xlabel('x')
        ylabel('y')
        title('Density','FontSize',15)
        subplot(3,2,4)
        pcolor(xx,yy1,(P));
        colormap("Jet")
        colorbar
        shading interp
        axis equal tight
        xlabel('x')
        ylabel('y')
        title('Pressure','FontSize',15)
        subplot(3,2,5)
        pcolor(xx,yy1,(e));
        colormap("Jet")
        colorbar
        shading interp
        axis equal tight
        xlabel('x')
        ylabel('y')
        title('Internal Energy','FontSize',15)
        subplot(3,2,6)
        pcolor(xx,yy1,(T));
        colormap("Jet")
        colorbar
        shading interp
        axis equal tight
        xlabel('x')
        ylabel('y')
        title('Temperature','FontSize',15)
        drawnow
        y=n;
%         frame=getframe(gcf);
%         writeVideo(vid,frame)
    end
    n=n+1;
   end
   machnumber=zeros(nx,ny);
   for i =1:nx
       for j=1:ny
           machnumber(i,j)=u(i,j)/(sqrt(gamma*R*T(i,j)));
       end
   end

   figure(2)
   nn=linspace(1,iter,iter);
   plot(nn,residual)
   xlabel('Number Of Time Steps')
   ylabel('U-velocity residual')
   title("Convergence plot of velocity for M="+M,'FontSize',15)

   figure(3)
       pcolor(xx,yy1,machnumber)
       colormap("Jet")
       colorbar
       shading interp
       axis equal tight
       xlabel('x')
       ylabel('y')
       title('Mach number','FontSize',15)
inde=floor(H/(deta*4))+1;
figure(4)
plot(yy1(:,inde),machnumber(:,inde));
xlabel('x')
ylabel('Machnumber')
title("Mach number variation along x-direction for M="+M,'FontSize',15)

figure(5)
plot(yy1(:,inde),P(:,inde));
xlabel('x')
ylabel('Pressure')
title("Pressure variation along x-direction for M="+M,'FontSize',15)

figure(6)
plot(yy1(:,inde),T(:,inde));
xlabel('x')
ylabel('Temperature')
title("Temperature variation along x-direction for M="+M,'FontSize',15)
   
end