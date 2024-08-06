function [u2,v2,T2,P2,rho2,dxi,deta,U1,e,etaeta,dydeta,J]=Supersonic_double_cone_gridtransformation_adiabatic(yy1,xx)
    nx=100;
    ny=100;
    %%%%%%sta=50;
    L=1*10^-6;
    H=1*10^-6;
    M=3;
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
    dt=2.35*10^-11;      %%%%%%%%%%%%%%time
    u2=zeros(nx,ny,1);
    v2=zeros(nx,ny,1);
    T2=zeros(nx,ny,1);
    P2=zeros(nx,ny,1);
    %Initial Condition
    u2(2:nx,2:ny-1,1)=Uinf;
    %%%%u(1:sta,1,1)=Uinf;%cha
    v2(2:nx,2:ny-1,1)=0;
    %%%%v(1:sta,1,1)=0;
    P2(2:nx,2:ny-1,1)=Pinf;
    T2(2:nx,2:ny-1,1)=Tinf;


    u2(:,1,:)=0;%wall cha

    v2(:,1,:)=0;%wall cha
   % Adiabatic wall condition
   for i=2:nx-1
       T2(i,1)=(4*T2(i,2)-T2(i,3)+(T2(i+1,1)-T2(i-1,1))*(deta/dxi)*(xiy(i,1)/etay(i,1)))/3;
   end

   c11=3*(xiy(nx,1)/(dxi)-etay(nx,1)/(deta));

  T2(nx,1)=((4*T2(nx-1,1)-T2(nx-2,1))*(xiy(nx,1))/(dxi) + (-4*T2(nx,2)+T2(nx,3))*(etay(nx,1)/(deta)))/c11;
   
   % T2(:,1)=Tinf;

    %far field u
    u2(1,:,:)=Uinf;
    u2(:,ny,:)=Uinf;
    %farfield v
    v2(1,:,:)=0;
    v2(:,ny,:)=0;
    %pressure field(inlet and farfield)
    P2(1,:,:)=Pinf;
    P2(:,ny,:)=Pinf;
    %farfield T
    T2(1,:,:)=Tinf;
    T2(:,ny,:)=Tinf;
    P2(:,1,1)=2*P2(:,2,1)-P2(:,3,1);%wall
    %leading edge u
    u2(1,1,:)=0;%cha
    v2(1,1,:)=0;
    P2(1,1,:)=Pinf;
    T2(1,1,:)=Tinf;

   
    rho2=zeros(nx,ny,1);
    for i=1:nx
        for j=1:ny
            rho2(i,j)=P2(i,j)/(R*T2(i,j));
        end
    end

    

   Enew2=zeros(4,nx,ny);

   Fnew2=zeros(4,nx,ny);

  

   
   n=0;
   y=0;
   figure(1)
   iter=3500;
    residual=zeros(iter,1);

    vid=VideoWriter('Wedge_adiabatic1','MPEG-4');
    vid.FrameRate=30;
    open(vid)

   while n<iter
       U1=prim2cons(rho2,u2,v2,T2,cv);
       Uinit=U1;
       U11=J.*squeeze(U1(1,:,:));
       U21=J.*squeeze(U1(2,:,:));
       U31=J.*squeeze(U1(3,:,:));
       U41=J.*squeeze(U1(4,:,:));
       U1(1,:,:)=U11;
       U1(2,:,:)=U21;
       U1(3,:,:)=U31;
       U1(4,:,:)=U41;
   
       Uinital=U1;
      
       
       mu2=sutherland(T2,mu_o);
       K2=(cp/Prandtl).*mu2;
       %updating E
       dudxi=ddx_bwd_xi(u2,dxi,deta,yy1,xx);
       dvdeta=ddy_central_eta(v2,dxi,deta,yy1,xx);
       dTdxi=ddx_bwd_xi(T2,dxi,deta,yy1,xx);
       dudeta=ddy_central_eta(u2,dxi,deta,yy1,xx);
       dvdxi=ddx_bwd_xi(v2,dxi,deta,yy1,xx);
       tau_xx_new=2*mu2.*(dudxi-1/3.*(dudxi+dvdeta));
       tau_xy_new=mu2.*(dudeta+dvdxi);
       for k=1:4
        for j=1:ny
            for i=1:nx
                if k==1
                    Enew2(k,i,j)=rho2(i,j)*u2(i,j);
                elseif k==2
                    Enew2(k,i,j)=rho2(i,j)*u2(i,j)*u2(i,j)+P2(i,j)-tau_xx_new(i,j);
                elseif k==3
                    Enew2(k,i,j)=rho2(i,j)*u2(i,j)*v2(i,j)-tau_xy_new(i,j);
                else
                    Enew2(k,i,j)=(rho2(i,j)*(cv*T2(i,j)+0.5*(u2(i,j)^2+v2(i,j)^2))+P2(i,j))*u2(i,j)-u2(i,j)*tau_xx_new(i,j)-v2(i,j)*tau_xy_new(i,j)-K2(i,j)*dTdxi(i,j);
                end
            end
        end
       end

       %updating F
       dudeta=ddy_bwd_eta(u2,dxi,deta,yy1,xx);
       dvdxi=ddx_central_xi(v2,dxi,deta,yy1,xx);
       dvdeta=ddy_bwd_eta(v2,dxi,deta,yy1,xx);
       dudxi=ddx_central_xi(u2,dxi,deta,yy1,xx);
       dTdeta=ddy_bwd_eta(T2,dxi,deta,yy1,xx);
   
       tau_xy_new=mu2.*(dudeta+dvdxi);
       tau_yy_new=2*mu2.*(dvdeta-1/3*(dudxi+dvdeta));

       for k=1:4
        for j=1:ny
            for i=1:nx
                if k==1
                    Fnew2(k,i,j)=rho2(i,j)*v2(i,j);
                elseif k==2
                    Fnew2(k,i,j)=rho2(i,j)*u2(i,j)*v2(i,j)-tau_xy_new(i,j);
                elseif k==3
                    Fnew2(k,i,j)=rho2(i,j)*v2(i,j)*v2(i,j)+P2(i,j)-tau_yy_new(i,j);  
                else
                    Fnew2(k,i,j)=(rho2(i,j)*(cv*T2(i,j)+0.5*(u2(i,j)^2+v2(i,j)^2))+P2(i,j))*v2(i,j)-u2(i,j)*tau_xy_new(i,j)-v2(i,j)*tau_yy_new(i,j)-K2(i,j)*dTdeta(i,j);
                end
            end
        end
       end

       E_eps(1,:,:)=J.*squeeze(Enew2(1,:,:)).*xix+J.*squeeze(Fnew2(1,:,:)).*xiy;
       E_eps(2,:,:)=J.*squeeze(Enew2(2,:,:)).*xix+J.*squeeze(Fnew2(2,:,:)).*xiy;
       E_eps(3,:,:)=J.*squeeze(Enew2(3,:,:)).*xix+J.*squeeze(Fnew2(3,:,:)).*xiy;
       E_eps(4,:,:)=J.*squeeze(Enew2(4,:,:)).*xix+J.*squeeze(Fnew2(4,:,:)).*xiy;

       F_eps(1,:,:)=J.*squeeze(Enew2(1,:,:)).*etax+J.*squeeze(Fnew2(1,:,:)).*etay;
       F_eps(2,:,:)=J.*squeeze(Enew2(2,:,:)).*etax+J.*squeeze(Fnew2(2,:,:)).*etay;
       F_eps(3,:,:)=J.*squeeze(Enew2(3,:,:)).*etax+J.*squeeze(Fnew2(3,:,:)).*etay;
       F_eps(4,:,:)=J.*squeeze(Enew2(4,:,:)).*etax+J.*squeeze(Fnew2(4,:,:)).*etay;
       
       
       
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
       
       Ubar=U1-dt.*(dEdx+dFdy);
      
       U11=squeeze(Ubar(1,:,:))./J;
       U21=squeeze(Ubar(2,:,:))./J;
       U31=squeeze(Ubar(3,:,:))./J;
       U41=squeeze(Ubar(4,:,:))./J;
       Ubar(1,:,:)=U11;
       Ubar(2,:,:)=U21;
       Ubar(3,:,:)=U31;
       Ubar(4,:,:)=U41;
    
   
       [rho2,u2,v2,T2,P2] = cons2prim(Ubar,R,cv);
       
       
      
    %outflow
    u2(nx,:,:)=2*u2(nx-1,:,:)-u2(nx-2,:,:);
    v2(nx,:,:)=2*v2(nx-1,:,:)-v2(nx-2,:,:);
    T2(nx,2:ny,:)=2*T2(nx-1,2:ny,:)-T2(nx-2,2:ny,:);
    P2(nx,:,:)=2*P2(nx-1,:,:)-P2(nx-2,:,:);

    for i=2:nx-1
       T2(i,1)=(4*T2(i,2)-T2(i,3)+(T2(i+1,1)-T2(i-1,1))*(deta/dxi)*(xiy(i,1)/etay(i,1)))/3;
   end

   c11=3*(xiy(nx,1)/(dxi)-etay(nx,1)/(deta));

  T2(nx,1)=((4*T2(nx-1,1)-T2(nx-2,1))*(xiy(nx,1))/(dxi) + (-4*T2(nx,2)+T2(nx,3))*(etay(nx,1)/(deta)))/c11;

   % T2(:,1)=Tinf;
    u2(:,1,:)=0;
    v2(:,1,:)=0;%wall cha
    %%T2(:,1,:)=Tinf;no wall temp const.
    P2(:,1,:)=2*P2(:,2,:)-P2(:,3,:);%wall pressure extrapolation
    %far field u
    u2(1,:,:)=Uinf;
    u2(:,ny,:)=Uinf;
    %farfield v
    v2(1,:,:)=0;
    v2(:,ny,:)=0;
     %pressure field(inlet and farfiel)
    P2(1,:,:)=Pinf;
    P2(:,ny,:)=Pinf;
    %farfield T
    T2(1,:,:)=Tinf;
    T2(:,ny,:)=Tinf;
    %leading edge u
    u2(1,1,:)=0;%cha
    v2(1,1,:)=0;
    P2(1,1,:)=Pinf;
    T2(1,1,:)=Tinf; 
        Utmp=prim2cons(rho2,u2,v2,T2,cv);
       U11=squeeze(Utmp(1,:,:)).*J;
       U21=squeeze(Utmp(2,:,:)).*J;
       U31=squeeze(Utmp(3,:,:)).*J;
       U41=squeeze(Utmp(4,:,:)).*J;
       Ubar(1,:,:)=U11;
       Ubar(2,:,:)=U21;
       Ubar(3,:,:)=U31;
       Ubar(4,:,:)=U41;

    for i=1:nx
        for j=1:ny
            rho2(i,j)=P2(i,j)/(R*T2(i,j));
        end
    end
       mu2=sutherland(T2,mu_o);
       K2=(cp/Prandtl).*mu2;
       dudxi=ddx_fwd_xi(u2,dxi,deta,yy1,xx);
       dvdeta=ddy_central_eta(v2,dxi,deta,yy1,xx);
       dTdxi=ddx_fwd_xi(T2,dxi,deta,yy1,xx);
       dudeta=ddy_central_eta(u2,dxi,deta,yy1,xx);
       dvdxi=ddx_fwd_xi(v2,dxi,deta,yy1,xx);
       tau_xx_new=2*mu2.*(dudxi-1/3.*(dudxi+dvdeta));
       tau_xy_new=mu2.*(dudeta+dvdxi);

       for k=1:4
        for j=1:ny
            for i=1:nx
                if k==1
                    Enew2(k,i,j)=rho2(i,j)*u2(i,j);
                elseif k==2
                    Enew2(k,i,j)=rho2(i,j)*u2(i,j)*u2(i,j)+P2(i,j)-tau_xx_new(i,j);
                elseif k==3
                    Enew2(k,i,j)=rho2(i,j)*u2(i,j)*v2(i,j)-tau_xy_new(i,j);
                else
                    Enew2(k,i,j)=(rho2(i,j)*(cv*T2(i,j)+0.5*(u2(i,j)^2+v2(i,j)^2))+P2(i,j))*u2(i,j)-u2(i,j)*tau_xx_new(i,j)-v2(i,j)*tau_xy_new(i,j)-K2(i,j)*dTdxi(i,j);
                end
            end
        end
       end

      %updating F
       dudeta=ddy_fwd_eta(u2,dxi,deta,yy1,xx);
       dvdxi=ddx_central_xi(v2,dxi,deta,yy1,xx);
       dvdeta=ddy_fwd_eta(v2,dxi,deta,yy1,xx);
       dudxi=ddx_central_xi(u2,dxi,deta,yy1,xx);
       dTdeta=ddy_fwd_eta(T2,dxi,deta,yy1,xx);
     
       tau_xy_new=mu2.*(dudeta+dvdxi);
       tau_yy_new=2*mu2.*(dvdeta-1/3*(dudxi+dvdeta)); 

       for k=1:4
        for j=1:ny
            for i=1:nx
                if k==1
                    Fnew2(k,i,j)=rho2(i,j)*v2(i,j);
                elseif k==2
                    Fnew2(k,i,j)=rho2(i,j)*u2(i,j)*v2(i,j)-tau_xy_new(i,j);
                elseif k==3
                    Fnew2(k,i,j)=rho2(i,j)*v2(i,j)*v2(i,j)+P2(i,j)-tau_yy_new(i,j);  
                else
                    Fnew2(k,i,j)=(rho2(i,j)*(cv*T2(i,j)+0.5*(u2(i,j)^2+v2(i,j)^2))+P2(i,j))*v2(i,j)-u2(i,j)*tau_xy_new(i,j)-v2(i,j)*tau_yy_new(i,j)-K2(i,j)*dTdeta(i,j);
                end
            end
        end
       end

       E_eps(1,:,:)=J.*squeeze(Enew2(1,:,:)).*xix+J.*squeeze(Fnew2(1,:,:)).*xiy;
       E_eps(2,:,:)=J.*squeeze(Enew2(2,:,:)).*xix+J.*squeeze(Fnew2(2,:,:)).*xiy;
       E_eps(3,:,:)=J.*squeeze(Enew2(3,:,:)).*xix+J.*squeeze(Fnew2(3,:,:)).*xiy;
       E_eps(4,:,:)=J.*squeeze(Enew2(4,:,:)).*xix+J.*squeeze(Fnew2(4,:,:)).*xiy;

       F_eps(1,:,:)=J.*squeeze(Enew2(1,:,:)).*etax+J.*squeeze(Fnew2(1,:,:)).*etay;
       F_eps(2,:,:)=J.*squeeze(Enew2(2,:,:)).*etax+J.*squeeze(Fnew2(2,:,:)).*etay;
       F_eps(3,:,:)=J.*squeeze(Enew2(3,:,:)).*etax+J.*squeeze(Fnew2(3,:,:)).*etay;
       F_eps(4,:,:)=J.*squeeze(Enew2(4,:,:)).*etax+J.*squeeze(Fnew2(4,:,:)).*etay;
       
        
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

       U1=0.5*((Ubar+Uinital)-dt*(dEdx)-dt*(dFdy));

       U11=squeeze(U1(1,:,:))./J;
       U21=squeeze(U1(2,:,:))./J;
       U31=squeeze(U1(3,:,:))./J;
       U41=squeeze(U1(4,:,:))./J;
       U1(1,:,:)=U11;
       U1(2,:,:)=U21;
       U1(3,:,:)=U31;
       U1(4,:,:)=U41;

       [rho2,u2,v2,T2,P2] = cons2prim(U1,R,cv);

     %outflow
    u2(nx,:,:)=2*u2(nx-1,:,:)-u2(nx-2,:,:);
    v2(nx,:,:)=2*v2(nx-1,:,:)-v2(nx-2,:,:);
    T2(nx,:,:)=2*T2(nx-1,:,:)-T2(nx-2,:,:);
    P2(nx,:,:)=2*P2(nx-1,:,:)-P2(nx-2,:,:);

    for i=2:nx-1
       T2(i,1)=(4*T2(i,2)-T2(i,3)+(T2(i+1,1)-T2(i-1,1))*(deta/dxi)*(xiy(i,1)/etay(i,1)))/3;
   end

   c11=3*(xiy(nx,1)/(dxi)-etay(nx,1)/(deta));

  T2(nx,1)=((4*T2(nx-1,1)-T2(nx-2,1))*(xiy(nx,1))/(dxi) + ((-4*T2(nx,2)+T2(nx,3))*etay(nx,1))/(deta))/c11;

%T2(:,1)=Tinf;
u2(:,1,:)=0;
    v2(:,1,:)=0;%wall cha
    %%T2(:,1,:)=Tinf;no wall temp const
    P2(:,1,:)=2*P2(:,2,:)-P2(:,3,:);%wall pressure extrapolation
    %far field u
    u2(1,:,:)=Uinf;
    u2(:,ny,:)=Uinf;
    %farfield v
    v2(1,:,:)=0;
    v2(:,ny,:)=0;
    %pressure field(inlet and farfiel)
    P2(1,:,:)=Pinf;
    P2(:,ny,:)=Pinf;
    %farfield T
    T2(1,:,:)=Tinf;
    T2(:,ny,:)=Tinf;
    %leading edge u
    u2(1,1,:)=0;%cha
    v2(1,1,:)=0;
    P2(1,1,:)=Pinf;
    T2(1,1,:)=Tinf;
    e=cv*T2;
    for i=1:nx
        for j=1:ny
            rho2(i,j)=P2(i,j)/(R*T2(i,j));
        end
    end

    U1=prim2cons(rho2,u2,v2,T2,cv);
    
    residual(n+1)=(norm(squeeze(Uinit(2,:,:))-squeeze(U1(2,:,:))));
    %Uinital=U;
    drawnow;

    % Creating Video

    if n==y+1000
        subplot(3,2,1)
        pcolor(xx,yy1,real(u2));
        colormap("Jet")
        colorbar
        shading interp
        axis equal tight
        xlabel('x')
        ylabel('y')
        title('u-velocity','FontSize',15)
        subplot(3,2,2)
        pcolor(xx,yy1,real(v2));
        colormap("Jet")
        colorbar
        shading interp
        axis equal tight
        xlabel('x')
        ylabel('y')
        title('v-velocity','FontSize',15)
        subplot(3,2,3)
        pcolor(xx,yy1,real(rho2));
        colormap("Jet")
        colorbar
        shading interp
        axis equal tight
        xlabel('x')
        ylabel('y')
        title('Density','FontSize',15)
        subplot(3,2,4)
        pcolor(xx,yy1,real(P2));
        colormap("Jet")
        colorbar
        shading interp
        axis equal tight
        xlabel('x')
        ylabel('y')
        title('Pressure','FontSize',15)
        subplot(3,2,5)
        pcolor(xx,yy1,real(e));
        colormap("Jet")
        colorbar
        shading interp
        axis equal tight
        xlabel('x')
        ylabel('y')
        title('Internal Energy','FontSize',15)
        subplot(3,2,6)
        pcolor(xx,yy1,real(T2));
        colormap("Jet")
        colorbar
        shading interp
        axis equal tight
        xlabel('x')
        ylabel('y')
        title('Temperature','FontSize',15)
        drawnow
        y=n;
        frame=getframe(gcf);
        writeVideo(vid,frame)

    end
    n=n+1;
   end
   machnumber1=zeros(nx,ny);
   for i =1:nx
       for j=1:ny
           machnumber1(i,j)=real(u2(i,j))/(sqrt(gamma*R*real(T2(i,j))));
       end
   end

   figure(2)
   nn=linspace(1,iter,iter);
   plot(nn,residual)
   xlabel('Number Of Time Steps')
   ylabel('U-velocity residual')
   title('Convergence plot of velocity','FontSize',15)

   figure(3)
       pcolor(xx,yy1,machnumber1)
       colormap("Jet")
       colorbar
       shading interp
       axis equal tight
       xlabel('x')
       ylabel('y')
       title('Mach number','FontSize',15)



inde=floor(H/(deta*2.1))+1;
figure(4)
plot(yy1(:,inde),machnumber1(:,inde));
xlabel('x')
ylabel('Machnumber')
title("Mach number variation along x-direction for M="+M,'FontSize',15)

figure(5)
plot(yy1(:,inde),real(P2(:,inde)));
xlabel('x')
ylabel('Pressure')
title("Pressure variation along x-direction for M="+M,'FontSize',15)

figure(6)
plot(yy1(:,inde),real(T2(:,inde)));
xlabel('x')
ylabel('Temperature')
title("Temperature variation along x-direction for M="+M,'FontSize',15)

figure(7)
plot(yy1(:,1),real(T2(:,1)));
xlabel('x')
ylabel('Temperature')
title("Temperature variation along the wall for M="+M,'FontSize',15)
end