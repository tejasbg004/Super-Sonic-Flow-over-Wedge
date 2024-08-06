%Numerical Schelerin Images.
kap=10;
Beta=0.8;
addpath('Derivative_Functions')
drhodx=ddx_fwd_xi(rho,dxi,deta,yy1,xx1);
drhody=ddy_fwd_eta(rho,dxi,deta,yy1,xx1);
drho=sqrt(drhodx.*drhodx+drhody.*drhody);
rhomax=max(drho,[],"all");
Sch1=Beta*exp(-kap*drho/rhomax);

% dUdx=ddx_fwd(u,dx);
% dUdy=ddy_fwd(u,dy);
% dU=sqrt(dUdx.*dUdx+dUdy.*dUdy);
% Umax=max(dU,[],"all");
% Sch2=Beta*exp(-kap*dU/Umax);
% 
% dVdx=ddx_fwd(v,dx);
% dVdy=ddy_fwd(v,dy);
% dV=sqrt(dVdx.*dVdx+dVdy.*dVdy);
% Vmax=max(dV,[],"all");
% Sch3=Beta*exp(-kap*dV/Vmax);
% 
% dPdx=ddx_fwd(P,dx);
% dPdy=ddy_fwd(P,dy);
% dP=sqrt(dPdx.*dPdx+dPdy.*dPdy);
% Pmax=max(dP,[],"all");
% Sch4=Beta*exp(-kap*dP/Pmax);
% 
% dedx=ddx_fwd(e,dx);
% dedy=ddy_fwd(e,dy);
% de=sqrt(dedx.*dedx+dedy.*dedy);
% emax=max(de,[],"all");
% Sch5=Beta*exp(-kap*de/emax);
% 
% dTdx=ddx_fwd(T,dx);
% dTdy=ddy_fwd(T,dy);
% dT=sqrt(dTdx.*dTdx+dTdy.*dTdy);
% Tmax=max(dT,[],"all");
% Sch6=Beta*exp(-kap*dT/Tmax);

subplot(1,1,1)
pcolor(xx1,yy1,Sch1)
hold on
plot([0,1.48413e-5],[0,1.07677e-5],color='red',LineWidth=3)
xlabel('X')
ylabel('Y')
title('Schlieren Image of Density ','FontSize',15)
colormap('gray')
shading interp
axis equal tight

% subplot(3,2,2)
% pcolor(xx,yy,Sch2)
% xlabel('X')
% ylabel('Y')
% title('Schlieren Image of U-velocity ','FontSize',15)
% colormap('gray')
% shading interp
% axis equal tight
% 
% subplot(3,2,3)
% pcolor(xx,yy,Sch3)
% xlabel('X')
% ylabel('Y')
% title('Schlieren Image of V-velocity ','FontSize',15)
% colormap('gray')
% shading interp
% axis equal tight
% 
% subplot(3,2,4)
% pcolor(xx,yy,Sch4)
% xlabel('X')
% ylabel('Y')
% title('Schlieren Image of Pressure','FontSize',15)
% colormap('gray')
% shading interp
% axis equal tight
% 
% subplot(3,2,5)
% pcolor(xx,yy,Sch5)
% xlabel('X')
% ylabel('Y')
% title('Schlieren Image of Internal Energy ','FontSize',15)
% colormap('gray')
% shading interp
% axis equal tight
% 
% subplot(3,2,6)
% pcolor(xx,yy,Sch6)
% xlabel('X')
% ylabel('Y')
% title('Schlieren Image of Temperature ','FontSize',15)
% colormap('gray')
% shading interp
% axis equal tight