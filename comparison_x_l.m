% Plots of P2 T2 at X/L=0.25 X/L= 0.25,0.75 and 0.75

Tinf=288.15;
Pinf=101300;
L=1.5*10^-5;
H=1.3*10^-5;
x=linspace(0,L,nx);
y=linspace(0,H,ny);
figure(1)
hold on 
   plot(P1(20,:,:)/Pinf,y/H);
   plot(P1(40,:,:)/Pinf,y/H);
   plot(P1(60,:,:)/Pinf,y/H);
   plot(P1(80,:,:)/Pinf,y/H);
   plot(P01(20,:,:)/Pinf,y/H,'--');
   plot(P01(40,:,:)/Pinf,y/H,'--');
   plot(P01(60,:,:)/Pinf,y/H,'--');
   plot(P01(80,:,:)/Pinf,y/H,'--');
   xlabel('P/Pinf')
   ylabel('y-Distance')
   title('Comparison of Normalized Pressure between  various x/L location','FontSize',15)
   legend('x/L=0.2','x/L=0.4','x/L=0.6','x/L=0.8','adiabatic x/L=0.2','adiabatic x/L=0.4','adiabatic x/L=0.6','adiabatic x/L=0.8')
hold off

figure(2)
hold on 
   plot(T1(20,:,:)/Tinf,y/H);
   plot(T1(40,:,:)/Tinf,y/H);
   plot(T1(60,:,:)/Tinf,y/H);
   plot(T1(80,:,:)/Tinf,y/H);
   plot(T01(20,:,:)/Tinf,y/H,'--');
   plot(T01(40,:,:)/Tinf,y/H,'--');
   plot(T01(60,:,:)/Tinf,y/H,'--');
   plot(T01(80,:,:)/Tinf,y/H,'--');
   xlabel('T/Tinf')
   ylabel('y-Distance')
   title('Comparison of Normalized Temperature between  various x/L location','FontSize',15)
   legend('x/L=0.2','x/L=0.4','x/L=0.6','x/L=0.8','adiabatic x/L=0.2','adiabatic x/L=0.4','adiabatic x/L=0.6','adiabatic x/L=0.8')
hold off
