%Comparison for two cases.
Tinf=288.15;
Pinf=101300;
x=linspace(0,0.00005,100);
y=linspace(0,0.00005,100);
figure(1)
hold on 
   plot(P2(20,:,:)/Pinf,y);
   plot(P(20,:,:)/Pinf,y);
   xlabel('P/Pinf at x/l=0.25')
   ylabel('y-Distance')
   title('Comparison of Normalized Pressure between Adiabatic and Constant wall Temperature at x/L=0.25 ','FontSize',15)
   legend('Adiabatic','Constant wall Temp')
hold off

figure(2)
    hold on
   plot(T2(20,:,:)/Tinf,y);
   plot(T(20,:,:)/Tinf,y);
   xlabel('T/tinf at x/l=0.25')
   ylabel('y-Distance')
   title('Comparison of Normalized Temperature between Adiabatic and Constant wall Temperature at x/L=0.25 ','FontSize',15)
   legend('Adiabatic','Constant wall Temp')
  hold off

  figure(3)
hold on 
   plot(P2(40,:,:)/Pinf,y);
   plot(P(40,:,:)/Pinf,y);
   xlabel('P/Pinf at x/l=0.5')
   ylabel('y-Distance')
   title('Comparison of Normalized Pressure between Adiabatic and Constant wall Temperature at x/L=0.5 ','FontSize',15)
   legend('Adiabatic','Constant wall Temp')
hold off

figure(4)
hold on 
   plot(T2(40,:,:)/Tinf,y);
   plot(T(40,:,:)/Tinf,y);
   xlabel('T/tinf at x/l=0.5')
   ylabel('y-Distance')
   title('Comparison of Normalized Temperature between Adiabatic and Constant wall Temperature at x/L=0.5 ','FontSize',15)
   legend('Adiabatic','Constant wall Temp')
  hold off

  figure(5)
hold on 
   plot(P2(60,:,:)/Pinf,y);
   plot(P(60,:,:)/Pinf,y);
   xlabel('P/Pinf at x/l=0.75')
   ylabel('y-Distance')
   title('Comparison of Normalized Pressure between Adiabatic and Constant wall Temperature at x/L=0.75 ','FontSize',15)
   legend('Adiabatic','Constant wall Temp')
hold off

figure(6)
hold on 
   plot(T2(60,:,:)/Tinf,y);
   plot(T(60,:,:)/Tinf,y);
   xlabel('T/tinf at x/l=0.75')
   ylabel('y-Distance')
   title('Comparison of Normalized Temperature between Adiabatic and Constant wall Temperature at x/L=0.75 ','FontSize',15)
   legend('Adiabatic','Constant wall Temp')
  hold off

  
  figure(7)
  hold on 
  plot(x,T2(:,1,:))
  plot(x,T(:,1,:))
  xlabel('x')
  ylabel('T(x)')
  title('Comparison of wall Temperature between Adiabatic and Constant wall Temperature as a function of x ','FontSize',15)
  legend('Adiabatic','Constant wall Temp')
  hold off