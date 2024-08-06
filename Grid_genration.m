nx=100;
ny=100;
L=1.5*10^-5;
H=1.3*10^-5;
x=linspace(0,L,nx);
y=linspace(0,H,ny);

[xx,yy]=ndgrid(x,y);
yy1=yy;
yy2=yy;
xx1=xx;


n=0.265;
j=1;

while n>0

    yy1(:,j)=n*xx(:,j)+y(j)-n*(x(1));
    n=n-0.0040;
    j=j+1;
end

nn=-0.2;

    xx1(2:end-1,1:3)=nn*yy1(2:end-1,1:3)+xx(2:end-1,1)-nn*yy1(2:end-1,1);

for i=2:nx
    m=((xx1(i,3)-xx1(i,ny))/(yy1(i,3)-yy1(i,ny)));
    xx1(i,3:ny)=m*yy1(i,3:ny)+(xx1(i,3))-m*yy1(i,3);
end

%{
coeff=0.073;
fst=(yy1(60,1));
lst=yy1(60,80);
b=fst-coeff;
alp=log((lst-b)/(coeff*lst));





yy1(61:80,:)=coeff*exp(alp*yy1(61:80,:))+b;
yy2(1:70,1)=0.4*xx(1:70,1)+y(1)-0.4*(x(1));

%}

%{
while n>0

    yy1(:,j)=n*xx(:,j)+y(j)-n*(x(1));
    n=n-0.0056;
    j=j+1;
end
%}
%{
nn=0.6;
l=1;
while nn>0
    yy1(20:45,l)=nn*xx(20:45,l)+y(l)-nn*(x(14));
    nn=nn-0.01;
    l=l+1;
end

coeff=2.75;
fst=(yy1(45,1));
disp(fst)
disp(yy1(20,1));
lst=yy1(45,80);
disp(lst)
alp=0.25374219;
bb=-2.544308586;

yy1(46:80,:)=coeff*exp(alp*yy1(46:80,:))+bb;
yy2(1:70,1)=0.4*xx(1:70,1)+y(1)-0.4*(x(1));

%}

u=zeros(nx,ny);
figure(1)
pcolor(xx1,yy1,u);
u=zeros(nx,ny);
figure(2)
pcolor(xx,yy1,u);

figure(3)
hold on
for n = 1:numel(x) %// loop over vertical lines
    plot([xx1(n,1) xx1(n,end)], [yy1(n,1) yy1(n,end)], 'k-'); %// change 'k-' to whatever you need
end
for n = 1:numel(y)%// loop over horizontal lines
    plot([xx1(1,n) xx1(end,n)], [yy1(1,n) yy1(end,n)], 'k-'); %// change 'k-' to whatever you need
end

figure(4)
hold on
for n = 1:numel(x) %// loop over vertical lines
    plot([xx(n,1) xx(n,end)], [yy1(n,1) yy1(n,end)], 'k-'); %// change 'k-' to whatever you need
end
for n = 1:numel(y)%// loop over horizontal lines
    plot([xx(1,n) xx(end,n)], [yy1(1,n) yy1(end,n)], 'k-'); %// change 'k-' to whatever you need
end