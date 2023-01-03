clc;
clear;
close all;
newData1 = load('-mat', 'Data.mat');
vars = fieldnames(newData1);
for i = 1:length(vars)
    assignin('base', vars{i}, newData1.(vars{i}));
end

a=find(time==813);
aa=find(time==814.86);  %Start of Flight
bb=find(time==816.68); 
cc=find(time==819);

input=xlsread("input1.xlsx","Sheet1");

m0dot0f=input(5,1); %kg/sec
s=input(3,1); %m^2
mAC=@(t)input(1,1)-m0dot0f*t;  %kg
g=9.81; %m/s^2
b=input(2,1); %m
MAC=input(4,1);   %m
 %kg/hr

phiB=(input(16,1))*pi/180;             % Incidence angle of Booster  (rad)
phiT=input(15,1)*pi/180;
phi_T=(phiT+phiB);                                   % Incidence angle of Booster's phiT=0;

 %kg*m^2
del0e=input(25,1)*pi/180;   %rad
Th0=input(9,1);
Th=@(t)Th0;%round((Th0.*t./0.05).*heaviside(t)+(Th0-0.07*Th0*t-(Th0.*t./0.05)).*heaviside(t-0.05)+(-Th0).*heaviside(t-1.82),2);  %N
AR=b^2/s;   %Airplane Aspect Ratio

osw=1.78*(1-0.045*AR^0.68)-0.64;
pho=1.225;   %kg/m^3

mB0=input(6,1);                  % mass of Booster and fuel - Empty (kg)
mB_fuel=input(7,1);              % mass of Booster's fuel (kg)
mdot_fB=input(8,1);              % rate of burning fuel (kg/s)

lB=input(13,1);                  % lenght of Booster  (m)

%position nose of booster from nose of Aircraft
CGB=[input(11,1); input(11,2)];   % C.G of Booster (m) from JATO nose

      
TR_Boost=[input(13,1); input(13,2)];

TR=[input(14,1); input(14,2)];



 CGa=[input(10,1);input(10,2)];
mB=@(t)(mB0+mB_fuel-mdot_fB.*t).*heaviside(t)+(-mB_fuel+mdot_fB.*t).*heaviside(t-input(26,1));
m=@(t)mAC(t)+mB(t);

lbcg=@(t)(input(27,1)*mB0+(27/1000+input(28,1)/2)*(mB_fuel-mdot_fB.*t))./mB(t);
xbcg=@(t)(input(29,1)+lbcg(t)*cosd(phi_T));
zbcg=@(t)(input(29,2)+lbcg(t)*sind(phi_T));

xcg=@(t)(xbcg(t)*mB(t)+mAC(t)*CGa(1,1))./m(t);
zcg=@(t)((zbcg(t)*mB(t)+mAC(t)*CGa(2,1))./m(t));

A=+tan(180-phiT); B=-1;
C=-tan(180-phiT)*TR_Boost(1,1)+TR_Boost(2,1);

TL1=@(t)(-xcg(t)+TR(1,1));
TL2=@(t)(-zcg(t)+TR(2,1));
TA0=@(t)atan(TL2(t)/TL1(t));
r0=@(t)sqrt(TL1(t)^2+TL2(t)^2);
 
TV1=@(t)(-xcg(t)+TR_Boost(1,1)); 
TV2=@(t)(-zcg(t)+TR_Boost(2,1));

d2=@(t)(A*-xcg(t)+B*zcg(t)+C)/sqrt(A^2+B^2);

r=@(t)sqrt(TV1(t)^2+TV2(t)^2);
TA=@(t)atan(TV2(t)/TV1(t));

%iT=@(t)phiT*pi/180-TA(t);
d1=@(t)TV1(t)*sin(phiT)+TV2(t)*cos(phiT);



cL=@(x)input(18,1)+input(19,1)*(atan(x(2)./x(1)))+input(20,1)*input(25,1);  %cL of Aircraft cL=@(x)0.0314+4.659*(atan(x(2)./x(1)))+0.2189*del0e;  %cL=@(x)0.06+5.01*atan(x(2)/x(1))+0.2189*del0e; 
cD=@(x)input(21,1)+1/(pi()*osw*AR)*cL(x)^2; %cD of Aircraft
cM=@(x)input(22,1)+input(23,1)*(atan(x(2)./x(1)))+input(24,1)*input(25,1); %cM of Aircraft cM=@(x)0.0703-0.903*(atan(x(2)./x(1)))-0.6572*del0e; cM=@(x)0-0.954*atan(x(2)/x(1))-0.6572*del0e;
alpha=@(x)atan(x(2)./x(1));

qD=@(x)0.5*pho*(x(1).^2+x(2).^2);  %Dynamic Preassure
Iyy=@(t)input(12,1)*m(t)/m(0); %Inertia Moment

 %total mass of Aircraft Plus Bosster and Fuel %320
 

Xdot=@(t,x)[-x(2).*x(3)-g.*sin(x(4))+(Th(t)./m(t)).*cos(phiT)+input(9,2)./m(t)+qD(x).*s./m(t).*(-cD(x)*cos(alpha(x))+cL(x).*sin(alpha(x)));  %794

x(1).*x(3)+g.*cos(x(4))-(Th(t)./m(t)).*sin(phiT)+qD(x).*s./m(t).*(-cD(x).*sin(alpha(x))-cL(x).*cos(alpha(x)));

qD(x).*s.*MAC.*cM(x)./Iyy(t)+Th(t)*r(t)*sin(phiT-TA(t))./Iyy(t)-input(9,2)*r0(t)*sin(TA0(t))./Iyy(t);     %794

x(3);

x(1).*cos(x(4))+x(2).*sin(x(4));

-x(1).*sin(x(4))+x(2).*cos(x(4));];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




x0=[0.000001,0,0,pi*input(17,1)/180,0,0];

[T,Y]=ode45(Xdot,[0,input(26,1)],x0);
y=Y(:,1);

figure("Name","u");
plot(T,Y(:,1));
hold on
%plot(T1,Y1(:,1));

figure("Name","w");
plot(T,Y(:,2));
hold on
%plot(T1,Y1(:,2));

figure("Name","q");
plot(T,180*Y(:,3)/pi);
hold on
plot(time(aa:bb)-814.86,Q(aa:bb),"--");

%plot(T1,Y1(:,3));


figure("Name","Theta");
plot(T,Y(:,4)*180/pi);
hold on
plot(time(aa:bb)-814.86,theta(aa:bb)+1.32,"--");
%plot(T1,Y1(:,4)*180/pi);

figure("Name","X");
plot(T,Y(:,5));
hold on
%plot(T1,Y1(:,5));
%plot(time(aa:bb)-814.86,Range(aa:bb)-11.751,"--");

figure("Name","Z");
plot(T,-Y(:,6));

hold on

%plot(T1,Y1(:,6));
%plot(time(aa:bb)-814.86,BAR_height(aa:bb)-650.577,"--");
%plot(time(aa:bb)-814.86,Height(aa:bb)-650.577-169.656,"--");

figure("Name","X-Z");
plot(Y(:,5),-Y(:,6));
hold on

%plot(Y1(:,5),-Y1(:,6));
%plot(Range(bb:cc)-11.751,BAR_height(bb:cc)-650.577,"--");

cL0=0.06+5.01*atan(Y(:,2)./Y(:,1))+0.2189*del0e; %cL of Aircraft
cD0=0.034+1/(pi()*osw*AR)*cL0.^2; %cD of Aircraft
qD0=0.5*pho*(Y(:,1).^2+Y(:,2).^2);
Th0=Th(T); m0=m(T);
wdot=Y(:,1).*Y(:,3)+g.*cos(Y(:,4))-(Th(T)./m(T)).*sin(phiB)+qD0.*s./m(T).*(-cD0.*sin(atan(Y(:,2)./Y(:,1)))-cL0.*cos(atan(Y(:,2)./Y(:,1))));


xdot=Y(:,1).*cos(Y(:,4))+Y(:,2).*sin(Y(:,4));
zdot=-Y(:,1).*sin(Y(:,4))+Y(:,2).*cos(Y(:,4));



figure("Name","xdot");
plot(T,xdot);  % velocity of X
hold on
%plot(T1,x1dot);
%plot(time(bb:cc)-815,NORTH_VEL(bb:cc),"--");

figure("Name","zdot");
plot(T,-zdot);
hold on

%plot(T1,z1dot);
%plot(time(aa:bb)-814.86,Ver_Speed(aa:bb),"--")

figure("Name","V");
plot(T,sqrt(xdot.^2+zdot.^2));
hold on
%plot(T1,sqrt(x1dot.^2+z1dot.^2));
%plot(time(bb:cc)-815,sqrt(Ver_Speed(bb:cc).^2+NORTH_VEL(bb:cc).^2),"--")


V1=sqrt(xdot.^2+zdot.^2);
%cf=fit(T,sqrt(xdot.^2+zdot.^2),"rat55");
%V=@(t)(cf.p1*t^5 + cf.p2*t^4 + cf.p3*t^3 + cf.p4*t^2 + cf.p5*t + cf.p6) / (t^5 + cf.q1*t^4 + cf.q2*t^3 + cf.q3*t^2 + cf.q4*t + cf.q5);

% figure
% Z1=Y(:,6);
% cf=fit(T,Z1,"rat55");
% syms t;
% Z=(cf.p1*t^5 + cf.p2*t^4 + cf.p3*t^3 + cf.p4*t^2 + cf.p5*t + cf.p6) / (t^5 + cf.q1*t^4 + cf.q2*t^3 + cf.q3*t^2 + cf.q4*t + cf.q5);
% Z=diff(Z);




% opt=ga(@optfun,2,[],[],[],[],[-1.6 330*9.81],[-1.2 400*9.81]);
% 
% XboosterFromNose=opt(1)
% Thrust=opt(2)/9.81


