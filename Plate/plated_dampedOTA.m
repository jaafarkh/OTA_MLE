% transverse response of rectangular plate with length a, width b and hight h
% using analytical solution of Vibration of Continuous Systems(S. S. Rao)
% see sec. 14.5 and example 14.1, p489-490
% assuming simply supported boundary conditions for all edges
% assuming point force applied at point (x0,y0) 
% assume zero initial conditions
%this will generate time time for order tracking OBMA
close all
clear
clc
%plate dimensions; a: width, b: height, h: thickness
a = 0.5;
b = 0.4;
h = 0.001;

%material properties
rho = 7850;%2770; 
nu = 0.3; %0.33;
E = 200e9; %71e9;

D = E*h^3/(12*(1-nu^2)); %flextural rigidity
A1mn = 2/sqrt(rho*h*a*b); %normal modes normalized amplitude
w1 = pi^2*sqrt(D/(rho*h)); %OMEGAmn pre-product

% F0: the amplitude of the applied constant force 
% (x0,y0): position of the force
% (x1,y1): point at which the response will be calculated
F0 = 160; %force amplitude
x0=0.125;
y0=0.075;
MaxmMod = 4; %maximum number of modes in x-direction (m)
MaxnMod = 3; %maximum number of modes in y-direction (n)
            %total modes = MaxMod * MaxMod

fs = 1280;
TotalTime = 40;
Np = TotalTime*fs;
Tstep = 1/fs;
W=@(x,y,m,n) A1mn*sin(m*pi*x/a)*sin(n*pi*y/b);
tt=Tstep*(0:Np-1)';
t0=0.009397; %duration of impulse force
ft = zeros(Np,1); %force as function of t
tachodat = zeros(length(tt),1); %angle data
ft(1:t0/Tstep) = F0;
freqHz=5;
theta=0;
for i=1:Np
    %theta = 2*pi*(5.0*tt(i)+55.0*tt(i)*tt(i)/(2*TotalTime));%direct integration
    ft(i) = 4*sin(theta)+4*sin(4*theta)+12*sin(2*pi*40*tt(i));%4,4,12 generate chirp excitation
    w = freqHz*2*pi;
    tachodat(i) = 5*sin(theta);
    theta = theta + w*Tstep;
    freqHz = 5 + 55 *(i/Np)^2;%70
end
plot(ft)
res=zeros(Np,1);
% return
Res2=zeros(length(tt),13);
zetamn = zeros (MaxmMod,MaxnMod);

alphad = 20; %Rayliegh damping coefficients
betad= 0.000;
disp('Calculating the natural frequencies')
freqmat=zeros(MaxmMod*MaxnMod,1);
i=1;
for  m = 1:MaxmMod        
        for n = 1:MaxnMod
           wmn = w1*((m/a)^2+(n/b)^2);
           zetamn(m,n) = alphad/(2*wmn)+betad*wmn/2; %Rayliegh damping
           
           fmn = wmn/(2*pi);
           fprintf('m = %d, n = %d, fmn = %8.4f Hz\n',m,n, fmn);
           freqmat(i) = fmn;
           i=i+1;
        end
end
freqmat=sort(freqmat) %sorted frequencies

i=1;
for y1=0.1:0.1:0.3
    for x1=0.1:0.1:0.4
        DuhAmn = zeros(MaxmMod,MaxnMod);
        DuhBmn = zeros(MaxmMod,MaxnMod);
        for k=1:Np
            t=tt(k);
            r=0;
            ri=0;
            for  m = 1:MaxmMod
                for n = 1:MaxnMod
                    wmn = w1*((m/a)^2+(n/b)^2);
                    wd = wmn*sqrt(1-zetamn(m,n)^2);
                    etamn = (1/wmn)*W(x0,y0,m,n) * exp(-zetamn(m,n)*wmn*t)*(DuhAmn(m,n)*sin(wd*t)-DuhBmn(m,n)*cos(wd*t));
                    % etamn is the response of the decoupled equations
                    r = r + W(x1,y1,m,n) * etamn; %sum up the decoupled responses to find the overall response
                    ri = ri + 0.01*sin(m*pi*x1/a)*sin(n*pi*y1/b)*cos(wmn*t); %due to initial conditions only
                    DuhAmn(m,n) = DuhAmn(m,n) + Tstep * ft(k)*exp(zetamn(m,n)*wmn*t)*cos(wd*t);
                    DuhBmn(m,n) = DuhBmn(m,n) + Tstep * ft(k)*exp(zetamn(m,n)*wmn*t)*sin(wd*t);
                end

            end
            res(k) = r;
        end
        figure();
        plot(tt,res);
        Res2(:,i) = res;
        disp([i x1 y1]);
        i = i+1;
    end
end
figure();
plot(tt,tachodat)
Res2(:,13) = tachodat;
pathFolderResults = [pwd,'\otaP\input\timeDataPlateOTA2.txt'];
%pathFolderResults = [pwd,'\otaP\input\timeDataPlateOTANoF.txt'];

writematrix(Res2,pathFolderResults,'Delimiter','tab')
