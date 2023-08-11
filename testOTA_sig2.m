%this code applies OTA for simulated data or data obtained from simulation
%of 2-DOF system
% Case-2: 5 to 60 Hz, duration: 30 sec, fs = 1280, bnr= rev
clear;
close all;
addpath("otaP\functions\");
caseNo = 2; %signal case for simulation, 1: constant-amp, 2: two-DOF
thetSrc = 0; %angular position source, 0: actual, 1: corrected by interpolation of speed, 2: uncorrected; 
resSrc = 0; %response source, 0: simulation, 
opt = 5; %1:constant DT ML, 2: constant Rev ML, 3: VKF, 4:: const DT Blough TVDFT, 5: const Rev Blough TVDFT
%For constant amp vibration, ML is always better
%for 2-DOF vibration
%constant DT with bsize 150 or less, ML is better than Blough especially
%cases: Constant DT = 150 point ML is better with theta correction whether
%using actual theta or interpolated theta
%constant Rev=4, ML is better whether using uncorrected theta or
%interplolated theta
%6: VK with overlapping 
rm=[150;150;190]; %VKF weighting factor, 
forder=0;  %VK filter order: 0,1,2, No. of poles = filter order + 1
corrEnable = 1; %1: enable omega correction using GD, 0: disabled
OCMEnabled = 1; % 0: OCM is diabled, 1: OCM is enabled in TVDFT
fs = 1280;
dwnsamp = 1;%downsampling factor of vibration data, 1: no downsampling
bsize = 150;%processing length for constant DT in samples
overlap = 0;%overlapping for constant DT
bnr = 4; %processing block in rev
nroverlp=0;

%define orders; 
% maximum order frequency MUST NOT exceed half sampling rate
% Orders must be in ascending order for Run-up test and descending order
% for coast-down test
% constant frequency items can be put in any location, better at the
% end of orders
% OrderType: 0: Order, 1: constant frequency
NTraces = 3;
TraceOrder= zeros(NTraces,1);
TraceType=zeros(NTraces,1);

TraceOrder(1)=1;
TraceType(1)=0;

TraceOrder(2)=2;
TraceType(2)=0;

TraceOrder(3)=40;
TraceType(3)=1;

TrigLev1=0;
HystV = 0.2;

TotalTime = 30;
Np = TotalTime*fs;
Res2 = zeros(Np,2);
ResCF = zeros(Np,1);
ft = zeros(Np,1);
reconssig = zeros(Np,1);

dt = 1/fs;
ts=dt*(0:Np-1)';
freqHz=5;
thetaM = zeros(Np,1);
f1 = zeros(Np,1);
f3 = TraceOrder(3)*ones(Np,1);
ang=0;
%% define the simulation for two-DOF system
K = 3.0*[30000 -16000;-16000 26000];
M =[1.2 0; 0 1.0];
C=[36 0; 0 30];
[phi,D]=eig(K,M); % modes shapes and eigenvalues, 
%               Note: matlab provides mass normalized mode shapes
wn=sqrt(D);
z1=(C(1,1)/M(1,1))/(2*wn(1,1));%alpha/2w
z2=(C(2,2)/M(2,2))/(2*wn(2,2));
wr=[wn(1,1) wn(2,2)];
zr=[z1 z2];
alphar=@(i,j,r,w) 1000*phi(i,r)*phi(j,r)/(wr(r)*wr(r)-w*w+2i*zr(r)*w*wr(r));
%% start signal generation
acta = zeros(Np,3);
actb = zeros(Np,3);
actamp = zeros(Np,3);
% f1=5+(55/(Np-1))*(0:Np-1)';
for i=1:Np
    ang = 2*pi*(5.0*ts(i)+55.0*ts(i)*ts(i)/(2*TotalTime));%direct integration of omega
    w = freqHz*2*pi;
    if caseNo == 1
        for j=1:3
         acta(i,j) = 0;
         actb(i,j) = 4;
        end
    else
      a11_1X=alphar(1,1,1,TraceOrder(1)*w)+alphar(1,1,2,TraceOrder(1)*w);
      acta(i,1)=100*real(a11_1X);
      actb(i,1)=100*imag(a11_1X); 
      a11_4X=alphar(1,1,1,TraceOrder(2)*w)+alphar(1,1,2,TraceOrder(2)*w);
      acta(i,2)=100*real(a11_4X);
      actb(i,2)=100*imag(a11_4X);       
      acta(i,3)=0;
      actb(i,3)=6; 
    end
    for j=1:3
         actamp(i,j) = sqrt(acta(i,j)^2+actb(i,j)^2);
    end
    Res2(i,1) = acta(i,1)*cos(TraceOrder(1)*ang) + actb(i,1)*sin(TraceOrder(1)*ang) ...
        + acta(i,2)*cos(TraceOrder(2)*ang) + actb(i,2)*sin(TraceOrder(2)*ang) ...
        + acta(i,3)*cos(2*pi*TraceOrder(3)*ts(i)) + actb(i,3)*sin(2*pi*TraceOrder(3)*ts(i));

    ResCF(i,1) = 150*sin(2*pi*TraceOrder(3)*ts(i));
    f1(i) = freqHz;
    
    Res2(i,2) = 5*sin(ang);
    thetaM(i) = ang;
%     ang = ang + w*dt;%digital integration of Omega
    freqHz = 5 + 55 *i/(Np-1); % = 5 + 55 t/T
end
% thetaM=2*pi*cumsum(f1,1)*dt;
% plot(env1x)
% return
f2 = f1*TraceOrder(2);
r = Res2(:,2); %last column is tacho signal

L = size(Res2,1);
firstdetected = 0;
theta = zeros(L,1);%extracted from tacho signal
thetaJ = zeros(L,1);%extracted from poly fitting
cpos = zeros(L,1);
OMG = zeros(L,1);
OMG0 = zeros(L,1); %from unfiltered tacho signal
pindex0=1;

nr=0;
Jterm=0;
stp=0; %starting or first edge index
endp=L; %last edge index
if r(1) > TrigLev1 
    PState=1;
else
    PState=0;
end
for i=2:L
    if r(i) > (TrigLev1 + HystV)
        if PState == 0
            if firstdetected > 0
                Count1 = i - pindex0 ;
                endp = i;
                nr = nr+1; %number of revoultions
                cpos(nr) = i;

                t2=  Count1;
                if nr > 1
                    w1= 2*pi*fs / t1; %assuming linear speed variation, calculate the slope
                    w2=2*pi*fs / t2;
                    a1 = (w2-w1)/t2;%dw/dt
                else
                    w1 = 2 * pi *fs / Count1;
                    w2 = w1;
                    a1 = 0;
                end
                tt=0;
                for j = 1:Count1
                    w = w1 + a1*j; %instantaneous speed
                    tt = tt + w/fs; %to ensure complete 2pi per cycle
                end
                Jterm = w1 - 2 * pi *fs / (Count1 + 2*randn(1,1));
                for j = 1:Count1
                    w = w1 + a1*j; %instantaneous speed
                   % OMG0(j+pindex0) = (w1+w2)/2; %good estimate
                    OMG0(j+pindex0) = w1+(w2-w1)*j/Count1; %good estimate
                    theta(j+pindex0) = theta(j-1+pindex0) + (w*dt)*2*pi/tt;
                   
                   % theta(j+pindex0) = theta(j-1+pindex0) + (w2*dt)+0.000002*i/L;
                end
                psi0 = theta(Count1+pindex0);
                
                t1= Count1;
                pindex0 = i;
            else
                stp = i;
                firstdetected = 1;
                pindex0 = i;
                theta(pindex0) = 0.05;
               
                nr = 0;                
                t1 = 0 ;
            end
        end
        PState = 1;
    end
    if r(i) < (TrigLev1 - HystV)
        PState = 0;
    end
end
rf=zeros(1,L);
rev=2;
for i = 1:L
    if i <= cpos(1)
        r(i)= 0;%r(i)/5;
    end
    if cpos(rev) == i %update filter parameters each rev
        i0 = cpos(rev-1);       
       
        ci=cpos(rev)-i0;
        for j=1:ci
            r(i0+j-1)=sin((j-1)*2*pi/ci); %this gives better results
        end
        rev = rev+1;
    end
end


fncr= zeros(nr+12,1);
fncrm= zeros(nr+12,1);
cposm= zeros(nr+12,1);
for i=2:nr
    fncr(i) = fs / (cpos(i) - cpos(i-1));%just to show the accuracy of the speed signal extracted      
end
for i=2:nr-1
    fncrm(i) = (fncr(i) + fncr(i+1))/2;       
end
for i=2:nr   
    cposm(i)=(cpos(i)+cpos(i-1))/2;
end
fncrm(nr)=fncr(nr);
figure()
%plot(cpos(2:rev),fncr(2:rev),'r');
hold on
%plot(cpos2(2:rev),fcr(2:rev),'b')
hold on
%yy = csaps(cpos(2:rev),fncr(2:rev),0,cpos(2:rev));
fspd=fit(dt*cpos(2:nr),fncrm(2:nr),'poly3','Robust','Bisquare') ;%poly3, try smoothingspline ,'Robust','Bisquare' ,'LAR'
%fspd=fit(dt*cposm(2:nr),fncr(2:nr),'poly3','Robust','Bisquare') ;%poly3, try smoothingspline ,'Robust','Bisquare' ,'LAR'
yy=fspd(dt*(0:L-1));
plot(dt*(0:L-1),yy)
hold on
plot(dt*(0:L-1),f1);
thetaJ=2*pi*cumsum(yy,1)*dt;
thetaJ=thetaJ-thetaJ(1);



fprintf('the number of revolutions:  %d \n',nr);
%  return
figure
theterr=thetaM-thetaJ;
plot(theterr(1:Np-1000))
    %return
if(thetSrc == 0)
    theta2=thetaM;%use actual
end
if(thetSrc == 1)
    theta2=thetaJ;%use interpolated
end
if(thetSrc == 2)
    theta2=theta;% use un-corrected theta
    %theta2= theta2+0.03*randn(Np,1);
end
%   return
%% begin OTA
No = 1;
if opt == 1 || opt == 4 % constant DT
       
    nblocks = 1 + fix(L/bsize);%initial number of blocks
    %FRF=zeros(nblocks,No,1); 
    FRFOrder = zeros(Np,NTraces);
    freqBandOrder = zeros(nblocks,NTraces);
    blktime=zeros(nblocks,1);
    nblocks = 0;
    st =  cpos(2);%starting point 
    Nst=st;
    while (st+bsize) <= cpos(nr)
        stp = st;
        endp = st+bsize-1;
        nblocks = nblocks + 1;           
        
        for i =1:No         
            x=Res2(stp:dwnsamp:endp,i);
            xtheta=theta2(stp:dwnsamp:endp);                       
           
            speedCorr=corrEnable;
            %[OA,OB] = itvdft(x,size(x,1),TraceType,TraceOrder,NTraces,xtheta,stp/dwnsamp,dwnsamp/fs,OCMEnabled);
            if opt == 1
                [OA,OB] = otamlkh(x,size(x,1),TraceType,TraceOrder,NTraces,xtheta,stp/dwnsamp,dwnsamp/fs,speedCorr);
            else
                [OA,OB] = blougtvdft(x,size(x,1),TraceType,TraceOrder,NTraces,xtheta,stp/dwnsamp,dwnsamp/fs,OCMEnabled);
            end
            
            %
            k=1;
            for j=1:NTraces
                if TraceType(j) >= 0 
                    FRFOrder((stp:endp),k) = OA(j) + 1i*OB(j);% Cn =0.5(an - i bn)  
                   k = k+1;
                end
            end
            
        end
        %evaluating frequency bands for tracked orders
        k=1;
        for j=1:NTraces
            if TraceType(j) == 0 %include only orders, discard const freq orders                
                freqBandOrder(nblocks,k) = TraceOrder(j)*(1.8*OMG(stp) + 0.20* OMG(endp))/(4*pi);
                k = k+1;
            end
        end
        st = st + bsize-overlap;
    end
    Nend = endp; %trim edges
    k=1;
    lastFreq = 0;
    freqBandOrder =  freqBandOrder(1:nblocks,:); %trim to include the processed blocks only
%     FRFOrder = FRFOrder(1:nblocks,:,:);
     blktime = dt*(0:Np-1);
     figure()
  
     plot(blktime(Nst:bsize:Nend),abs(FRFOrder(Nst:bsize:Nend,1)));ylim([0 max(abs(FRFOrder(:,1,1)))])
     hold on
     plot(blktime(Nst:bsize:Nend),actamp(Nst:bsize:Nend,1));ylim([0 1.1*max(abs(actamp(Nst:Nend,1)))])
     figure()
     plot(blktime(Nst:bsize:Nend),abs(FRFOrder(Nst:bsize:Nend,2)));ylim([0 max(abs(FRFOrder(:,2)))])
     hold on
     plot(blktime(Nst:bsize:Nend),actamp(Nst:bsize:Nend,2));ylim([0 1.1*max(abs(actamp(Nst:Nend,2)))])
     if(NTraces > 2)
       figure()
       plot(blktime(Nst:bsize:Nend),abs(FRFOrder(Nst:bsize:Nend,3)));ylim([0 max(abs(FRFOrder(:,3)))])
     end
     figure
%      bsize=2*bsize;
     endplot=fix(Nend/bsize);
     plot(blktime(Nst:bsize:Nend),abs(FRFOrder(Nst:bsize:Nend,1)),'-bs','LineWidth',0.5,'MarkerEdgeColor','b', 'MarkerSize',3,'MarkerIndices',1:6:endplot);ylim([0 1.1*max(abs(FRFOrder(Nst:Nend,1,:)))]);
     hold on
     plot(blktime(Nst:bsize:Nend),abs(FRFOrder(Nst:bsize:Nend,2)),'-rd','LineWidth',0.5,'MarkerEdgeColor','r', 'MarkerSize',3,'MarkerIndices',1:6:endplot)
     plot(blktime(Nst:bsize:Nend),abs(FRFOrder(Nst:bsize:Nend,3)),'-go','LineWidth',0.5,'MarkerEdgeColor','g', 'MarkerSize',3,'MarkerIndices',1:6:endplot)
     plot(blktime(Nst:Nend),abs(actamp(Nst:Nend,1)),':b','LineWidth',1.2);ylim([0 1.1*max((actamp(Nst:Nend,1)))]);        
     plot(blktime(Nst:Nend),abs(actamp(Nst:Nend,2)),':r','LineWidth',1.2);
     plot(blktime(Nst:Nend),abs(actamp(Nst:Nend,3)),':g','LineWidth',1.2);
        
%      plot(blktime,abs(FRFOrder(:,1,1)));ylim([0 max(abs(FRFOrder(:,1,1)))]);
%      hold on
%      plot(blktime,abs(FRFOrder(:,1,2)))
%      plot(blktime,abs(FRFOrder(:,1,3)))
%      ylabel({'Amp.'}, 'FontSize', 10);
    set(gcf,'position',[200,200,400,300])
    ylabel('Amp. (mm)');xlabel('time (s)')
    legend('o1=1X','o2=2X','o3=40Hz','actual');
    
    sumamp2=sum(actamp(Nst:Nend,:).^2,1);
    prda=100*sqrt(sum((acta(Nst:Nend,:)-real(FRFOrder(Nst:Nend,:))).^2,1)./sumamp2);
    prdb=100*sqrt(sum((actb(Nst:Nend,:)-imag(FRFOrder(Nst:Nend,:))).^2,1)./sumamp2);
    RMSEa=sqrt(sum((acta(Nst:Nend,:)-real(FRFOrder(Nst:Nend,:))).^2,1)/(Nend-Nst))
    RMSEb=sqrt(sum((actb(Nst:Nend,:)-imag(FRFOrder(Nst:Nend,:))).^2,1)/(Nend-Nst))
    MPEA=100*max(abs(actamp(Nst:Nend,:)-abs(FRFOrder(Nst:Nend,:))))./max(actamp(Nst:Nend,:))
   
%     for i=1:Np
%          for j=1:NTraces
%              if TraceType(j) == 0
%                reconssig(i) = reconssig(i) + real(FRFOrder(i,j))*cos(TraceOrder(j)*theta2(i))+imag(FRFOrder(i,j))*sin(TraceOrder(j)*theta2(i));
%              else
%                reconssig(i) = reconssig(i) + real(FRFOrder(i,j))*cos(TraceOrder(j)*2*pi*(i-1)*dt)+imag(FRFOrder(i,j))*sin(TraceOrder(j)*2*pi*(i-1)*dt);
%              end
%         end        
%     end
%     figure
%     plot((reconssig(Nst:Nend)));
%     title('Reconstructed signal', 'FontSize', 11, 'Color', 'k')
%     ylabel({'Recons'}, 'FontSize', 10);
%     figure
%     plot(Res2(Nst:Nend,1)-reconssig(Nst:Nend));%
%     title('Differnce actual - (Reconstructed)', 'FontSize', 11, 'Color', 'k')
%     ylabel({'Diff.'}, 'FontSize', 10);
    fprintf('constant DT, the number of blocks is:  %d \n',nblocks);       
    
end

if opt == 2 || opt == 5 % constant number of rev     
    nblocks =1 + fix(nr/bnr);
    %FRF=zeros(nblocks,No,1); 
    FRFOrder = zeros(Np,NTraces);
    freqBandOrder = zeros(nblocks,NTraces);    
    nblocks = 0;
    st =  2;%starting cycle 
    Nst = cpos(st);
    tindex=zeros(nblocks,1);
    while (st+bnr) <= nr

        stp = cpos(st);
        endp = cpos(st+bnr);
        bs = endp - stp+1;
        nblocks = nblocks + 1;    
        tindex(nblocks) =fix((stp+endp)/2);
        
        for i =1:No
            x=Res2(stp:dwnsamp:endp,i);
            xtheta=theta2(stp:dwnsamp:endp);
            %[OA,OB] = itvdft(x,size(x,1),TraceType,TraceOrder,NTraces,xtheta,stp/dwnsamp,dwnsamp/fs,OCMEnabled);
            
            speedCorr=corrEnable;
            if opt == 2
                [OA,OB] = otamlkh(x,size(x,1),TraceType,TraceOrder,NTraces,xtheta,stp/dwnsamp,dwnsamp/fs,speedCorr);
            else
                [OA,OB] = blougtvdft(x,size(x,1),TraceType,TraceOrder,NTraces,xtheta,stp/dwnsamp,dwnsamp/fs,OCMEnabled);
            end
                      
            k=1;
            for j=1:NTraces
                if TraceType(j) >= 0 
                    FRFOrder((stp:endp),k) = OA(j) + 1i*OB(j);% Cn =0.5(an - i bn)  
                   k = k+1;
                end
            end
            
        end
        %evaluating frequency bands for tracked orders
        k=1;
        for j=1:NTraces
            if TraceType(j) == 0 %include only orders, discard const freq orders                
                freqBandOrder(nblocks,k) = TraceOrder(j)*(1.8*OMG(stp) + 0.20* OMG(endp))/(4*pi);
                k = k+1;
            end
        end
        st = st + bnr-nroverlp;
    end
    Nend = endp; %trim edges
    tindex=tindex(1:nblocks);
    k=1;
    lastFreq = 0;
    freqBandOrder =  freqBandOrder(1:nblocks,:); %trim to include the processed blocks only
   
     blktime = dt*(0:Np-1);
     figure()
     
     plot(blktime(Nst:Nend),abs(FRFOrder(Nst:Nend,1)));ylim([0 max(abs(FRFOrder(:,1)))])
     hold on
     plot(blktime(Nst:bsize:Nend),actamp(Nst:bsize:Nend,1));ylim([0 1.1*max(abs(actamp(Nst:Nend,1)))])
     figure()
     plot(blktime(Nst:Nend),abs(FRFOrder(Nst:Nend,2)));ylim([0 max(abs(FRFOrder(:,2)))])
     hold on
     plot(blktime(Nst:Nend),actamp(Nst:Nend,2));ylim([0 1.1*max(abs(actamp(Nst:Nend,2)))])
     if(NTraces > 2)
       figure()
       plot(blktime(tindex),abs(FRFOrder(tindex,3)));ylim([0 max(abs(FRFOrder(:,3)))])
     end
     figure
%      bsize=2*bsize;
     endplot=nblocks;     
     plot(blktime(tindex),abs(FRFOrder(tindex,1)),'-bs','LineWidth',0.5,'MarkerEdgeColor','b', 'MarkerSize',3,'MarkerIndices',1:6:endplot);ylim([0 1.1*max(abs(FRFOrder(Nst:Nend,1,:)))]);
     hold on
     plot(blktime(tindex),abs(FRFOrder(tindex,2)),'-rd','LineWidth',0.5,'MarkerEdgeColor','r', 'MarkerSize',3,'MarkerIndices',1:6:endplot)
     plot(blktime(tindex),abs(FRFOrder(tindex,3)),'-go','LineWidth',0.5,'MarkerEdgeColor','g', 'MarkerSize',3,'MarkerIndices',1:6:endplot)
     plot(blktime(Nst:Nend),abs(actamp(Nst:Nend,1)),':b','LineWidth',1.2);ylim([0 1.1*max((actamp(Nst:Nend,1)))]);        
     plot(blktime(Nst:Nend),abs(actamp(Nst:Nend,2)),':r','LineWidth',1.2);
     plot(blktime(Nst:Nend),abs(actamp(Nst:Nend,3)),':g','LineWidth',1.2);
        
%      plot(blktime,abs(FRFOrder(:,1,1)));ylim([0 max(abs(FRFOrder(:,1,1)))]);
%      hold on
%      plot(blktime,abs(FRFOrder(:,1,2)))
%      plot(blktime,abs(FRFOrder(:,1,3)))
%      ylabel({'Amp.'}, 'FontSize', 10);
    set(gcf,'position',[200,200,400,300])
    ylabel('Amp. (mm)');xlabel('time (s)')
    legend('o1=1X','o2=2X','o3=40Hz','actual');
    
    sumamp2=sum(actamp(Nst:Nend,:).^2,1);
    prda=100*sqrt(sum((acta(Nst:Nend,:)-real(FRFOrder(Nst:Nend,:))).^2,1)./sumamp2);
    prdb=100*sqrt(sum((actb(Nst:Nend,:)-imag(FRFOrder(Nst:Nend,:))).^2,1)./sumamp2);
    RMSEa=sqrt(sum((acta(Nst:Nend,:)-real(FRFOrder(Nst:Nend,:))).^2,1)/(Nend-Nst))
    RMSEb=sqrt(sum((actb(Nst:Nend,:)-imag(FRFOrder(Nst:Nend,:))).^2,1)/(Nend-Nst))
    MPEA=100*max(abs(actamp(Nst:Nend,:)-abs(FRFOrder(Nst:Nend,:))))./max(actamp(Nst:Nend,:))
   
    
    fprintf('constant Number of revolutions, the number of blocks is:  %d \n',nblocks);       
    
end

if opt == 3 % vkf
%      Np=2000;
    FRFOrder = zeros(Np,No,NTraces);
      
    Nst=600;Nend = Np-200; %trim edges
    blktime = dt*(0:Np-1); 
    if NTraces == 2
        if(thetSrc == 0)
            fp=[f1 f2];%use actual
        end
        if(thetSrc == 1)
            fp=[yy TraceOrder(2)*yy]; % use corrected frequencies
        end
        if(thetSrc == 2)
            fp=[OMG0/(2*pi) TraceOrder(2)*OMG0/(2*pi)]; % use un-corrected frequencies
        end
        fp=fp(1:Np,:);
        xe=Res2(1:dwnsamp:Np,1);
        [xm,bwm] = vkmmy(xe,fp,fs,rm,forder);
        FRFOrder(:,1,1) = xm(:,1);
        FRFOrder(:,1,2) = xm(:,2);       
      
        figure()
        plot(blktime,abs(FRFOrder(Nst:Nend,1,1)));ylim([0 max(abs(FRFOrder(Nst:Nend,1,1)))])
        figure()
        plot(blktime,abs(FRFOrder(Nst:Nend,1,2)));ylim([0 max(abs(FRFOrder(Nst:Nend,1,2)))])
        
    else
        if(thetSrc == 0)
            fp=[f1 f2 f3];%use actual
        end
        if(thetSrc == 1)
            fp=[yy TraceOrder(2)*yy f3]; % use corrected frequencies
        end
        if(thetSrc == 2)
            fp=[OMG0/(2*pi) TraceOrder(2)*OMG0/(2*pi) f3]; % use un-corrected frequencies
        end
        fp=fp(1:Np,:);
        xe=Res2(1:dwnsamp:Np,1);
        [xm,bwm] = vkmmy(xe,fp,fs,rm,forder);
        
        FRFOrder(:,1,1) = xm(:,1);
        FRFOrder(:,1,2) = xm(:,2);
        FRFOrder(:,1,3) = xm(:,3);
       
        figure()
        plot(blktime(Nst:Nend),abs(FRFOrder(Nst:Nend,1,1)));%ylim([0 max(abs(FRFOrder(Nst:Nend,1,1)))])
        hold on
        plot(blktime(Nst:Nend),actamp(Nst:Nend,1));ylim([0 max(abs(actamp(Nst:Nend,1)))])
        figure()
        plot(blktime(Nst:Nend),abs(FRFOrder(Nst:Nend,1,2)));ylim([0 max(abs(FRFOrder(Nst:Nend,1,2)))])
        figure()
        plot(blktime(Nst:Nend),abs(FRFOrder(Nst:Nend,1,3)));ylim([0 max(abs(FRFOrder(Nst:Nend,1,3)))]);
        figure
        mytheta = exp(1i*2*pi*cumsum(fp,1)*dt);        
%         for i=1:Np
%             reconssig(i) = xm(i,1)*mytheta(i,1) + xm(i,2)*mytheta(i,2) + xm(i,3)*mytheta(i,3);
%         end
%         plot(imag(reconssig));
%         title('Imag. of Reconstructed signal', 'FontSize', 11, 'Color', 'k')
%         ylabel({'Imag.'}, 'FontSize', 10);
%         figure
%         plot(real(reconssig));
%         title('Real of Reconstructed signal', 'FontSize', 11, 'Color', 'k')
%         ylabel({'Real'}, 'FontSize', 10);
%         figure
%         plot(xe-real(reconssig));
%         title('Differnce actual - real(Reconstructed)', 'FontSize', 11, 'Color', 'k')
%         ylabel({'Diff.'}, 'FontSize', 10);

        figure        
        endplot=fix(Nend/bsize);
        plot(blktime(Nst:bsize:Nend),abs(FRFOrder(Nst:bsize:Nend,1)),'-bs','LineWidth',0.5,'MarkerEdgeColor','b', 'MarkerSize',3,'MarkerIndices',1:6:endplot);ylim([0 max(1.1*abs(FRFOrder(Nst:bsize:Nend,1,1)))]);
        hold on
        plot(blktime(Nst:bsize:Nend),abs(FRFOrder(Nst:bsize:Nend,2)),'-rd','LineWidth',0.5,'MarkerEdgeColor','r', 'MarkerSize',3,'MarkerIndices',1:6:endplot)
        plot(blktime(Nst:bsize:Nend),abs(FRFOrder(Nst:bsize:Nend,3)),'-go','LineWidth',0.5,'MarkerEdgeColor','g', 'MarkerSize',3,'MarkerIndices',1:6:endplot)
        plot(blktime(Nst:Nend),abs(actamp(Nst:Nend,1)),':b','LineWidth',1.2);ylim([0 1.1*max((actamp(Nst:Nend,1)))]);        
        plot(blktime(Nst:Nend),abs(actamp(Nst:Nend,2)),':r','LineWidth',1.2);
        plot(blktime(Nst:Nend),abs(actamp(Nst:Nend,3)),':g','LineWidth',1.2);
        
       ylabel('Amp. (mm)');xlabel('time (s)')
       set(gcf,'position',[200,200,400,300])
       legend('o1=1X','o2=2X','o3=40Hz','actual');

       RMSEa=sqrt(sum((acta(Nst:Nend,:)-real(FRFOrder(Nst:Nend,:))).^2,1)/(Nend-Nst))
       RMSEb=sqrt(sum((actb(Nst:Nend,:)+imag(FRFOrder(Nst:Nend,:))).^2,1)/(Nend-Nst))
       MPEA=100*max(abs(actamp(Nst:Nend,:)-abs(FRFOrder(Nst:Nend,:))))./max(actamp(Nst:Nend,:))

    end
  
end
if opt == 6  % VK with overlapping 
    
    if(thetSrc == 0)
        fp=[f1 f2 f3];%use actual
    end
    if(thetSrc == 1)
        fp=[OMG/(2*pi) 4*OMG/(2*pi) f3]; % use corrected frequencies
    end
    if(thetSrc == 2)
        fp=[OMG0/(2*pi) 4*OMG0/(2*pi) f3]; % use un-corrected frequencies
    end
    nblocks = 1 + fix(L/bsize);%initial number of blocks
    %FRF=zeros(nblocks,No,1); 
    FRFOrder = zeros(Np,No,NTraces);
    
    blktime=zeros(nblocks,1);
    nblocks = 0;
    st =  1;%starting point 
    k=1;
    while (st+bsize) <= Np/1
        stp = st;
        endp = st+bsize-1;
        nblocks = nblocks + 1;   
        tp =dt*(stp+endp)/2;% dt*stp;
        blktime(nblocks) = tp;
        xe=Res2(stp:dwnsamp:endp,1);
        fpb=fp(stp:dwnsamp:endp,:);
       
        
         if k < 2
             x0=[];
             [xm,bwm] = vkmmy(xe,fpb,fs,rm,forder);
             Lb = bsize-0.5*overlap;
             FRFOrder(stp:stp+Lb-1,1,1) = xm(1:Lb,1);
             FRFOrder(stp:stp+Lb-1,1,2) = xm(1:Lb,2);
             FRFOrder(stp:stp+Lb-1,1,3) = xm(1:Lb,3);
         else
             x0=xm(bsize-overlap:bsize-overlap/2-1,:);
             [xm,bwm] = vkmmy(xe,fpb,fs,rm,forder);
             Lb = bsize-overlap;
             stpb = overlap/2;
             stp2= stp+stpb;
             FRFOrder(stp2:stp2+Lb-1,1,1) = xm(stpb:stpb+Lb-1,1);
             FRFOrder(stp2:stp2+Lb-1,1,2) = xm(stpb:stpb+Lb-1,2);
             FRFOrder(stp2:stp2+Lb-1,1,3) = xm(stpb:stpb+Lb-1,3);
         end
                   
        st = st + (bsize-overlap);
        k = k+1;
    end
    Nst=1;Nend = endp-overlap/2; %trim edges
    figure()
    plot(abs(FRFOrder(Nst:Nend,1,1)));ylim([0 max(abs(FRFOrder(Nst:Nend,1,1)))])
    figure()
    plot(abs(FRFOrder(Nst:Nend,1,2)));ylim([0 max(abs(FRFOrder(Nst:Nend,1,2)))])
    figure()
    plot(abs(FRFOrder(Nst:Nend,1,3)));ylim([0 max(abs(FRFOrder(Nst:Nend,1,3)))]);

    
end

theta=theta(1:end-1000);
thetaM=thetaM(1:end-1000);
% figure();
% plot(thetaM,theta);
max(thetaM-theta);

