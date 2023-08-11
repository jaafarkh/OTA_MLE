function [OA,OB] = otamlkhRev(x,N,TraceType,TraceOrder,NTraces,Theta,StP,deltaT,showerr,cpos,nrev)
%ITVDFT performs improved time variant discrete Fourier transform
% usage: itvdft(x,N,TraceType,TraceOrder,NTraces,Theta,StP,deltaT,OCMEnabled)
% INPUTS:
% x: input signal (starting from 1 index)
% Theta : angular data (1-based index)
% TraceType,TraceOrder,NTrace: order information
% NTraces: number of traces
% TraceType: 0 : Order, 1: fixed frequency
% TraceOrder: order multplier or frequency
% OCMEnable: 1: OCM enabled, 0:disabled
% StP: processing time start point (used for fixed frequency traces)
% deltaT = 1 /(Sample Rate)
% OUTPUTS:
% OA, OB: vectors containing order amplitude and phase (Real and Imaginary
% data)
z=100*hann(N);
x=x.*z;
OA = zeros(NTraces,1);
OB = zeros(NTraces,1);
Phi = zeros(N,NTraces);
B = zeros(N,NTraces*2);
%prepare angular data
NOrder = NTraces;

for j=1:NTraces
     if TraceType(j) == 0 %order
            Phi(:,j) = TraceOrder(j)*Theta(:);
     end
     if TraceType(j) == 1 %fixed frequency in Hz
            Phi(:,j) = 2 * pi * TraceOrder(j) * (StP + (0:N-1)) * deltaT;
     end
end
% for k=1:N   
%     for i = 0:NOrder - 1
%         B(k,i * 2 + 1) = cos(Phi(k,i+1));
%         B(k,i * 2 + 2) = sin(Phi(k,i+1));
%     end   
% end
for i = 0:NOrder - 1
    B(:,i * 2 + 1) = cos(Phi(:,i+1)).*z;
    B(:,i * 2 + 2) = sin(Phi(:,i+1)).*z;
end
AC=pinv(B)*x;
% AC=pinv(B'*B)*B'*x;%yet another method
NT=40;
e=zeros(NT,nrev);
cgr=zeros(NT,nrev);
if(showerr > 0)
    err=x-B*AC;
    g=0.05;
    cg = 0.0;
    e(1)=sum(err.^2);
    %e(1)=max(abs(err));
    for tr=2:NT
        cg = cg+g;
        cgr(tr)=cg;
        for j=1:NTraces
            if TraceType(j) == 0 %order
                Phi(:,j) = TraceOrder(j)*(Theta(:)+cg*deltaT*(0:N-1)');%*deltaT*(0:N-1)'
            end
            if TraceType(j) == 1 %fixed frequency in Hz
                Phi(:,j) = 2 * pi * TraceOrder(j) * (StP + (0:N-1)) * deltaT;
            end
        end
        
        for i = 0:NOrder - 1
            B(:,i * 2 + 1) = cos(Phi(:,i+1)).*z;
            B(:,i * 2 + 2) = sin(Phi(:,i+1)).*z;
        end
       
        AC=pinv(B)*x;
        err=x-B*AC;
        e(tr)=sum(err.^2);
        %e(tr)=max(abs(err));
        if(e(tr)>e(tr-1))
            %g=-0.5*g; %*(e(tr)-e(tr-1))
            g=-0.5*g;
           % tr
        else
            g=1.5*g;
        end
        
    end
    
    [~, minIndex] = min(e);
    %e
    cg=cgr(minIndex);
    for j=1:NTraces
        if TraceType(j) == 0 %order
            Phi(:,j) = TraceOrder(j)*(Theta(:)+cg*deltaT*(0:N-1)');
        end
        if TraceType(j) == 1 %fixed frequency in Hz
            Phi(:,j) = 2 * pi * TraceOrder(j) * (StP + (0:N-1)) * deltaT;
        end
    end
    for i = 0:NOrder - 1
        B(:,i * 2 + 1) = cos(Phi(:,i+1)).*z;
        B(:,i * 2 + 2) = sin(Phi(:,i+1)).*z;
    end
    AC=pinv(B)*x;
end

 for i = 0:NOrder - 1
     OA(i + 1) = AC(2 * i + 1);
     OB(i + 1) = AC(2 * i + 2);
 end
 

end


