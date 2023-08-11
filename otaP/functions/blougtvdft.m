function [OA,OB] = blougtvdft(x,N,TraceType,TraceOrder,NTraces,Theta,StP,deltaT,OCMEnabled)
%ITVDFT performs Blough time variant discrete Fourier transform
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

Phi = zeros(N,NTraces);
%prepare angular data
for j=1:NTraces
     if TraceType(j) == 0 %order
            Phi(:,j) = TraceOrder(j)*Theta(:);
     end
     if TraceType(j) == 1 %fixed frequency in Hz
            Phi(:,j) = 2 * pi * TraceOrder(j) * (StP + (0:N-1)) * deltaT;
     end
end
z=1000*hann(N); %to avoid small values
% z=1000*ones(N,1);
NOrder = NTraces;

A = zeros(NOrder,NOrder); 
B = zeros(NOrder,1); 
x=x.*z;
for i = 1:NOrder 
        B(i,1) = (1/N)*sum(exp(-1i*Phi(:,i)).*x);   
        g1=Phi(:,i);
        for j=1:NOrder
            g2=Phi(:,j);
            A(i,j) = (1/N)*sum(exp(1i*(g2-g1)).*z) ;
        end
end

if OCMEnabled > 0
    AC = pinv(A)*B;
    OA=2*real(AC);
    OB=-2*imag(AC);
else
    OA=2*real(B);
    OB=-2*imag(B);
end

end


