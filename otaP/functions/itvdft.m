function [OA,OB] = itvdft(x,N,TraceType,TraceOrder,NTraces,Theta,StP,deltaT,OCMEnabled)
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
OA = zeros(NTraces,1);
OB = zeros(NTraces,1);
Phi = zeros(N,NTraces);
%prepare angular data
for j=1:NTraces
    for i=1:N
        if TraceType(j) == 0 %order
            Phi(i,j) = TraceOrder(j)*Theta(i);
        end
         if TraceType(j) == 1 %fixed frequency in Hz
            Phi(i,j) = 2 * pi * TraceOrder(j) * (StP + i -1) * deltaT;
        end
    end
end
z=1000; %to avoid small values
NOrder = NTraces;
NO2 = 2*NOrder;
A = zeros(NO2,NO2); g=zeros(NOrder,1);
B = zeros(NO2,1); Tf= zeros(NO2,1);
for k=1:N
    for i=1:NOrder
        g(i) = Phi(k,i);
    end
    for i = 0:NOrder - 1
        Tf(i * 2 + 1) = cos(g(i + 1));
        Tf(i * 2 + 2) = sin(g(i + 1));
    end
    for i = 0:NOrder - 1
        B(i * 2 + 1) = B(i * 2 + 1) + x(k) * z * cos(g(i + 1));
        B(i * 2 + 2) = B(i * 2 + 2) + x(k) * z * sin(g(i + 1));
    end
    for i = 1:NO2
        for m = i:NO2
            A(i, m) = A(i, m) + Tf(i) * z * Tf(m);
        end
    end
    for i = 2:NO2
        for m = 1 : i - 1
            A(i, m) = A(m, i);
        end
    end
    
    if OCMEnabled > 0
        AC = pinv(A)*B;
        for i = 0:NOrder - 1
            OA(i + 1) = AC(2 * i + 1);
            OB(i + 1) = AC(2 * i + 2);
        end
    else
        for i = 0:NOrder - 1
            OA(i + 1) = 2*B(2 * i + 1)/N/z;
            OB(i + 1) = 2*B(2 * i + 2)/N/z;
        end

    end
end

end


