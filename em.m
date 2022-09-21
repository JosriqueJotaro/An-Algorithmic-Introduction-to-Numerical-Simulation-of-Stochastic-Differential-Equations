%EM  Euler-Maruyama method on linear SDE
%
% SDE is  dX = lambda*X dt + mu*X dW,   X(0) = Xzero,
%      where lambda = 2, mu = 1 and Xzero = 1.
%
% Discretized Brownian path over [0,1] has dt = 2^(-8).
% Euler-Maruyama uses timestep R*dt.

rng(100,'v5normal');

% problem parameters
lambda = 2; mu = 1; Xzero = 1;    
T = 1; N = 2^8; dt = T/N;       

% Brownian increments
dW = sqrt(dt)*randn(1,N);       

% discretized Brownian path 
W = cumsum(dW);                   

Xtrue = Xzero*exp((lambda-0.5*mu^2)*(dt:dt:T)+mu*W); 
plot(0:dt:T,[Xzero,Xtrue],'m-'), hold on

% L EM steps of size Dt = R*dt
R = 4; Dt = R*dt; L = N/R;       

% preallocate for efficiency
Xem = zeros(1,L);                 
Xtemp = Xzero;

for j = 1:L
    Winc = sum(dW(R*(j-1)+1:R*j)); 
    Xtemp = Xtemp + Dt*lambda*Xtemp + mu*Xtemp*Winc;
    Xem(j) = Xtemp;
end

plot(0:Dt:T,[Xzero,Xem],'r--*'), hold off
xlabel('t','FontSize',12)
ylabel('X','FontSize',16,'Rotation',0,'HorizontalAlignment','right')

% suppress code analyzer message
%#ok<*NOPTS>

emerr = abs(Xem(end)-Xtrue(end))