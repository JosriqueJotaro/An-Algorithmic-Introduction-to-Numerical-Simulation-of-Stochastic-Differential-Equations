%EMSTRONG  Test strong convergence of Euler-Maruyama
%
% Solves    dX = lambda*X dt + mu*X dW,   X(0) = Xzero,
%         where lambda = 2, mu = 1 and Xzer0 = 1.
%
% Discretized Brownian path over [0,1] has dt = 2^(-9).
% E-M uses 5 different timesteps: 16dt, 8dt, 4dt, 2dt, dt.
% Examine strong convergence at T=1:  E | X_L - X(T) |.

rng(100,'v5normal');

% problem parameters
lambda = 2; mu = 1; Xzero = 1;    
T = 1; N = 2^9; dt = T/N;    

% number of paths sampled
M = 1000;                         

% preallocate array
Xerr = zeros(M,5);                

% sample over discrete Brownian paths
for s = 1:M                       
    % Brownian increments
    dW = sqrt(dt)*randn(1,N);    
    
    % discrete Brownian path 
    W = cumsum(dW);               
    Xtrue = Xzero*exp((lambda-0.5*mu^2)+mu*W(end));
    
    for p = 1:5                       
        % L Euler steps of size Dt = R*dt
        R = 2^(p-1); Dt = R*dt; L = N/R;     
        Xtemp = Xzero;
        for j = 1:L
             Winc = sum(dW(R*(j-1)+1:R*j));
             Xtemp = Xtemp + Dt*lambda*Xtemp + mu*Xtemp*Winc;
        end
        
        % store the error at t = 1
        Xerr(s,p) = abs(Xtemp - Xtrue);      
    end
end

Dtvals = dt*(2.^(0:4));           

% top LH picture
subplot(221)                                 
loglog(Dtvals,mean(Xerr),'b*-'), hold on

% reference slope of 1/2 
loglog(Dtvals,(Dtvals.^(.5)),'r--'), hold off 
axis([1e-3 1e-1 1e-4 1])
xlabel('\Delta t'), ylabel('Sample average of | X(T) - X_L |')
title('emstrong.m','FontSize',10)

% suppress code analyzer message
%#ok<*NOPTS>

%%%% Least squares fit of error = C * Dt^q %%%%
A = [ones(5,1), log(Dtvals)']; rhs = log(mean(Xerr)');
sol = A\rhs; q = sol(2)
resid = norm(A*sol - rhs)