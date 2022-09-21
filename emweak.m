% EMWEAK  Test weak convergence of Euler-Maruyama
% Listing 7
%
% Solves    dX = lambda*X dt + mu*X dW,   X(0) = Xzero,
%        where lambda = 2, mu = 0.1 and Xzer0 = 1.
%
% E-M uses 5 different timesteps: 2^(p-10),  p = 1,2,3,4,5.
% Examine weak convergence at T=1:   | E (X_L) - E (X(T)) |.
%
% Different paths are used for each E-M timestep.
% Code is vectorized over paths.
%
% Uncommenting the line indicated below gives the weak E-M method.

rng(100,'v5normal');                     
% problem parameters
lambda = 2; mu = 0.1; Xzero = 1; T = 1;  

% number of paths sampled
M = 50000;                               

% preallocate arrays
Xem = zeros(5,1);          

% take various Euler timesteps
for p = 1:5                 
    % L Euler steps of size Dt 
    Dt = 2^(p-10); L = T/Dt;          
    Xtemp = Xzero*ones(M,1);            
    for j = 1:L
        Winc = sqrt(Dt)*randn(M,1);
        
        % Winc = sqrt(Dt)*sign(randn(M,1)); %% use for weak E-M %%
        Xtemp = Xtemp + Dt*lambda*Xtemp + mu*Xtemp.*Winc;
    end
    Xem(p) = mean(Xtemp);
end
Xerr = abs(Xem - exp(lambda));

Dtvals = 2.^((1:5)-10);          

% top RH picture
subplot(222)                             
loglog(Dtvals,Xerr,'b*-'), hold on

% reference slope of 1
loglog(Dtvals,Dtvals,'r--'), hold off    
axis([1e-3 1e-1 1e-4 1])
xlabel('\Delta t'), ylabel('| E(X(T)) - Sample average of X_L |')
title('emweak.m','FontSize',10)

% suppress code analyzer message
%#ok<*NOPTS>

%%%% Least squares fit of error = C * dt^q %%%%
A = [ones(p,1), log(Dtvals)']; rhs = log(Xerr);
sol = A\rhs; q = sol(2)
resid = norm(A*sol - rhs)
