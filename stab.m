% STAB  Mean-square and asymptotic stability test for E-M
% Listing 9
%
% SDE is  dX = lambda*X dt + mu*X dW,   X(0) = Xzero,
%      where lambda and mu are constants and Xzero = 1.

rng(100,'v5normal');
T = 20; M = 50000; Xzero = 1;          

% linetypes for plot
ltype = {'b-','r--','m-.'};             

subplot(211)  %%%%%%%%%%%% Mean Square %%%%%%%%%%%%%

% problem parameters
lambda = -3; mu = sqrt(3);              
for k = 1:3
    Dt = 2^(1-k);                      
    N = T/Dt;
    Xms = zeros(1,N); Xtemp = Xzero*ones(M,1);
    for j = 1:N
           Winc = sqrt(Dt)*randn(M,1);  
           Xtemp = Xtemp + Dt*lambda*Xtemp + mu*Xtemp.*Winc;
           
           % mean-square estimate
           Xms(j) = mean(Xtemp.^2);     
    end
    semilogy(0:Dt:T,[Xzero,Xms],ltype{k},'Linewidth',2), hold on
end
legend('\Delta t = 1','\Delta t = 1/2','\Delta t = 1/4')
title('Mean-Square: \lambda = -3, \mu = \surd 3','FontSize',16)
ylabel('E[X^2]','FontSize',12), axis([0,T,1e-20,1e+20]), hold off

subplot(212)  %%%%% Asymptotic: a single path %%%%%%%

% problem parameters
T = 500;
lambda = 0.5; mu = sqrt(6);             
for k = 1:3
    Dt = 2^(1-k);                     
    N = T/Dt;
    Xemabs = zeros(1,N); Xtemp = Xzero;
    for j = 1:N
           Winc = sqrt(Dt)*randn;      
           Xtemp = Xtemp + Dt*lambda*Xtemp + mu*Xtemp*Winc;
           Xemabs(j) = abs(Xtemp);
    end
    semilogy(0:Dt:T,[Xzero,Xemabs],ltype{k},'Linewidth',2), hold on
end
legend('\Delta t = 1','\Delta t = 1/2','\Delta t = 1/4')
title('Single Path: \lambda = 1/2, \mu = \surd 6','FontSize',16)
ylabel('|X|','FontSize',12), axis([0,T,1e-50,1e+100]), hold off
