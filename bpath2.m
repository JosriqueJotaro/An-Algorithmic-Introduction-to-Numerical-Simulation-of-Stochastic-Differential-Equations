% BPATH2  Brownian path simulation: vectorized
% Listing 2

% set the state of randn
rng(100,'v5normal');
T = 1; N = 500; dt = T/N;

% increments
dW = sqrt(dt)*randn(1,N);

% cumulative sum
W = cumsum(dW);             

% plot W against t
plot(0:dt:T,[0,W],'r-')   
xlabel('t','FontSize',16)
ylabel('W(t)','FontSize',16,'Rotation',0)