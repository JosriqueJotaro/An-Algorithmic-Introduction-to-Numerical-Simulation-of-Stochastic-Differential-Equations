% BPATH1  Brownian path simulation
% Listing 1

% set the state of randn
rng(100,'v5normal');

% static parameters
T = 1; N = 500; dt = T/N;

% preallocate arrays for efficiency
dW = zeros(1,N);
W = zeros(1,N);

% first approximation outside the loop since W(0) = 0 is not allowed
dW(1) = sqrt(dt)*randn;
W(1) = dW(1);
for j = 2:N
    % general increment
    dW(j) = sqrt(dt)*randn;   
    W(j) = W(j-1) + dW(j); 
end

% plot W against t
plot(0:dt:T,[0,W],'r-')
xlabel('t','FontSize',16) 
ylabel('W(t)','FontSize',16,'Rotation',0)