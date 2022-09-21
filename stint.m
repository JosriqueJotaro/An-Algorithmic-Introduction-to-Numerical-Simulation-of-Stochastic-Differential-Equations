%STINT  Approximate stochastic integrals
%
% Ito and Stratonovich integrals of W dW

% set the state of randn
rng(100,'v5normal');                      
T = 1; N = 500; dt = T/N;

% increments
dW = sqrt(dt)*randn(1,N);               

% cumulative sum
W = cumsum(dW);                         

% suppress code analyzer message
%#ok<*NOPTS>

ito = sum([0,W(1:end-1)].*dW)
strat = sum((0.5*([0,W(1:end-1)]+W) + 0.5*sqrt(dt)*randn(1,N)).*dW) 

itoerr = abs(ito - 0.5*(W(end)^2-T))
straterr = abs(strat - 0.5*W(end)^2)