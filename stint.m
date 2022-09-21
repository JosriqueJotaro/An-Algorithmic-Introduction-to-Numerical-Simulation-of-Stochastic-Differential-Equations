% STINT  Approximate stochastic integrals
% Listing 4
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


% the Ito integral refers to the limiting case of the stochastic left hand
% Riemann sum, whereas the Stratonovich integral refers to the limiting
% case of the stochastic midpoint Riemann sum. The purpose of this
% calculation is to demonstrate that as N approaches infinity (limiting
% case), that the Ito and Stratonovich Integals do not converge as one may
% suspect.
ito = sum([0,W(1:end-1)].*dW)
strat = sum((0.5*([0,W(1:end-1)]+W) + 0.5*sqrt(dt)*randn(1,N)).*dW) 

itoerr = abs(ito - 0.5*(W(end)^2-T))
straterr = abs(strat - 0.5*W(end)^2)