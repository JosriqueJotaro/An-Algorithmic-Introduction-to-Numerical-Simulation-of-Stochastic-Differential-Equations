%BPATH3  Function along a Brownian path

% set the state of randn
rng(100,'v5normal');
T = 1; N = 500; dt = T/N; t = dt:dt:1;

% M paths simultaneously
M = 1000;

% increments
dW = sqrt(dt)*randn(M,N);

% cumulative sum
W = cumsum(dW,2);
U = exp(repmat(t,[M 1]) + 0.5*W);
Umean = mean(U);

% plot mean over M paths
plot([0,t],[1,Umean],'b-'), hold on

% plot 5 individual paths
plot([0,t],[ones(5,1),U(1:5,:)],'r--'), hold off
xlabel('t','FontSize',16)
ylabel('U(t)','FontSize',16,'Rotation',0,'HorizontalAlignment','right')
legend('mean of 1000 paths','5 individual paths','location','northwest')

% suppress code analyzer message
%#ok<*NOPTS>

% sample error
averr = norm((Umean - exp(9*t/8)),'inf')