

alphaR = @(V) 30/(1 + exp((V+230)/(25)));
alphaN = @(V) 0.09/(1 + exp((V+100)/(9)));
betaR = @(V) .15/(1 + exp(-(V+120)/(20)));
betaN = @(V) 0.00035*exp(0.07*(V+25));
n_inf = @(A,B) A./(A+B);
tau_n = @(A,B) 1./(A+B); 


V = linspace(-100,70,1e3);
tau_N = 0*V;
tau_R = 0*V;
N_inf = 0*V;
R_inf = 0*V;

for i = 1:length(V)
	N_inf(i) = n_inf(alphaN(V(i)),betaN(V(i)));
	R_inf(i) = n_inf(alphaR(V(i)),betaR(V(i)));

	tau_N(i) = tau_n(alphaN(V(i)),betaN(V(i)));
	tau_R(i) = tau_n(alphaR(V(i)),betaR(V(i)));
end

figure, hold on
subplot(1,2,1); hold on
plot(V,N_inf,'k')
plot(V,R_inf,'r')
plot(V,R_inf.*N_inf,'k--')


subplot(1,2,2); hold on
plot(V,tau_N,'k')
plot(V,tau_R,'r')
set(gca,'YScale','log')