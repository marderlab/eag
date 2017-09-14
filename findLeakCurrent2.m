% findLeakCurrent
% this function segments a I, V curve into a leak current and everything else
% using heuristics based on Ohm's Law 

function [g_leak, E_leak, g_leak_err, E_leak_err] = findLeakCurrent(I,V)

% fit lines progressively and estimate slopes
r2 = NaN*V;
slopes = NaN*V;

for j = 2:length(V)
	x = V(1:j);
	y = I(1:j);
	[ff, gof] = fit(y(:),x(:),'poly1');
	r2(j) = gof.rsquare;
	slopes(j) = 1/ff.p1;

end

% calculate the slopes using [-80 to -60 range]
x = V(1:5);
y = I(1:5);
[ff, gof] = fit(y(:),x(:),'poly1');

E_leak = ff.p2;
g_leak = 1/ff.p1;

ci = confint(ff);
E_leak_err = diff(ci(:,2))/2;
g_leak_err = abs(diff(1./ci(:,1)))/2;


