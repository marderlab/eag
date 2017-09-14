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

% calculate the slopes using r^2 vs V
dynamic_cutoff = find(r2(V<0)>.9,1,'last');
x = V(1:dynamic_cutoff);
y = I(1:dynamic_cutoff);

V(dynamic_cutoff)

[ff, gof] = fit(y(:),x(:),'poly1');

E_leak = ff.p2;
g_leak = 1/ff.p1;
if dynamic_cutoff > 2
	ci = confint(ff);
	E_leak_err = diff(ci(:,2))/2;
	g_leak_err = abs(diff(1./ci(:,1)))/2;
else
	E_leak_err = Inf;
	g_leak_err = Inf;
end


