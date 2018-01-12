
%% Characterizing Potassium Channels in simulated data 
% First, before I move on to real data and the EAG channels, I set up a simple simulation of a cell with only potassium channels, and perform a virtual voltage clamp experiment. In the following figure, (a) shows how the asymptotic injected current varies as a function of the holding voltage. Now we can convert these currents into conductances using:
%
% $$ g(t)=\frac{I(t)}{V(t)-E} $$
% 

%%
% However, we still don't know how conductance depends on the gating variable. Specifically, what is the exponent $p$? In (b-e), I attempt to find the exponent by plotting various fractional powers of the conductance, and fitting Boltzman functions (red) to the data. It looks like an exponent of 4 looks like the best fit. This can be confirmed in (f), where I am plotting the goodness of fit vs. the exponent. Notice the large minimum corresponding to $p = 4$. Thus, we have determined both the exponent $p$, and how the asymptotic values of the gating variable depend on voltage. 

% warning off
% n = neuron;
% warning on

% n.parameters.gNa=  0;
% n.parameters.gCaT=  0;
% n.parameters.gCaS=  0;
% n.parameters.gLeak=  0;
% n.parameters.gA=  0;
% n.parameters.gKCa=  0;
% n.parameters.gH=  0;

% n.implementation = 'SGS/MATLAB';
% n.t_end = 100;

% time = n.dt:n.dt:n.t_end;

% V_clamp = 0*(time) - 50;
% n.voltage_clamp = V_clamp;
% [~,~,N] = n.integrate;


% V = linspace(-80,80,33);
% I_inf = 0*V;

% all_I = zeros(length(time),length(V));

% for i = 1:length(V)
% 	V_clamp = 0*(time) - 80;
% 	V_clamp(1e3:end) = V(i);
% 	n.voltage_clamp = V_clamp;
% 	[~, ~, N] = n.integrate;
% 	all_I(:,i) = N(:,14);
% end

% mI = all_I(end,:);

% figure('outerposition',[200 200 1501 902],'PaperUnits','points','PaperSize',[1501 902]); hold on
% subplot(2,3,1); hold on
% plot(V,mI,'k+')
% xlabel('Voltage (mV)')
% ylabel('Injected current (nA)')


% % now convert currents to conductances
% K_rev_pot = -80;
% g = mI./(V - K_rev_pot);


boltz = @(A,B,x) (1./(1 + exp((x+A)./B)));
tauX = @(A,B,D,E,x) (A - B./(1+exp((x+D)./E)));
% mKinf = @(V) boltz(12.3,-11.8,V);
% taumK = @(V) tauX(7.2,6.4,28.3,-19.2,V);

% % now plot the conductances, with various roots 
% p = [1:8];
% r2 = NaN*p;
% for i = 1:length(p)
% 	this_g = g.^(1/p(i));
% 	this_g = this_g/this_g(end);

% 	ff = fit(vectorise(V(2:end)),vectorise(this_g(2:end)),boltz,'StartPoint',[12.3, -11.8]);

% 	if i < 5
% 		subplot(2,3,1+p(i)); hold on
% 		plot(V,this_g,'k+')
% 		ylabel(['g^{1/' oval(p(i)) '} (norm)'])
% 		xlabel('Voltage (mV)')
% 		plot(V,ff(V),'r')
% 	end

% 	r2(i) = rsquare(ff(V(2:end)),this_g(2:end));

% end

% g = mI./(V - K_rev_pot);
% this_g = g.^(1/4);
% this_g = this_g/this_g(end);

% ff = fit(vectorise(V(2:end)),vectorise(this_g(2:end)),boltz,'StartPoint',[12.3, -11.8]);

% subplot(2,3,6); hold on
% plot(p,1-r2)
% xlabel('p')
% ylabel('1 - r^2')
% set(gca,'YScale','log')

% prettyFig();

% labelFigure('x_offset',-.01,'y_offset',-.0,'font_size',24)

% if being_published
% 	snapnow
% 	delete(gcf)
% end

%% Estimating the timescale of the conductance: simulations
% Now I attempt to recover how the timescale of action of the gating variable depends on voltage from this synthetic voltage clamp data. (a) show the raw time traces of the injected currents. Different colors indicate different holding potentials. (b) shows the same data, by raised to the inverse exponent to determine the time course of the gating variable. Note that in (b), I've normalized the data and also that the time series look exponential, unlike the sigmoidal time series in (a). I fit exponentials to all the time series in (b), and show the timescales estimated this way in (c) (black crosses). Also superimposed are the ground truth (the actual dependence of timescale on voltage, blue) and the recovered fit (red). It's reasonably good for most of the range. 

% % convert current traces into conductance traces
% % then convert conductance into gating variable
% g = all_I;
% n = all_I;
% for i = 1:length(V)
% 	g(:,i) = all_I(:,i)./(V(i) - K_rev_pot);
% 	n(:,i) = g(:,i).^(1/4);
% end
% n = n/max(max(n));

% % fit exponentials to this

% c = parula(34);

% figure('outerposition',[0 0 1501 500],'PaperUnits','points','PaperSize',[1501 500]); hold on

% subplot(1,3,1); hold on
% c = parula(34);
% for i = 1:length(V)
% 	plot(time-50,all_I(:,i),'Color',c(i,:))
% end
% set(gca,'XLim',[0 10])
% xlabel('Time (ms)')
% ylabel('Injected Current (nA)')

% subplot(1,3,2); hold on

% for i = 1:length(V)
% 	plot(time-50,n(:,i),'Color',c(i,:))
% end
% set(gca,'XLim',[0 10])
% xlabel('Time (ms)')
% ylabel('Estimated n (norm)')

% tau_n = NaN*V;
% tau_r2 =  NaN*V;
% f = fittype('1 - exp(-x./tau)');

% time = time(:);
% for i = 1:length(V)
% 	y = n(1e3:end,i);
% 	y = y/max(y);
% 	try
% 		[ft,gof] = fit(time(1:1e3+1),y(:),f,'StartPoint',1,'Upper',100,'Lower',0);
% 		tau_n(i) = ft.tau;
% 		tau_r2(i) = gof.rsquare;
% 	catch
% 	end
% end

% subplot(1,3,3); hold on
% plot(V,tau_n,'k+')

% plot(V,taumK(V),'b')

% ft = fit(vectorise(V(5:end)),vectorise(tau_n(5:end)),tauX,'StartPoint',[7,6,28,-20],'Upper',[10,10,100,10],'Lower',[-1, -1, -1, -30]);

% plot(V,ft(V),'r')
% xlabel('Voltage (mV)')
% ylabel('\tau_{m} (ms)')

% legend({'Voltage Clamp data','Ground truth','Fit'})

% prettyFig();

% labelFigure('x_offset',-.01,'y_offset',-.01,'font_size',28)

% if being_published
% 	snapnow
% 	delete(gcf)
% end

%% Simulations: confounding effects of leak currents
% What happens when leak currents also exist in the cell? The following figure shows the same analysis, but with a small leak conductance in the cell together with the potassium conductance. (a) shows the resulting current-voltage curve. It looks very similar to the old I-V curve, but a more careful examination reveals some key differences: that it is actually negative for very low membrane potentials, and that it looks linear at the very beginning. (b) shows this more clearly, which is the same curve, but zoomed into the region of interest. Note that it is strikingly linear. The membrane potential at which it crosses 0 current is indicated, and it turns out that it corresponds to the reversal potential of the leak current (-50 mV). So we have already characterized one key parameter of the leak current. 

%%
% Since leak currents are Ohmic, the region of the I-V curve where leak currents dominate should be a straight line. To determine this, I plot $r^2$ of linear fits starting from the lowest holding potential, as a function of membrane potential (c). We can clearly see that for membrane potentials < -20 mV, leak currents dominate, since $r^2$ is close to 1 here, suggesting perfect linearity. In (d), I plot the slope of these linear fits as a function of holding potential. The slope has units of $\mu S$. I've also indicated the true leak conductance with a red line, and we can see that the calculated slope is good agreement with the ground truth for membrane potentials < -20 mV. We have thus characterized both the reversal potential and the absolute strength of the leak conductance, despite it co-existing with a voltage-gated Potassium current. 

% warning off
% n = neuron;
% warning on

% n.parameters.gNa=  0;
% n.parameters.gCaT=  0;
% n.parameters.gCaS=  0;
% n.parameters.gLeak=  100;
% n.parameters.gA=  0;
% n.parameters.gKCa=  0;
% n.parameters.gH=  0;

% n.implementation = 'SGS/MATLAB';
% n.t_end = 100;

% time = n.dt:n.dt:n.t_end;

% V_clamp = 0*(time) - 50;
% n.voltage_clamp = V_clamp;
% [~,~,N] = n.integrate;


% V = linspace(-80,80,33);
% I_inf = 0*V;

% all_I = zeros(length(time),length(V));

% for i = 1:length(V)
% 	V_clamp = 0*(time) - 80;
% 	V_clamp(1e3:end) = V(i);
% 	n.voltage_clamp = V_clamp;
% 	[~, ~, N] = n.integrate;
% 	all_I(:,i) = N(:,14);
% end

% mI = all_I(end,:);

% figure('outerposition',[200 200 1401 1001],'PaperUnits','points','PaperSize',[1401 1001]); hold on
% subplot(2,2,1); hold on
% plot(V,mI,'k+-')
% xlabel('Membrane potential (mV)')
% ylabel('Injected current (nA)')

% subplot(2,2,2); hold on
% plot(V,mI,'k+-')
% plot([-50 -50],[-1e3 1e3],'r:')
% plot([-500 50],[0 0],'k:')
% xlabel('Membrane potential (mV)')
% ylabel('Injected current (nA)')
% set(gca,'XLim',[-80 -10],'YLim',[-200 400])

% % fit lines progressively and estimate slopes
% r2 = NaN*V;
% slopes = NaN*V;

% for i = 2:length(V)
% 	x = V(1:i);
% 	y = mI(1:i);
% 	[ff, gof] = fit(x(:),y(:),'poly1');
% 	r2(i) = gof.rsquare;
% 	slopes(i) = ff.p1;
% end

% subplot(2,2,3); hold on
% plot(V,r2,'k+')
% xlabel('Membrane potential (mV)')
% ylabel('r^2 of Ohmic fit')

% subplot(2,2,4); hold on
% plot(V,slopes,'k+')
% xlabel('Membrane potential (mV)')
% ylabel('Slope of Ohmic fit (\muS)')

% plot([-90 90],[6.3 6.3],'r:')
% set(gca,'YScale','log','XLim',[-80 80])

% prettyFig();

% labelFigure('x_offset',-.01,'y_offset',-.0,'font_size',28)

% if being_published
% 	snapnow
% 	delete(gcf)
% end
