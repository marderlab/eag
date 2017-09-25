% this script estimates m_inf and tau_inf values for the EAG current
% from experimental data 



pHeader;

%%
% In this document, I attempt to characterize EAG channels in the Hodgkin-Huxley formalism, where the current of a channel $i$ is given by
%
% $$ I_{i}=g_{i}m^{p}h(V-E_{i}) $$ 
%

%%
% where $g$ is the conductance, $m$ and $h$ are gating variables, and $E_{i}$ is the reversal potential. 

%%
% When cells are voltage clamped, the current injected into the cell to keep it at the clamped voltage is, by definition, equal to the membrane currents. If there exists only one conductance in the cell, then the injected currents are the currents that flow through the conductance of interest. 

%% Characterizing Potassium Channels in simulated data 
% First, before I move on to real data and the EAG channels, I set up a simple simulation of a cell with only potassium channels, and perform a virtual voltage clamp experiment. In the following figure, (a) shows how the asymptotic injected current varies as a function of the holding voltage. Now we can convert these currents into conductances using:
%
% $$ g(t)=\frac{I(t)}{V(t)-E} $$
% 

%%
% However, we still don't know how conductance depends on the gating variable. Specifically, what is the exponent $p$? In (b-e), I attempt to find the exponent by plotting various fractional powers of the conductance, and fitting Boltzman functions (red) to the data. It looks like an exponent of 4 looks like the best fit. This can be confirmed in (f), where I am plotting the goodness of fit vs. the exponent. Notice the large minimum corresponding to $p = 4$. Thus, we have determined both the exponent $p$, and how the asymptotic values of the gating variable depend on voltage. 

warning off
n = neuron;
warning on

n.parameters.gNa=  0;
n.parameters.gCaT=  0;
n.parameters.gCaS=  0;
n.parameters.gLeak=  0;
n.parameters.gA=  0;
n.parameters.gKCa=  0;
n.parameters.gH=  0;

n.implementation = 'SGS/MATLAB';
n.t_end = 100;

time = n.dt:n.dt:n.t_end;

V_clamp = 0*(time) - 50;
n.voltage_clamp = V_clamp;
[~,~,N] = n.integrate;


V = linspace(-80,80,33);
I_inf = 0*V;

all_I = zeros(length(time),length(V));

for i = 1:length(V)
	V_clamp = 0*(time) - 80;
	V_clamp(1e3:end) = V(i);
	n.voltage_clamp = V_clamp;
	[~, ~, N] = n.integrate;
	all_I(:,i) = N(:,14);
end

mI = all_I(end,:);

figure('outerposition',[200 200 1501 902],'PaperUnits','points','PaperSize',[1501 902]); hold on
subplot(2,3,1); hold on
plot(V,mI,'k+')
xlabel('Voltage (mV)')
ylabel('Injected current (nA)')


% now convert currents to conductances
K_rev_pot = -80;
g = mI./(V - K_rev_pot);


boltz = @(A,B,x) (1./(1 + exp((x+A)./B)));
tauX = @(A,B,D,E,x) (A - B./(1+exp((x+D)./E)));
mKinf = @(V) boltz(12.3,-11.8,V);
taumK = @(V) tauX(7.2,6.4,28.3,-19.2,V);

% now plot the conductances, with various roots 
p = [1:8];
r2 = NaN*p;
for i = 1:length(p)
	this_g = g.^(1/p(i));
	this_g = this_g/this_g(end);

	ff = fit(vectorise(V(2:end)),vectorise(this_g(2:end)),boltz,'StartPoint',[12.3, -11.8]);

	if i < 5
		subplot(2,3,1+p(i)); hold on
		plot(V,this_g,'k+')
		ylabel(['g^{1/' oval(p(i)) '} (norm)'])
		xlabel('Voltage (mV)')
		plot(V,ff(V),'r')
	end

	r2(i) = rsquare(ff(V(2:end)),this_g(2:end));

end

g = mI./(V - K_rev_pot);
this_g = g.^(1/4);
this_g = this_g/this_g(end);

ff = fit(vectorise(V(2:end)),vectorise(this_g(2:end)),boltz,'StartPoint',[12.3, -11.8]);

subplot(2,3,6); hold on
plot(p,1-r2)
xlabel('p')
ylabel('1 - r^2')
set(gca,'YScale','log')

prettyFig();

labelFigure('x_offset',-.01,'y_offset',-.0,'font_size',24)

if being_published
	snapnow
	delete(gcf)
end

%% Estimating the timescale of the conductance: simulations
% Now I attempt to recover how the timescale of action of the gating variable depends on voltage from this synthetic voltage clamp data. (a) show the raw time traces of the injected currents. Different colors indicate different holding potentials. (b) shows the same data, by raised to the inverse exponent to determine the time course of the gating variable. Note that in (b), I've normalized the data and also that the time series look exponential, unlike the sigmoidal time series in (a). I fit exponentials to all the time series in (b), and show the timescales estimated this way in (c) (black crosses). Also superimposed are the ground truth (the actual dependence of timescale on voltage, blue) and the recovered fit (red). It's reasonably good for most of the range. 

% convert current traces into conductance traces
% then convert conductance into gating variable
g = all_I;
n = all_I;
for i = 1:length(V)
	g(:,i) = all_I(:,i)./(V(i) - K_rev_pot);
	n(:,i) = g(:,i).^(1/4);
end
n = n/max(max(n));

% fit exponentials to this

c = parula(34);

figure('outerposition',[0 0 1501 500],'PaperUnits','points','PaperSize',[1501 500]); hold on

subplot(1,3,1); hold on
c = parula(34);
for i = 1:length(V)
	plot(time-50,all_I(:,i),'Color',c(i,:))
end
set(gca,'XLim',[0 10])
xlabel('Time (ms)')
ylabel('Injected Current (nA)')

subplot(1,3,2); hold on

for i = 1:length(V)
	plot(time-50,n(:,i),'Color',c(i,:))
end
set(gca,'XLim',[0 10])
xlabel('Time (ms)')
ylabel('Estimated n (norm)')

tau_n = NaN*V;
tau_r2 =  NaN*V;
f = fittype('1 - exp(-x./tau)');

time = time(:);
for i = 1:length(V)
	y = n(1e3:end,i);
	y = y/max(y);
	try
		[ft,gof] = fit(time(1:1e3+1),y(:),f,'StartPoint',1,'Upper',100,'Lower',0);
		tau_n(i) = ft.tau;
		tau_r2(i) = gof.rsquare;
	catch
	end
end

subplot(1,3,3); hold on
plot(V,tau_n,'k+')

plot(V,taumK(V),'b')

ft = fit(vectorise(V(5:end)),vectorise(tau_n(5:end)),tauX,'StartPoint',[7,6,28,-20],'Upper',[10,10,100,10],'Lower',[-1, -1, -1, -30]);

plot(V,ft(V),'r')
xlabel('Voltage (mV)')
ylabel('\tau_{m} (ms)')

legend({'Voltage Clamp data','Ground truth','Fit'})

prettyFig();

labelFigure('x_offset',-.01,'y_offset',-.01,'font_size',28)

if being_published
	snapnow
	delete(gcf)
end

%% Simulations: confounding effects of leak currents
% What happens when leak currents also exist in the cell? The following figure shows the same analysis, but with a small leak conductance in the cell together with the potassium conductance. (a) shows the resulting current-voltage curve. It looks very similar to the old I-V curve, but a more careful examination reveals some key differences: that it is actually negative for very low membrane potentials, and that it looks linear at the very beginning. (b) shows this more clearly, which is the same curve, but zoomed into the region of interest. Note that it is strikingly linear. The membrane potential at which it crosses 0 current is indicated, and it turns out that it corresponds to the reversal potential of the leak current (-50 mV). So we have already characterized one key parameter of the leak current. 

%%
% Since leak currents are Ohmic, the region of the I-V curve where leak currents dominate should be a straight line. To determine this, I plot $r^2$ of linear fits starting from the lowest holding potential, as a function of membrane potential (c). We can clearly see that for membrane potentials < -20 mV, leak currents dominate, since $r^2$ is close to 1 here, suggesting perfect linearity. In (d), I plot the slope of these linear fits as a function of holding potential. The slope has units of $\mu S$. I've also indicated the true leak conductance with a red line, and we can see that the calculated slope is good agreement with the ground truth for membrane potentials < -20 mV. We have thus characterized both the reversal potential and the absolute strength of the leak conductance, despite it co-existing with a voltage-gated Potassium current. 

warning off
n = neuron;
warning on

n.parameters.gNa=  0;
n.parameters.gCaT=  0;
n.parameters.gCaS=  0;
n.parameters.gLeak=  100;
n.parameters.gA=  0;
n.parameters.gKCa=  0;
n.parameters.gH=  0;

n.implementation = 'SGS/MATLAB';
n.t_end = 100;

time = n.dt:n.dt:n.t_end;

V_clamp = 0*(time) - 50;
n.voltage_clamp = V_clamp;
[~,~,N] = n.integrate;


V = linspace(-80,80,33);
I_inf = 0*V;

all_I = zeros(length(time),length(V));

for i = 1:length(V)
	V_clamp = 0*(time) - 80;
	V_clamp(1e3:end) = V(i);
	n.voltage_clamp = V_clamp;
	[~, ~, N] = n.integrate;
	all_I(:,i) = N(:,14);
end

mI = all_I(end,:);

figure('outerposition',[200 200 1401 1001],'PaperUnits','points','PaperSize',[1401 1001]); hold on
subplot(2,2,1); hold on
plot(V,mI,'k+-')
xlabel('Membrane potential (mV)')
ylabel('Injected current (nA)')

subplot(2,2,2); hold on
plot(V,mI,'k+-')
plot([-50 -50],[-1e3 1e3],'r:')
plot([-500 50],[0 0],'k:')
xlabel('Membrane potential (mV)')
ylabel('Injected current (nA)')
set(gca,'XLim',[-80 -10],'YLim',[-200 400])

% fit lines progressively and estimate slopes
r2 = NaN*V;
slopes = NaN*V;

for i = 2:length(V)
	x = V(1:i);
	y = mI(1:i);
	[ff, gof] = fit(x(:),y(:),'poly1');
	r2(i) = gof.rsquare;
	slopes(i) = ff.p1;
end

subplot(2,2,3); hold on
plot(V,r2,'k+')
xlabel('Membrane potential (mV)')
ylabel('r^2 of Ohmic fit')

subplot(2,2,4); hold on
plot(V,slopes,'k+')
xlabel('Membrane potential (mV)')
ylabel('Slope of Ohmic fit (\muS)')

plot([-90 90],[6.3 6.3],'r:')
set(gca,'YScale','log','XLim',[-80 80])

prettyFig();

labelFigure('x_offset',-.01,'y_offset',-.0,'font_size',28)

if being_published
	snapnow
	delete(gcf)
end


%% EAG channels: estimating voltage dependence of gating variable 
% To model EAG channels, I need to understand how these channels vary their activity as a function of both voltage and Calcium concentrations. These channels don't inactivate, so their behavior can be described by a single gating variable ($m$). In Hodgkin-Huxley terms, we want to determine how the steady state value of $m$ depends on the voltage and Calcium, and how the timescale of m changes with voltage and Calcium. 

%%
% First, I estimate the voltage dependence on the gating variable by looking at currents in voltage clamp in HEK cells. In this data, EAG channels were inserted into HEK cells, and currents recorded from them in voltage clamp. The experiment was repeated with HEK cells without EAG.  


% load all the data and convert it into a nicer format
path_name = '/Users/srinivas/data/calcium-microdomain/HEK';
load(joinPath(path_name,'WTData.mat'))

time = 1:length(WTFile{1,1}{1});

all_I = zeros(length(time),0);
Ca_levels = zeros(0);
cell_id = zeros(0);
all_V = zeros(0);
genotype = zeros(0); % 2 = WT EAG, 1 = No EAG, 3 = Mutated EAG

V_space = linspace(-80,80,33);
Ca_space = [13e-9 100e-9 10e-6 100e-6 500e-6]; % M
c = 1;

for i = 1:size(WTFile,1)
	for j = 1:size(WTFile,2)
		if ~isempty(WTFile{i,j})
			thisI = WTFile{i,j}{1};
			for k = 1:33
				% assemble
				all_I(:,c) = thisI(:,k);
				Ca_levels(c) = Ca_space(j);
				cell_id(c) = i;
				all_V(c) = V_space(k);
				genotype(c) = 2; 
				c = c + 1;
			end
		end
	end
end

load(joinPath(path_name,'NoEAGData.mat'))

for i = 1:size(NoEAGFile,1)
	for j = 1:size(NoEAGFile,2)
		if ~isempty(NoEAGFile{i,j})
			thisI = NoEAGFile{i,j}{1};
			for k = 1:33
				% assemble
				all_I(:,c) = thisI(:,k);
				Ca_levels(c) = Ca_space(j);
				cell_id(c) = i;
				all_V(c) = V_space(k);
				genotype(c) = 1; 
				c = c + 1;
			end
		end
	end
end

load(joinPath(path_name,'MutData.mat'))

for i = 1:size(MutFile,1)
	for j = 1:size(MutFile,2)
		if ~isempty(MutFile{i,j})
			thisI = MutFile{i,j}{1};
			for k = 1:33
				% assemble
				all_I(:,c) = thisI(:,k);
				Ca_levels(c) = Ca_space(j);
				cell_id(c) = i;
				all_V(c) = V_space(k);
				genotype(c) = 3; 
				c = c + 1;
			end
		end
	end
end


% compute asymptotic currents for all traces
I_inf = mean(all_I(4e3:5e3,:));

%% A principled way to determine leak conductances from noisy data
% First, I look at the currents of each neuron and see if I can determine leak currents in these cells that I can clearly separate from the EAG currents. Leak currents are Ohmic, and are therefore manifest as linear curves of I vs. V. (a) show sthe I-V curves for all the HEK cells expressing EAG at 13 nM Calcium. Each line is a different cell. (b) shows the $r^2$ of the linear fit to these curves, as a function of membrane potential, similar to the curves I showed in the simulated data. Only 4 curves are shown for clarity. (c) shows the slopes of these fits, as a function of membrane potential, for each cell. 

%%
% In the simulations, there was no problem identifying the region where leak currents dominated, as there was no noise, and it was perfectly linear. In real data, deviations from linearity can arise due to either voltage gated channels, or due to random noise. How do we tell the difference? What subset of the data should I use to estimate the leak conductance? There isn't one perfect answer: in (d-e), I choose different regions to illustrate what happens to my estimates of leak conductance (on the Y-axis) and the leak reversal potential (X-axis). Each cross is a single neuron, and the error bars indicate 95% confidence intervals. A more principled way to do this is to pick the range of V based on the curves in (b) (basically, you find the point where the $r^2$ curve drops < .9). This is what I did in (f), and this leads to estimates of both leak conductance and reversal potentials with much lower errors. 

figure('outerposition',[0 0 1602 1201],'PaperUnits','points','PaperSize',[1602 1201]); hold on

subplot(2,3,1); hold on
for i = 1:max(cell_id)
	this_I = I_inf(cell_id == i & genotype == 2 & Ca_levels == Ca_space(1));
	plot(V_space,this_I,'+-')
end
set(gca,'XLim',[-80 0],'YLim',[-2000 4000])
xlabel('Membrane potential (mV)')
ylabel('Current (pA)')

E_leak = NaN(max(cell_id),3);
E_leak_err = NaN(max(cell_id),3);
g_leak = NaN(max(cell_id),3);
g_leak_err = NaN(max(cell_id),3);

cutoff = [5 10 13];

for i = 1:max(cell_id)
	% fit lines progressively and estimate slopes
	r2 = NaN*V_space;
	slopes = NaN*V_space;

	this_I = I_inf(cell_id == i & genotype == 2 & Ca_levels == Ca_space(1));

	for j = 2:length(V_space)
		x = V_space(1:j);
		y = this_I(1:j);
		[ff, gof] = fit(y(:),x(:),'poly1');
		r2(j) = gof.rsquare;
		slopes(j) = 1/ff.p1;

		for k = 1:length(cutoff)-1
			if j == cutoff(k)
				E_leak(i,k) = ff.p2;
				ci = confint(ff);
				E_leak_err(i,k) = diff(ci(:,2))/2;
				g_leak(i,k) = slopes(j);
				g_leak_err(i,k) = abs(diff(1./ci(:,1)))/2;
			end
		end
	end

	% also calculate it dynamically 
	dynamic_cutoff = find(r2(V_space<0)>.9,1,'last');
	x = V_space(1:dynamic_cutoff);
	y = this_I(1:dynamic_cutoff);
	[ff, gof] = fit(y(:),x(:),'poly1');

	E_leak(i,3) = ff.p2;
	g_leak(i,3) = 1/ff.p1;
	if dynamic_cutoff > 2
		ci = confint(ff);
		E_leak_err(i,3) = diff(ci(:,2))/2;
		g_leak_err(i,3) = abs(diff(1./ci(:,1)))/2;
	end

	subplot(2,3,2); hold on
	if ismember(i, [1 3 9 12])
		plot(V_space,r2,'+-')
	else
		plot(NaN,r2,'+-') % preserves color order
	end

	subplot(2,3,3); hold on
	plot(V_space,slopes,'+-')
	

end

subplot(2,3,2); hold on
xlabel('Membrane potential (mV)')
ylabel('r^2 of Ohmic fit')

subplot(2,3,3); hold on
xlabel('Membrane potential (mV)')
ylabel('Slope of Ohmic fit (nS)')
set(gca,'YScale','linear','XLim',[-80 80],'YLim',[0 130])

for k = 1:length(cutoff)
	subplot(2,3,3 + k); hold on
	for i = 1:max(cell_id)
		errorbar(E_leak(i,k),g_leak(i,k),g_leak_err(i,k),g_leak_err(i,k),E_leak_err(i,k),E_leak_err(i,k),'o')
	end
	set(gca,'XLim',[-80 80],'YLim',[0 130])
	xlabel('E_{leak} (mV)')
	ylabel('g_{leak} (nS)')
	title(['-80 to ' oval(V_space(cutoff(k))) ' mV'])
end

title('Dynamic V range')

prettyFig();

labelFigure('x_offset',-.01,'y_offset',-.0,'font_size',28)




%%
% However, I later realized that this algorithm tends to systematically overestimate the leak conductances. I default back to using the region from -80 mV to -60 mV to estimating the leak currents. Now, I apply this method to every trace in the dataset, and estimate leak currents for every neuron. 

g_leak = NaN(max(cell_id),3,5);
g_leak_err = NaN(max(cell_id),3,5);
E_leak = NaN(max(cell_id),3,5);
E_leak_err = NaN(max(cell_id),3,5);


for i = 1:max(genotype)
	for j = 1:length(Ca_space)
		for k = 1:max(cell_id)
			this_I = I_inf(cell_id == k & genotype == i & Ca_levels == Ca_space(j));
			if ~isempty(this_I)
				[g_leak(k,i,j), E_leak(k,i,j), g_leak_err(k,i,j), E_leak_err(k,i,j)] = findLeakCurrent2(this_I,V_space);
			end
		end
	end
end



%%
% In this figure, I look at how the leak currents of individual cells vary with Calcium concentration. Each figure corresponds to data from a single cell. The top row shows HEK cells w/o EAG, the second row HEK cells with WT EAG, and the third row shows data from HEK cells with mutant EAG. 

figure('outerposition',[0 0 1100 1111],'PaperUnits','points','PaperSize',[1100 1111]); hold on

c = parula(6); % for Calcium colors

for i = 1:max(genotype)

	gL = squeeze(g_leak(:,i,:));
	gL_err = squeeze(g_leak_err(:,i,:));
	EL = squeeze(E_leak(:,i,:));
	EL_err = squeeze(E_leak_err(:,i,:));

	for j = 1:3
		subplot(3,3,(i-1)*3 + j); hold on
		for k = 1:5
			errorbar(EL(j,k),gL(j,k),gL_err(j,k),gL_err(j,k),EL_err(j,k),EL_err(j,k),'o','Color',c(k,:))
		end
		set(gca,'XLim',[-100 10],'YLim',[0 20])
		xlabel('E_{leak} (mV)')
		ylabel('g_{leak} (nS)')

		switch i 
		case 1
			title(['No EAG #' oval(j)])
		case 2
			title(['WT EAG #' oval(j)])
		case 3
			title(['Mut EAG #' oval(j)])
		end
	end
end

prettyFig();

if being_published
	snapnow
	delete(gcf)
end


%%
% Now I compare all cells across genotypes and group by Calcium concentration. In the following figure, each plot groups all cells together by the intracellular Calcium concentration. HEK cells without EAG are shown in black, HEK cells with WT EAG are shown in blue, and HEK cells with mutant EAG are shown in red. Each dot is a single cell, and the error bars are 95% confidence intervals. There doesn't seem to be a very strong pattern here.  

L = {'13nM','100nM','10uM','100uM','500uM'};

figure('outerposition',[0 0 1201 901],'PaperUnits','points','PaperSize',[1201 901]); hold on
for i = 1:5
	subplot(2,3,i); hold on
	l1 = errorbar(E_leak(:,1,i),g_leak(:,1,i),g_leak_err(:,1,i),g_leak_err(:,1,i),E_leak_err(:,1,i),E_leak_err(:,1,i),'o','Color','k');
	l2 = errorbar(E_leak(:,2,i),g_leak(:,2,i),g_leak_err(:,2,i),g_leak_err(:,2,i),E_leak_err(:,2,i),E_leak_err(:,2,i),'o','Color','b');
	l3 = errorbar(E_leak(:,3,i),g_leak(:,3,i),g_leak_err(:,3,i),g_leak_err(:,3,i),E_leak_err(:,3,i),E_leak_err(:,3,i),'o','Color','r');

	title(['Ca = ' L{i}])
	xlabel('E_{leak} (mV)')
	ylabel('g_{leak} (nS)')
	set(gca,'XLim',[-100 80],'YLim',[0 30])
end

l = legend([l1,l2,l3],{'No EAG','WT EAG','Mut EAG'});
l.Position = [0.6561 0.3116 0.0883 0.0785];
prettyFig();

if being_published
	snapnow
	delete(gcf)
end


% subtract the estimated leak currents from the raw current traces 
fixed_I = all_I;
good_cells = g_leak*0;
for i = 1:3
	for j = 1:5
		for k = 1:max(cell_id)
			these_currents = all_I(:,cell_id == k & genotype == i & Ca_levels == Ca_space(j));

			if size(these_currents,2) == 0
				continue
			end

			fixed_currents = these_currents;

			m = g_leak(k,i,j);
			c = -E_leak(k,i,j)*g_leak(k,i,j);

			for l = 1:33
				estimated_leak_I = m*(V_space(l)) + c;
				fixed_currents(:,l) = these_currents(:,l) - estimated_leak_I;
			end

			allowed_err = Inf;

			if g_leak_err(k,i,j)/g_leak(k,i,j) < allowed_err & E_leak_err(k,i,j)/E_leak(k,i,j) < allowed_err & g_leak(k,i,j) > 0
				fixed_I(:,cell_id == k & genotype == i & Ca_levels == Ca_space(j)) = fixed_currents;

				good_cells(k,i,j) = 1;

			end
		end
	end
end


old_I  = all_I;
all_I = fixed_I;

% convert all_I into conductances
g = all_I;
E_K = -106.17; % mV
for i = 1:size(all_I,2)
	g(:,i) = all_I(:,i)./(all_V(i) - E_K);
end

% compute asymptotic conductances for all traces
g_inf = mean(g(4e3:5e3,:));

%%
% In the following figure, I plot all the current traces in current clamp mode for the lowest Calcium case (13 nM). For now, I neglect the Calcium dependence and focus on the voltage dependence. The top row shows data from HEK cells without the EAG channel, and the bottom row shows data from HEK cells with the EAG channel. (a) shows the current traces for a single HEK cell without EAG channels at 13nM Calcium. (b) is the same data, but the currents are converted into conductances. (c) Shows the conductances for all the cells in the dataset. (d-f) is in the same format, but for cells that do express EAG. Currents here have been leak subtracted on a per-cell, per-condition basis as previously discussed. 

figure('outerposition',[0 0 1501 901],'PaperUnits','points','PaperSize',[1501 901]); hold on

% first show the cells with no EAG
these_traces = find(cell_id == 1 & Ca_levels == 13e-9 & genotype == 1);

c = parula(34);
time = ((1:length(all_I)) - 51)/10;

subplot(2,3,1); hold on
for i = 1:length(these_traces)
	idx = these_traces(i);
	plot(time,all_I(:,idx),'Color',c(find(V_space == all_V(idx)),:))
end
xlabel('Time (ms)')
ylabel('Current')
set(gca,'YLim',[0 2000])
title('w/o EAG, 1 cell')

subplot(2,3,2); hold on
for i = 1:length(these_traces)
	idx = these_traces(i);
	plot(time,g(:,idx),'Color',c(find(V_space == all_V(idx)),:))
end
xlabel('Time (ms)')
ylabel('Conductance')
set(gca,'YLim',[0 100])
title('w/o EAG, 1 cell')

% now show the conductance-voltage curves for all cells 
subplot(2,3,3); hold on
for i = 1:max(cell_id)
	y = g_inf(cell_id == i & Ca_levels == 13e-9 & genotype == 1);
	if ~isempty(y)
		plot(V_space,y)
	end
end
xlabel('Membrane potential (mV)')
ylabel('Conductance')
title('w/o EAG, all cells')
set(gca,'YLim',[0 100])

% now show the cells with EAG
these_traces = find(cell_id == 1 & Ca_levels == 13e-9 & genotype == 2);

c = parula(34);
time = ((1:length(all_I)) - 51)/10;

subplot(2,3,4); hold on
for i = 1:length(these_traces)
	idx = these_traces(i);
	plot(time,all_I(:,idx),'Color',c(find(V_space == all_V(idx)),:))
end
xlabel('Time (ms)')
ylabel('Current')
set(gca,'YLim',[0 2000])
title('with EAG, 1 cell')

subplot(2,3,5); hold on
for i = 1:length(these_traces)
	idx = these_traces(i);
	plot(time,g(:,idx),'Color',c(find(V_space == all_V(idx)),:))
end
xlabel('Time (ms)')
ylabel('Conductance')
set(gca,'YLim',[0 100])
title('with EAG, 1 cell')

% now show the conductance-voltage curves for all cells 

subplot(2,3,6); hold on
for i = 1:max(cell_id)
	y = g_inf(cell_id == i & Ca_levels == 13e-9 & genotype == 2);
	if ~isempty(y)
		plot(V_space,y,'DisplayName',oval(i))
	end
end
xlabel('Membrane potential (mV)')
ylabel('Conductance')
title('with EAG, cells')
set(gca,'YLim',[0 100])

prettyFig();

labelFigure('x_offset',-.01,'y_offset',-.0,'font_size',28)

if being_published
	snapnow
	delete(gcf)
end


%%
% Notice that there is one massive outlier in the WT EAG cells -- for whatever reason, this cell has much bigger conductances for every membrane potential. I'm going to remove this from the dataset. 

rm_this = find(cell_id < 3 & genotype == 2); % 1 & 2 are bad
genotype(rm_this) = [];
all_I(:,rm_this) = [];
Ca_levels(rm_this) = [];
cell_id(rm_this) = [];
all_V(rm_this) = [];
g(:,rm_this) = [];
g_inf(rm_this) = [];
	 

%%
% Now I average the conductances of the cells expressing EAG and those that don't, and subtract the latter from the former to get a conductance curve for just EAG. This lets me estimate the gating function of the EAG channel, exactly like I did in the simulated data. In the following figure, (a-e) show estimated asymptotic gating variables as a function of membrane potential. (a-e) differ in the assumed exponent $p$. In (a-e), data is shown as black crosses, and the red lines are Boltzmann fits to the data. In (f), I plot the goodness-of-fit vs. the exponent $p$, and see that there is a clear minimum at $p = 2$, suggesting that this is the exponent. 

temp = g_inf(Ca_levels == 13e-9 & genotype == 2);
g_with_EAG = mean(reshape(temp,33,length(temp)/33),2);
temp = g_inf(Ca_levels == 13e-9 & genotype == 1);
g_wo_EAG = mean(reshape(temp,33,length(temp)/33),2);
g_EAG = g_with_EAG - g_wo_EAG;

boltz = @(A,B,x) (1./(1 + exp((x+A)./B)));

figure('outerposition',[0 0 1301 997],'PaperUnits','points','PaperSize',[1301 997]); hold on

% now plot the conductances, with various roots 
p = [1:8];
r2 = NaN*p;
for i = 1:length(p)
	this_g = g_EAG;
	this_g = this_g - min(this_g);
	this_g = real(this_g.^(1/p(i)));
	this_g = this_g/this_g(end);

	ff = fit(vectorise(V_space(2:end)),vectorise(this_g(2:end)),boltz,'StartPoint',[50, -20]);

	if i < 6
		subplot(2,3,p(i)); hold on
		plot(V_space,this_g,'k+')
		ylabel(['g^{1/' oval(p(i)) '} (norm)'])
		xlabel('Voltage (mV)')
		plot(V_space,ff(V_space),'r')
	end

	r2(i) = rsquare(ff(V_space(2:end)),this_g(2:end));

end

% use p = 2 because I think that's correct, and it's close to the minimum
this_g = real(g_EAG.^(1/2));
this_g = this_g/this_g(end);

mInfEAG = fit(vectorise(V_space(2:end)),vectorise(this_g(2:end)),boltz,'StartPoint',[12.3, -11.8]);

subplot(2,3,6); hold on
plot(p,1-r2)
xlabel('p')
ylabel('1 - r^2')

suptitle('WT EAG channels -- gating functions')
prettyFig();

labelFigure('x_offset',-.01,'y_offset',-.01,'font_size',28)

if being_published
	snapnow
	delete(gcf)
end


%% Mutated EAG Channels: gating functions
% Now I do the same thing for the mutated EAG channels, and characterise their gating functions. Once again, I get $p = 2$, and the gating function looks very similar. 


temp = g_inf(Ca_levels == 13e-9 & genotype == 3);
g_with_mut_EAG = mean(reshape(temp,33,length(temp)/33),2);
temp = g_inf(Ca_levels == 13e-9 & genotype == 1);
g_wo_EAG = mean(reshape(temp,33,length(temp)/33),2);
g_mut_EAG = g_with_mut_EAG - g_wo_EAG;
g_mut_EAG(g_mut_EAG<0) = 0;
boltz = @(A,B,x) (1./(1 + exp((x+A)./B)));

figure('outerposition',[0 0 1301 997],'PaperUnits','points','PaperSize',[1301 997]); hold on

% now plot the conductances, with various roots 
p = [1:8];
r2 = NaN*p;
for i = 1:length(p)
	this_g = g_mut_EAG.^(1/p(i));
	this_g = this_g/this_g(end);

	ff = fit(vectorise(V_space(2:end)),vectorise(this_g(2:end)),boltz,'StartPoint',[12.3, -11.8]);

	if i < 6
		subplot(2,3,p(i)); hold on
		plot(V_space,this_g,'k+')
		ylabel(['g^{1/' oval(p(i)) '} (norm)'])
		xlabel('Voltage (mV)')
		plot(V_space,ff(V_space),'r')
	end

	r2(i) = rsquare(ff(V_space(2:end)),this_g(2:end));

end


this_g = g_mut_EAG.^(1/2);
this_g = this_g/this_g(end);

mInfMutEAG = fit(vectorise(V_space(2:end)),vectorise(this_g(2:end)),boltz,'StartPoint',[12.3, -11.8]);

subplot(2,3,6); hold on
plot(p,1-r2)
xlabel('p')
ylabel('1 - r^2')

suptitle('Mutated EAG channels -- gating functions')
prettyFig();

labelFigure('x_offset',-.01,'y_offset',-.01,'font_size',28)

if being_published
	snapnow
	delete(gcf)
end

%%
% Now I compare the gating functions for the WT and the mutated EAG channels. They look almost identical. This makes sense, because the mutation in the EAG channel should affect Calcium binding, not activation. 

figure('outerposition',[200 200 601 601],'PaperUnits','points','PaperSize',[1200 601]); hold on

plot(V_space,mInfMutEAG(V_space),'r')
plot(V_space,mInfEAG(V_space),'b')
legend({'Mutated EAG','WT EAG'},'Location','southeast')
xlabel('Membrane potential (mV)')
ylabel('m_{Inf}')

prettyFig();

if being_published
	snapnow
	delete(gcf)
end

% save this for later use
save('EAG_activation_functions.mat','mInfMutEAG','mInfEAG');


%% Timescale dependence on voltage
% In this section, I look at how the timescale of activation of the EAG current depends on the voltage. First, I construct time series of EAG conductances by subtracted cell-averaged currents in HEK cells without EAG from cell-averaged currents in HEK cells with EAG, and then converting those corrects into conductances. These conductances are shown below: (a) shows these conductances for WT EAG channels, and (b) shows the same conductances for mutant EAG channels. Colors indicate holding voltage (blue = more hyperpolarized). Note that these raw curves look quite similar. 


% make conductance time series for the WT and mut EAG
temp = g(:,Ca_levels == 13e-9 & genotype == 1);
temp = reshape(temp,size(temp,1),33,size(temp,2)/33);
g_without_EAG =  squeeze(mean(temp,3));

temp = g(:,Ca_levels == 13e-9 & genotype == 2);
temp = reshape(temp,size(temp,1),33,size(temp,2)/33);
g_with_WT_EAG =  squeeze(mean(temp,3));

temp = g(:,Ca_levels == 13e-9 & genotype == 3);
temp = reshape(temp,size(temp,1),33,size(temp,2)/33);
g_with_mut_EAG =  squeeze(mean(temp,3));

g_EAG = g_with_WT_EAG - g_without_EAG;
g_mut_EAG = g_with_mut_EAG - g_without_EAG;

% normalize these conductances 
for i = 1:length(V_space)
	g_EAG(:,i) = g_EAG(:,i)/mean(g_EAG(4900:5000,i));
	g_mut_EAG(:,i) = g_mut_EAG(:,i)/mean(g_mut_EAG(4900:5000,i));
end

% show these
figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1); hold on
for i = 10:5:length(V_space)
	y = g_EAG(:,i);
	y(1:100) = NaN;
	plot(time,y,'Color',c(i,:))
end
set(gca,'YLim',[0 1.1])
xlabel('Time (ms)')
ylabel('Conductance (norm)')
title('WT EAG')

subplot(1,2,2); hold on
for i = 10:5:length(V_space)
	y = g_mut_EAG(:,i);
	y(1:100) = NaN;
	plot(time,y,'Color',c(i,:))
end
set(gca,'YLim',[0 1.1])
xlabel('Time (ms)')
ylabel('Conductance (norm)')
title('Mutant EAG')
prettyFig();

labelFigure('x_offset',-.01,'y_offset',-.04,'font_size',28)

if being_published
	snapnow
	delete(gcf)
end

%%
% Now, I fit exponentials to these time series, after taking their inverse power for a integers up to 8. I then compute the mean goodness of fit for all these fits, and plot that as a function of the exponent $p$. Note that the best fit is for $p=2$, which is what I use throughout to model these channels. 

% fit a timescale to every time trace, for every integer power
p = 1:8;
tau_mut = NaN(length(V_space),length(p));
tau_wt = NaN(length(V_space),length(p));
tau_r2_mut = NaN(length(V_space),length(p));
tau_r2_wt = NaN(length(V_space),length(p));
f = fittype('1 - exp(-x./tau)');


for i = 1:length(V_space)
	for j = 1:length(p)
		% first do WT
		y = g_EAG(101:end-1,i);
		y = y.^(1/p(j));
		x = time(101:end-1);
		if isreal(y)
			[ff,gof] = fit(x(:),y(:),f,'StartPoint',10,'Lower',0);
			tau_wt(i,j) = ff.tau;
			tau_r2_wt(i,j) = gof.rsquare;
		else
			tau_r2_wt(i,j) = -1;
		end

		% do mutant
		y = g_mut_EAG(101:end-1,i);
		y = y.^(1/p(j));
		x = time(101:end-1);
		if isreal(y)
			[ff,gof] = fit(x(:),y(:),f,'StartPoint',10,'Lower',0);
			tau_mut(i,j) = ff.tau;
			tau_r2_mut(i,j) = gof.rsquare;
		else
			tau_r2_wt(i,j) = -1;
		end
	end
end


figure('outerposition',[200 200 601 601],'PaperUnits','points','PaperSize',[1200 601]); hold on
for i = 1:length(p)
	l1 = plot(i,mean(tau_r2_wt(tau_r2_wt(:,i) > 0,i)),'k+');
end
set(gca,'XLim',[0 8.5],'YLim',[0 1])
xlabel('p')
ylabel('<r^2> of exponential fit')

for i = 1:length(p)
	l2 = plot(i,mean(tau_r2_mut(tau_r2_mut(:,i) > 0,i)),'r+');
end
legend([l1,l2],{'WT EAG','Mutant EAG'})

prettyFig();

if being_published
	snapnow
	delete(gcf)
end

%%
% Now I plot the timescales from this fit as a function of membrane potential, keeping only timescales where the goodness of exponential fit was > 0.8. I do this for various values of $p$, for both the mutant (red) and the WT (black) EAG channels. In each plot below, the black circles are the timescales of exponential fits to time series of estimated gating variable for the WT EAG channels, and the red crosses are timescales of exponential fits to to time series of estimated gating variable for the mutant EAG channels. The dashed lines are functional fits to the crosses or the circles, using standard functions as in Liu et al. or Prinz et al. Note that these channels appear to be much slower than the Potassium channels I simulated at the very beginning.

r2_cutoff = .5;
tau_r2_mut(isnan(tau_r2_mut)) = -1;
figure('outerposition',[0 0 1001 901],'PaperUnits','points','PaperSize',[1001 901]); hold on
for i = 1:4
	subplot(2,2,i); hold on
	y = tau_wt(:,i);
	x = V_space;
	x(tau_r2_wt(:,i) < r2_cutoff) = [];
	y(tau_r2_wt(:,i) < r2_cutoff) = [];
	plot(x,y,'ko')

	% fit functions to this
	ff = fit(x(:),y(:),tauX,'StartPoint',[7 6 30 -20]);
	plot(V_space,ff(V_space),'k--')

	tauEAG = ff;

	y = tau_mut(:,i);
	x = V_space;
	x(tau_r2_mut(:,i) < r2_cutoff) = [];
	y(tau_r2_mut(:,i) < r2_cutoff) = [];
	plot(x,y,'r+')

	% fit functions to this
	ff = fit(x(:),y(:),tauX,'StartPoint',[7 6 30 -20]);
	plot(V_space,ff(V_space),'r--')

	title(['p = ' oval(i)])

	set(gca,'YLim',[0 200])
	xlabel('Membrane potential (mV)')
	ylabel('\tau_m (ms)')

	if i == 1
		clear l
		l(1) = plot(NaN,NaN,'k');
		l(2) = plot(NaN,NaN,'r');
		legend(l,{'WT EAG','Mut EAG'},'location','northeast')
	end

	% save these for later use
	if i == 2
		save('EAG_activation_functions.mat','tauEAG','-append');
	end
end

prettyFig();

labelFigure('x_offset',-.01,'y_offset',-.0,'font_size',28)

if being_published
	snapnow
	delete(gcf)
end


%% Calcium dependence of gating variable
% In this section, I characterize how these channels depend on intracellular Calcium concentration. In the following figure, I plot conductance of WT and mutant channels as a function of membrane potential, for various intracellular Calcium concentrations. Note that at $500\mu M$, the mutant channels have a greatly decreased conductance compared to the WT channels, but they're not that different at other Calcium concentrations. 

% the trick here is to only use data from cells that have every Calcium level.
% it is an unfair comparison if some cells have some calcium levels, 
% but not all 

all_ca_levels = unique(Ca_levels(find(genotype == 1))); % because control doesn't have all cal levels

good_ctrl_cells = [];
good_wt_cells = [];
good_mut_cells = [];

for i = 1:max(cell_id)
	if all(ismember(all_ca_levels,unique(Ca_levels(find(cell_id == i & genotype == 1)))))
		good_ctrl_cells = [good_ctrl_cells; i];
	end
	if all(ismember(all_ca_levels,unique(Ca_levels(find(cell_id == i & genotype == 2)))))
		good_wt_cells = [good_wt_cells; i];
	end
	if all(ismember(all_ca_levels,unique(Ca_levels(find(cell_id == i & genotype == 3)))))
		good_mut_cells = [good_mut_cells; i];
	end
end


g_inf_0 = zeros(length(V_space),length(Ca_space));
g_inf_wt = zeros(length(V_space),length(Ca_space));
g_inf_mut = zeros(length(V_space),length(Ca_space));
for i = 1:length(Ca_space)
	temp = g_inf(ismember(cell_id,good_ctrl_cells) & Ca_levels == Ca_space(i) & genotype == 1);
	temp = mean(reshape(temp,33,length(temp)/33),2);
	g_inf_0(:,i) = temp;

	temp = g_inf(ismember(cell_id,good_wt_cells) & Ca_levels == Ca_space(i) & genotype == 2);
	temp = mean(reshape(temp,33,length(temp)/33),2);
	g_inf_wt(:,i) = temp;

	temp = g_inf(ismember(cell_id,good_mut_cells) & Ca_levels == Ca_space(i) & genotype == 3);
	temp = mean(reshape(temp,33,length(temp)/33),2);
	g_inf_mut(:,i) = temp;
end

g_inf_wt = g_inf_wt - g_inf_0;
g_inf_mut = g_inf_mut - g_inf_0;

figure('outerposition',[20 200 1501 500],'PaperUnits','points','PaperSize',[1501 500]); hold on
subplot(1,3,1); hold on
c = parula(6);
clear l
for i = 1:length(Ca_space)
	l(i) = plot(V_space,g_inf_wt(:,i),'Color',c(i,:));
end
L = {'13nM','100nM','10uM','100uM','500uM'};
legend(l,L,'location','northwest')
title('WT EAG')
xlabel('Membrane potential (mV)')
ylabel('Conductance')

subplot(1,3,2); hold on
c = parula(6);
for i = 1:length(Ca_space)
	plot(V_space,g_inf_mut(:,i),'Color',c(i,:))
end
title('Mutated EAG')
xlabel('Membrane potential (mV)')
ylabel('Conductance')

subplot(1,3,3); hold on
c = parula(6);
for i = 1:length(Ca_space)
	y_wt = g_inf_wt(:,i)/max(max(g_inf_wt));
	y_mut = g_inf_mut(:,i)/max(max(g_inf_mut));
	plot(V_space,y_wt,'Color',c(i,:))
	plot(V_space,y_mut,'--','Color',c(i,:))
end
title('WT vs Mut (dashed)')
xlabel('Membrane potential (mV)')
ylabel('Conductance')



prettyFig();

if being_published
	snapnow
	delete(gcf)
end

%%
% Now I plot the conductance as a function of the Calcium concentration, for various holding potentials, for both the WT and mutant EAG channels. Note that the WT EAG channels do not have a monotonic relationship with Calcium. For some (unknown) reason, the conductance *increases* from 13 nM to 100 nM, and then decreases. The mutated EAG channels, on the other hand, have a monotonically decreasing conductance that depends on Calcium concentration. 

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1); hold on
c = parula(34);
for i = 1:3:length(V_space)
	x = Ca_space([1:3 5]);
	y = g_inf_wt(i,[1:3 5]);
	plot(x,y,'-+','Color',c(i,:))
end
set(gca,'XScale','log','XTick',Ca_space(([1:3 5])),'XTickLabel',L([1:3 5]))
xlabel('[Ca^2^+] (M)')
ylabel('Conductance')
title('WT EAG')

subplot(1,2,2); hold on
for i = 1:3:length(V_space)
	x = Ca_space([1:3 5]);
	y = g_inf_mut(i,[1:3 5]);
	plot(x,y,'-+','Color',c(i,:))
end
set(gca,'XScale','log','XTick',Ca_space(([1:3 5])),'XTickLabel',L([1:3 5]))
xlabel('[Ca^2^+] (M)')
ylabel('Conductance')
title('Mut EAG')

prettyFig();

if being_published
	snapnow
	delete(gcf)
end


%%
% Now, I fit some simple functions to the conductances at the highest holding potential (because that is the easiest to interpret). Using these functions, I can capture the fact that mutated EAG channels are less inhibited by $500 \mu M$ Calcium than WT EAG channels, but I neglect the non-monotonic nature of the WT EAG channels, and the corresponding relative decrease in WT EAG channels compared to mutated channels at very low Calcium concentrations. Note that the effect of the mutation is quite subtle: a small rightward shift (loss in sensitivity to Calcium) together with a small offset (persistent activity even at very high Calcium).

f = fittype('A.*(1-B)./(x+A) + B');
M = max([g_inf_mut(:); g_inf_wt(:)]);

figure('outerposition',[0 0 1501 500],'PaperUnits','points','PaperSize',[1501 500]); hold on
subplot(1,3,1); hold on

c = parula(length(V_space)+1);

all_x = []; all_y = [];
for i = 33:1:length(V_space)
	x = Ca_space([1:3 5]);
	y = g_inf_wt(i,[1:3 5]);
	y = y/M;
	y = y.^2; % assuming p = 2
	plot(x,y,'k+')
	all_x = [all_x(:); x(:)];
	all_y = [all_y(:); y(:)];
end

set(gca,'XScale','log','XTick',Ca_space(([1:3 5])),'XTickLabel',L([1:3 5]),'YLim',[0 1.5])
xlabel('[Ca^2^+] (M)')
ylabel('m_{Inf}')
title('WT EAG')

ff = fit(all_x(:),all_y(:),f,'StartPoint',[1 .1],'Lower',[0 0],'Upper',[1 0]);
temp = logspace(log10(min(Ca_space)),log10(max(Ca_space)),1e3);
plot(temp,ff(temp),'k')

Ca_dep = ff;
save('EAG_activation_functions.mat','Ca_dep','-append');

subplot(1,3,3); hold on
plot(temp,ff(temp),'k')

subplot(1,3,1); hold on
t =  ['$\frac{' oval((ff.A)*1e6) '\mu M}{[Ca^{2+}]+ ' oval((ff.A)*1e6) ' \mu M}+ ' oval(ff.B) '$'];
h = text(1e-7,1.1,t,'interpreter','latex','FontSize',24,'Color','r');

subplot(1,3,2); hold on
all_x = []; all_y = [];
for i = 33:1:length(V_space)
	x = Ca_space([1:3 5]);
	y = g_inf_mut(i,[1:3 5]);
	y = y/M;
	y = y.^2; % assuming p = 2
	plot(x,y,'r+')
	all_x = [all_x(:); x(:)];
	all_y = [all_y(:); y(:)];
end

set(gca,'XScale','log','XTick',Ca_space(([1:3 5])),'XTickLabel',L([1:3 5]),'YLim',[0 1.5]')
xlabel('[Ca^2^+] (M)')
ylabel('m_{Inf}')
title('Mut EAG')

ff = fit(all_x(:),all_y(:),f,'StartPoint',[1 .1],'Lower',[0 min(all_y)]);
temp = logspace(log10(min(Ca_space)),log10(max(Ca_space)),1e3);
plot(temp,ff(temp),'r')

Ca_dep_mut = ff;
save('EAG_activation_functions.mat','Ca_dep_mut','-append');

t =  ['$\frac{' oval((ff.A)*1e6) '\mu M}{[Ca^{2+}]+ ' oval((ff.A)*1e6) ' \mu M}+ ' oval(ff.B) '$'];
h = text(1e-7,1.1,t,'interpreter','latex','FontSize',24,'Color','r');

subplot(1,3,3); hold on
plot(temp,ff(temp),'r')
legend({'WT','Mutated'})
set(gca,'XScale','log','XTick',Ca_space(([1:3 5])),'XTickLabel',L([1:3 5]),'YLim',[0 1.5]')

prettyFig();

if being_published
	snapnow
	delete(gcf)
end




%% Version Info
%
pFooter;














