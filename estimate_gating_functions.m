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

n.set('g_Na',0)
n.set('g_CaT',0)
n.set('g_CaS',0)
n.set('g_leak',0)
n.set('g_A',0)
n.set('g_KCa',0)
n.set('g_H',0)

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


%% Estimating voltage dependence of gating variable 
% To model EAG channels, I need to understand how these channels vary their activity as a function of both voltage and Calcium concentrations. These channels don't inactivate, so their behavior can be described by a single gating variable ($m$). In Hodgkin-Huxley terms, we want to determine how the steady state value of $m$ depends on the voltage and Calcium, and how the timescale of m changes with voltage and Calcium. 

%%
% First, I estimate the voltage dependence on the gating variable by looking at currents in voltage clamp in HEK cells.In this data, EAG channels were inserted into HEK cells, and currents recorded from them in voltage clamp. The experiment was repeated with HEK cells without EAG, so subtracting one from the other should give us "pure" EAG currents. 


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

% remove baseline from every current trace
for i = 1:size(all_I,2)
	all_I(:,i) = all_I(:,i) - mean(all_I(1:50,i));
end


% convert all_I into conductances
g = all_I;
E_K = -106.17; % mV
for i = 1:size(all_I,2)
	g(:,i) = all_I(:,i)./(all_V(i) - E_K);
end

% compute asymptotic conductances for all traces
g_inf = mean(g(4e3:5e3,:));

%%
% In the following figure, I plot all the current traces in current clamp mode for the lowest Calcium case (13 nM). For now, I neglect the Calcium dependence and focus on the voltage dependence. (a) shows the current traces for a single cell that expressed EAG channels at 13nM Calcium. (b) is the same, but from a cell that doesn't express EAG channels. Note that the amplitudes and the kinetics of the currents are different. Different colors in (a-b) are different holding potentials. (c) and (d) show the I-V relationships for all the data at 13 nM Calcium. Each line corresponds to a different cell. 

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

subplot(2,3,2); hold on
for i = 1:length(these_traces)
	idx = these_traces(i);
	plot(time,g(:,idx),'Color',c(find(V_space == all_V(idx)),:))
end
xlabel('Time (ms)')
ylabel('Conductance')
set(gca,'YLim',[0 100])

% now show the conductance-voltage curves for all cells 
clear ax1 ax2
ax1 = subplot(2,3,3); hold on
for i = 1:max(cell_id)
	y = g_inf(cell_id == i & Ca_levels == 13e-9 & genotype == 1);
	if ~isempty(y)
		plot(ax1,V_space,y)
	end
end
xlabel(ax1,'Membrane potential (mV)')
ylabel(ax1,'Conductance')
title(ax1,'All cells')
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

subplot(2,3,5); hold on
for i = 1:length(these_traces)
	idx = these_traces(i);
	plot(time,g(:,idx),'Color',c(find(V_space == all_V(idx)),:))
end
xlabel('Time (ms)')
ylabel('Conductance')
set(gca,'YLim',[0 100])

% now show the conductance-voltage curves for all cells 
clear ax1 ax2
ax1 = subplot(2,3,6); hold on
for i = 1:max(cell_id)
	y = g_inf(cell_id == i & Ca_levels == 13e-9 & genotype == 2);
	if ~isempty(y)
		plot(ax1,V_space,y)
	end
end
xlabel(ax1,'Membrane potential (mV)')
ylabel(ax1,'Conductance')
title(ax1,'All cells')
set(gca,'YLim',[0 100])

prettyFig();


if being_published
	snapnow
	delete(gcf)
end

%%
% Now I average the conductances of the cells expressing EAG and those that don't, and subtract the latter from the former to get a conductance curve for just EAG. 

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
	this_g = g_EAG.^(1/p(i));
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

% use p = 4 because I think that's correct, and it's close to the minimum
this_g = g_EAG.^(1/4);
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


return


%% Timescale dependence on voltage
% In this section, I look at how the timescale of activation of the EAG current depends on the voltage. To estimate this, I fit exponentials to the rising phases of the EAG currents, and plot these timescales as a function of voltage. (a) shows the time traces of the mean EAG current as a function of voltage (colors). For each trace, I fit an exponential, and retain fits only when the r^2 of the fit is > 0.9. (b) shows the estimated timescales vs. voltage, together with a function I fit to this. I can't estimate the timescales at very low voltages, so we have to extrapolate there. 


V = linspace(-80,80,33);
t = linspace(0,500,5000);
I_WT = NaN(length(V),12,5000);
I_control = NaN(length(V),12,5000);

% first do WT. 
load(joinPath(path_name,'WTData.mat'))

for i = 1:12
	for j = 1
		if ~isempty(WTFile{i,j})
			temp = WTFile{i,j}{1}; % lowest calcium conc.
			for k = 1:33
				I_WT(k,i,:) = temp(52:end,k);
			end
		end
	end
end

% average over all neurons 
I_WT = squeeze(mean(I_WT,2))';

clear temp
% do the controls
load(joinPath(path_name,'NoEAGData.mat'))

for i = 1:5
	for j = 1
		if ~isempty(NoEAGFile{i,j})
			temp = NoEAGFile{i,j}{1}; % lowest calcium conc.
			for k = 1:33
				I_control(k,i,:) = temp(52:end,k);
			end
		end
	end
end

% average over all neurons 
I_control = squeeze(nanmean(I_control,2))';

% subtract
I = I_WT - I_control;

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1); hold on
c = parula(34);
for i = 1:33
	plot(t,I(:,i),'Color',c(i,:))
end
set(gca,'XLim',[0 100])
ylabel('Mean EAG current')
xlabel('Time (ms)')
set(gca,'YLim',[0 4000])

% fit a timescale to each time trace
tau = NaN*V;
r2 = NaN*V;
f = fittype('1 - exp(-x./tau)');


for i = 1:length(V)
	y = I(:,i);
	y = y/mean(y(2e3:2.5e3));
	[~,idx] = min(abs(y));
	y = y(idx:2e3);
	x = t(idx:2e3);
	try
		[ff,gof] = fit(x(:),y(:),f,'StartPoint',1);
		tau(i) = ff.tau;
		r2(i) = gof.rsquare;
	catch
	end
end

c = colorbar();
title(c,'V (mV)')
caxis([-80 80])

% and fit a function to this 
subplot(1,2,2); hold on
f = fittype('A - B./(1+exp((x+D)./E))');

r2_cutoff = .9;

ff = fit(V(r2>r2_cutoff)',tau(r2>r2_cutoff)',f,'StartPoint',[1 1 1 1]);
plot(V,ff(V),'r')
plot(V(r2>r2_cutoff),tau(r2>r2_cutoff),'k+')
legend({'Fit','Data'})

xlabel('Voltage (mV)')
ylabel('Timescale (ms)')
set(gca,'YLim',[0 400])


prettyFig()

t =  '$-0.67+\frac{26500}{1+\exp\left(\frac{V+283}{46.03}\right)}$';
h = text(20,250,t,'interpreter','latex','FontSize',24);

labelFigure('x_offset',-.01,'y_offset',-.01,'font_size',24)

if being_published	
	snapnow	
	delete(gcf)
end


%% Calcium dependence of gating variable
% In this section, I look at how the steady state currents (a proxy for the gating variable) vary with Calcium concentration in HEK cells. First, I plot the steady state currents vs. Calcium concnetration for each cell, both for the HEK cells that express EAG (top row), and for the HEK cells that don't (bottom row). In each plot, I plot the raw current vs. the Calcium concentration. The different lines correspond to different holding potentials (colors). 

figure('outerposition',[0 0 1444 901],'PaperUnits','points','PaperSize',[1444 901]); hold on
for i = 1:4
	ax(i) = subplot(2,4,i); hold on
	set(ax(i),'XScale','log','YScale','log')
end

V = linspace(-80,80,33);
Imax_WT = NaN(length(V),12,5);
Imax_Control = NaN(length(V),5,5);
Ca = [13e-9 100e-9 10e-6 100e-6 500e-6]; % M

% first do WT. 
load(joinPath(path_name,'WTData.mat'))

for i = 1:12
	for j = 1:5
		if ~isempty(WTFile{i,j})
			I = WTFile{i,j}{1}; % lowest calcium conc.
			temp = mean(I(4e3:5e3,:));
			Imax_WT(:,i,j) = temp; 
		end
	end
end

c = parula(34);

for i = 1:4
	temp = squeeze(Imax_WT(:,i,:));
	for j = 1:3:33
		plot(ax(i),Ca,temp(j,:),'-+','Color',c(j,:))
		xlabel(ax(i),'Ca^2^+ (M)')
		ylabel(ax(i),'Current')
		set(ax(i),'YLim',[0 18000],'XTick',logspace(-9,-3,7),'YLim',[1 1e5],'YTick',logspace(1,5,5))
		title(ax(i),['Cell ' oval(i) ' (+EAG)'])
	end
end


% average across all cells for which we have all Calcium levels 
Imax_WT = squeeze(nanmean(Imax_WT(:,1:4,:),2));


% now get the control 
load(joinPath(path_name,'NoEAGData.mat'))

clear ax
for i = 1:3
	ax(i) = subplot(2,3,i+3); hold on
	set(ax(i),'XScale','log','YScale','log')
end

for i = 1:5
	for j = 1:5
		if ~isempty(NoEAGFile{i,j})
			I = NoEAGFile{i,j}{1}; 
			temp = mean(I(4e3:5e3,:));
			Imax_Control(:,i,j) = temp; 
		end
	end
end

c = parula(34);

for i = 1:3
	temp = squeeze(Imax_Control(:,i,:));
	for j = 1:3:33
		plot(ax(i),Ca,temp(j,:),'-+','Color',c(j,:))
		xlabel(ax(i),'Ca^2^+ (M)')
		ylabel(ax(i),'Current')
		set(ax(i),'YLim',[0 18000],'XTick',logspace(-9,-3,7),'YLim',[1 1e5],'YTick',logspace(1,5,5))
		title(ax(i),['Cell ' oval(i) ' (-EAG)'])
	end
end

h = colorbar;
caxis([-80 80])

title(h,'V (mV)')

prettyFig();

if being_published
	snapnow
	delete(gcf)
end


% average across all cells
Imax_Control = squeeze(nanmean(Imax_Control(:,1:3,:),2));


%%
% Now I average over all cells, and compare the curves between cells that have EAG and cells that don't. (a) shows the mean EAG currents as a function of the Calcium concentration -- they seem to decrease with increasing Calcium concentration. In (b), I plot the highest curve (the curve @ 80mV) and fit a simple function to this. I choose a Hill-like function, and it seems to a decent job at predicting how Calcium inhibits EAG currents. 

Imax = Imax_WT - Imax_Control;

% since the fourth column is empty, skip it
Imax(:,4) = [];

figure('outerposition',[0 0 1400 701],'PaperUnits','points','PaperSize',[1400 701]); hold on
subplot(1,2,1); hold on
for i = 1:5:33
	plot(Ca([1 2 3 5]),Imax(i,:),'-+','Color',c(i,:))
end
set(gca,'YLim',[-500 7000],'XTick',logspace(-9,-3,7),'YScale','linear','XScale','log')
xlabel('Ca^2^+ (M)')
ylabel('Mean EAG Current')

h = colorbar;
caxis([-80 80])

title(h,'V (mV)')

subplot(1,2,2); hold on
i = 31;
x = Ca([1 2 3 5]);
y = Imax(i,:);
y = y/max(y);

f = fittype('A./(x+A)');
ff = fit(x(:),y(:),f,'StartPoint',[1]);
xx = logspace(-9,1,1e4);
plot(xx,ff(xx),'r')

set(gca,'YLim',[0 1],'XTick',logspace(-9,-3,7),'XLim',[1e-9 1e-3],'XScale','log')
for i = 31
	y = Imax(i,:);
	y = y/max(y);
	plot(x,y,'+','Color',c(i,:),'MarkerSize',24)
end
ylabel('m_{EAG}')
xlabel('Ca^2^+ (M)')
legend({'Fit','Data @ 80mV'})

prettyFig();

t =  '$\frac{6\mu M}{[Ca^{2+}]+6\mu M}$';
h = text(1e-8,.5,t,'interpreter','latex','FontSize',24);


labelFigure('x_offset',-.05,'y_offset',-.01,'font_size',24)

if being_published
	snapnow
	delete(gcf)
end


%% How do activity timescales depend on Calcium? 
% In this section, I check to see if the timescale of EAG current activation depends on the Calcium concentration. 


V = linspace(-80,80,33);
t = linspace(0,500,5000);
I_WT = NaN(length(V),12,5000,j);
I_control = NaN(length(V),12,5000,j);

% first do WT. 
load(joinPath(path_name,'WTData.mat'))

for i = 1:12
	for j = 1:5
		if ~isempty(WTFile{i,j})
			temp = WTFile{i,j}{1}; % lowest calcium conc.
			for k = 1:33
				I_WT(k,i,:,j) = temp(52:end,k);
			end
		end
	end
end

% average over all neurons 
I_WT = squeeze(nanmean(I_WT,2));

clear temp
% do the controls
load(joinPath(path_name,'NoEAGData.mat'))

for i = 1:5
	for j = 1:5
		if ~isempty(NoEAGFile{i,j})
			temp = NoEAGFile{i,j}{1}; % lowest calcium conc.
			for k = 1:33
				I_control(k,i,:,j) = temp(52:end,k);
			end
		end
	end
end

% average over all neurons 
I_control = squeeze(nanmean(I_control,2));

% subtract
I = I_WT - I_control;

figure('outerposition',[0 0 1666 902],'PaperUnits','points','PaperSize',[1666 902]); hold on
clear ax ax2
% plot time traces 
example_V = 20:20:80;
for i = 1:4
	ax(i) = subplot(2,4,i); hold on
	ax2(i) = subplot(2,4,i+4); hold on
	xlabel(ax(i),'Time (ms)')
	xlabel(ax2(i),'Time (ms)')
	ylabel(ax(i),'Mean EAG current')
	ylabel(ax2(i),'Mean EAG current (norm)')
	set(ax(i),'YLim',[-500 6000])
	set(ax2(i),'YLim',[0 1.1])
end

c = parula(6);
L = {'13nM','100nM','10uM','100uM','500uM'};
clear l
for i = 1:4
	thisI = squeeze(I(find(V == example_V(i)),:,:));
	for j = 1:5
		y = thisI(:,j);
		l(j) = plot(ax(i),t,y,'Color',c(j,:));
		y = y/mean(y(4e3:4.5e3));
		plot(ax2(i),t,y,'Color',c(j,:))
		title(ax(i),['V = ' oval(example_V(i)) ' mV'])


	end
	if i == 1
		legend(l,L);
	end
end


prettyFig()

if being_published	
	snapnow	
	delete(gcf)
end

%%
% These traces don't look very different in their kinetics, so I will neglect the effect of Calcium on the kinetics. 

%% Comparison to Potassium Channels
% In this section, I compare the gating functions of the EAG channels that I estimate here to "normal" Potassium channels as in Liu et al. to get a sense of how different these channels are, and if they make sense. 



% define functions
boltz = @(V,A,B) (1./(1 + exp((V+A)./B)));
tauX = @(V,A,B,D,E) (A - B./(1+exp((V+D)./E)));

% m_inf
mNainf = @(V) (boltz(V,25.5,-5.29));
mKinf = @(V) boltz(V,12.3,-11.8);
mEAGinf = @(V,Ca) boltz(V,-22.28,-23.39).*(6./(6+Ca));


% tau_m
taumNa = @(V) tauX(V,1.32,1.26,120.0,-25.0);
taumK = @(V) tauX(V,7.2,6.4,28.3,-19.2);
taumEAG = @(V) tauX(V,-0.67,-26500,283,46.03);

V = linspace(-80,80,1e3);

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on

subplot(1,2,1); hold on
plot(V,mNainf(V))
plot(V,mKinf(V))
plot(V,mEAGinf(V,0.05))
ylabel('z_{Inf}')
xlabel('Membrane potential (mV)')

subplot(1,2,2); hold on
plot(V,taumNa(V))
plot(V,taumK(V))
plot(V,taumEAG(V))
ylabel('\tau (ms)')
xlabel('Membrane potential (mV)')
set(gca,'YScale','log')
legend({'NaV','K','EAG'})

prettyFig();

if being_published
	snapnow
	delete(gcf)
end

%% Version Info
%
pFooter;














