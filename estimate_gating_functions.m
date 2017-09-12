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
% In the following figure, I plot all the current traces in current clamp mode for the lowest Calcium case (13 nM). For now, I neglect the Calcium dependence and focus on the voltage dependence. The top row shows data from HEK cells without the EAG channel, and the bottom row shows data from HEK cells with the EAG channel. (a) shows the current traces for a single HEK cell without EAG channels at 13nM Calcium. (b) is the same data, but the currents are converted into conductances. (c) Shows the conductances for all the cells in the dataset. (d-f) is in the same format, but for cells that do express EAG. 

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
title(ax1,'w/o EAG, all cells')
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
clear ax1 ax2
ax1 = subplot(2,3,6); hold on
for i = 1:max(cell_id)
	y = g_inf(cell_id == i & Ca_levels == 13e-9 & genotype == 2);
	if ~isempty(y)
		plot(ax1,V_space,y,'DisplayName',oval(i))
	end
end
xlabel(ax1,'Membrane potential (mV)')
ylabel(ax1,'Conductance')
title(ax1,'with EAG, cells')
set(gca,'YLim',[0 100])

prettyFig();

 labelFigure('x_offset',-.01,'y_offset',-.0,'font_size',28)

if being_published
	snapnow
	delete(gcf)
end



%%
% Notice that there is one massive outlier in the WT EAG cells -- for whatever reason, this cell has much bigger conductances for every membrane potential. I'm going to remove this from the dataset. 

rm_this = find(cell_id < 3 & genotype == 2); % both 1 and 2 are bad
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

% use p = 2 because I think that's correct, and it's close to the minimum
this_g = g_EAG.^(1/2);
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

r2_cutoff = .8;
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














