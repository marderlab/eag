% this script estimates m_inf and tau_inf values for the EAG current
% from experimental data 



pHeader;


%% Estimating voltage dependence of gating variable 
% First, I estimate the voltage dependence on the gating variable by looking at currents in voltage clamp in HEK cells. 


path_name = '/Users/srinivas/data/calcium-microdomain/HEK';


%%
% In the following figure, I plot all the current traces in current clamp mode for the lowest Calcium case (13 nM). For now, I neglect the Calcium dependence and focus on the voltage dependence.  

figure('outerposition',[0 0 900 901],'PaperUnits','points','PaperSize',[900 901]); hold on
for i = 1:4
	ax(i) = subplot(2,2,i); hold on
end




V = linspace(-80,80,33);
Imax_WT = NaN(length(V),12,5);
Imax_Control = NaN(length(V),5,5);

% first do WT. 
load(joinPath(path_name,'WTData.mat'))

for i = 1:12
	for j = 1:5
		if ~isempty(WTFile{i,j})
			I = WTFile{i,j}{1}; % lowest calcium conc.
			if i == 1 & j == 1
				c = parula(34);
				t = linspace(0,500,5000);
				for k = 1:33
					plot(ax(1),t,I(52:end,k),'Color',c(k,:))
				end
			end
			temp = mean(I(4e3:5e3,:));
			Imax_WT(:,i,j) = temp; 
		end
	end
end
set(ax(1),'YLim',[0 2000])
xlabel(ax(1),'Time (ms)')
ylabel(ax(1),'Current')
title(ax(1),'HEK cells w EAG')

% show I-V curves for all cells 
for i = 1:12
	plot(ax(3),V,squeeze(Imax_WT(:,i,1)))
end
xlabel(ax(3),'Voltage (mV)')
ylabel(ax(3),'Current')
set(ax(3),'YLim',[-1000 20000])

% average across all cells
Imax_WT = squeeze(nanmean(Imax_WT,2));



% now get the control 
load(joinPath(path_name,'NoEAGData.mat'))

for i = 1:5
	for j = 1:5
		if ~isempty(NoEAGFile{i,j})
			I = NoEAGFile{i,j}{1}; 
			if i == 1 & j == 1
				c = parula(34);
				t = linspace(0,500,5000);
				for k = 1:33
					plot(ax(2),t,I(52:end,k),'Color',c(k,:))
				end
			end
			temp = mean(I(4e3:5e3,:));
			Imax_Control(:,i,j) = temp; 
		end
	end
end

title(ax(2),'HEK cells w/o EAG')
set(ax(2),'YLim',[0 2000])

% show I-V curves for all cells 
for i = 1:5
	plot(ax(4),V,squeeze(Imax_WT(:,i,1)))
end
xlabel(ax(4),'Voltage (mV)')
ylabel(ax(4),'Current')
set(ax(4),'YLim',[-1000 20000])

% average across all cells
Imax_Control = squeeze(nanmean(Imax_Control,2));

prettyFig();

labelFigure('x_offset',-.05,'y_offset',-.01,'font_size',24)

if being_published
	snapnow
	delete(gcf)
end


%%
% Now, I subtract the no EAG curves from the EAG curves to estimate the I-V relationship of the EAG conductance alone. 


figure('outerposition',[0 0 1501 500],'PaperUnits','points','PaperSize',[1501 500]); hold on

% subtract no EAG data from WT to get only EAG data
Imax = Imax_WT - Imax_Control;

subplot(1,3,1); hold on
plot(V,Imax_Control(:,1),'k')
plot(V,Imax_WT(:,1),'r')
xlabel('Voltage (mV)')
legend({'-EAG','+EAG'},'Location','northwest')
ylabel('Measured current')

subplot(1,3,2); hold on
plot(V,Imax(:,1),'k')
xlabel('Voltage (mV)')
ylabel('Estimated EAG current')

% normalize -- just for the lowest calcium case 
Imax = Imax - min(Imax(:,1));
Imax = Imax/max(Imax(:,1));

subplot(1,3,3); hold on
y = Imax(:,1);
plot(V,y,'k+')
f = fittype('1./(1  + exp((x+A)./(B)))');
ff = fit(V(:),y,f,'StartPoint',[-1 -10]);
xlabel('Voltage (mV)')
plot(V,ff(V),'r')
legend({'Data','Fit'},'Location','northwest')
ylabel('m_{Infinity}','interpreter','tex')

t =  '$\frac{1}{1+\exp\left(\frac{V-22.28}{-23.39}\right)}$';


prettyFig();

h = text(-40,.7,t,'interpreter','latex','FontSize',24);

if being_published
	snapnow
	delete(gcf)
end


%% Timescale dependence on voltage
% In this section, I look at how the timescale of activation of the EAG current depends on the voltage. To estimate this, I fit exponentials to the rising phases of the EAG currents, and plot these timescales as a function of voltage. 


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

if being_published	
	snapnow	
	delete(gcf)
end


%% Calcium dependence of gating variable
% In this section, I look at how the steady state currents (a proxy for the gating variable) vary with Calcium concentration in HEK cells. First, I plot the steady state currents vs. Calcium concnetration for each cell, both for the HEK cells that express EAG (top row), and for the HEK cells that don't (bottom row). 

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
% Now I average over all cells, and compare the curves between cells that have EAG and cells that don't. 

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

%% Version Info
%
pFooter;














