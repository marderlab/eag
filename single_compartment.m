% first make some single-compartment models

% parameters
s = synTerm(false,false);
p = s.parameters;
eag_gbar = 300;

x = xolotl;
x.cleanup;

x.dt = s.dt;
x.t_end = s.t_end;
T = (s.dt:s.dt:s.t_end);
a = find(T>700,1,'first');
z = find(T>703,1,'first');
V_clamp = 0*(x.dt:x.dt:x.t_end) - 50;
V_clamp(a:z) = 70;

x.addCompartment('C0',-70,p.Ca_in,s.Cm,s.A,s.A,p.phi,s.Ca_out,p.Ca_in,s.tau_Ca,0);
x.V_clamp = V_clamp; % only the first compartment can be clamped 


x.addCompartment('Term',-70,p.Ca_in,s.Cm,s.A,s.A,p.phi,s.Ca_out,p.Ca_in,s.tau_Ca,0);

x.addConductance('Term','aldrich/DmNaV',500,s.E_Na);
x.addConductance('Term','Shaker',60,s.E_K);
x.addConductance('Term','EAGwt',eag_gbar,s.E_K);
x.addConductance('Term','EAGmut',eag_gbar,s.E_K);
x.addConductance('Term','EAGes',0,s.E_K);
x.addConductance('Term','hara/Cac',p.gCa2,20);
x.addConductance('Term','Leak',p.gLeak,s.E_leak)

x.addSynapse('Elec','C0','Term',0.06); 

x.transpile;
x.compile;

x.skip_hash_check = true;
x.closed_loop = false;


figure('outerposition',[0 0 1000 901],'PaperUnits','points','PaperSize',[1000 901]); hold on


subplot(2,2,1); hold on

% show the AP waveforms for null, WT, mutant and ES

c = [0 0 0; lines(3)];

% null
x.Term.EAGwt.gbar = 0;
x.Term.EAGmut.gbar = 0;
x.Term.EAGes.gbar = 0;
V = x.integrate;
plot(T-700,V(:,2),'k')
ylabel('V_m (mV)')
xlabel('Time (ms)')

% WT
x.Term.EAGwt.gbar = eag_gbar;
x.Term.EAGmut.gbar = 0;
x.Term.EAGes.gbar = 0;
V = x.integrate;
plot(T-700,V(:,2),'Color',c(2,:))

% mutant
x.Term.EAGmut.gbar = eag_gbar;
x.Term.EAGwt.gbar = 0;
x.Term.EAGes.gbar = 0;
V = x.integrate;
plot(T-700,V(:,2),'Color',c(3,:))

% % extra-sensitive
% x.Term.EAGmut.gbar = 0;
% x.Term.EAGwt.gbar = 0;
% x.Term.EAGes.gbar = eag_gbar;
% V = x.integrate;
% plot(T-700,V(:,2),'Color',c(4,:))

set(gca,'XLim',[-10 100])
legend({'EAG null','EAG wt','EAG mutant','EAG high sensitive'})
title('Ca_{ext} = 1500uM')


% now vary the external calcium 
all_Ca = logspace(log10(15),log10(1500),20);
all_ap_width = NaN(length(all_Ca),3);
all_time_to_peak = NaN(length(all_Ca),3);


for i = 1:length(all_Ca)

	x.Term.Ca_out = all_Ca(i);

	% WT
	x.Term.EAGwt.gbar = eag_gbar;
	x.Term.EAGmut.gbar = 0;
	x.Term.EAGes.gbar = 0;
	
	V = x.integrate; V = V(:,2);

	[m,idx] = max(V);
	if m > 0
		all_time_to_peak(i,1) = idx*x.dt - 700;
	end

	a = find(V>0,1,'first');
	z = find(V>0,1,'last');
	if z == length(V)
		z = Inf;
	end
	if ~isempty(z) && ~isempty(a)
		all_ap_width(i,1) = (z-a)*x.dt;
	end

	% mutant
	x.Term.EAGwt.gbar = 0;
	x.Term.EAGmut.gbar = eag_gbar;
	x.Term.EAGes.gbar = 0;
	V = x.integrate; V = V(:,2);

	[m,idx] = max(V);
	if m > 0
		all_time_to_peak(i,2) = idx*x.dt - 700;
	end

	a = find(V>0,1,'first');
	z = find(V>0,1,'last');
	if z == length(V)
		z = Inf;
	end
	if ~isempty(z) && ~isempty(a)
		all_ap_width(i,2) = (z-a)*x.dt;
	end

	% extra-sensitive
	x.Term.EAGwt.gbar = 0;
	x.Term.EAGmut.gbar = 0;
	x.Term.EAGes.gbar = eag_gbar;
	V = x.integrate; V = V(:,2);

	[m,idx] = max(V);
	if m > 0
		all_time_to_peak(i,3) = idx*x.dt - 700;
	end

	a = find(V>0,1,'first');
	z = find(V>0,1,'last');
	if z == length(V)
		z = Inf;
	end
	if ~isempty(z) && ~isempty(a)
		all_ap_width(i,3) = (z-a)*x.dt;
	end




end

subplot(2,2,2); hold on
for i = 1:2
	plot(all_Ca,all_ap_width(:,i),'o','Color',c(i+1,:))
end
legend({'EAG wt','EAG mutant','EAG hi sensitive'})
ylabel('AP width at 0mV (ms)')
xlabel('Ca_{ext} (uM)')
set(gca,'XScale','log','YLim',[0 2.5])



% now vary the gbar of eag and show effect of AP
x.Term.Ca_out = 1500;

all_eag_gbar = linspace(0,600,20);
all_ap_width = NaN(length(all_eag_gbar),2);
all_time_to_peak = NaN(length(all_eag_gbar),2);

for i = 1:length(all_eag_gbar)


	% WT
	x.Term.EAGwt.gbar = all_eag_gbar(i);
	x.Term.EAGmut.gbar = 0;
	V = x.integrate; V = V(:,2);

	[m,idx] = max(V);
	if m > 0
		all_time_to_peak(i,1) = idx*x.dt - 700;
	end

	a = find(V>0,1,'first');
	z = find(V>0,1,'last');
	if z == length(V)
		z = Inf;
	end
	if ~isempty(z) && ~isempty(a)
		all_ap_width(i,1) = (z-a)*x.dt;
	end

	% mutant
	x.Term.EAGwt.gbar = 0;
	x.Term.EAGmut.gbar = all_eag_gbar(i);
	V = x.integrate; V = V(:,2);

	[m,idx] = max(V);
	if m > 0
		all_time_to_peak(i,2) = idx*x.dt - 700;
	end

	a = find(V>0,1,'first');
	z = find(V>0,1,'last');
	if z == length(V)
		z = Inf;
	end
	if ~isempty(z) && ~isempty(a)
		all_ap_width(i,2) = (z-a)*x.dt;
	end


end

subplot(2,2,3); hold on
plot(all_eag_gbar,all_time_to_peak(:,1),'bo-')
plot(all_eag_gbar,all_time_to_peak(:,2),'ro-')
legend({'EAG wt','EAG mutant'})
ylabel('Time to peak (ms)')
xlabel('g_{EAG} (uS/mm^2)')
set(gca,'YLim',[0 4])

subplot(2,2,4); hold on
plot(all_eag_gbar,all_ap_width(:,1),'bo-')
plot(all_eag_gbar,all_ap_width(:,2),'ro-')
legend({'EAG wt','EAG mutant'})
ylabel('AP width at 0mV (ms)')
xlabel('g_{EAG} (uS/mm^2)')
set(gca,'YLim',[0 25])

prettyFig();

