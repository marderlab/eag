% makes fig 

s = synTerm(false,false);

[x,x2] =  makeEAGXolotl(s);


[V_wt,Ca_wt] = x.integrate;
[V_mut,Ca_mut] = x2.integrate;
time = x.dt:x.dt:x.t_end;
time = time - 700;

x.skip_hash_check = true;
x2.skip_hash_check = true;


% wire up some parameters together
for i = 1:4
	x.synapses(i).gbar = @() evalin('base','params(1)');
	x2.synapses(i).gbar = @() evalin('base','params(1)');
end

x.C1.EAGwt.gbar = @() evalin('base','params(2)');
x.C2.EAGwt.gbar = @() evalin('base','params(2)');
x2.C1.EAGmut.gbar = @() evalin('base','params(2)');
x2.C2.EAGmut.gbar = @() evalin('base','params(2)');

x.C1.Leak.gbar = @() evalin('base','params(3)');
x.C2.Leak.gbar = @() evalin('base','params(3)');
x2.C1.Leak.gbar = @() evalin('base','params(3)');
x2.C2.Leak.gbar = @() evalin('base','params(3)');

x.C2.Cac.gbar = @() evalin('base','params(4)');
x2.C2.Leak.gbar = @() evalin('base','params(4)');

x.C1.Shaker.gbar = @() evalin('base','params(5)');
x.C2.Shaker.gbar = @() evalin('base','params(5)');
x2.C1.Shaker.gbar = @() evalin('base','params(5)');
x2.C2.Shaker.gbar = @() evalin('base','params(5)');

x.C1.DmNaV.gbar = @() evalin('base','params(6)');
x.C2.DmNaV.gbar = @() evalin('base','params(6)');
x2.C1.DmNaV.gbar = @() evalin('base','params(6)');
x2.C2.DmNaV.gbar = @() evalin('base','params(6)');

x.C1.phi = @() evalin('base','params(7)');
x.C2.phi = @() evalin('base','params(7)');
x2.C1.phi = @() evalin('base','params(7)');
x2.C2.phi = @() evalin('base','params(7)');


param_names = {'g_{internal} (nS)','g_{EAG} (\muS/mm^2)','g_{Leak} (\muS/mm^2)','g_{Ca} (\muS/mm^2)','g_{K} (\muS/mm^2)','g_{Na} (\muS/mm^2)','\phi'};
lb = [.01, 0,   0,  0,   0,   0,   0];
ub = [1 , 1000, 10, 40, 100, 2e3, 1000];


for i = 1:7
	params = [.1, 600, 1, 20, 50, 1e3, 500];

	figure('outerposition',[3 3 1200 900],'PaperUnits','points','PaperSize',[1200 600]); hold on
	for j = 1:4
		ax(j) = subplot(2,2,j); hold on
	end

	this_param_values = linspace(lb(i),ub(i),20);

	time = x.dt:x.dt:x.t_end;
	time = time - 700;

	area_under_spike_1_wt = NaN*this_param_values;
	area_under_spike_2_wt = NaN*this_param_values;

	area_under_spike_1_mut = NaN*this_param_values;
	area_under_spike_2_mut = NaN*this_param_values;

	ca_width_wt = NaN*this_param_values;
	ca_width_mut = NaN*this_param_values;

	ca_peak_wt = NaN*this_param_values;
	ca_peak_mut = NaN*this_param_values;

	for j = 1:length(this_param_values)
		params(i) = this_param_values(j);

		[V_wt,Ca_wt] = x.integrate;
		[V_mut,Ca_mut] = x2.integrate;

		% do compartment 1
		V = V_wt(:,2);
		V = V - V(time==0);
		area_under_spike_1_wt(j) = sum(abs(V(find(time==0):find(time==100))))*x.dt;
		
		% do compartment 2
		V = V_wt(:,3);
		V = V - V(time==0);
		area_under_spike_2_wt(j) = sum(abs(V(find(time==0):find(time==100))))*x.dt;

		% find Ca
		Ca = Ca_wt(:,3);
		ca_max = max(Ca);
		a = find(Ca>ca_max/2,1,'first');
		z = find(Ca>ca_max/2,1,'last');
		ca_width_wt(j) = (z-a)*x.dt;
		ca_peak_wt(j) = ca_max;

		% do compartment 1
		V = V_mut(:,2);
		V = V - V(time==0);
		area_under_spike_1_mut(j) = sum(abs(V(find(time==0):find(time==100))))*x.dt;
		
		% do compartment 2
		V = V_mut(:,3);
		V = V - V(time==0);
		area_under_spike_2_mut(j) = sum(abs(V(find(time==0):find(time==100))))*x.dt;

		% find Ca
		Ca = Ca_mut(:,3);
		ca_max = max(Ca);
		a = find(Ca>ca_max/2,1,'first');
		z = find(Ca>ca_max/2,1,'last');
		ca_width_mut(j) = (z-a)*x.dt;
		ca_peak_mut(j) = ca_max;

	end

	plot(ax(1),this_param_values,area_under_spike_1_wt,'k+-')
	plot(ax(1),this_param_values,area_under_spike_1_mut,'r+-')

	plot(ax(2),this_param_values,area_under_spike_2_wt,'ko-')
	plot(ax(2),this_param_values,area_under_spike_2_mut,'ro-')

	plot(ax(3),this_param_values,ca_width_wt,'ko-')
	plot(ax(3),this_param_values,ca_width_mut,'ro-')

	plot(ax(4),this_param_values,ca_peak_wt,'ko-')
	plot(ax(4),this_param_values,ca_peak_mut,'ro-')

	for j = 1:length(ax)
		xlabel(ax(j),param_names{i})
	end
	set(ax(3),'YLim',[0 200])
	ylabel(ax(3),'Width of Calcium response (ms)')

	ylabel(ax(1),'Area under spike (mV ms)')
	ylabel(ax(2),'Area under spike (mV ms)')

	ylabel(ax(4),'Peak calcium (\muM)')

	prettyFig();
	


end



return

axes(ax(1));
o = imread('./pics/cartoon.png');
imagesc(o);
axis ij
axis image
axis off
ax(1).Position = [.05 .52 .5 .4];

axes(ax(2))
plot(time,V_wt(:,3),'k')
plot(time,V_mut(:,3),'r')
xlabel('Time (ms)')
set(gca,'XLim',[-50 250],'YLim',[-90 80])
ylabel('V (mV)')
legend({'WT','Mutant'})

axes(ax(3))
plot(time,V_wt(:,3),'k')
plot(time,V_wt(:,2),'k:')
xlabel('Time (ms)')
set(gca,'XLim',[-50 250],'YLim',[-90 80])
ylabel('V (mV)')
legend({'In calcium microdomain','Outside calcium microdomian'})

axes(ax(4))
plot(time,Ca_wt(:,3),'k')
plot(time,Ca_mut(:,3),'r')
xlabel('Time (ms)')
set(gca,'XLim',[-50 250])
ylabel('[Ca^2^+] (\muM)')


prettyFig();

labelFigure('x_offset',-.01,'y_offset',.01,'font_size',28)
ax(4).Box = 'off';