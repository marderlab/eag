% makes fig 

s = synTerm(false,false);
[x,x2] =  makeEAGXolotl(s);

[V_wt,Ca_wt] = x.integrate;
[V_mut,Ca_mut] = x2.integrate;
time = x.dt:x.dt:x.t_end;
time = time - 700;

figure('outerposition',[3 3 1200 900],'PaperUnits','points','PaperSize',[1200 600]); hold on
for i = 1:4
	ax(i) = subplot(2,2,i); hold on
end
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