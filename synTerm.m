% interactive model for studying effects of EAG
% close to and far away from calcium
% microdomains in a simple model of the pre-
% synaptic terminal in Drosophila NMJs
% written by Srinivas Gorur-Shandilya (http://srinivas.gs)
% This code is free software, and its use is governed by the 
% GNU GPL v3. 

classdef synTerm < handle 
	properties 
		handles
		parameters
            units
		lb
		ub
		dt = 50e-3;
		t_end = 1e3;
		xolotl_obj@xolotl

        % define some fixed parameters 
        Cm = 10; % uS/mm^2
        A = 1.25e-5; % mm^2
        Ca_out = 1500; % uM
        tau_Ca = 50; % ms
        E_Na = 72.8;
        E_K = -106;
        E_leak = -70;

	end % end props

    methods 

        function self = synTerm(make_fig)

            if nargin < 1
                make_fig = true;
            end

            p.Ca_in = 23e-3;
            lb.Ca_in = 23e-3;
            ub.Ca_in = 23e-2;
            units.Ca_in = 'uM';

            p.phi = 500;
            lb.phi = .1;
            ub.phi = 1e3;
            units.phi = '';

            p.gNa1 = 1e3;
            lb.gNa1 = 0;
            ub.gNa1 = 1e3;
            units.gNa1 = 'uS/mm^2';

            p.gNa2 = 1e3;
            lb.gNa2 = 0;
            ub.gNa2 = 1e3;
            units.gNa2 = 'uS/mm^2';

            p.gK1 = 50;
            lb.gK1 = 0;
            ub.gK1 = 1e3;
            units.gK1 = 'uS/mm^2';

            p.gK2 = 50;
            lb.gK2 = 0;
            ub.gK2 = 1e3;
            units.gK2 = 'uS/mm^2';

            p.gCa2 = 20;
            lb.gCa2 = 0;
            ub.gCa2 = 1e3;
            units.gCa2 = 'uS/mm^2';

            p.gLeak = 1;
            lb.gLeak = 0;
            ub.gLeak = 1e3;
            units.gLeak = 'uS/mm^2';

            p.gEAG1 = 600;
            lb.gEAG1 = 0;
            ub.gEAG1 = 1e3;
            units.gEAG1 = 'uS/mm^2';

            p.gEAG2 = 600;
            lb.gEAG2 = 0;
            ub.gEAG2 = 1e3;
            units.gEAG2 = 'uS/mm^2';

            p.g1 = .1;
            lb.g1 = 0;
            ub.g1 = 10;
            units.g1 = 'uS';

            p.g2 = .1;
            lb.g2 = 0;
            ub.g2 = 10;
            units.g2 = 'uS';

            p.g12 = .1;
            lb.g12 = 0;
            ub.g12 = 10;
            units.g12 = 'uS';

            % check if there is something in the cache
            temp = cache('syn_term_params');
            if ~isempty(temp)
              disp('Reading cached values...')
              % replace parameters with cached parameters
              p = temp;
              % pick some reasonable lower and upper bounds
              f = fieldnames(temp);
              for i = 1:length(f)
                  if  lb.(f{i}) > temp.(f{i})
                      lb.(f{i}) = temp.(f{i});
                  end
                  if  ub.(f{i}) < temp.(f{i})
                      ub.(f{i}) = temp.(f{i});
                  end
              end
            end

            self.parameters = p;
            self.lb = lb;
            self.ub = ub;
            self.units = units;


             % create a window to show the voltage traces 
            if make_fig
                self.handles.time_series_fig = figure('position',[100 250 1701 1100],'NumberTitle','off','IntegerHandle','off','Name',['manipulate[synapticTerminal]']);
                for i = 1:6
                  self.handles.ax(i) = subplot(3,2,i); hold on
                  xlabel(self.handles.ax(i),'Time (ms)')
                end

                ylabel(self.handles.ax(1),'V (mV)')
                ylabel(self.handles.ax(2),'Ca (\mu M)')
                ylabel(self.handles.ax(3),'I (nA)')
                ylabel(self.handles.ax(4),'I (nA)')
                title(self.handles.ax(4),'Currents in Calcium microdomain')
                title(self.handles.ax(3),'Currents in control compartment')

                ylabel(self.handles.ax(5),'Gating variable')
                ylabel(self.handles.ax(6),'Gating variable')

                linkaxes(self.handles.ax,'x')
                linkaxes(self.handles.ax(3:4),'y')

                set(self.handles.ax(1),'XLim',[0 self.t_end])
                set(self.handles.ax(1),'YLim',[-100 100])

                self.configureXolotl;

                % spawn a puppeteer object
                pp = puppeteer(self.parameters,self.lb,self.ub,units);

                % configure 
                attachFigure(pp,self.handles.time_series_fig)
                pp.callback_function = @self.updatePlots;

                self.makePlotHandles;
            end

        end %end constructor 

        function makePlotHandles(self)

            % create vectors to store simulations
            T = (self.dt:self.dt:self.t_end);
            V = NaN*(T);
            C = NaN*(T);

            self.handles.V_trace(1,1) = plot(self.handles.ax(1),T,V,'k','LineWidth',2);
            self.handles.V_trace(2,1) = plot(self.handles.ax(1),T,V,'r','LineWidth',2);
            self.handles.V_trace(3,1) = plot(self.handles.ax(1),T,V,'b','LineWidth',2);
            self.handles.V_trace(1,2) = plot(self.handles.ax(1),T,V,'k','LineWidth',2);
            self.handles.V_trace(2,2) = plot(self.handles.ax(1),T,V,'r','LineWidth',2);
            self.handles.V_trace(3,2) = plot(self.handles.ax(1),T,V,'b','LineWidth',2);
            for i = 1:3
                self.handles.V_trace(i,2).LineStyle = ':';
            end
            legend(self.handles.V_trace(:,1),{'Incoming AP','V far from VGCC','V close to VGCC'})

            self.handles.C_trace(1) = plot(self.handles.ax(2),T,V,'b','LineWidth',2);
            self.handles.C_trace(2) = plot(self.handles.ax(2),T,V,'b','LineWidth',2,'LineStyle',':');

            % also plot the internal and external calcium concentrations in dotted lines 
            self.handles.Ca_out = plot(self.handles.ax(2),T,T*0 + self.Ca_out,'k:');
            self.handles.Ca_in = plot(self.handles.ax(2),T,T*0 + self.parameters.Ca_in,'k:');
            set(self.handles.ax(2),'YLim',[self.parameters.Ca_in/2 self.Ca_out*2],'YScale','log')


            % currents
            self.handles.Itrace(1) = plot(self.handles.ax(3),T,V,'LineWidth',2);
            self.handles.Itrace(2) = plot(self.handles.ax(3),T,V,'LineWidth',2);
            self.handles.Itrace(3) = plot(self.handles.ax(3),T,V,'LineWidth',2);
            self.handles.Itrace(4) = plot(self.handles.ax(4),T,V,'LineWidth',2);
            self.handles.Itrace(5) = plot(self.handles.ax(4),T,V,'LineWidth',2);
            self.handles.Itrace(6) = plot(self.handles.ax(4),T,V,'LineWidth',2);
            self.handles.Itrace(7) = plot(self.handles.ax(4),T,V,'LineWidth',2);

            legend(self.handles.Itrace(1:3),{'gNa','gK','gEAG'})
            legend(self.handles.Itrace(4:7),{'gNa','gK','gEAG','gCa'})

            % gating variables
            self.handles.m1(1) = plot(self.handles.ax(5),T,V,'LineWidth',2);
            self.handles.m1(2) = plot(self.handles.ax(5),T,V,'LineWidth',2);
            legend(self.handles.m1,{'m_{K}','m_{EAG}'})

            self.handles.m2(1) = plot(self.handles.ax(6),T,V,'LineWidth',2);
            self.handles.m2(2) = plot(self.handles.ax(6),T,V,'LineWidth',2);
            legend(self.handles.m2,{'m_{K}','m_{EAG}'})

            figure(self.handles.time_series_fig)
            prettyFig();

            self.handles.V_trace(1).YData = self.xolotl_obj(1).V_clamp;
        end

		function updatePlots(self,parameters)
			self.parameters = parameters;
            self.handles.Ca_in.YData = parameters.Ca_in + 0*self.handles.Ca_in.YData;
            self.handles.ax(2).YLim(1) = parameters.Ca_in/2;

			% update xolotl objects
            for i = 1:2
      			self.xolotl_obj(i).C1.Ca_in = parameters.Ca_in;
      			self.xolotl_obj(i).C2.Ca_in = parameters.Ca_in;
      			self.xolotl_obj(i).C2.phi = parameters.phi;

      			% update conductances
      			self.xolotl_obj(i).C1.DmNaV.gbar = parameters.gNa1;
      			self.xolotl_obj(i).C2.DmNaV.gbar = parameters.gNa2;
      			self.xolotl_obj(i).C1.Shaker.gbar = parameters.gK1;
      			self.xolotl_obj(i).C2.Shaker.gbar = parameters.gK2;
      			self.xolotl_obj(i).C2.Cac.gbar = parameters.gCa2;
      			self.xolotl_obj(i).C2.Leak.gbar = parameters.gLeak;
      			self.xolotl_obj(i).C1.Leak.gbar = parameters.gLeak;
      			self.xolotl_obj(i).synapses(1).gbar = parameters.g1;
      			self.xolotl_obj(i).synapses(2).gbar = parameters.g2;
      			self.xolotl_obj(i).synapses(3).gbar = parameters.g12;
      			self.xolotl_obj(i).synapses(4).gbar = parameters.g12;
            end

            self.xolotl_obj(1).C1.EAGwt.gbar = parameters.gEAG1;
            self.xolotl_obj(1).C2.EAGwt.gbar = parameters.gEAG2;
            self.xolotl_obj(2).C1.EAGmut.gbar = parameters.gEAG1;
            self.xolotl_obj(2).C2.EAGmut.gbar = parameters.gEAG2;

			% integrate the WT model
			[V,Ca,I_clamp,C] = self.xolotl_obj(1).integrate;

			% extract data and update plots
			self.handles.V_trace(2,1).YData = V(:,2);
			self.handles.V_trace(3,1).YData = V(:,3);
			self.handles.C_trace(1).YData = Ca(:,3);

            % update gating variables
            self.handles.m1(1).YData = C(:,3);
            self.handles.m1(2).YData = C(:,5);
            self.handles.m2(1).YData = C(:,11);
            self.handles.m2(2).YData = C(:,13);

            % compute currents (actually this isn't correct---didn't update for the new drosophila channels)
            self.handles.Itrace(1).YData = (C(:,1).^3).*C(:,2).*parameters.gNa1.*(V(:,2) - self.E_Na); % I_Na
            self.handles.Itrace(2).YData = (C(:,3).^4).*parameters.gK1.*(V(:,2) - self.E_K); % I_K1
            self.handles.Itrace(3).YData = (C(:,5).^2).*parameters.gEAG1.*(V(:,2) - self.E_K); % I_EAG1

            self.handles.Itrace(4).YData = (C(:,9).^3).*C(:,10).*parameters.gNa2.*(V(:,3) - self.E_Na); % I_Na
            self.handles.Itrace(5).YData = (C(:,11).^4).*parameters.gK2.*(V(:,3) - self.E_K); % I_K1
            self.handles.Itrace(6).YData = (C(:,13).^2).*parameters.gEAG2.*(V(:,3) - self.E_K); % I_EAG1
            E_Ca = Ca(:,6);
            self.handles.Itrace(7).YData = (C(:,15).^3).*C(:,16).*parameters.gCa2.*(V(:,3) - E_Ca); % I_Ca

            my = nanmax(abs([self.handles.Itrace.YData]))*1.1;
            set(self.handles.ax(3),'YLim',[-my my])
            set(self.handles.ax(4),'YLim',[-my my])


            % integrate the mutant 
            [V,Ca,I_clamp,C] = self.xolotl_obj(2).integrate;
            self.handles.V_trace(2,2).YData = V(:,2);
            self.handles.V_trace(3,2).YData = V(:,3);

            self.handles.C_trace(2).YData = Ca(:,3);


		end

        function delete(self)
              % cache the parameters for later use 
              cache('syn_term_params',self.parameters)
        end

        function configureXolotl(self)
              [x,x2] =  makeEAGXolotl(self);
              self.xolotl_obj = [x; x2];
        end        

	end % end methods 

end % end classdef


