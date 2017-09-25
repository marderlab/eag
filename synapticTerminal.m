% synapticTerminal.m

classdef synapticTerminal < handle & matlab.mixin.CustomDisplay

    properties
        dt = 100e-3 % ms
        t_end =  5000 % ms

        initial_conditions = [-50 -50 0.05 0 0 0 0 0 0 0 0];

        % these propoerties are for manipulate
        handles % for GUI elements
        % order is as follows:

        parameters
        lb
        ub
        units

        E_Na = 72.8;
        E_K = -106;
        E_leak = -50;

        V_drive

        implementation

    end

    properties (Access = protected)
        legal_parameter_names
        available_implementations = {'SGS/MATLAB'};
    end


    methods (Access = protected)
        function displayScalarObject(self)
            url = 'https://gitlab.com/marderlab/neuron-network-models/';
            fprintf(['<a href="' url '">neuron</a> object with the following properties:\n\n']);


            fprintf('Parameter         Value\n')
            fprintf('------------     --------\n')
            f = fieldnames(self.parameters);
            for i = 1:length(f)
                padding_string_size = 17 - length(f{i});
                padding_string = repmat(' ',1,padding_string_size);
                v = [oval(self.parameters.(f{i})) ' ' self.units.(f{i})];
                this_string = [f{i} padding_string v '\n'];
                fprintf(this_string)
            end

            fprintf('\n')
            disp('Underling implementation: ')
            fprintf('-----------------------\n')
            disp(self.implementation)

            fprintf('\n')
            disp('Version (git hash): ')
            fprintf('-----------------------\n')
            disp(gitHash((which(mfilename))))

        end % end displayScalarObject
   end % end protected methods

    methods
        function self = synapticTerminal()
    

            % define parameters 
            p.C_m = 10;
            lb.C_m = 1;
            ub.C_m = 20;
            units.C_m = 'nF/mm^2';

            p.A1 = 1.25e-5;
            lb.A1 = 1e-6;
            ub.A1 = 1e-3;
            units.A1 = 'mm^2';

            p.A2 = 1.25e-5;
            lb.A2 = 1e-6;
            ub.A2 = 1e-3;
            units.A2 = 'mm^2';

            p.Ca_out = 1500;
            lb.Ca_out = 0.03;
            ub.Ca_out = 5e3;
            units.Ca_out = 'uM';

            p.Ca_in = 250;
            lb.Ca_in = 100;
            ub.Ca_in = 500;
            units.Ca_in = 'uM';

            p.tau_Ca = 50;
            lb.tau_Ca = 1;
            ub.tau_Ca = 200;
            units.tau_Ca = 'ms';

            p.f = 1.496;
            lb.f = .1;
            ub.f = 10;
            units.f = 'uM/nA';

            p.gNa1 = 1830;
            lb.gNa1 = 0;
            ub.gNa1 = 1e4;
            units.gNa1 = 'uS/mm^2';

            p.gNa2 = 1830;
            lb.gNa2 = 0;
            ub.gNa2 = 1e4;
            units.gNa2 = 'uS/mm^2';

            p.gK1 = 610;
            lb.gK1 = 0;
            ub.gK1 = 1e4;
            units.gK1 = 'uS/mm^2';

            p.gK2 = 610;
            lb.gK2 = 0;
            ub.gK2 = 1e4;
            units.gK2 = 'uS/mm^2';

            p.gCa2 = 610;
            lb.gCa2 = 0;
            ub.gCa2 = 1e4;
            units.gCa2 = 'uS/mm^2';

            p.gLeak = .99;
            lb.gLeak = 0;
            ub.gLeak = 1e4;
            units.gLeak = 'uS/mm^2';

            p.gEAG1 = .99;
            lb.gEAG1 = 0;
            ub.gEAG1 = 1e4;
            units.gEAG1 = 'uS/mm^2';

            p.gEAG2 = .99;
            lb.gEAG2 = 0;
            ub.gEAG2 = 1e4;
            units.gEAG2 = 'uS/mm^2';

            p.g1 = 23;
            lb.g1 = 0;
            ub.g1 = 1e4;
            units.g1 = 'uS';

            p.g2 = 23;
            lb.g2 = 0;
            ub.g2 = 1e4;
            units.g2 = 'uS';

            p.g12 = 23;
            lb.g12 = 0;
            ub.g12 = 1e4;
            units.g12 = 'uS';



            self.implementation = self.available_implementations{1};

            % check if there is something in the cache
            temp = cache('syn_term_params');
            if ~isempty(temp)
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

        end % end constructor



        function set(self,param,value)
            % check if parameter exists 
            if ~any(strcmp(self.legal_parameter_names,param)) && ~any(strcmp(self.parameter_names,param))
                error('Unknown parameter.')
            end
            idx1 = find(strcmp(self.parameter_names,param));
            idx2 = find(strcmp(self.legal_parameter_names,param));
            idx = unique([idx1 idx2]);
            assert(length(idx) == 1, 'Could not determine parameter to set')
            self.parameters(idx) = value;
        end

        function set.V_drive(self,value)
            time = self.dt:self.dt:self.t_end;
            assert(length(value) == length(time),'V_drive does not match the timing info set in dt and t_end')
            self.V_drive = value;
        end

        function [V1, V2, Ca, currents,gv] = integrate(self)
            switch self.implementation
            case 'SGS/MATLAB'
                [V1, V2, Ca, currents,gv] = self.integrate_SGS_MATLAB;
            case 'SGS/Julia'
                [V2, Ca, n] = self.integrate_SGS_JULIA;
            otherwise
                error('Unknown implementation')
            end
                   
        end

        function delete(self)
            % cache the parameters for later use 
            cache('syn_term_params',self.parameters)
        end

        function [V1, V2, Ca2, currents, gating_variables] = integrate_SGS_MATLAB(self)
            % pure MATLAB implementation of this model 

            % assemble parameters
            P = self.parameters;


            % define functions
            boltz = @(V,A,B) (1./(1 + exp((V+A)./B)));
            tauX = @(V,A,B,D,E) (A - B./(1+exp((V+D)./E)));

            % m_inf
            mNainf = @(V) (boltz(V,25.5,-5.29));
            mKinf = @(V) boltz(V,12.3,-11.8);
            mCaTinf = @(V) boltz(V,27.1,-7.2);
            mEAGinf = @(V,Ca) boltz(V,23.12,-16.94)*(9.29e-2/(9.29e-3+Ca));

            % h_inf 
            hNainf = @(V) (boltz(V,48.9,5.18));
            hCaTinf = @(V) boltz(V,32.1,5.5);

            % tau_m
            taumNa = @(V) tauX(V,1.32,1.26,120.0,-25.0);
            taumK = @(V) tauX(V,7.2,6.4,28.3,-19.2);
            taumCaT = @(V) tauX(V,21.7,21.3,68.1,-20.5);
            taumEAG = @(V) tauX(V,5497,5500,251.5,-51.5);
            % let's just use a simpler tau for now
            % taumEAG = taumK;

            % tau_h
            tauhNa = @(V) (0.67/(1+exp((V+62.9)/-10.0)))*(1.5 + 1/(1+exp((V+34.9)/3.6)));
            tauhCaT = @(V) tauX(V,105.,89.8,55.,-16.9);

            nsteps = floor(self.t_end/self.dt);


            % pre-allocate arrays
            V1 = zeros(nsteps,1);
            mNa1 = mNainf(self.E_leak);
            hNa1 = hNainf(self.E_leak);
            mK1 = mKinf(self.E_leak);
            mEAG1 = mEAGinf(self.E_leak,0.05);


            V2 = zeros(nsteps,1);
            Ca2 = zeros(nsteps,1);
            gNa1 = zeros(nsteps,1);
            gK1 = zeros(nsteps,1);
            gEAG1 = zeros(nsteps,1);
            gNa2 = zeros(nsteps,1);
            gK2 = zeros(nsteps,1);
            gEAG2 = zeros(nsteps,1);
            gCa2 = zeros(nsteps,1);
            E_Ca = zeros(nsteps,1);

            I1 = zeros(nsteps,1);
            I2 = zeros(nsteps,1);


            mNa2 = mNainf(self.E_leak);
            hNa2 = hNainf(self.E_leak);
            mK2 = mKinf(self.E_leak);
            mEAG2 = mEAGinf(self.E_leak,0.05);
            mCaT2 = mCaTinf(self.E_leak);
            hCaT2 = hCaTinf(self.E_leak);

            % set some initial conditions 
            V1(1) = self.E_leak;
            V2(1) = self.E_leak;
            Ca2(1) = self.parameters.Ca_in;

            dt = self.dt;
            R_by_zF = 500.0*(8.6174e-005);
            T = 10 + 273.15;
            RT_by_zF = R_by_zF*T;
            f = P.f;
            fA = f*P.A2;
            exp_dt_by_tau_Ca = exp(-dt/P.tau_Ca);

            for i = 2:nsteps

                V_ext = self.V_drive(i-1);

                % compartment 1  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

                minf_Na = mNainf(V1(i-1));
                hinf_Na = hNainf(V1(i-1));
                minf_K  = mKinf(V1(i-1));
                minf_EAG = mEAGinf(V1(i-1),self.parameters.Ca_in);

                % integrate Na channel
                mNa1 = minf_Na + (mNa1 - minf_Na)*exp(-dt/taumNa(V1(i-1)));
                hNa1 = hinf_Na + (hNa1 - hinf_Na)*exp(-dt/tauhNa(V1(i-1)));

                % integrate K and EAG channels
                mK1 = minf_K + (mK1 - minf_K)*exp(-dt/taumK(V1(i-1)));
                mEAG1 = minf_EAG + (mEAG1 - minf_EAG)*exp(-dt/taumEAG(V1(i-1)));
                
                % compute effective conductances 
                gNa1(i) = P.gNa1*(mNa1^3)*hNa1; 
                gK1(i) = P.gK1*(mK1^4);
                gEAG1(i) = P.gEAG1*(mEAG1^2); 

                sigma_g1 = gNa1(i) + gK1(i) + gEAG1(i) + P.gLeak; 

                % compute currents from electric coupling 
                if ~isnan(V2(i-1))
                    I1(i) = (P.g1*(V_ext - V1(i-1)) + P.g12*(V2(i-1) - V1(i-1)))*1e-6; % in nA
                else
                    I1(i) = (P.g1*(V_ext - V1(i-1)))*1e-6;
                end

                V_inf1 = (gNa1(i)*self.E_Na + gK1(i)*self.E_K + gEAG1(i)*self.E_K + P.gLeak*self.E_leak + I1(i)/P.A1)/sigma_g1;

                tau_v1 = P.C_m/(sigma_g1); % unit correction

                V1(i) = V_inf1 + (V1(i-1) - V_inf1)*exp(-dt/tau_v1);


                % compartment 2  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

                % compute Calcium reversal potential
                E_Ca(i) = RT_by_zF*log(P.Ca_out/Ca2(i-1)); 

                minf_Na = mNainf(V2(i-1));
                hinf_Na = hNainf(V2(i-1));
                minf_K  = mKinf(V2(i-1));
                minf_EAG = mEAGinf(V2(i-1),Ca2(i-1));
                minf_CaT = mCaTinf(V2(i-1));
                hinf_CaT = hCaTinf(V2(i-1));

                % integrate Na channel
                mNa2 = minf_Na + (mNa2 - minf_Na)*exp(-dt/taumNa(V2(i-1)));
                hNa2 = hinf_Na + (hNa2 - hinf_Na)*exp(-dt/tauhNa(V2(i-1)));

                % integrate K and EAG channels
                mK2 = minf_K + (mK2 - minf_K)*exp(-dt/taumK(V2(i-1)));
                mEAG2 = minf_EAG + (mEAG2 - minf_EAG)*exp(-dt/taumEAG(V2(i-1)));

                % integrate Calcium channel
                mCaT2 = minf_Na + (mCaT2 - minf_Na)*exp(-dt/taumCaT(V2(i-1)));
                hCaT2 = hinf_Na + (hCaT2 - hinf_Na)*exp(-dt/tauhCaT(V2(i-1)));
                
                % compute effective conductances 
                gNa2(i) = P.gNa2*(mNa2^3)*hNa2; 
                gK2(i) = P.gK2*(mK2^4);
                gEAG2(i) = P.gEAG2*(mEAG2^2); % assuming p =2
                gCa2(i) = P.gCa2*(mCaT2^3)*hCaT2; 

                sigma_g2 = gNa2(i) + gK2(i) + gEAG2(i) + gCa2(i) + P.gLeak; 

                % compute currents from electric coupling 
                I2(i) = (P.g2*(V_ext - V2(i-1)) + P.g12*(V1(i-1) - V2(i-1)))*1e-6; % in nA

                V_inf2 = (gNa2(i)*self.E_Na + gK2(i)*self.E_K + gCa2(i)*E_Ca(i) + gEAG2(i)*self.E_K + P.gLeak*self.E_leak +  I2(i)/P.A2)/sigma_g2;
                tau_v2 = P.C_m/(sigma_g2); % unit correction

                V2(i) = V_inf2 + (V2(i-1) - V_inf2)*exp(-dt/tau_v2);


                % update Calcium 
                ca_current = gCa2(i)*(V2(i-1) - E_Ca(i));
                cinf = P.Ca_in - fA*ca_current; 
                Ca2(i) = cinf + (Ca2(i-1) - cinf)*exp_dt_by_tau_Ca; 

            end


            currents.gNa1 = gNa1.*(V1 - self.E_Na).*P.A2;
            currents.gK1 = gK1.*(V1 - self.E_K).*P.A2;
            currents.gEAG1 = gEAG1.*(V1 - self.E_K).*P.A2;
            currents.gNa2 = gNa2.*(V2 - self.E_Na).*P.A2;
            currents.gK2 = gK2.*(V2 - self.E_K).*P.A2;
            currents.gEAG2 = gEAG2.*(V2 - self.E_K).*P.A2;
            currents.gCa2 = gCa2.*(V2 - E_Ca).*P.A2;


            % also convert the conductances back into gating variables so we can examine what the channels are actually doing 
            gating_variables.EAG1 = (gEAG1./self.parameters.gEAG1).^(1/2);
            gating_variables.EAG2 = (gEAG2./self.parameters.gEAG2).^(1/2);
            gating_variables.K1 = (gK1./self.parameters.gK1).^(1/4);
            gating_variables.K2 = (gK2./self.parameters.gK2).^(1/4);

            currents.I1 = I1;
            currents.I2 = I2;

        end

        function manipulate(self)
            % create a window to show the voltage traces 
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
            title(self.handles.ax(3),'Currents in Compartment #1')

            ylabel(self.handles.ax(5),'Gating variable')
            ylabel(self.handles.ax(6),'Gating variable')

            linkaxes(self.handles.ax,'x')
            linkaxes(self.handles.ax(3:4),'y')

            set(self.handles.ax(1),'XLim',[0 self.t_end])
            set(self.handles.ax(1),'YLim',[-80 80])
            

            % spawn a puppeteer object
            pp = puppeteer(self.parameters,self.lb,self.ub,self.units);

            % configure 
            attachFigure(pp,self.handles.time_series_fig)
            pp.callback_function = @self.updatePlots;

            dt =self.dt;
            t_end = self.t_end;

            % create vectors to store simulations
            T = (dt:dt:t_end);
            V = NaN*(dt:dt:t_end);
            C = NaN*(dt:dt:t_end);

            self.handles.V_trace(1) = plot(self.handles.ax(1),T,V,'k','LineWidth',2);
            self.handles.V_trace(2) = plot(self.handles.ax(1),T,V,'r','LineWidth',2);
            self.handles.V_trace(3) = plot(self.handles.ax(1),T,V,'b','LineWidth',2);
            legend(self.handles.V_trace,{'Incoming AP','V far from VGCC','V close to VGCC'})

            self.handles.C_trace = plot(self.handles.ax(2),T,V,'b','LineWidth',2);
            % also plot the internal and external calcium concentrations in dotted lines 
            self.handles.Ca_out = plot(self.handles.ax(2),T,T*0 + self.parameters.Ca_out,'k:');
            self.handles.Ca_in = plot(self.handles.ax(2),T,T*0 + self.parameters.Ca_in,'k:');
            set(self.handles.ax(2),'YLim',[self.parameters.Ca_in/2 self.parameters.Ca_out*2],'YScale','log')


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
        

            self.handles.V_trace(1).YData = self.V_drive;

 
        end % end manipulate 

        function updatePlots(self,parameters)

            self.parameters = parameters;

            self.handles.Ca_in.YData = self.handles.Ca_in.YData*0 + self.parameters.Ca_in;
            self.handles.Ca_out.YData = self.handles.Ca_out.YData*0 + self.parameters.Ca_out;

            [V1,V2,Ca2,currents,gv] = self.integrate;

            try
                self.handles.m1(1).YData = gv.K1;
                self.handles.m1(2).YData = gv.EAG1;

                self.handles.m2(1).YData = gv.K2;
                self.handles.m2(2).YData = gv.EAG2;

                self.handles.V_trace(2).YData = V1;
                self.handles.V_trace(3).YData = V2;
                self.handles.C_trace.YData = Ca2;
                
                self.handles.Itrace(1).YData = currents.gNa1;
                self.handles.Itrace(2).YData = currents.gK1;
                self.handles.Itrace(3).YData = currents.gEAG1;

                self.handles.Itrace(4).YData = currents.gNa2;
                self.handles.Itrace(5).YData = currents.gK2;
                self.handles.Itrace(6).YData = currents.gEAG2;
                self.handles.Itrace(7).YData = currents.gCa2;

                set(self.handles.ax(2),'YLim',[self.parameters.Ca_in/2 self.parameters.Ca_out*2],'YScale','log')

                Y = max(max(abs(struct2mat(currents))));
                self.handles.ax(3).YLim = [-Y Y];
            catch
                warning('Bad integraton')
            end
        end


    end % end methods 

end % end classdef 



