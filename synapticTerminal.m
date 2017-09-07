% synapticTerminal.m

classdef synapticTerminal < handle & matlab.mixin.CustomDisplay

    properties
        dt = 100e-3 % ms
        t_end =  5000 % ms

        initial_conditions = [-50 -50 0.05 0 0 0 0 0 0 0 0];

        % these propoerties are for manipulate
        handles % for GUI elements
        % order is as follows:
        %     cm  A extCa Ca0  tauCa f   conductances(8)
        parameter_names = {'C_{m}','A1','A2','Ca_{o}','Ca_{i}','tau_{Ca}','f'  ,'g_{Na}1','g_{Na}2','g_{Ca}','g_{K}1','g_{K}2','g_1','g_2','g_1_2','g_{leak}','g_{EAG}1','g_{EAG}2'}
        parameters      = [10      .06   .06   3000     0.05      200      15    183       183        2.3      30       30     700   700     300    .0993,    10, 10]
        lb              = [0       0     0     1        0.01      1        1     0         0          0        0        0       0     0       0      0, 0, 0];
        ub              = [100     1     1     1e4      100       1000     100   300       300        300      300      300     999   999     300    300, 300, 300];
        parameter_units = {'nF/mm^2','mm^2','mm^2','uM','uM',    'ms','uM/nA','mS/cm^2','mS/cm^2','mS/cm^2','mS/cm^2','mS/cm^2','mS/cm^2','mS/cm^2','mS/cm^2','mS/cm^2','mS/cm^2','mS/cm^2'}


        E_Na = 50;
        E_K = -80;
        E_leak = -50;

        V_drive

        implementation

    end

    properties (Access = protected)
        legal_parameter_names
        available_implementations = {'SGS/Julia','SGS/MATLAB'};
    end


    methods (Access = protected)
        function displayScalarObject(self)
            url = 'https://gitlab.com/marderlab/neuron-network-models/';
            fprintf(['<a href="' url '">synapticTerminal</a> object with the following properties:\n\n']);


            fprintf('Parameter         Value\n')
            fprintf('------------     --------\n')
            for i = 1:length(self.parameter_names)
                padding_string_size = 17 - length(self.parameter_names{i});
                padding_string = repmat(' ',1,padding_string_size);
                v = [oval(self.parameters(i)) ' ' self.parameter_units{i}];
                this_string = [self.parameter_names{i} padding_string v '\n'];
                fprintf(this_string)
            end

            fprintf('\n')


        end % end displayScalarObject
   end % end protected methods

    methods
        function self = synapticTerminal()
            self.constructLegalParameterNames;

            self.implementation = self.available_implementations{2};

            % % add my Julia code 
            % p = fileparts(which(mfilename));
            % jl.include(joinPath(p,'param_types.jl')) 
            % jl.include(joinPath(p,'liu_gating_functions.jl'))
            % jl.include(joinPath(p,'synapticTerminal.jl'))
            % jl.include(joinPath(p,'neuronMex.jl'))
        end % end constructor

        function constructLegalParameterNames(self)

            % construct legal parameter names from parameter names
            self.legal_parameter_names = self.parameter_names;
            for i = 1:length(self.parameter_names)
                this_name = self.parameter_names{i};
                this_name = strrep(this_name,'{','');
                this_name = strrep(this_name,'}','');
                this_name = strrep(this_name,'/','');
                this_name = strrep(this_name,'\','');
                this_name = strrep(this_name,'.','');
                self.legal_parameter_names{i} = this_name;
            end

        end

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

        function [V, Ca, n] = integrate(self)
            switch self.implementation
            case 'SGS/MATLAB'
                [V, Ca, n] = self.integrate_SGS_MATLAB;
            case 'SGS/Julia'
                [V, Ca, n] = self.integrate_SGS_JULIA;
            otherwise
                error('Unknown implementation')
            end
                   
        end

        function [V1, V2, Ca2] = integrate_SGS_MATLAB(self)
            % pure MATLAB implementation of this model 

             % assemble parameters
            for i = 1:length(self.parameter_names)
                feval(@()assignin('caller',self.legal_parameter_names{i}, self.parameters(i)));
            end


            % define functions
            boltz = @(V,A,B) (1/(1 + exp((V+A)/B)));
            tauX = @(V,A,B,D,E) (A - B/(1+exp((V+D)/E)));

            % m_inf
            mNainf = @(V) (boltz(V,25.5,-5.29));
            mKinf = @(V) boltz(V,12.3,-11.8);
            mCaTinf = @(V) boltz(V,27.1,-7.2);
            mEAGinf = @(V,Ca) boltz(V,-22.28,-23.39).*(6./(6+Ca));


            % h_inf 
            hNainf = @(V) (boltz(V,48.9,5.18));
            hCaTinf = @(V) boltz(V,32.1,5.5);


            % tau_m
            taumNa = @(V) tauX(V,1.32,1.26,120.0,-25.0);
            taumK = @(V) tauX(V,7.2,6.4,28.3,-19.2);
            taumCaT = @(V) tauX(V,21.7,21.3,68.1,-20.5);
            taumEAG = @(V) tauX(V,-0.67,-26500,283,46.03);

            % tau_h
            tauhNa = @(V) (0.67/(1+exp((V+62.9)/-10.0)))*(1.5 + 1/(1+exp((V+34.9)/3.6)));
            tauhCaT = @(V) tauX(V,105.,89.8,55.,-16.9);

            nsteps = floor(self.t_end/self.dt);


            % pre-allocate arrays
            V1 = zeros(nsteps,1);
            mNa1 = 0;
            hNa1 = 0;
            mK1 = 0;
            mEAG1 = 0;


            V2 = zeros(nsteps,1);
            Ca2 = zeros(nsteps,1);
            mNa2 = 0;
            hNa2 = 0;
            mK2 = 0;
            mEAG2 = 0;
            mCaT2 = 0;
            hCaT2 = 0;

            % set some initial conditions 
            Ca2(1) = 0.05;

            dt = self.dt;
            R_by_zF = 500.0*(8.6174e-005);
            T = 10 + 273.15;
            RT_by_zF = R_by_zF*T;
            f = 14.96;
            fA = f*A2;
            exp_dt_by_tau_Ca = exp(-dt/tau_Ca);

            for i = 2:nsteps

                V_ext = self.V_drive(i-1);

                % compartment 1  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

                minf_Na = mNainf(V1(i-1));
                hinf_Na = hNainf(V1(i-1));
                minf_K  = mKinf(V1(i-1));
                minf_EAG = mEAGinf(V1(i-1),0.05);

                % integrate Na channel
                mNa1 = minf_Na + (mNa1 - minf_Na)*exp(-dt/taumNa(V1(i-1)));
                hNa1 = hinf_Na + (hNa1 - hinf_Na)*exp(-dt/tauhNa(V1(i-1)));

                % integrate K and EAG channels
                mK1 = minf_K + (mK1 - minf_K)*exp(-dt/taumK(V1(i-1)));
                mEAG1 = minf_EAG + (mEAG1 - minf_EAG)*exp(-dt/taumEAG(V1(i-1)));
                
                % compute effective conductances 
                gNa1 = g_Na1*(mNa1^3)*hNa1; 
                gK1 = g_K1*(mK1^4);
                gEAG1 = g_EAG1*(mEAG1^4); %% How do we know this exponent?

                sigma_g1 = gNa1 + gK1 + gEAG1 + g_leak; 

                % compute currents from electric coupling 
                I1 = (g_1*(V_ext - V1(i-1)) + g_1_2*(V2(i-1) - V1(i-1)))*1e-6; % in nA

                V_inf1 = (gNa1*self.E_Na + gK1*self.E_K + gEAG1*self.E_K + g_leak*self.E_leak + I1/A1)/sigma_g1;
                tau_v1 = C_m/(sigma_g1*10); % unit correction

                V1(i) = V_inf1 + (V1(i-1) - V_inf1)*exp(-dt/tau_v1);

                % compartment 2  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

                % compute Calcium reversal potential
                E_Ca = RT_by_zF*log(Ca_o/Ca2(i-1)); 

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
                gNa2 = g_Na2*(mNa2^3)*hNa2; 
                gK2 = g_K2*(mK2^4);
                gEAG2 = g_EAG2*(mEAG2^4); %% How do we know this exponent?
                gCa2 = g_Ca*(mCaT2^3)*hCaT2; 

                sigma_g2 = gNa2 + gK2 + gEAG2 + gCa2 + g_leak; 

                % compute currents from electric coupling 
                I2 = (g_2*(V_ext - V2(i-1)) + g_1_2*(V1(i-1) - V2(i-1)))*1e-6; % in nA

                V_inf2 = (gNa2*self.E_Na + gK2*self.E_K + gCa2*E_Ca + gEAG2*self.E_K + g_leak*self.E_leak +  I2/A2)/sigma_g2;
                tau_v2 = C_m/(sigma_g2*10); % unit correction

                V2(i) = V_inf2 + (V2(i-1) - V_inf2)*exp(-dt/tau_v2);

                % update Calcium 
                ca_current = gCa2*(V2(i-1) - E_Ca);
                cinf = Ca_i - fA*ca_current; 
                Ca2(i) = cinf + (Ca2(i-1) - cinf)*exp_dt_by_tau_Ca; 

            end
        end



        function N = integrate_SGS_JULIA(self)

            % assemble things to pass into Julia
            parameters = self.parameters(:);
            reversal_potentials  = [self.E_Na; self.E_K; self.E_leak];
            initial_conditions = self.initial_conditions(:);
            time_params = [self.dt; self.t_end];
            S = self.V_drive(:); 


            % run using Julia
            N = jl.call('intSynTerm',parameters,reversal_potentials, initial_conditions, time_params',S);


            % update initial conditions 
            %self.initial_conditions = N(end,:);
        

        end % end integrate Julia code


    

        function quitManipulateCallback(self,~,~)
            % destroy every object in self.handles
            d = fieldnames(self.handles);
            for i = 1:length(d)
                try
                    delete(self.handles.(d{i}))
                    self.handles = rmfield(self.handles,d{i});
                catch
                end
            end

        end % end quitManipulateCallback 


        function checkBounds(self)
            assert(~any(self.lb >= self.ub),'At least one lower bound is greater than a upper bound');
            assert(all(self.parameters <= self.ub) && all(self.parameters >= self.lb),'At least one parameter out of bounds');
            max_conductances = self.parameters(~cellfun(@isempty,cellfun(@(x) strfind(x,'g_'),self.parameter_names,'UniformOutput',false)));
            assert(min(max_conductances) >= 0,'Conductances cannot be < 0');
        end

        function manipulate(self)
            % check if a manipulate control window is already open. otherwise, create it
            make_gui = true;
            if isfield(self.handles,'manipulate_control')
                if isvalid(self.handles.manipulate_control)
                    make_gui = false;
                end
            end

            if make_gui
                Height = 1100;
                self.handles.manipulate_control = figure('position',[1000 250 1400 Height], 'Toolbar','none','Menubar','none','NumberTitle','off','IntegerHandle','off','CloseRequestFcn',@self.quitManipulateCallback,'Name',['manipulate[neuron]']);

                % add 4 axes
                dummy = NaN*(self.dt:self.dt:self.t_end);
                for i = 1:4
                    self.handles.ax(i) = axes;
                    self.handles.plot_handles(i) = plot(dummy,dummy,'k');
                    self.handles.ax(i).XLim = [0 self.t_end];
                end 
                self.handles.ax(1).Position = [.05 .55 .3 .4];
                self.handles.ax(2).Position = [.4 .55 .3 .4];
                self.handles.ax(3).Position = [.05 .05 .3 .4];
                self.handles.ax(4).Position = [.4 .05 .3 .4];
  
                % fix the Y-lim for all the voltage axes
                self.handles.ax(1).YLim = [-70 70];
                self.handles.ax(2).YLim = [-70 70];
                self.handles.ax(3).YLim = [-70 70];

                % plot the V_drive onto the first axes
                for i = 1:4
                    self.handles.plot_handles(i).XData = self.dt:self.dt:self.t_end;
                end
                self.handles.plot_handles(1).YData = (self.V_drive);


                % draw for the first time
                f = self.parameter_names;
                pvec = self.parameters;

                % make sure the bounds are OK
                checkBounds(self);
        
                Height = 1010;
                nspacing = 50;
                x_offset = self.t_end*2.5;
                for i = 1:length(f)
                    self.handles.control(i) = uicontrol(self.handles.manipulate_control,'Position',[1070 Height-i*nspacing 230 20],'Style', 'slider','FontSize',12,'Callback',@self.sliderCallback,'Min',self.lb(i),'Max',self.ub(i),'Value',pvec(i));
       
                    try    % R2013b and older
                       addlistener(self.handles.control(i),'ActionEvent',@self.sliderCallback);
                    catch  % R2014a and newer
                       addlistener(self.handles.control(i),'ContinuousValueChange',@self.sliderCallback);
                    end

                    % hat tip: http://undocumentedmatlab.com/blog/continuous-slider-callback
                    thisstring = [f{i} '= ',mat2str(self.parameters(i))];
                    self.handles.controllabel(i) = text(self.handles.ax(1),1,1,thisstring,'FontSize',20);
                    self.handles.lbcontrol(i) = uicontrol(self.handles.manipulate_control,'Position',[1305 Height-i*nspacing+3 40 20],'style','edit','String',mat2str(self.lb(i)),'Callback',@self.resetSliderBounds);
                    self.handles.ubcontrol(i) = uicontrol(self.handles.manipulate_control,'Position',[1350 Height-i*nspacing+3 40 20],'style','edit','String',mat2str(self.ub(i)),'Callback',@self.resetSliderBounds);


                end
                for i = 1:length(f)
                    self.handles.controllabel(i).Position = [x_offset 55-(i-1)*16];
                end


                % non stop integration
                while isfield(self.handles,'ax')
                    if isvalid(self.handles.ax)
                        [V1, V2, Ca] = self.integrate;
                        self.handles.plot_handles(2).YData = V1;
                        self.handles.plot_handles(3).YData = V2;
                        self.handles.plot_handles(4).YData = Ca;
                        drawnow limitrate
         
                    end
                end


            end % end if make-gui
        end % end manipulate 



        function sliderCallback(self,src,~)
            % evaluate
            this_param = find(self.handles.control == src);
            f = self.parameter_names;

            self.parameters(this_param) = src.Value;
            thisstring = [f{this_param} '= ',mat2str(src.Value)];

            % update the values shown in text
            self.handles.controllabel(this_param).String = thisstring;

        end % end sliderCallback

        function resetSliderBounds(self,src,~)
            checkStringNum(self,src.String);
            if any(self.handles.lbcontrol == src)
                % some lower bound being changed
                this_param = find(self.handles.lbcontrol == src);
                new_bound = str2double(src.String);
                self.lb(this_param) = new_bound;
                
                if self.handles.control(this_param).Value < new_bound
                    self.handles.control(this_param).Value = new_bound;
                    self.parameters.(self.parameter_names{this_param}) = new_bound;
                end
                checkBounds(self);
                self.handles.control(this_param).Min = new_bound;
            elseif any(self.handles.ubcontrol == src)
                % some upper bound being changed
                this_param = find(self.handles.ubcontrol == src);
                new_bound = str2double(src.String);
                self.ub(this_param) = new_bound;
                
                if self.handles.control(this_param).Value > new_bound
                    self.handles.control(this_param).Value = new_bound;
                    self.parameters.(self.parameter_names{this_param}) = new_bound;
                end
                checkBounds(self);
                self.handles.control(this_param).Max = new_bound;
            else
                error('error 142')
            end
        end % end resetSliderBounds

        function checkStringNum(self,value)
            assert(~isnan(str2double(value)),'Enter a real number')
        end % end checkStringNum

    end % end methods 

end % end classdef 



