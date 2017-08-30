% synapticTerminal.m

classdef synapticTerminal < handle & matlab.mixin.CustomDisplay

    properties
        dt = 50e-3 % ms
        t_end =  5000 % ms

        initial_conditions = [-50 -50 0.05 0 0 0 0 0 0 0 0];

        % these propoerties are for manipulate
        handles % for GUI elements
        % order is as follows:
        %     cm  A extCa Ca0  tauCa f   conductances(8)
        parameter_names = {'C_{m}','A1','A2','Ca_{o}','Ca_{i}','tau_{Ca}','f'  ,'g_{Na}1','g_{Na}2','g_{Ca}','g_{K}1','g_{K}2','g_1','g_2','g_1_2','g_{leak}'}
        parameters      = [10      .06   .06   3000     0.05      200      15    183       183        2.3      2.7     24.6     98    61     1.01    .0993]
        lb              = [0       0     0     1        0.01      1        1     0         0          0        0        0       0     0       0      0];
        ub              = [100     1     1     1e4      100       1000     100   300       300        300      300      300     300   300     300    300];
        parameter_units = {'nF/mm^2','mm^2','mm^2','uM','uM',    'ms','uM/nA','mS/cm^2','mS/cm^2','mS/cm^2','mS/cm^2','mS/cm^2','mS/cm^2','mS/cm^2','mS/cm^2','mS/cm^2'}


        E_Na = 50;
        E_K = -80;
        E_leak = -50;

        V_drive

    end

    properties (Access = protected)
        legal_parameter_names
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

            % add my Julia code 
            p = fileparts(which(mfilename));
            jl.include(joinPath(p,'param_types.jl')) 
            jl.include(joinPath(p,'liu_gating_functions.jl'))
            jl.include(joinPath(p,'synapticTerminal.jl'))
            jl.include(joinPath(p,'neuronMex.jl'))
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
            idx = [idx1 idx2];
            assert(length(idx) == 1, 'Could not determine parameter to set')
            self.parameters(idx) = value;
        end

        function set.V_drive(self,value)
            time = self.dt:self.dt:self.t_end;
            assert(length(value) == length(time),'V_drive does not match the timing info set in dt and t_end')
            self.V_drive = value;
        end


        function N = integrate(self)

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
                Height = 900;
                self.handles.manipulate_control = figure('position',[1000 250 1400 Height], 'Toolbar','none','Menubar','none','NumberTitle','off','IntegerHandle','off','CloseRequestFcn',@self.quitManipulateCallback,'Name',['manipulate[neuron]']);

                % add an axes to show the plots
                self.handles.ax = axes;
                self.handles.ax.Position = [.05 .05 .67 .4];
                xlabel(self.handles.ax,'Time (s)');
                ylabel(self.handles.ax,'V_{m} (mV)');
                self.handles.V_trace = plot(self.handles.ax,NaN,NaN,'k','LineWidth',1.5);
                self.handles.ax.YLim = [-70 70];

                self.handles.ax2 = axes;
                self.handles.ax2.Position = [.05 .5 .67 .4];
                self.handles.C_trace = plot(self.handles.ax2,NaN,NaN,'k','LineWidth',1.5);

                % draw for the first time
                % draw for the first time
                f = self.parameter_names;
                pvec = self.parameters;

                % make sure the bounds are OK
                checkBounds(self);
        
                Height = 830;
                nspacing = 58.5;
                for i = 1:length(f)
                    self.handles.control(i) = uicontrol(self.handles.manipulate_control,'Position',[1070 Height-i*nspacing 230 20],'Style', 'slider','FontSize',12,'Callback',@self.sliderCallback,'Min',self.lb(i),'Max',self.ub(i),'Value',pvec(i));
       
                    try    % R2013b and older
                       addlistener(self.handles.control(i),'ActionEvent',@self.sliderCallback);
                    catch  % R2014a and newer
                       addlistener(self.handles.control(i),'ContinuousValueChange',@self.sliderCallback);
                    end

                    % hat tip: http://undocumentedmatlab.com/blog/continuous-slider-callback
                    thisstring = [f{i} '= ',mat2str(self.parameters(i))];
                        
                    self.handles.controllabel(i) = text(self.handles.ax,1150,70,thisstring,'FontSize',20);
                    self.handles.lbcontrol(i) = uicontrol(self.handles.manipulate_control,'Position',[1305 Height-i*nspacing+3 40 20],'style','edit','String',mat2str(self.lb(i)),'Callback',@self.resetSliderBounds);
                    self.handles.ubcontrol(i) = uicontrol(self.handles.manipulate_control,'Position',[1350 Height-i*nspacing+3 40 20],'style','edit','String',mat2str(self.ub(i)),'Callback',@self.resetSliderBounds);


                end
                for i = 1:length(f)
                    self.handles.controllabel(i).Position = [1150 230-(i-1)*23];
                end

                % while the window is valid, simulate in 2-second blocks and show it
                window_size = self.t_end/self.dt;
                self.handles.ax.XLim = [0 window_size];
                self.handles.ax2.XLim = [0 window_size];


    
                T = self.dt:self.dt:self.t_end;
                self.handles.V_trace.XData = T;

                while isfield(self.handles,'ax')
                    if isvalid(self.handles.ax)

                        N = self.integrate;
                        
        
                        self.handles.V_trace.YData = N(:,1);


                        drawnow limitrate
                    end
                end


            end % end if make-gui
        end % end manipulate 

        function sliderCallback(self,src,~)
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



