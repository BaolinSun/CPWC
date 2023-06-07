classdef fresnel < handle
%fresnel   fresnel definition
%
%   See also PULSE, BEAM, PHANTOM, PROBE

%   authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
%   $Date: 2017/02/24 $

    %% public properties
    properties  (Access = public)
        phantom             % phantom class
        pulse               % pulse class
        probe               % probe class
        sequence            % collection of wave classes
        sampling_frequency  % sampling frequency [Hz]
        PRF                 % pulse repetition frequency [Hz]
    end
    
    %% dependent properties
    properties  (Dependent)   
        N_elements         % number of elements in the probe
        N_points           % number of points in the phantom
        N_waves            % number of waves
        N_events           % number of events (waves*frames)
        N_frames           % number of frames
    end
    
    %% private properties
    properties  (Access = private)   
        version='v1.0.7';  % fresnel version
    end
    
    %% constructor
    methods (Access = public)
        function h=fresnel()
            %fresnel   Constructor of fresnel class
            %
            %   Syntax:
            %   h = fresnel()
            %
            %   See also BEAM, PHANTOM, PROBE, PULSE                      
            
        end
    end
    
    %% set methods
    methods  
        function out_dataset=go(h)
            disp(sprintf('USTB''s Fresnel impulse response simulator (%s)',h.version));
            disp('---------------------------------------------------------------');
            
            %% checking we have all we need
            assert(numel(h.probe)>0,'The PROBE parameter is not set.');
            assert(numel(h.phantom)>0,'The PHANTOM parameter is not set.');
            assert(numel(h.pulse)>0,'The PULSE parameter is not set.');
            assert(numel(h.sequence)>0,'The SEQUENCE parameter is not set.');
            assert(numel(h.sampling_frequency)>0,'The SAMPLING_FREQUENCY parameter is not set.');
            
            % phantom dimensions
            if(length(h.phantom)>1)
                assert(mod(h.N_events,h.N_waves)==0,'The number of waves and phantoms does not macth. length(phantom)=length(waves)*frames');
            end
            
            % checking number of elements
            assert(h.probe.N_elements==h.sequence(1).N_elements,'Mismatch in the number of elements in probe and the size of delay and apodization vectors in beam');
           
            %% unwrapping the signal
            focusing_delay=zeros(h.N_elements,1,h.N_waves);
            apodization=zeros(h.N_elements,1,h.N_waves);
            for n_w=1:h.N_waves 
                focusing_delay(:,1,n_w)=h.sequence(n_w).delay_values;
                apodization(:,1,n_w)=h.sequence(n_w).apodization_values;
            end
            
            %% reference sound speed
            c0=h.phantom(1).sound_speed; % choose the first speed of sound as reference
            
            %% minimum distance for including geometric dispersion
            delta0=4*pi*0.1e-3;
            
            %% time vector
            max_range=0;
            min_range=Inf;
            for n_e=1:h.N_events
                for n_p=1:h.N_points
                    max_range=max([max_range; sqrt(sum((ones(h.N_elements,1)*h.phantom(n_e).points(n_p,1:3)-h.probe.geometry(:,1:3)).^2,2))]);
                    min_range=min([min_range; sqrt(sum((ones(h.N_elements,1)*h.phantom(n_e).points(n_p,1:3)-h.probe.geometry(:,1:3)).^2,2))]);
                end
            end
            time_1w=(min_range/c0-8/h.pulse.center_frequency/h.pulse.fractional_bandwidth+min(focusing_delay(:))):(1/h.sampling_frequency):(max_range/c0+4/h.pulse.center_frequency/h.pulse.fractional_bandwidth + max(focusing_delay(:)));                                                  % time vector [s]
            time_2w=(2*min_range/c0-8/h.pulse.center_frequency/h.pulse.fractional_bandwidth+min(focusing_delay(:))):(1/h.sampling_frequency):(2*max_range/c0+4/h.pulse.center_frequency/h.pulse.fractional_bandwidth+ max(focusing_delay(:)));                                               % time vector [s]
            N_samples=numel(time_2w);                                                                              % number of time samples

            % save the data into a CHANNEL_DATA structure
            out_dataset=uff.channel_data();
            out_dataset.probe=h.probe();
            out_dataset.pulse=h.pulse();
            out_dataset.phantom=h.phantom();
            out_dataset.sequence=h.sequence();
            out_dataset.sampling_frequency=h.sampling_frequency();
            out_dataset.sound_speed=c0;
            wave_delays=false;
            
            % check if there wave delays are imposed in the sequence
            % definition
            for n_w=1:h.N_waves
                if abs(h.sequence(n_w).delay)>0
                    wave_delays=true;
                    continue;
                end
            end
            if wave_delays
                out_dataset.initial_time=0;
            else
                out_dataset.initial_time=time_2w(1);
            end
            
            out_dataset.data=zeros(N_samples,h.N_elements,h.N_waves,h.N_frames);
            out_dataset.PRF=h.PRF;
                  
            % the frame loop
            tools.workbar();
            N=h.N_points*h.N_waves*h.N_frames;
            for n_f=1:h.N_frames
                
                % the wave loop
                for n_w=1:h.N_waves 
                
                    % support for static phantoms
                    if(h.N_events==1)
                       current_phantom=h.phantom;
                    else
                       current_phantom=h.phantom((n_f-1)*h.N_waves+n_w);
                    end

                    % wavenumber
                    k=2*pi*h.pulse.center_frequency/current_phantom.sound_speed;  

                    %% points loop
                    for n_p=1:h.N_points
                        % progress bar
                        n=(n_p+h.N_points*(n_w-1)+h.N_points*h.N_waves*(n_f-1));
                        if 1%mod(n,round(N/100))==1
                            tools.workbar(n/N,sprintf('Fresnel simulator (%s)',h.version),'USTB');
                        end
                        
                        % computing geometry relations to the point
                        distance  = sqrt(sum((h.probe.geometry(:,1:3)-ones(h.N_elements,1)*current_phantom.points(n_p,1:3)).^2,2));
                        theta     = atan2(current_phantom.x(n_p)-h.probe.x, current_phantom.z(n_p)-h.probe.z)-h.probe.theta;
                        phi       = asin((current_phantom.y(n_p)-h.probe.y)./distance)-h.probe.phi;

                        % directivity between probe and the point
                        directivity = sinc(k*h.probe.width/2/pi.*tan(theta)).*sinc(k*h.probe.height/2/pi.*tan(phi)./cos(theta));
                        
                        % delay between probe and the point
                        propagation_delay = distance/current_phantom.sound_speed;%

                        % attenuation (absorption & geometrical dispersion)
                        attenuation = 10.^(-current_phantom.alpha*(distance/1e-2)*(h.pulse.center_frequency/1e6)).*directivity.*delta0./(4*pi*distance);
                        
                        % computing the transmit signal 
                        delayed_time=ones(h.N_elements,1)*time_1w-(propagation_delay+h.sequence(n_w).delay_values)*ones(1,numel(time_1w));                
                        transmit_signal=sum(bsxfun(@times,h.pulse.signal(delayed_time),apodization(:,:,n_w).*attenuation),1);  

                        % computing the receive signal
                        delayed_time=ones(h.N_elements,1)*time_2w-propagation_delay*ones(1,N_samples);                
                        if wave_delays
                            delayed_time=delayed_time+h.sequence(n_w).delay-time_2w(1);
                        end
                        out_dataset.data(:,:,n_w,n_f)=out_dataset.data(:,:,n_w,n_f)+bsxfun(@times,interp1(time_1w,transmit_signal,delayed_time,'linear',0),attenuation).';
                        
                        % computing second order scattering
%                         extra_distance = sqrt(sum((bsxfun(@minus,current_phantom.points([1:n_p-1 n_p+1:h.N_points],1:3),current_phantom.points(n_p,1:3))).^2,2));
%                         extra_delay = extra_distance/current_phantom.sound_speed;
%                         extra_attenuation= 10.^(-current_phantom.alpha*(extra_distance/1e-2)*(h.pulse.center_frequency/1e6)).*delta0./(4*pi*extra_distance);
%                         for nnp=1:length(extra_distance)
%                             h.reverb(:,:,n_w,n_f)=h.reverb(:,:,n_w,n_f)+bsxfun(@times,interp1(time_1w,transmit_signal,delayed_time-extra_delay(nnp),'linear',0),attenuation.*extra_attenuation(nnp)).';
%                         end
                    end                   
                end
            end
            tools.workbar(1);
        end
    end
    
    %% set methods
    methods  
        function h=set.phantom(h,in_phantom)
            assert(isa(in_phantom,'uff.phantom'), 'The phantom is not a PHANTOM class. Check HELP PHANTOM.');
            h.phantom=in_phantom;
        end
        function h=set.pulse(h,in_pulse)
            assert(isa(in_pulse,'uff.pulse'), 'The pulse is not a PULSE class. Check HELP PULSE.');
            h.pulse=in_pulse;
        end
        function h=set.probe(h,in_probe)
            assert(isa(in_probe,'uff.probe'), 'The probe is not a PROBE class. Check HELP PROBE.');
            h.probe=in_probe;
        end
        function h=set.sequence(h,in_sequence)
            assert(isa(in_sequence,'uff.wave'), 'The sequence members are not a WAVE class. Check HELP WAVE.');
            h.sequence=in_sequence;
        end
        function h=set.sampling_frequency(h,in_sampling_frequency)
            assert(numel(in_sampling_frequency)==1, 'The sampling frequency must be a escalar');
            h.sampling_frequency=in_sampling_frequency;
        end       
        function h=set.PRF(h,in_PRF)
            assert(numel(in_PRF)==1, 'The PRF must be a escalar');
            h.PRF=in_PRF;
        end       
    end
    
    %% get methods
    methods  
        function value=get.N_elements(h)
            value=h.probe.N_elements;
        end
        function value=get.N_points(h)
            value=h.phantom.N_points;
        end
        function value=get.N_waves(h)
            value=numel(h.sequence);
        end
        function value=get.N_events(h)
            value=length(h.phantom);
        end
        function value=get.N_frames(h)
            value=ceil(h.N_events/h.N_waves);
        end
    end
    
end