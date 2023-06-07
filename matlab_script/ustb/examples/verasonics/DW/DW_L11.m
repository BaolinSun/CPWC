%% NB! This is a preliminary example, we will add DW from Verasonics to USTB soon :)
%
%
% Notice: 
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility 
%   for its use
%
% File name: SetUpL11_4vFlashAngles.m - Example of plane wave imaing with 
%                                       steering angle transmits
% Description: 
%   Sequence programming file for L11-4v Linear array, using plane wave
%   transmits with multiple steering angles. All 128 transmit and receive
%   channels are active for each acquisition. Processing is asynchronous
%   with respect to acquisition.
%
% Last update:
% 11/10/2015 - modified for SW 3.0

clear all
P.startDepth = 5;   % Acquisition depth in wavelengths
P.endDepth = 192;   % This should preferrably be a multiple of 128 samples.

% DW sequence
frames = 1;
PRF=6250;
alpha_max= 45;
radius = 20e-3;
na = 15;      % Set na = number of angles.
if (na > 1), 
    dtheta = (alpha_max*pi/180)/(na-1); 
    P.startAngle = -alpha_max*pi/180/2; 
else
    dtheta = 0; 
    P.startAngle=0; 
end 

% Define system parameters.
Resource.Parameters.numTransmit = 128;      % number of transmit channels.
Resource.Parameters.numRcvChannels = 128;    % number of receive channels.
Resource.Parameters.speedOfSound = 1540;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = 1;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.

% Specify Trans structure array.
Trans.name = 'L11-4v';
Trans.units = 'mm'; % Explicit declaration avoids warning message when selected by default
Trans = computeTrans(Trans);  % L11-4v transducer is 'known' transducer so we can use computeTrans.
Trans.maxHighVoltage = 50;  % set maximum high voltage limit for pulser supply.

% Specify PData structure array.
PData(1).PDelta = [Trans.spacing, 0, 0.5];
PData(1).Size(1) = ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3)); % startDepth, endDepth and pdelta set PData(1).Size.
PData(1).Size(2) = ceil((Trans.numelements*Trans.spacing)/PData(1).PDelta(1));
PData(1).Size(3) = 1;      % single image page
PData(1).Origin = [-Trans.spacing*(Trans.numelements-1)/2,0,P.startDepth]; % x,y,z of upper lft crnr.
% No PData.Region specified, so a default Region for the entire PData array will be created by computeRegions.

% Specify Media point targets at 0 0 20e-3 for verification purposes
lambda = Resource.Parameters.speedOfSound / (Trans.frequency*1e6);
Media.MP=[0, 0, 20e-3/lambda, 1;];

% Specify Resources.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = na*4096*2; % this size allows for maximum range
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = frames;    % 30 frames stored in RcvBuffer.
Resource.InterBuffer(1).numFrames = 1;   % one intermediate buffer needed.
Resource.ImageBuffer(1).numFrames = 10;
Resource.DisplayWindow(1).Title = 'L11-4vFlashAngles';
Resource.DisplayWindow(1).pdelta = 0.35;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)];   % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).numFrames = frames;
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow(1).Colormap = gray(256);

% Specify Transmit waveform structure.  
TW(1).type = 'parametric';
TW(1).Parameters = [Trans.frequency,.67,2,1];   % A, B, C, D

% Specify TX structure array.
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'Apod', kaiser(Resource.Parameters.numTransmit,1)', ...
                   'focus', -radius/lambda, ...
                   'Steer', [0.0,0.0], ...
                   'Delay', zeros(1,Trans.numelements)), 1, na);
% - Set event specific TX attributes.
if fix(na/2) == na/2       % if na even
    P.startAngle = (-(fix(na/2) - 1) - 0.5)*dtheta;
else
    P.startAngle = -fix(na/2)*dtheta;
end
for n = 1:na   % na transmit events
    TX(n).Steer = [(P.startAngle+(n-1)*dtheta),0.0];
    TX(n).Delay = computeTXDelays(TX(n));
end

% Specify TGC Waveform structure.
TGC.CntrlPts = [0,297,424,515,627,764,871,1000];
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

% Specify Receive structure arrays. 
% - We need na Receives for every frame.
maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
Receive = repmat(struct('Apod', ones(1,Trans.numelements), ...
                        'startDepth', P.startDepth, ...
                        'endDepth', maxAcqLength,...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', 'NS200BW', ...
                        'mode', 0, ...
                        'callMediaFunc', 0), 1, na*Resource.RcvBuffer(1).numFrames);
                    
% - Set event specific Receive attributes for each frame.
for i = 1:Resource.RcvBuffer(1).numFrames
    Receive(na*(i-1)+1).callMediaFunc = 1;
    for j = 1:na
        Receive(na*(i-1)+j).framenum = i;
        Receive(na*(i-1)+j).acqNum = j;
    end
end

% Specify Recon structure arrays.
% - We need one Recon structures which will be used for each frame. 
Recon = struct('senscutoff', 0.6, ...
               'pdatanum', 1, ...
               'rcvBufFrame',-1, ...
               'IntBufDest', [1,1], ...
               'ImgBufDest', [1,-1], ...
               'RINums', 1:na);

% Define ReconInfo structures.
% We need na ReconInfo structures for na steering angles.
ReconInfo = repmat(struct('mode', 'accumIQ', ...  % default is to accumulate IQ data.
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'regionnum', 1), 1, na);
% - Set specific ReconInfo attributes.
if na>1
    ReconInfo(1).mode = 'replaceIQ'; % replace IQ data
    for j = 1:na  % For each row in the column
        ReconInfo(j).txnum = j;
        ReconInfo(j).rcvnum = j;
    end
    ReconInfo(na).mode = 'accumIQ_replaceIntensity'; % accum and detect
else
    ReconInfo(1).mode = 'replaceIntensity';
end

% Specify Process structure array.
pers = 20;
Process(1).classname = 'Image';
Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum',1,...   % number of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'pdatanum',1,...    % number of PData structure to use
                         'pgain',1.0,...            % pgain is image processing gain
                         'reject',2,...      % reject level 
                         'persistMethod','simple',...
                         'persistLevel',pers,...
                         'interpMethod','4pt',...  %method of interp. (1=4pt)
                         'grainRemoval','none',...
                         'processMethod','none',...
                         'averageMethod','none',...
                         'compressMethod','power',...
                         'compressFactor',40,...
                         'mappingMethod','full',...
                         'display',1,...      % display image after processing
                         'displayWindow',1};

% Specify SeqControl structure arrays.
SeqControl(1).command = 'jump'; % jump back to start
SeqControl(1).argument = 1;
SeqControl(2).command = 'timeToNextAcq';  % time between synthetic aperture acquisitions
SeqControl(2).argument = 1/PRF*1e6;   
SeqControl(3).command = 'timeToNextAcq';  % time between frames
SeqControl(3).argument = 20000 - (na-1)*160;  % 20 msec
SeqControl(4).command = 'returnToMatlab';
nsc = 5; % nsc is count of SeqControl objects

% Specify Event structure arrays.
n = 1;
for i = 1:Resource.RcvBuffer(1).numFrames
    for j = 1:na                      % Acquire frame
        Event(n).info = 'Full aperture.';
        Event(n).tx = j;   % use next TX structure.
        Event(n).rcv = na*(i-1)+j;   
        Event(n).recon = 0;      % no reconstruction.
        Event(n).process = 0;    % no processing
        Event(n).seqControl = 2;
        n = n+1;
    end
    Event(n-1).seqControl = [3,nsc]; % modify last acquisition Event's seqControl
      SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
      nsc = nsc+1;

    Event(n).info = 'recon and process'; 
    Event(n).tx = 0;         % no transmit
    Event(n).rcv = 0;        % no rcv
    Event(n).recon = 1;      % reconstruction
    Event(n).process = 1;    % process
    Event(n).seqControl = 0;
    if floor(i/5) == i/5     % Exit to Matlab every 5th frame 
        Event(n).seqControl = 4; % return to Matlab
    else
        Event(n).seqControl = 0;
    end
    n = n+1;
end

Event(n).info = 'Jump back';
Event(n).tx = 0;        % no TX
Event(n).rcv = 0;       % no Rcv
Event(n).recon = 0;     % no Recon
Event(n).process = 0; 
Event(n).seqControl = 1;

% User specified UI Control Elements
% - Sensitivity Cutoff
UI(1).Control =  {'UserB7','Style','VsSlider','Label','Sens. Cutoff',...
                  'SliderMinMaxVal',[0,1.0,Recon(1).senscutoff],...
                  'SliderStep',[0.025,0.1],'ValueFormat','%1.3f'};
UI(1).Callback = text2cell('%SensCutoffCallback');

% - Range Change
MinMaxVal = [64,300,P.endDepth]; % default unit is wavelength
AxesUnit = 'wls';
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm');
        AxesUnit = 'mm';
        MinMaxVal = MinMaxVal * (Resource.Parameters.speedOfSound/1000/Trans.frequency);
    end
end
UI(2).Control = {'UserA1','Style','VsSlider','Label',['Range (',AxesUnit,')'],...
                 'SliderMinMaxVal',MinMaxVal,'SliderStep',[0.1,0.2],'ValueFormat','%3.0f'};
UI(2).Callback = text2cell('%RangeChangeCallback');

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 5;

% Save all the structures to a .mat file.
filename ='a';
save(filename);

% call VSX
VSX;

%% converting to UFF
disp('Storing data into UFF format');

%% Reading data
data=zeros(Receive(1).endSample, Resource.Parameters.numRcvChannels, na, Resource.RcvBuffer(1).numFrames);
probe_geometry=Trans.ElementPos(1:Resource.Parameters.numRcvChannels,1:3)*1e-3; % in      

%% calculate offset 
offset_distance=(TW.peak)*lambda;   % in [m]
if strcmp(Trans.units,'mm')
    offset_distance=offset_distance+2*Trans.lensCorrection*1e-3;
elseif strcmp(Trans.units,'wavelengths')
    offset_distance=offset_distance+2*Trans.lensCorrection*lambda;
end
offset_time=offset_distance/Resource.Parameters.speedOfSound;   % in [s]

%% wave delay
figure;
angles = P.startAngle:dtheta:-P.startAngle;
c0 = Resource.Parameters.speedOfSound;
t0_1=zeros(1,na);
for n_angle=1:na
    source = [TX(n_angle).focus*lambda*sin(angles(n_angle)) 0 TX(n_angle).focus*lambda*cos(angles(n_angle))]; 
    dst = sqrt(sum(bsxfun(@plus,source,-probe_geometry).^2,2));
    dst0 = sqrt(sum(bsxfun(@plus,source,[0 0 0]).^2,2));
    
    if TX(n_angle).focus>0
        delay=max(dst)/c0-dst/c0;
        t0_1(n_angle)=-(dst0-max(dst))/c0;
    else
        delay=dst/c0-min(dst)/c0;
        t0_1(n_angle)=(dst0-min(dst))/c0;
    end
    plot(probe_geometry(:,1),TX(n_angle).Delay*lambda/Resource.Parameters.speedOfSound-t0_1(n_angle),'b-'); grid on; hold on;
    plot(probe_geometry(:,1),delay-t0_1(n_angle),'r--'); grid on; hold on;
    plot(0,0,'bo');
    %pause();
end

%% convert data
n=1;
Fs=4*Trans.frequency*1e6;

t_out=0:(1/Fs):((Receive(1).endSample-1)/Fs);
for n_frame = 1:Resource.RcvBuffer(1).numFrames
    for n_angle = 1:na
        % compute time vector for this line
        t_ini=2*Receive(n).startDepth*lambda/c0;
        t_end=2*Receive(n).endDepth*lambda/c0;
        no_t=(Receive(n).endSample-Receive(n).startSample+1);
        t_in=linspace(t_ini,t_end,no_t)-offset_time-t0_1(n_angle);
        
        % read data
        data(:,:,n_angle,n_frame)=interp1(t_in,double(RcvData{1}(Receive(n).startSample:Receive(n).endSample,:,n_frame)),t_out,'linear',0);
        n=n+1;
        
        % check delays
        check=1;
        if check
            t00=-(20/Fs):(0.5/Fs):(20/Fs);
            z0=20e-3;
            x0=0;
            source = [TX(n_angle).focus*lambda*sin(angles(n_angle)) 0 TX(n_angle).focus*lambda*cos(angles(n_angle))]; 
            delay=  sqrt(z0^2+(probe_geometry(:,1)-x0).^2)/c0 + sign(TX(n_angle).focus).*( sqrt(sum((source - [0 0 0]).^2))/c0 - sqrt(sum((source - [x0 0 z0]).^2))/c0);
            delayeddata=zeros(128,length(t00));
            for nch=1:128
                delayeddata(nch,:)=interp1(t_out-delay(nch),data(:,nch,n_angle,n_frame),t00);
            end
            
            figure(102); hold off;
            pcolor(1:128,t00,delayeddata.'); shading flat; colormap gray; colorbar; hold on;
            plot(1:128,zeros(1,128),'r--');
            title(n_angle);
            drawnow;
            %pause();
        end
    end
end

% probe
probe = uff.linear_array();
probe.pitch = Trans.spacingMm*1e-3;
probe.N = Trans.numelements;
probe.element_width = Trans.elementWidth*1e-3;
probe.element_height = 5e-3;
h_fig = probe.plot(); hold on;

% sequence
sequence = repmat(uff.wave(),[1 na]);
for n_angle = 1:na
    sequence(n_angle).source.distance = TX(n_angle).focus*lambda;
    sequence(n_angle).source.azimuth = angles(n_angle);
    sequence(n_angle).source.plot(h_fig);
end

% Create channel_data object
channel_data = uff.channel_data();
channel_data.data=data;
channel_data.sampling_frequency = Fs;
channel_data.initial_time = t_out(1);
channel_data.sound_speed = Resource.Parameters.speedOfSound;
channel_data.sequence = sequence;
channel_data.probe = probe;

channel_data.name = 'DW dataset recorded on Verasonics Vantage 256 and L11 probe';
channel_data.version = 'v1.0';
channel_data.PRF = PRF;

%% SCAN
sca=uff.linear_scan();
sca.x_axis=lambda*(PData(1).Origin(1)+(0:PData(1).Size(2)-1)*PData(1).PDelta(1)).';
sca.z_axis=lambda*(PData(1).Origin(3)+(0:PData(1).Size(1)-1)*PData(1).PDelta(3)).';
 
%% BEAMFORMER
bmf=beamformer();
bmf.channel_data=channel_data;
bmf.scan=sca;

bmf.receive_apodization.window=uff.window.flat;
bmf.receive_apodization.f_number=2*tan(Recon.senscutoff);
bmf.receive_apodization.apex.distance=Inf;

bmf.transmit_apodization.window=uff.window.flat;
bmf.transmit_apodization.f_number=2*tan(Recon.senscutoff);
bmf.transmit_apodization.apex.distance=Inf;

% beamforming
b_data=bmf.go({process.das process.coherent_compounding});

%% show
b_data.plot();

%% write channel_data to path
% uff_filename = 'DW_simulation_L11';
% uff_file=uff([ustb_path '/data/' uff_filename],'write');
% uff_file.write(channel_data,'channel_data');
% uff_file.write(b_data,'beamformed_data');

%% Verasonics vs USTB beamforming
vb=abs(ImgData{1}(:,:,1,1));
ub=abs(reshape(b_data.data(:,1,1,1),[sca.N_z_axis sca.N_x_axis 1 1]));
figure;
subplot(1,2,1)
imagesc(20*log10(vb/max(vb(:)))); caxis([-60 0]); colormap gray; colorbar; title('Verasonics')
subplot(1,2,2)
imagesc(20*log10(ub/max(ub(:)))); caxis([-60 0]); colormap gray; colorbar; title('USTB')

return

% **** Callback routines to be converted by text2cell function. ****
%SensCutoffCallback - Sensitivity cutoff change
ReconL = evalin('base', 'Recon');
for i = 1:size(ReconL,2)
    ReconL(i).senscutoff = UIValue;
end
assignin('base','Recon',ReconL);
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'Recon'};
assignin('base','Control', Control);
return
%SensCutoffCallback

%RangeChangeCallback - Range change
simMode = evalin('base','Resource.Parameters.simulateMode');
% No range change if in simulate mode 2.
if simMode == 2
    set(hObject,'Value',evalin('base','P.endDepth'));
    return
end
Trans = evalin('base','Trans');
Resource = evalin('base','Resource');
scaleToWvl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);

P = evalin('base','P');
P.endDepth = UIValue;
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm');
        P.endDepth = UIValue*scaleToWvl;    
    end
end
assignin('base','P',P);

evalin('base','PData(1).Size(1) = ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3));');
evalin('base','PData(1).Region = computeRegions(PData(1));');
evalin('base','Resource.DisplayWindow(1).Position(4) = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);');
Receive = evalin('base', 'Receive');
maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
for i = 1:size(Receive,2)
    Receive(i).endDepth = maxAcqLength;
end
assignin('base','Receive',Receive);
evalin('base','TGC.rangeMax = P.endDepth;');
evalin('base','TGC.Waveform = computeTGCWaveform(TGC);');
evalin('base','if VDAS==1, Result = loadTgcWaveform(1); end');
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'PData','InterBuffer','ImageBuffer','DisplayWindow','Receive','Recon'};
assignin('base','Control', Control);
assignin('base', 'action', 'displayChange');
return
%RangeChangeCallback
