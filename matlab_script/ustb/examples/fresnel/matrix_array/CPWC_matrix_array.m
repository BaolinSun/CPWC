%% CPWC simulation on a matrix array with the USTB built-in Fresnel simulator
%
% In this example we show how to use the built-in fresnel simulator in USTB
% to generate a Coherent Plane-Wave Compounding (CPWC) dataset using a 
% 2D matrix array probe and illustrate how it can be beamformed with USTB.
%
% This tutorial assumes familiarity with the contents of the 
% <./CPWC_linear_array.html 'CPWC simulation with the USTB built-in Fresnel 
% simulator'> tutorial. Please feel free to refer back to that for more 
% details.
% 
% by Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no> and Arun 
% Asokan Nair <anair8@jhu.edu> 29.03.2017_

%%
% 
% Clear the memory of any lingering settings and data, close all 
% previously opened plots and define the F-Number we will be using.

clear all;
close all;

F_number=1.7;

%% Phantom
%
% We start off defining an appropriate *phantom* structure to image. 
% Our phantom here is a single point scatterer. USTB's implementation 
% of *phantom* comes with a *plot* method to visualize the phantom for free!

pha=uff.phantom();
pha.sound_speed=1540;               % speed of sound [m/s]
pha.points=[0, 0, 20e-3, 1];...     % point scatterer position [m]
fig_handle=pha.plot();             
             
%% Probe
%
% The next UFF structure we look at is *probe*. It contains information 
% about the probe's geometry. USTB's implementation of *probe* comes with a 
% *plot* method too. When combined with the previous figure we can see the
% position of the probe respect to the phantom.
% Here, we shall set the probe to be a *matrix_array* type.

prb=uff.matrix_array();
prb.N_x=16;                 % number of elements in azimuthal direction
prb.N_y=16;                 % number of elements in elevational direction
prb.pitch_x=600e-6;         % probe pitch in azimuth [m]
prb.pitch_y=600e-6;         % probe pitch in elevation [m]
prb.element_width=570e-6;   % element width [m]
prb.element_height=570e-6;  % element height [m]
prb.plot(fig_handle);

%% Pulse
% 
% We then define the pulse-echo signal which is done here using the 
% *fresnel* simulator's *pulse* structure. We could also use 
% <http://field-ii.dk/ 'Field II'> for a more accurate model.

pul=uff.pulse();
pul.center_frequency=3e6;         % transducer frequency [MHz]
pul.fractional_bandwidth=0.6;     % fractional bandwidth [unitless]
pul.plot([],'2-way pulse');

%% Sequence generation
%
% Now, we shall generate our sequence! Keep in mind that the *fresnel* simulator
% takes the same sequence definition as the USTB beamformer. In UFF and
% USTB a sequence is defined as a collection of *wave* structures. 
% 
% For our example here, we define a sequence of plane waves at various 
% transmit angles (as defined in tx_angles). Apodization is set to none.

alpha_max=1/2/F_number;                     % Maximum (in absolute value) plane wave angle [rad]
N=15;                                       % Number of plane waves
tx_angles=linspace(-alpha_max,alpha_max,N); % Transmit angles for plane waves [rad]
seq=uff.wave();
for n=1:N
    seq(n)=uff.wave();
    seq(n).probe=prb;
    
    seq(n).source.azimuth=tx_angles(n);
    seq(n).source.distance=Inf;
    
    seq(n).sound_speed=pha.sound_speed;
    
    % show source
    fig_handle=seq(n).source.plot(fig_handle);
end

%% The Fresnel simulator
%
% Finally, we launch the built-in simulator. The simulator takes in a
% *phantom*, *pulse*, *probe* and a sequence of *wave* structures along 
% with the desired sampling frequency, and returns a *channel_data* UFF 
% structure.

sim=fresnel();

% setting input data 
sim.phantom=pha;                % phantom
sim.pulse=pul;                  % transmitted pulse
sim.probe=prb;                  % probe
sim.sequence=seq;               % beam sequence
sim.sampling_frequency=41.6e6;  % sampling frequency [Hz]

% we launch the simulation
channel_data=sim.go();
 
%% Scan
%
% The scan area is defines as a collection of pixels spanning our region of 
% interest. For our example here, we use the *linear_3D_scan* option, 
% which is defined with three components: the x-dimension (azimuthal) range, 
% the z-dimension (depth) range and the y-dimension (elevational) range. 
% *scan* too has a useful *plot* method it can call.

scan=uff.linear_3D_scan('radial_axis',linspace(-4e-3,4e-3,200).','axial_axis',linspace(18e-3,22e-3,100).','roll',0);
scan.plot(fig_handle,'Scenario');    % show mesh
 
%% Pipeline
%
% With *channel_data* and a *scan* we have all we need to produce an
% ultrasound image. We now use a USTB structure *beamformer*, that takes an
% *apodization* structure in addition to the *channel_data* and *scan*.

mid=midprocess.das();
mid.dimension = dimension.both;
mid.channel_data=channel_data;
mid.scan=scan;

mid.receive_apodization.window=uff.window.tukey50;
mid.receive_apodization.f_number=F_number;

mid.transmit_apodization.window=uff.window.tukey50;
mid.transmit_apodization.f_number=F_number;

b_data=mid.go();
b_data.plot();

%% 
%
% The *pipeline* structure allows you to implement different beamformers 
% by combination of multiple built-in *processes*. By changing the *process*
% chain other beamforming sequences can be implemented. It returns yet 
% another *UFF* structure: *beamformed_data*.

% To achieve the goal of this example, we use delay-and-sum (implemented in 
% the *das_matlab()* process) as well as coherent compounding.

