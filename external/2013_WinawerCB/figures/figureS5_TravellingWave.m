%% Script to plot Figure S6, Traveling Wave Simulation
%
% Winawer, Kay, Foster, Parvizi, Wandell
% *Asynchronous broadband signals are the principal source of the BOLD
% response in human visual cortex*
% _Current Biology, 2013_
%
% This figure shows the predicted effect of a travelling of cortical
% activity on measured ECoG signals, which is the upper right panel in
% Figure 6S.
%
% The plot measures the amount of signal attenuation If a periodic cortical signal originates from a point on cortex
% and spreads at a fixed speed across cortex, how much signal attenuation
% and how much phase delay will there be in the aggregate signal measured
% by an electrode, as compared to the measurement from a coherent signal
% (no traveling wave)?
%
% Copyright Jonathan Winawer, 2013


%%
% Traveling wave assumed to be 0.3 meter per second
%  *Robustness of Traveling Waves in Ongoing Activity of Visual Cortex*
%  Nauhaus, Busse, Ringach, Carandini, _JNS 2012_
%  http://www.jneurosci.org/content/32/9/3088.full


%% Simulate travelling wave

% Propagation speed of traveling wave across cortex (~0.3 m/s in macaque)
speed = 0.01:0.01:1; % meters / second, or mm/ms

% Radius of antenna (electrode is ~1.2 mm radius)
radii = .5:.5:3;

% loop across radii and compute attenuated signal for each TW speed
amp = zeros(length(radii), length(speed));
ph  = zeros(length(radii), length(speed));

for r = 1:length(radii)
    
    % antenna function of electrode (mm, radius)
    radius = radii(r);
    
    % positions under electrode
    [x, y] = meshgrid(linspace(-radius,radius,101), linspace(-radius,radius,101));
    distance = sqrt(x.^2 + y.^2);
    distance = distance(:);
    
    % pointwise signal delay, in ms
    delay   = distance * (1./speed);
    
    % periodic signal is 15 Hz. period is 1000/15 ms
    period  = 1000/15;
    
    % pointwise phase lag, in radians
    phaseLag   = 2*pi*delay/period;
    
    % spatial integration under electrode (assume a complex sum of harmonics)
    ch1  = mean(cos(phaseLag));
    ch2  = mean(cos(phaseLag-pi/2));
    amp(r,:)  = sqrt(ch1.^2 + ch2.^2);
    ph(r,:)   = atan(ch2./ch1);
 
end

%% Plot 

% Set up figure window
fH = figure; clf; 
set(fH, 'Color', 'w'); 
set(gca, 'FontSize', 20, 'ColorOrder', jet(size(amp,1))); hold all;
xlabel('Speed (mm/ms)')
ylabel('Relative amplitude')

% Plot the signal amplitude as a function of traveling wave speed and
% antenna radius
plot(speed, amp, '-', 'LineWidth', 2); 

% Legend
h = legend({num2str(radii', '%3.1f')}, 'Location', 'Best', 'FontSize', 16);
v = get(h,'title');
set(v,'string','Antenna radius (mm)', 'FontSize', 16);
   
% Draw a line to indicate speed of TW from macaque literature 
plot([.3 .3], [0 1], 'k--')
text(0.31, 0.2, 'Speed of TW in macaque', 'FontSize', 14)

%% save
savepth = fullfile(ecogPRFrootPath, 'scratch');
hgexport(fH, fullfile(savepth, 'FigureS5_travelingWave.eps'))


%% Analysis 2 
% The vector plots below (not included in the publication) show the
% predicted modulation in phase and amplitude arising from a traveling wave
% of cortical activity assuming a periodic input at 15 Hz, and an antenna
% function of 1.15 mm.

%% (a) Assume no traveling wave  

% antenna function of electrode (mm, radius)
radius  = 1.15;        

% periodic signal is 15 Hz. period is 1000/15 ms
period  = 1000/15; 

% positions under electrode
[x, y] = meshgrid(linspace(-radius,radius,1001), linspace(-radius,radius,1001));
distance = sqrt(x.^2 + y.^2);

delay    = 0; 

% pointwise phase lag, in radians
phaseLag   = 2*pi*delay/period;

% spatial integration under electrode (assume a complex sum of harmonics)
ch1      = mean(cos(phaseLag(:)));
ch2      = mean(cos(phaseLag(:)-pi/2));
amp      = sqrt(ch1^2 + ch2^2);
ph       = atan(ch2/ch1);

fprintf('\nCoherent signal (no traveling wave).\n')
fprintf('Amplitude:%5.3f\nPhase:%5.3f\n',amp, ph)
fprintf('Peak shift:  %5.3f ms\n',period*ph/(2*pi))

fH = figure; set(fH, 'Color', 'w');
c = compass(ch1, ch2, 'g'); 
set(c, 'LineWidth', 2);
hold on

%% (b) Assume traveling wave of 0.3 meters per second 
%     = 300 mm per 1000 ms

% propagation speed across cortex (ms/mm)
speed = 0.3; % meters / second, or mm/ms
      
% pointwise signal delay, in ms
delay   = distance/speed; 

% pointwise phase lag, in radians
phaseLag   = 2*pi*delay/period;

% spatial integration under electrode (assume a complex sum of harmonics)
ch1  = mean(cos(phaseLag(:)));
ch2  = mean(cos(phaseLag(:)-pi/2));
amp  = sqrt(ch1^2 + ch2^2);
ph   = atan(ch2/ch1);
fprintf('\nTraveling wave %3.1f ms/mm.\n', 1/speed)
fprintf('Amplitude:%5.3f\nPhase:%5.3f\n',amp, ph)
fprintf('Peak shift:  %5.3f ms\n',period*ph/(2*pi))


c = compass(ch1, ch2, 'r'); 
set(c, 'LineWidth', 2);

%% (c) Assume slow traveling wave (sanity check)

speed    = 0.01; % m/s (or mm/ms)
delay    = distance/speed; 
phaseLag = 2*pi*delay/period;
ch1      = mean(cos(phaseLag(:)));
ch2      = mean(cos(phaseLag(:)-pi/2));
amp      = sqrt(ch1^2 + ch2^2);
ph       = atan(ch2/ch1);
fprintf('\nTraveling wave %4.1f ms/mm.\n', 1/speed)
fprintf('Amplitude:%5.3f\nPhase:%5.3f\n',amp, ph)
fprintf('Peak shift:  %5.3f ms\n',period*ph/(2*pi))

c = compass(ch1, ch2, 'b'); 
set(c, 'LineWidth', 2);

legend('No traveling wave', 'TW = 0.3 m/s',  'TW = 0.01 m/s')
