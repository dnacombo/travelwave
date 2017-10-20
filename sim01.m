%% Simulating travelling waves
% In this script, we will try to simulate traveling waves in a model head
% and see what methods we can use to recover the simulated activity.
%
% Our rationale is the following:
% # Create a sensor array and head model
% # Simulate two dipoles at distance in this head
% # Apply a LCMV beamformer analysis to recover the filters associated with
% the two dipoles
% # Recover the signal at each location.
% 
% This script is based on http://www.fieldtriptoolbox.org/example/compute_forward_simulated_data_and_apply_a_beamformer_scan

addpath(fullfile(fileparts(mfilename('fullpath')),'misc'));
rm_frompath('[Ff]ield[Tt]rip')
addpath(fullfile(fileparts(mfilename('fullpath')),'fieldtrip'));
ft_defaults% add fieldtrip to the path and run ft_defaults

%% create a gradiometer array and head volume
% magnetometers at 12cm distance from the origin
[pnt, tri] = icosahedron162;
pnt = pnt(pnt(:,3)>=0,:);
grad.coilpos = 12*pnt;
grad.coilori = pnt;
for i=1:length(pnt)
  grad.label{i} = sprintf('chan%03d', i);
end

% spherical volume conductor with 10cm radius
hm      = [];
hm.r    = 10;
hm.o    = [0 0 0];

%% Simulate dipoles
% sufficiently far appart and pointing to different axis directions

cfg = [];
cfg.headmodel = hm;
cfg.grad = grad;

cfg.triallength = 10;
cfg.fsample = 500;

cfg.dip.pos = [0 4 4
               0 -4 4];    
cfg.dip.mom = [0 1 0
               0 -1 0]';

cfg.dip.amplitude   = [1 1];
cfg.dip.phase       = [0 pi/2];
cfg.dip.frequency   = [10 10];

nsamples = round(cfg.triallength*cfg.fsample);
time                = (0:(nsamples-1))/cfg.fsample;
% adding signal here to be able to plot it later
for i=1:size(cfg.dip.pos,1)
    cfg.dip.signal(i,:) = cos(cfg.dip.frequency(i)*time*2*pi + cfg.dip.phase(i)) * cfg.dip.amplitude(i);
end

cfg.relnoise = 10;

cfg.ntrials = 20;
cfg.randomseed = 123;
data = ft_dipolesimulation(cfg);

dip = cfg.dip;

%% compute the data covariance matrix for the LCMV
% will capture the activity of the simulated dipoles
cfg = [];
cfg.covariance = 'yes';
timelock = ft_timelockanalysis(cfg, data);
%% compute a FFT for the topoplot
cfg = [];
cfg.method = 'mtmfft';
cfg.output = 'pow';
cfg.foi    = 10;
cfg.taper  = 'hanning';
freq = ft_freqanalysis(cfg,data);

%% do the beamformer source reconstuction on a 2 cm grid
cfg = [];
cfg.headmodel = hm;
cfg.grad = grad;
cfg.grad.unit='cm'; %error otherwise
cfg.grid.resolution = 2;
cfg.method = 'lcmv';
cfg.keepfilter = 'yes'; % we keep filters to be able to recover time courses
source = ft_sourceanalysis(cfg, timelock);

%% plot head model
figure(589);clf
ft_plot_vol(hm,'vertexcolor','none','facecolor','none','edgecolor','k','facealpha',.5);
ft_plot_mesh(source.pos(source.inside,:));
ft_plot_topo3d(grad.coilpos,freq.powspctrm,'facealpha',.5,'contourstyle','black','isolines',linspace(min(freq.powspctrm),max(freq.powspctrm),5));
ft_plot_sens(grad);
hold on
k = dsearchn(source.pos,dip.pos);
quiver3(source.pos(k,1),source.pos(k,2),source.pos(k,3),dip.mom(1,:)',dip.mom(2,:)',dip.mom(3,:)','r','linewidth',2)
view(120,30)
rotate3d on

%% compute time courses and measure phase shifts between pairs of dipoles
tc1 = source.avg.filter{k(1)} * timelock.avg;
% tc1 = cellfun(@(x)source.avg.filter{k(1)} * x,data.trial,'uniformoutput',0);
tc2 = source.avg.filter{k(2)} * timelock.avg;

fig(13,333);clf
toi = timepts([0 1],timelock.time);
for i = 1:3
    subplot(1,3,i)
    h = plot(timelock.time(toi),tc1(i,toi));
    hold on;
    plot(timelock.time(toi),dip.mom(i,1) * dip.signal(1,toi),'color',get(h,'color'),'linestyle',':','linewidth',2);
    h = plot(timelock.time(toi),tc2(i,toi));
    plot(timelock.time(toi),dip.mom(i,2) * dip.signal(2,toi),'color',get(h,'color'),'linestyle',':','linewidth',2);
    ylim([-1.2 1.2])
end

%% compute the neural activity index, i.e. projected power divided by
% projected noise
cfg = [];
cfg.powmethod = 'none'; % keep the power as estimated from the data covariance, i.e. the induced power
source = ft_sourcedescriptives(cfg, source);
source.avg.nai=source.avg.pow./source.avg.noise;  %neural activity index calculation

source = rmfield(source,'time');

cfg = [];
cfg.method = 'ortho';
cfg.funparameter = 'nai';
cfg.funcolorlim = [1.5 2];  % the voxel in the center of the volume conductor messes up the autoscaling
ft_sourceplot(cfg, source);



