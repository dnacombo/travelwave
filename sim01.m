% This example script shows you how to create some simulated channel-level
% MEG data with a single dipole at a specified location in the head.
% Subsequently it does a beamformer source reconstruction to localize that
% source.
add_fieldtrip

% create a gradiometer array with magnetometers at 12cm distance from the origin
[pnt, tri] = icosahedron162;
pnt = pnt(pnt(:,3)>=0,:);
grad.coilpos = 12*pnt;
grad.coilori = pnt;
for i=1:length(pnt)
  grad.label{i} = sprintf('chan%03d', i);
end

% create a spherical volume conductor with 10cm radius
vol.r = 10;
vol.o = [0 0 0];

%% note that beamformer scanning will be done with a 1cm grid, so you should
% not put the dipole on a position that will not be covered by a grid
% location later
cfg = [];
cfg.vol = vol;
cfg.grad = grad;
cfg.triallength = 10;
cfg.fsample = 500;
cfg.dip.pos = [0 0 4
    1 0 4];    
cfg.dip.mom = [1 0 0
    1 0 0]';
cfg.dip.phase = [0 pi/2];
cfg.dip.frequency = [10 10];
cfg.relnoise = 10;
cfg.ntrials = 20;
%   cfg.dip.pos     = [Rx Ry Rz] (size Nx3)
%   cfg.dip.mom     = [Qx Qy Qz] (size 3xN)
%   cfg.dip.frequency    in Hz
%   cfg.dip.phase        in radians
%   cfg.dip.amplitude    per dipole
%   cfg.ntrials          number of trials
%   cfg.triallength      time in seconds
%   cfg.fsample          sampling frequency in Hz
%   cfg.relnoise    = add noise with level relative to simulated signal
%   cfg.absnoise    = add noise with absolute level
%   cfg.randomseed  = 'yes' or a number or vector with the seed value (default = 'yes')
data = ft_dipolesimulation(cfg);
dip = cfg.dip;

%% compute the data covariance matrix, which will capture the activity of
% the simulated dipole
cfg = [];
cfg.covariance = 'yes';
timelock = ft_timelockanalysis(cfg, data);

%% do the beamformer source reconstuction on a 1 cm grid
cfg = [];
cfg.vol = vol;
cfg.grad = grad;
cfg.grad.unit='cm'; %error otherwise
cfg.grid.resolution = 1;
cfg.method = 'lcmv';
cfg.lcmv.projectnoise='yes'; %needed for neural activity index
cfg.keepfilter = 'yes';
source = ft_sourceanalysis(cfg, timelock);

%% plot head model
figure(589);clf
ft_plot_vol(vol,'vertexcolor','none','facecolor','none','edgecolor','k','facealpha',.5);
ft_plot_mesh(source.pos(source.inside,:));
hold on
k = dsearchn(source.pos,dip.pos);
for i = 1:size(dip.pos,1)
    plot3(source.pos(k(i),1),source.pos(k(i),2),source.pos(k(i),3),'r.','markersize',10)
end
rotate3d on

%% compute time courses and measure phase shifts between pairs of dipoles
tc1 = source.avg.filter{k(1)} * timelock.avg;
% tc1 = cellfun(@(x)source.avg.filter{k(1)} * x,data.trial,'uniformoutput',0);
tc2 = source.avg.filter{k(2)} * timelock.avg;

figure(333);clf
plot(tc1(1,:));
hold on;
plot(tc2(1,:));


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



