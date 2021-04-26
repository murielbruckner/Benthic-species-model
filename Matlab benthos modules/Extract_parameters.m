%% Code to extract hydro-morphodynamic data
% This code reads the output data from the Delft3D run and computes habitat
% parameters of the species. These are then saved as matrices in a
% structure d3dparameters.mat, which is saved at the end of the
% morphological year. Based on these matrices, the species distribution is
% computed.


% extract water levels
WL = vs_get(NFS,'map-series','S1','quiet'); % Water level data at zeta points for all time-steps
waterdepth = cell(numel(WL),1); % preallocate matrix
flood = zeros(Ndim, Mdim); % preallocate matrix

% compute bed elevations for inundation period
if mor ==1
    depth                   = vs_get(NFS,'map-sed-series','DPS','quiet'); % bed topography with morphological changes
    depth_begin             = depth{1}; % Bed level matrix at begin of timestep trimmed to fit grid
    dts = length(depth);
else
    % in case of no morphological computations
    depth                   = vs_get(NFS,'map-const','DPS0','quiet'); % bed topography (without morphology) at zeta points trimmed to fit grid
end
% save data in postprocessing matrix
d3dparameters.depth(1).PerYear(ets,1)= {depth};


% extract mud fraction
mud = vs_get(NFS,'map-sed-series','FRAC','quiet');
mud = mud{dts, 1};
mud_fract = mud(:,:,2); % extract silt fraction
d3dparameters.mudfract(1).PerYear(ets,1)= {mud_fract};

% extract salinity
sal=vs_get(NFS,'map-series', 'R1','quiet');
S_all = zeros(Ndim, Mdim); % salinity

for dts=1:numel(sal) % determine water depth of all saved hydrol ts
    SAL = sal{dts,1};
    S1 = SAL(:,:,1, 1); % structure of matrix: M,N, layers, salinity
    S_all = S_all+S1; % compute mean salinity by addition of all time-steps
end
S_all = S_all/dts; % mean

% save in post-processing matrix
d3dparameters.salinityC(1).PerYear(ets,1) = {S_all};


% calculate water depth from water levels and bathymetry
for dts=1:numel(WL) % determine water depth of all saved hydrol ts
    if mor==1
        waterdepth{dts,1}= depth{dts,1}+WL{dts,1}; % sum depth (+) and water level (-)
        % write flooded days to matrix
    else
        waterdepth{dts,1}= depth+WL{dts,1}; % sum depth (+) and water level (-)
    end
    b=zeros(Ndim, Mdim); % preallocate matrix
    fl=find(waterdepth{dts,1}>fl_dr); % water depth has to be higher than fl-drying threshold
    if dts==1
        flood(fl)=1; % for flooded cells set =1
    else
        b(fl)=1;
        flood=flood+b; % sum of all time-steps flooded
    end
end
flooding_current=flood./dts; % mean

% save in post-processing matrix
d3dparameters.Flooding(1).PerYear(ets,1)={flooding_current};


% save in post-processing matrix
waterdepth=struct2mat(waterdepth,2);
waterdepth_mean = mean(waterdepth,3);
d3dparameters.waterdepth(1).PerYear(ets,1)= {waterdepth_mean};


% velocities
U1                      = vs_get(NFS,'map-series','U1','quiet'); % extract U velocity in U point
V1                      = vs_get(NFS,'map-series','V1','quiet'); % extract V velocity in V point
U1Mat                   = struct2mat(U1, 2);         % putting U1 from d3d-output in 3D matrix
V1Mat                   = struct2mat(V1, 2);         % putting V1 from d3d-output in 3D matrix
xy_velocity_new         = sqrt(U1Mat.^2 + V1Mat.^2); % calculate combined velocity in U and V direction using Pythagoras for each cell for both time-steps (ts)

% compute 90% quantile of velocity
velocity_max = quantile(xy_velocity_new,0.9,3);      % find 0.9 quantile nd write to 2D-matrix

% save in post-processing matrix
d3dparameters.VelocityMax(1).PerYear(ets,1)={velocity_max};

% save output for postprocessing at the end of year
if ets==t_eco_year
    
    savefile = ['d3d',num2str(ets)];
    savefile = [directory, 'results_', num2str(year),'/', savefile];
    save(savefile, 'd3dparameters');
end