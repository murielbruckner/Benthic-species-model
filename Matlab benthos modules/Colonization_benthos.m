%% Compute bioturbation
% Author: M. Bruckner 21/4/2021
% computations include extraction of hydro-morphodynamic parameters,
% species parameters and computation of species distributions based on
% inundation period and mud content and optionally flow velocity and
% salinity. The latter two can be turned on in the text-files by adding
% minimum and maximum thresholds. When not needed, replace the minimum thresholds with -99,
% this will cut the loop of the computations. Competition is computed at
% the end and saved as tce and ero files, which are exported to the main
% model.


% open grid-file and tce-files
GRID = wlgrid('read',[directory,'work/',ID1,'.grd']);
Taucrit = (wldep('read',[directory,'work/',ID1,'.tce'],GRID))'; % read taucrit data to create pre-allocated matrices

% preallocation to reset matrices/variables
type ={};
Tau_final =zeros(size(Taucrit));
Ero_final =zeros(size(Taucrit));

% Extraction of environmental parameters from Delft3D files
Extract_parameters


% open physical pressures from d3d matrix
inundation = d3dparameters.Flooding(1).PerYear(ets,1);
velo_max   = d3dparameters.VelocityMax(1).PerYear(ets,1);
mud_fract  = d3dparameters.mudfract(1).PerYear(ets,1);
S_all      = d3dparameters.salinityC(1).PerYear(ets,1);

%% read bioturbator parameters from text-file

% loop to open benthos text-files to count number of species
for nn = 1:20 % maximum 20 species
    matFilename = sprintf('Bioturbator%d.txt', nn);
    Check = exist([directory, 'work/',matFilename],'file'); % if file is present value =2, else zero
    if Check ==2
        no_bio = nn; % save number of veg types in seperate vector
    else
        nn=21;
    end
end
clear Check


% open each bioturbator-file
for nb=1:no_bio % over several bioturbators
    
    % reset distribution matrices of bioturbator
    growth_matrix1=zeros(size(mud_fract));
    growth_matrix2=zeros(size(mud_fract));
    growth_matrix3=zeros(size(mud_fract));
    
    % read text-file to extract 1 growth parameters, 2 habitat parameters
    % 1
    FID = fopen([directory, '/work/', 'Bioturbator', num2str(nb), '.txt']);
    data_bio = textscan(FID, '%2.1f%5.4f%7.5f%5.4f%7.5f%4.1f%4.1f%4.1f%7.3f%7.3f%7.3f%7.3f%7.3f', 'HeaderLines', 4);
    fclose(FID);
    % 2
    FID = fopen([directory, '/work/', 'Bioturbator', num2str(nb), '.txt']);
    data_mort = textscan(FID, '%4.2f%4.2f%4.2f%4.2f%4.1f%4.1f%4.3f%4.1f%4.1f%4.1f', 'HeaderLines', 6);
    fclose(FID);
    
    % extract growth parameters
    ets_growth1 = data_bio{1,7}; % ETS where growth begins
    ets_growth2 = data_bio{1,8}; % ETS where growth ends
    
    tau_sed = data_bio{1,4};  % abiotic tau crit of sediment
    ero_sed = data_bio{1,5};  % abiotic erosion parameter of sediment
    
    
    % check if in growth period
    if ets >= ets_growth1 && ets <= ets_growth2
        
        % turbator parameters from txt-file
        tau_bio = data_bio{1,2};  % new tau crit for colonized cells
        ero_bio = data_bio{1,3};  % new erosion parameter for colonized cells
        rand_bio = data_bio{1,9};  % factor of random colonization (if ~= 1 than randomness included)
        
        
        % species parameters
        growth_rate = data_mort{1,1};  % if growth within cells is included
        max_biomass = data_mort{1,4};  % maximum relative biomass (here =1)
        grazing = data_mort{1,5};      % grazing on biofilms if phyto > 0
        
        % habitat parameters
        habitat_in1 = data_mort{1,2};  % inundation threshold small
        habitat_in2 = data_mort{1,3};  % inundation threshold large
        
        mud_perc1 = data_bio{1,11};    % minumim mud percentage in bed
        mud_perc2 = data_bio{1,12};    % maximum mud percentage in bed
        
        habitat_sal1 = data_mort{1,6};  % velocity threshold small
        habitat_sal2 = data_mort{1,7};  % velocity threshold large
        
        habitat_vel1 = data_mort{1,8};  % velocity threshold small
        habitat_vel2 = data_mort{1,9};  % velocity threshold large
        
        % calculation of turbator distribution
        
        % distribution based on inundation period
        growth_matrix1 = growth_distr_bioturb(inundation{1,1},habitat_in1,habitat_in2);
        growth_matrix1(inundation{1,1}==0) = 0; % remove supratidal cells from computations
        
        % distribution based on salinity (if included)
        if habitat_sal1==-99
            growth_matrix2 = ones(size(inundation{1,1}));
        else
            growth_matrix2 = growth_distr_bioturb(S_all{1,1},habitat_sal1,habitat_sal2);
        end
        
        % distribution based on mud fraction
        growth_matrix3 = growth_distr_bioturb(mud_fract{1,1},mud_perc1,mud_perc2);
        
        % compute minimum fraction of inundation and salinity and mud
        a=min(growth_matrix1,growth_matrix2);
        growth_all=min(a,growth_matrix3);
        
        % distribution based on flow velocity
        if habitat_vel1~= -99
            growth_matrix4 = growth_distr_bioturb(velo_max{1,1},habitat_vel1,habitat_vel2);
            
            % compute minimum based on inundation, salinity and mud and
            % velocity
            growth_all=min(growth_all,growth_matrix4);
            type{nb}.mort(7) = {growth_matrix4};  % save for post-processing
        end
        
        % include random colonization based on factor defined in
        % text-files
        bioturb = randsample(find(growth_all>0),round(length(find(growth_all>0))/rand_bio)); % random colonization if rand_bio ~= 1
        
        % save distribution in matrix
        bio_distribution = zeros(size(growth_matrix1));
        bio_distribution(bioturb)=growth_all(bioturb);
        bio_distribution(bio_distribution>1)=1;
        
        
        
        % calculation growth if defined
        bio_distribution(bio_distribution>0) = bio_distribution(bio_distribution>0)*(1+growth_rate);
        
        % save all matrices for post-processing
        type{nb}.mort(1) = {growth_matrix1};
        type{nb}.mort(2) = {growth_matrix2};
        type{nb}.mort(3) = {growth_matrix3};
        type{nb}.mort(4) = {bio_distribution};
        
        
        %% computations of new tau crit and erosion parameter
        
        % compute relative bioturbation linearly from biomass
        % tau crit
        tau_diff = tau_bio-tau_sed;
        m = tau_diff/max_biomass; % slope function
        bioturb_rate = tau_sed + bio_distribution.*m; % find values for taucrit
        
        Taucrit = bioturb_rate; % new tau critical
        clear bioturb_rate
        
        % erosion parameter
        ero_diff = ero_bio-ero_sed;
        m = ero_diff/max_biomass; % slope function
        bioturb_rate = ero_sed + bio_distribution.*m; % find values for erosion parameter
        
        EroPar=bioturb_rate; % new erosion parameter
        clear bioturb_rate
        
        
        % save keyword for grazing per species to compute in MPB
        % computations
        GRAZE{nb}=grazing;
        
    else % if in winter season
        bio_distr11 = zeros(size(inundation{1,1})); % clear matrix when not in growth period
        bio_distr22 = zeros(size(inundation{1,1})); % clear matrix when not in growth period
        bio_distribution = zeros(size(inundation{1,1})); % clear matrix when not in growth period
        Taucrit = ones(size(GRID.X)+1)'*tau_sed;   % reset critical bed shear stress matrix
        EroPar  = ones(size(GRID.X)+1)'*ero_sed;   % reset erosion parameter matrix
    end
    
    % save tau crit and erosion parameter per type for post-processing
    type{nb}.mort(5) = {Taucrit};
    type{nb}.mort(6) = {EroPar};
    
    % remove zeros from matrices to format for Delft3D computations
    Taucrit(Taucrit==0)=nan;
    EroPar(EroPar==0)=nan;
    % save matrices in structure with all species
    TAUCRIT{nb} = Taucrit;
    EROPAR{nb} = EroPar;
    
end



% for several bioturbators find all cells where both species are present

try
    competition = find(cell2mat(type{1}.mort(4))==0 & cell2mat(type{2}.mort(4))>0);
catch
    competition = 0;
end


% loop over all species where the second species is the recessive one
for nb = 1:no_bio
    
    % species one is the dominant one
    if nb==1
        Tau_final = TAUCRIT{nb};
        Ero_final = EROPAR{nb};
        
    else % recessive one
        
        % check if in growth period
        if ets >= ets_growth1 && ets <= ets_growth2
            
            
            tau = TAUCRIT{nb};
            ero = EROPAR{nb};
            
            % replace empty cells with bioturbator 2
            Tau_final(competition) = tau(competition); % replace empty cells wtih bioturbator 2
            Ero_final(competition) = ero(competition);

        end
        
    end
    % postprocessing
    type{nb}.mort(9) = {Tau_final};
    type{nb}.mort(10) = {Ero_final};
    
end



% saving vecotr with positions to test working
savefile = ['bio_type_',num2str(ets)];
savefile = [directory, 'results_', num2str(year),'/', savefile];
save(savefile, 'type');

% write file to work-folder which is used in computations
wldep('write',[directory,'work/',ID1,'.tce'],'',Tau_final');
wldep('write',[directory,'work/',ID1,'.ero'],'',Ero_final');


