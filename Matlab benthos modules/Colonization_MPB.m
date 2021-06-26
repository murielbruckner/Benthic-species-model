%% module for microphytobenthos computations
% Author: Muriel Bruckner
% Final version 23/4/2021
% The module computes microphytobenthos settling and expansion for growth
% period specified in txt-file. Data extraction from d3d-output in case no
% bioturbation in simulation


GRID = wlgrid('read',[directory,'work\',ID1,'.grd']);

% read MPB data from txt-files
FID = fopen([directory, '\work\', 'Phyto1', '.txt']);
data_phyto = textscan(FID, '%7.3f%7.3f%7.3f%7.3f%7.3f%7.3f%7.3f', 'HeaderLines', 3);
fclose(FID);
tau_sed = data_phyto{1,2};  % tau crit of sediment
P ={}; % allocate matrix to save phyto fractions

% if no bioturbators present in simulation then first import of d3d-output
% data
if bioturbation==0
    Taucrit = ones(size(GRID.X)+1)'*tau_sed ;  % make initial matrix with abiotic value
    
    % extract D3D output
    Extract_parameters
    
else
    % open grid and critical bed shear stress data computed from colonization_benthos
    Taucrit = (wldep('read',[directory,'work\',ID1,'.tce'],GRID))';
    Taucrit(isnan(Taucrit==1))=0;  % remove nan's
end
%% compute microphytobenthos distribution

% if within growth period of MPB

if ets >=data_phyto{1,3} && ets <=data_phyto{1,4}
    
    
    % turbator parameters from txt-file
    tau_bio = data_phyto{1,1};      % new tau crit for colonized cells
    habitat_in1 = data_phyto{1,5};  % inundation threshold min
    habitat_in2 = data_phyto{1,6};  % inundation threshold max
    mud_th = data_phyto{1,7};       % mud threshold (if mud smaller then no habitat)
    
    % physical pressures from d3d matrix
    inundation = d3dparameters.Flooding(1).PerYear(ets,1);
    inundation = cell2mat(inundation);
    PHYTO=zeros(size(inundation)); % reset PHYTO matrix
    
    % mud content from d3d matrix
    mud_fract=d3dparameters.mudfract(1).PerYear(ets,1);
    mud_fract = cell2mat(mud_fract);
    
    % calculation of MPB distribution as absence/presence
    phyto_pres = find(inundation>habitat_in1 & inundation<habitat_in2 & mud_fract>mud_th); % find suitable habitat for MPB
    phyto_mat = Taucrit;                % load taucrit distribution
    phyto_mat(phyto_pres) = tau_bio;    % add biotic value in colonized cells
    PHYTO(phyto_pres) = 1;              % save presence in matrix for postprocessing
    
    % grazing computations
    if bioturbation==1
        
        % loop over bioturbators
        for nb=1:no_bio
            graze=GRAZE{nb};
            if graze==1

                if grazing_comp ==1 % mean tau
                    graze_mat=cell2mat(type{nb}.mort(5));                   % find tau crit of bioturbator
                    cells = find(graze_mat<tau_sed & phyto_mat>tau_sed);    %  find the cells that contain both MPB and bioturbator 
                    phyto_mat(cells)=(phyto_mat(cells)+graze_mat(cells))/2; % mean value of MPB by the bioturbation of grazer                    
                elseif grazing_comp == 2 % linear computation of tau
                    graze_mat=tau_sed-cell2mat(type{nb}.mort(5));        % substract sed value to get difference between abiotc-biotic case
                    cells = find(graze_mat>0 & phyto_mat>tau_sed);      % find the cells that contain both MPB and bioturbator                    
                    bioturb_fract = cell2mat(type{nb}.mort(4));          % extract bioturbator fractions
                    delta_tauMPB = (tau_bio-tau_sed)*(PHYTO-bioturb_fract); % variation in tau induced by the remaining fraction of MPB
                    phyto_mat(cells) = tau_sed + delta_tauMPB(cells) ;  % final values after grazing
                end
                
                % save the remaining MPB fraction
                rel_mat = PHYTO-cell2mat(type{nb}.mort(4));  % new MPB fractions (always >0)
                PHYTO(cells) = rel_mat(cells);               % add grazed MPB fractions
                
            end
        end

    end
    
    % save phyto distribution from ETS
    P{1}.distr(1) = {PHYTO};
    Taucrit = phyto_mat; % new value for new colonization
    
    % clear temporal matrix for clearing memory space
    clear inundation data_bio
    
else % if in winter season
    
    % reset phyto distribution
    P{1}.distr(1) = {zeros(size(GRID.X)+1)'};
    
end


% saving MPB fractions for output
savefile = ['P',num2str(ets)];
savefile = [directory, 'results_', num2str(year),'\', savefile];
save(savefile, 'P');
d3dparameters.taucrit(1).PerYear(ets,1) = {Taucrit};

% write file to work-folder and save for postprocessing
Taucrit(Taucrit==0)=nan;
wldep('write',[directory,'work\',ID1,'.tce'],'',Taucrit');
