%% Benthos model coupled to hydromorphological model Delft3D
% This module defines the modules to be called. 
% After calling all relevant modules, D3D is started through the batch-file
% and the trim-data is saved in the results-folder.

%% Preprocessing
% set up D3D parameters from mdf-file
inid3d

% set up working directory and run initial run to create trim file with initial parameters
ini_work

%% Start simulation
% Start year loop
 for year = 1:years

     
    % create result folders per year for storage of output
    mkdir(strcat(directory,'results_', num2str(year), '\'));
    
        
    %% Start loop over benthos computations
    
    % loop over ecological time-steps (ETS)
     for ets = 1:t_eco_year 

        % handle Delft3D administration by reading and adjusting Delft3D
        % in- and output files such as time-steps
        d3dadmin  
       
        
        % computation of bioturbator distribution
        if bioturbation==1
            Colonization_benthos
        end
        
        
        % computation of microphytobenthos distribution
       if phyto == 1
            Colonization_MPB
       end
       
       
       %% Run Delft3D and save output
        run_line =  [directory, 'work\', 'run_flow2d3d_parallel.bat'];
        cd([directory, 'work']);
        system(run_line);

        
       % save d3d data from ref scenario
        if ets~=1 && bioturbation==0 && phyto==0
            Extract_parameters
        end
        
        % copy results to result folder for analysis
        
        % save one full results-file per year
        if ets==t_eco_year
            
            copyfile([directory, 'work\trim-', ID1, '.def'], [directory, 'results_', num2str(year), '/trim-', ID1, '_', num2str(ets),'.def']);
            copyfile([directory, 'work\trim-', ID1, '.dat'], [directory, 'results_', num2str(year), '/trim-', ID1, '_', num2str(ets),'.dat']);
           
        end
        
        
        % to reduce output sizes average Delft3D results file of each
        % coupling
%         average_trim([directory, 'work\trim-', ID1, '.dat'],[directory, 'results_', num2str(year), '/trim-avg', ID1, '_', num2str(ets)]);

     end % end loop over ecological timesteps
        
 
    % happy message if no crash over the entire simulation
    if year == years
        display('Yeah! You made it!!!');
    end

end % end year-loop
 



