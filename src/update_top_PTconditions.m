global pw_top T_top cl_top
global pw_surf
global g dnw_top 
global T_time surface_T N_surface_T
global sea_level_time  sea_level N_sea_level
global water_depth_simul water_depth_local
global curr_time
global ROW

% water_depth_right could be read from a real sea level curve
% T_surf could be read from the measured air temperature history
% T_seafloor could be a function of water depth
relative_sl = interp(sea_level_time,sea_level,N_sea_level,curr_time,1);

water_depth_local = relative_sl - (-water_depth_simul);
T_groundsurface = interp(T_time,surface_T,N_surface_T,curr_time,1);
% Relative sea level is lower than the surface of the simulated
% location, the simulated location is exposed to Arctic air
for j=1:ROW
    if water_depth_local(1,j) <= 0 % Exposed to Arctic air
        T_top(1,j) = T_groundsurface;
        pw_top(1,j) = pw_surf;    % size = 1 x ROW
        cl_top(1,j) = 0.0;

    else % Flooded by sea water Chuvilin et al. (2022)
        T_top(1,j) = (water_depth_local(1,j)<=2)*0 +...
            (water_depth_local(1,j)>2)*(water_depth_local(1,j)<=30)*((water_depth_local(1,j)-2)*(-1.3/28))+...
            (water_depth_local(1,j)>30)*(-1.3);
        pw_top(1,j) = pw_surf + (water_depth_local(1,j)*dnw_top(1,j)*g);    % size = 1 x ROW
        cl_top(1,j) = 0.032;
    end
end


%% ----------------- Correct for initial sea level
if curr_time < 0.1
    
    relative_sl = 0;
    water_depth_local_temp = relative_sl - (-water_depth_simul);

    for j=1:ROW
        if water_depth_local_temp(1,j) <= 0 % Exposed to Arctic air
            pw_top(1,j) = pw_surf;    % size = 1 x ROW
            cl_top(1,j) = 0.0;
            T_top(1,j) = 0.0;
        else % Flooded by sea water Chuvilin et al. (2022)
            
            pw_top(1,j) = pw_surf + (water_depth_local_temp(1,j)*dnw_top(1,j)*g);    % size = 1 x ROW
            cl_top(1,j) = 0.032;
        end
    end

end
