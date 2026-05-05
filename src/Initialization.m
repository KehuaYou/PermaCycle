
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global COL ROW
global INDC2 INDC3
global model_top dpth dpth_groundsurface dx dy dz  vb 
global pw 
global T T_0
global sw sw_0 sg sg_0 si si_0
global cl cl_0 cm cm_0
global dnw dnw_0 dns dni g
global phi   phi_0
global kx 
global krw krg sgr  swr wn gn 
global vsw  vsg  viscosity_g viscosity_w
global pcgw  
global lambda lambda_g  lambda_w lambda_s lambda_i conductivity_index
global Cp_w Cp_s Cp_i
global qw qg qt cl_bottom
global pw_surf 
global pw_top T_top cm_top cl_top
global eps
global pw_ini T_ini pw_scale T_scale cl_scale
global dpth_below_gs
global Msalt dnw_top 
global Ax Swr_freeze
global T_time surface_T N_surface_T
global sea_level_time  sea_level N_sea_level
global water_depth_simul curr_time 
global total_grid_num
global p_salinity T_salinity salinity
global Mm Mh Mw N_hydration INDC1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%    Numerical Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ------------ numerical method  ----------------------
eps=1E-6;   % change in parameter when calculate derivative and the Jacobian matrix

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%    Methane Hydrate Phase Boundary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Setting up P,T,S diagram
load methane_hydrate_phase_boundary.mat
T_salinity=(-20:0.01:40);
salinity=(0:0.5:22)/100;

%% Some constant values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%    Constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ---------- density ------------
g=9.81;     % gravitational constant
dni=917;    % ice density at 0 degree c
dns=2750;   % density of solid grain
Msalt=0.05844;  % molecular weight of salt
Mm=0.016;   % molecular weight of methane 
Mh=0.1195;  % molecular weight of methane hydrate
Mw=0.018;
N_hydration=5.75; % hydration number
%% ---------- relative perm parameters ------------
sgr=0.02;  % residual gas saturation 
swr=0.15;       % residual water saturation 
wn=3.8;         % constant used in equations to calculate relative water permeability
gn=2.2;         % constant used in equations to calculate relative gas permeability

viscosity_g = 2.0e-5;   % dynamic viscosity of gas
viscosity_w = 1.31e-3;  % dynamic viscosity of water
%% ------------- Soil Water Characteristic Curve --------------------
Ax = 2.0;           % determine the slope of the curve, the larger Ax, the steep the curve is
Swr_freeze = 0.14;   % residual water saturation that cannot be frozen

%% ---------- capillary entry pressure ------------
Pd_entry=0.027839e6; % Capillary entry pressure in Pa

%% ------------ thermal constant ----------------------
Cp_w=4200;       % J/K*kg - heat capacity of water
Cp_s=735;        % J/K*kg - heat capacity of quartz,  Farouki OT 1981
Cp_i=2009;      %J/K*kg - heat capacity of ice 

lambda_g=0.03281; % W/m*K - thermal conductivity of gas
lambda_w=0.58; % W/m*K - thermal conductivity of water
lambda_s=3.9;  % W/m*K - thermal conductivity of solid grains
lambda_i=2.2;   % W/m*K - thermal conductivity of ice

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%    Geometery & Discretization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
water_depth_simul = 19.0;   % current day water depth of simulation
ROW=1;                      % number of horizontal grids
COL=50;                     % number of vertical grids

dx = 20.*ones (COL,ROW);    % vertical grid size
dy = 10.*ones(COL,ROW);     % Horizontal grid size
dz=1*ones(COL,ROW);         % Horizontal grid size
total_grid_num = COL;       % Total number of grid cells

model_top=dx(1,1)*0.5;        % The highest model grid, below its ground surface

% Vertical Depth
dpth_temp=zeros(COL,1);
dpth_temp(1)=model_top;
for i=1:COL-1
    dpth_temp(i+1)=dpth_temp(i)+(dx(i,1)+dx(i+1,1))*0.5;
end
dpth=dpth_temp;

% Volume of each grid cell
vb = dx.*dy.*dz;

% depth of groundsurface
dpth_groundsurface = 0; 
dpth_below_gs = dpth - dpth_groundsurface; 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%    Boundary conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ---------------- Top Boundary ---------------------

% Read in surface temperature and sea level curve -------
surface_temperature_history = readmatrix('surface_T.csv');
T_time = surface_temperature_history(:,1);
surface_T = surface_temperature_history(:,2);
N_surface_T = length(surface_T);

sea_level_curve = readmatrix('sealevel.xlsx');

sea_level_time = sea_level_curve(:,1);
sea_level = sea_level_curve(:,2);
N_sea_level = length(sea_level);


% Initialize top boundary conditions
pw_top = zeros(1,ROW);  % Pressure at the seafloor/groundsurface
T_top = zeros(1,ROW);   % Temperature at the seafloor/ground surface
cl_top = zeros(1,ROW);  % Salinity at the seafloor/ground surface
cm_top = zeros(1,ROW);  % Dissolved methane concentration at the seafloor/ground surface

dnw_top = 1028.1.*ones(1,ROW);     % water density at seafloor/ground surface
pw_surf = 0.1e6;    % Pressure at ground surface

% Calcualte the conditions at the top
curr_time = 0;  % current simulation time

update_top_PTconditions;

% ---------------- Bottom Boundary ---------------------      
qw = 0;                    % Bottom water flux in kg/m2/yr
qg=zeros(1,ROW);           % Bottom gas flux in kg/m2/yr 2e-6 m/s
qt=0.07 * ones(1,ROW);     % Geothermal heat flux in W/m2, Hornbach (2020), Figure 5 of Phrampus B.J. (2014)
cl_bottom = 0.032 * ones(1,ROW);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Initialization Conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% -----------------------initialize saturation--------------------------------
si_0 = zeros(COL,ROW);  % ice saturation
sg_0 = zeros(COL,ROW);  % free gas saturation
sw_0=1-si_0-sg_0;       % water saturation

sw=sw_0;
si=si_0;
sg=sg_0; 

%% --------------------------initial salinity-------------------------------
cl_0 = 0.032.*ones(COL,ROW);
cl=cl_0;

%% ----------------------- Porosity --------------------------------
% Best-fitting curve of measurements at Miline Point, North Slope of Alaska
% (Lee, 2005).
e_fold = 2226.0;            % e fold depth of porosity
phi_seafloor = 0.532;       % porosity at seafloor or ground surface
phi_0 = (phi_seafloor.*exp(-dpth_below_gs./e_fold)).*(dpth_below_gs>0) +...
    1.*(dpth_below_gs<=0);
phi=phi_0;

%% -----------------------Intrinsic permeability--------------------------------
kx=1e-15.*ones(COL,ROW);    % intrinsic permeability of sediment

calc_average_perm;

%% ----------------- initial  bulk thermal conductivity -------------------
conductivity_index = 2;  % parallel: 0; Series:1; Random: 2
lambda = zeros(COL,ROW); 
for i=1:COL
    for j=1:ROW
        lambda(i,j)= bulk_thermal_conductivity(phi(i,j), sw(i,j), si(i,j), sg(i,j));
    end
end

%% -----------------------initial temeprature-----------
T_grad1= qt(1,1) / lambda(COL,1); 
T=zeros(COL,ROW);
for j=1:ROW
    for i=1:COL
        T(i,j) = T_top(1,j) + T_grad1.*dpth_below_gs(i,j);
    end
end
T_0 = T;

%% ----------------------initial pressure ---------------------------
pw = zeros(COL,ROW);
for j=1:ROW
    pw(:,j) = pw_top(1,j) + dnw_top(1,j).*g.*dpth_below_gs(:,j); %initial hydrostatic pore pressure
end

%% ----------------------initial capillary pressure---------------------------
pcgw = Pd_entry.*ones(COL,ROW);  

%% ---------------initial dissolved methane concentration, density, viscosity and porosity------------------------------
vsw = viscosity_w .* ones(COL,ROW);  % water viscosity
vsg = viscosity_g .* ones(COL,ROW); % gas viscosity
for i=1:COL
    for j=1:ROW
        cm_0(i,j)=0;  % dissolved methane concentration
        dnw_0(i,j)=dnw_top(1,j);  % brine density
    end
end
cm=cm_0;
dnw=dnw_0;

%% ---------initial index (INDC2, INDC3) values----------------------
INDC1=zeros(COL,ROW);
for i=1:COL
    for j=1:ROW
        % INDC1=1, methane hydrate stable, INDC1=0, free gas stable or
        % three-phase coexist
        if cl(i,j) <interpolation2(p_salinity,T_salinity,salinity,pw(i,j)/1e6,T(i,j))
            INDC1(i,j)=1;
        else
            INDC1(i,j)=0;
        end
        
        % INDC2=1, one phase present; INDC2=3, two phases present; INDC2=5,
        % three phases present. This set excludes gas hydrate
        INDC2(i,j)=1;
        if (sw(i,j)>0)+(sg(i,j)>0)+(si(i,j)>0) ==2
            INDC2(i,j)=3;
        end
        if (sw(i,j)>0)+(sg(i,j)>0)+(si(i,j)>0) ==3
            INDC2(i,j)=5;
        end

        if T(i,j) > -cl(i,j)*(164.49*cl(i,j)+49.462)
            INDC3(i,j)=1;
        else
            INDC3(i,j)=4;
        end
        
    end
end

%% ------------ initial relative permeability -----------------------
for i=1:COL
    for j=1:ROW
        [krw(i,j),krg(i,j)] = relative_perm(sw(i,j),si(i,j),sg(i,j));
    end
end

%% --------initial flow and transport operator --------
% flow operator, diffusive operator and heat conduction operator
for i=1:COL
    for j=1:ROW
        calc_trans(i,j);
    end
end

%% ---------------------- diffusion coefficient ----------------------------
Dsaltx=zeros(COL,ROW);
calc_diffusion_coefficient; 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  Numerical Computation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------- For non-dimensionlization -------------------------
pw_ini=pw;
T_ini=T;
pw_scale=1025*g*500;
T_scale=20;
cl_scale=0.03;

%% ------------- Display Initial Condition --------------------------
j_show = 1;
disp('depth    |    pw    |    T    |    phi    |    sw    |    si    |    sg    |    cl    |    cm    |    INDC2    |    INDC3    ')
display([dpth(:,j_show)./1e3 pw(:,j_show)./1e6  T(:,j_show)   phi(:,j_show) sw(:,j_show)  si(:,j_show) sg(:,j_show)  cl(:,j_show)   cm(:,j_show)   INDC2(:,j_show)   INDC3(:,j_show)] )

