% Plot the simulation results: temperature, pressure, salinity, ice
% saturation,ice and methane hydrate stability zones

clc
clear all

%% ====================== File name =================================
identifier = 't';
start = length(identifier)+1;
file_names = [identifier '*.mat'];
files = dir(file_names);
loc = zeros(length(files),1);
for i=1:length(files)
    loc(i) = str2double(files(i).name(start:end-7));
end
[B_mat , ix ] = sort(loc);

%% =========== Creat variables to store the results ===========
time_snaps = (0:1:401);
NN=length(time_snaps);
kkk=1;

load('t0kyr.mat')
time_collect = zeros(1,NN);  % in kyrs
pressure = zeros(COL+1,NN);    % in MPa
temperature = zeros(COL+1,NN);  % in oC
salt_concentration = zeros(COL+1,NN);     % in mass fraction
ice_saturation = zeros(COL+1,NN);  % dimensionless

permafrost_top = zeros(1,NN);  % in meters
permafrost_base = zeros(1,NN);  % in meters
hydrate_top = zeros(1,NN);     % in meters
hydrate_base = 1001.*ones(1,NN);     % in meters

%% =========== Obtain the results from saved files ===========
for ind = time_snaps

    work_progress = num2str(kkk/NN*100);
    disp([work_progress,'%'])

    if ind<max(B_mat)

        if max(B_mat==ind)
            this_file = files(ix(B_mat==ind)).name;
        else
            this_file = files(ix((abs(B_mat-ind))==(min(abs(B_mat-ind))))).name;
        end
        load(this_file);

        %******************************************************************************
        if kkk==1
            time_collect(1,kkk)=0;
        else
            time_collect(1,kkk)=time(timestep)/(86400*365*1e3);
        end


        pressure(2:COL+1,kkk) = pw./1e6;
        temperature(2:COL+1,kkk)= T;
        salt_concentration(2:COL+1,kkk)= cl;
        ice_saturation(2:COL+1,kkk)= si;

        pressure(1,kkk) = pw_top ./1e6;
        temperature(1,kkk) = T_top;
        salt_concentration(1,kkk) = cl_top;
        ice_saturation(1,kkk) = 0;

        for i=1:COL-1
            if si(i,1)==0 && si(i+1,1)>0 && permafrost_top(kkk) == 0
                permafrost_top(kkk)=dpth(i+1);
            end
            if si(i,1)>0 && si(i+1,1)==0
                permafrost_base(kkk)=dpth(i);
            end

            if INDC1(i,1)==0 && INDC1(i+1,1)==1
                hydrate_top(kkk)=dpth(i+1,1);
            end
            if INDC1(i,1)==1 && INDC1(i+1,1)==0
                hydrate_base(kkk)=dpth(i,1);
            end


            if hydrate_top(kkk) == 0
                hydrate_base(kkk) = NaN;
                hydrate_top(kkk) = NaN;
            end

        end

        kkk=kkk+1;
        %******************************************************************************

    end
end

%% ================ Display Some Key Input Parameters =====================

disp('Input Parameters:')
disp('Solid Thermal Conductivity | Geothermal Heat Flux (mW/m2) | Water Depth (m) | E Fold (km) | Minimum Unfrozen Water Saturation')
display([lambda_s, qt*1e3, water_depth_simul, e_fold/1e3, Swr_freeze])


%% =========== PLot the time evolution of variables =======================

%------------------- Time and Depth ---------------------------
dpth_plot=zeros(COL+1,NN);
dpth_plot(:,1) = [0 dpth'];
for i=2:NN
    dpth_plot(:,i)=dpth_plot(:,1);
end

time_collect_plot=zeros(COL+1,NN);
for j=1:COL+1
    time_collect_plot(j,:) = time_collect;
end

fig1=figure;
set(fig1,'Position',[20 50 2000 1000], 'Color','w')

%--------------------- Plot Temperature ----------------------------
subplot(2,2,1)

h=pcolor(time_collect_plot-400,dpth_plot,temperature);
set(gca,'Ydir','reverse','Tickdir','out')
set(h,'EdgeColor','none')
clim([-20 20])
colormap jet
colorbar

set(gca,'xlim',[-400 2])
set(gca,'ylim',[0 1000])
set(gca,'YTick',0:200:1000,'FontSize', 20)
set(gca,'XTick',-400:100:2,'FontSize', 20)
ylabel('Depth, mbsf', 'FontSize', 20)
xlabel('Time,kyr','FontSize', 20)
title('Temperature, (^oC)','FontSize', 20)


%--------------------- Plot Pore Pressure -------------------------------
subplot(2,2,2)

h=pcolor(time_collect_plot-400,dpth_plot,pressure);
set(gca,'Ydir','reverse','Tickdir','out')
set(h,'EdgeColor','none')
clim([0 15])
colormap jet
colorbar

set(gca,'xlim',[-400 2])
set(gca,'ylim',[0 1000])
set(gca,'YTick',0:200:1000,'FontSize', 20)
set(gca,'XTick',-400:100:2,'FontSize', 20)
ylabel('Depth, mbsf', 'FontSize', 20)
xlabel('Time,kyr','FontSize', 20)
title('Pore pessure, (MPa)','FontSize', 20)

%-------------------- Plot Salinity --------------------------------
subplot(2,2,3)
h=pcolor(time_collect_plot-400,dpth_plot,salt_concentration.*100);
set(gca,'Ydir','reverse','Tickdir','out')
set(h,'EdgeColor','none')
clim([0 20])
colormap jet
colorbar

set(gca,'xlim',[-400 2])
set(gca,'ylim',[0 1000])
set(gca,'YTick',0:200:1000,'FontSize', 20)
set(gca,'XTick',-400:100:2,'FontSize', 20)
ylabel('Depth, mbsf', 'FontSize', 20)
xlabel('Time,kyr','FontSize', 20)
title('Salinity, (wt.%)','FontSize', 20)

%--------------------- Plot Ice Saturation -------------------------------
subplot(2,2,4)
h=pcolor(time_collect_plot-400,dpth_plot,ice_saturation);
set(gca,'Ydir','reverse')
set(h,'EdgeColor','none')
clim([0 1.0])
colormap jet
cbh = colorbar;
cbh.Ticks = (0:0.2:1.0);
cbh.TickLabels = num2cell((0:0.2:1.0)) ;

set(gca,'xlim',[-400 2])
set(gca,'ylim',[0 1000])
set(gca,'YTick',0:200:1000,'FontSize', 20)
set(gca,'XTick',-400:100:2,'FontSize', 20)
ylabel('Depth, mbsf', 'FontSize', 20)
xlabel('Time,kyr','FontSize', 20)
title('Ice saturation, (-)','FontSize', 20)


%% ***************** Stability Zone plot *****************
fig2=figure;
set(fig2,'Position',[80 80 1500 800], 'Color','w')

plot(time_collect-400,permafrost_top,'b','LineWidth',2)
hold on
plot(time_collect-400,permafrost_base,'-.b','LineWidth',2)
plot(time_collect-400,hydrate_top,'g','LineWidth',2)
plot(time_collect-400,hydrate_base,'-.g','LineWidth',2)
set(gca,'FontSize',18,'LineWidth',1,'Ydir','reverse','Tickdir', 'out')
set(gca,'xlim',[-400 2])
set(gca,'ylim',[0 1000])
set(gca,'YTick',0:200:1000,'FontSize', 22)
set(gca,'XTick',-400:100:2,'FontSize', 22)
ylabel('Depth, mbsf', 'FontSize', 22)
xlabel('Time,kyr','FontSize', 22)
title('Ice & Methane Hydrate Stability Zones','FontSize', 22)
legend('permafrost top', 'permafrost base', 'hydrate top', 'hydrate base','FontSize',20,'Location','SouthWest')