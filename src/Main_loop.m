clc
clear all
global INDC2 INDC3
global ROW COL 
global pw T T_0
global sw sw_0 sg sg_0 si si_0
global cl cl_0 cm cm_0
global dnw dnw_0 dng dng_0
global phi phi_0
global krw krg
global pcgw
global lambda
global C_stable C_labile
global C_stable_old C_labile_old
global dt
global N timestep
global iteration back_track 
global curr_time
global go_back

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%      Load  Initialization File             %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------------ Option-1: run from initialization file ----------------
timestep=1;   % the start timestep

% ------------- Option-2: run from saved .mat file -------------------
% load('t311.91kyr.mat')              % load initial input from .mat file
% timestep=timestep+1;              % the start timestep=current timestep+1
% string_transient=('t311.91kyr.mat');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%      time discretizaton            %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt_ini=86400*365*20;           % Initial time step
dt_cut_inner=0.5;               % when cut time step, dt=dt*dt_cut
dt_cut_outter=0.1;              % when we go back to previous saved timestep, dt_ini=dt_ini*dt_cut_outter
dt_minimum=86400*365*0.01;                   % the smallest timestep is 0.1 year
N=2000000;                       % total number of timestep we will run
N_save=50;                      % we save our results every N_save timestep

N_display=N;                   % we print results every N_display timestep
j_show = 1;
i_show = 50;

if timestep==1
    time=zeros(N+1,1);          % initializing time for each step
    time(1)=0;                  % we start at time=zero
    % save_transient_results;
end


%% ------------ Option-1: run from initialization file ----------------
Initialization;                     % load initial file
string = ('t0kyr.mat');             % save initial input
save(string)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%   Main Subroutine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_iter=3;        % maximum number of iterations for new-raphson to converge
back_track = 1;
go_back=0;          % index of primariy variable shift or not
break_flag=0;

while timestep<=N
    disp(timestep)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%       Save the data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pw_old=pw;
    T_0=T;
    cm_0=cm;
    cl_0=cl;
    sw_0=sw;
    sg_0=sg;
    si_0=si;
    dnw_0=dnw;
    dng_0=dng;
    phi_0=phi;
    krw_old=krw;
    krg_old=krg;
    pcgw_old=pcgw;
    lambda_old=lambda;
    C_labile_old=C_labile;
    C_stable_old=C_stable;
    INDC2_old=INDC2;
    INDC3_old=INDC3;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % at new time step, we return to original dt
    max_iter=3;
    back_track=1;
    dt=dt_ini;
    time(timestep+1)=time(timestep)+dt;

    if time(timestep+1)>401*(86400*365*1e3)
        string = ['t',num2str(time(timestep+1)/(86400*365*1e3)),'kyr','.mat'];
        string_transient=string;
        save(string_transient)
        timestep = N+1;
        break
    end

    iteration=0;
    cycle=0;
    loop=1;
    while loop==1
        %*************** Set index values ****************************
        go_back=0;
        nan_flag=0;
        t_flag=0;
        iteration=iteration+1;

        %*************** Top Boundary Condition ****************************
        curr_time = time(timestep+1)/(86400*365*1e3);
        update_top_PTconditions;

        %************** Solve the PDEs *************************************
       
        newton_raphson;

        %************** Check if solution blows up ************************
        for j=1:ROW
            for i=1:COL
                if isnan(pw(i,j)) || pw(i,j)<0 || cl(i,j)<0 || cl(i,j)>1
                    nan_flag=1;
                    go_back=go_back+1;
                    iteration=max_iter+1;   % get out of the loop
                    t_flag=1;               % decrease the timestep
                end
            end
        end
        %************** Update Parameters **********************************
        if nan_flag==0
            update_secondary_variables;
        end

        %************** Newton Raphson Iteration reaches maximum **********
        if iteration > max_iter

            %************** Primary Variable Shift *******************
            if nan_flag==0
                cycle=cycle+1;
                primary_variables_switch;
                if go_back>0
                    max_iter=7;         % maximum number of iterations for new-raphson to converge
                    back_track =0.2^(1/max_iter);
                end
            else
                max_iter=9;         % maximum number of iterations for new-raphson to converge
                back_track =0.1^(1/max_iter);
            end

            %************** We fall into infite cycles *******************
            if cycle > 4 && nan_flag==0
                go_back = 0;
            end

            %************** Re-iteration or Go to next Time step **********
            if dt>=dt_minimum

                % ------- We need to cut timestep and rerun current timestep
                if t_flag==1
                    dt=dt_cut_inner*dt;
                    disp('We cut timestep size to:')
                    display(dt)
                    time(timestep+1)=time(timestep)+dt;
                    INDC2=INDC2_old;
                    INDC3=INDC3_old;

                    t_flag=0;
                    go_back=go_back+1;
                    cycle=0;
                end

                % -------- We need to rerun current timestep because of timestep cut or phases change
                if go_back>0
                    pw=pw_old;
                    T=T_0;
                    cm=cm_0;
                    cl=cl_0;
                    sw=sw_0;
                    sg=sg_0;
                    si=si_0;
                    phi=phi_0;
                    krw=krw_old;
                    krg=krg_old;
                    dnw=dnw_0;
                    dng=dng_0;
                    pcgw=pcgw_old;
                    lambda=lambda_old;
                    C_labile=C_labile_old;
                    C_stable=C_stable_old;

                    loop=1;
                    iteration=0;

                    % ------ Solution converged, go to next timestep
                else
                    loop=0;

                end % end of if go_back==1

                % ------We need to rerun from a previous saved results
            else
                disp('Timestep size is less than the minimum')
                break
            end % end of if dt>dt_timinum

        end % end of if iteration>max_iter

    end % end of while loop==1

    if loop==0

        %************** Clean Solution ****************
        clean_solutions;

        %************** Check the phase zone ****************
        phase_zone_update;

        %************** Calcualte diffusion coefficient ****************
        calc_diffusion_coefficient;

        %*********************** Save Results *******************
        if ~mod(timestep,N_save)
            if break_flag==1
                break_flag=0;
                dt_ini=dt_ini/dt_cut_outter;
            end
            string = ['t',num2str(time(timestep+1)/(86400*365*1e3)),'kyr','.mat'];
            string_transient=string;
            save(string_transient)
        end
        %*********************** Display Results *******************
        if ~mod(timestep,N_display)
            disp('Top Boundary Condition:')
            display([water_depth_local, pw_top/1e6, T_top])
            disp('depth    |    pw    |    T    |    sw    |    si    |   sg   |  cl  |  cm  |  INDC2  |  INDC3  |  C_stable')
            display([dpth(1:i_show,j_show)./1e3 pw(1:i_show,j_show)./1e6  T(1:i_show,j_show)   sw(1:i_show,j_show)  si(1:i_show,j_show) sg(1:i_show,j_show)  cl(1:i_show,j_show)   cm(1:i_show,j_show)   INDC2(1:i_show,j_show)   INDC3(1:i_show,j_show) C_stable(1:i_show,j_show)])
        end

        %*********************** Update timestep *******************
        timestep=timestep+1;
    end

    %*********************** Go Back to Previously Saved Step *************
    if dt<dt_minimum
        if timestep<=N_save
            load('../output/t0kyr.mat');
            timestep=1;
        else
            load(string_transient);
        end
        disp('we have gone back to previous saved timestep')
        dt_ini=dt_ini*dt_cut_outter;
        break_flag=1;

        if timestep>1
            timestep=timestep+1;
        end
    end

end












