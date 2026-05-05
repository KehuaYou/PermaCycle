function Rt=calc_Rt(i,j)
global pw pg  pcgw  dpth vb
global twx twx1 ttx 
global sw_0 dnw dnw_0 si_0  T T_0
global phi phi_0 dt
global INDC2
global dns
global Cp_w Cp_s
global qt dy dz
global T_top pw_top dnw_top dpth_groundsurface
global dni Cp_i
global dx kx krw vsw g lambda
global COL

%---------- some parameters ------------
T0=273.15;
Li=334000;                      % latent heat for ice change to water
pg=pw+pcgw;

%---------- some simulation control -----------
adv_opt=1;   % if it =0, means no heat advection by water flow
T_opt=1;    % T_opt=1, only temperature change, T_opt=2, full derivative of accumulation term

%% Initialize Residual
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rt=0;

%% Flow and Transport
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if i==1
    if adv_opt==1
        if twx(i,j)*(pw(i+1,j)-pw(i,j)) > twx1(i,j)*(dpth(i+1,j)-dpth(i,j))
            Rt=Rt + twx(i,j)*(pw(i+1,j)-pw(i,j))*dnw(i+1,j)*Cp_w*(T0+T(i+1,j));
            Rt=Rt - twx1(i,j)*(dpth(i+1,j)-dpth(i,j))*dnw(i+1,j)*Cp_w*(T0+T(i+1,j));
        else
            Rt=Rt + twx(i,j)*(pw(i+1,j)-pw(i,j))*dnw(i,j)*Cp_w*(T0+T(i,j));
            Rt=Rt - twx1(i,j)*(dpth(i+1,j)-dpth(i,j))*dnw(i,j)*Cp_w*(T0+T(i,j));
        end
    end
    Rt=Rt + ttx(i,j)*(T(i+1,j)-T(i,j));
    
elseif i==COL
    if adv_opt==1
        if twx(i-1,j)*(pw(i,j)-pw(i-1,j)) > twx1(i-1,j)*(dpth(i,j)-dpth(i-1,j))
            Rt=Rt - twx(i-1,j)*(pw(i,j)-pw(i-1,j))*dnw(i,j)*Cp_w*(T0+T(i,j));
            Rt=Rt + twx1(i-1,j)*(dpth(i,j)-dpth(i-1,j))*dnw(i,j)*Cp_w*(T0+T(i,j));
        else
            Rt=Rt - twx(i-1,j)*(pw(i,j)-pw(i-1,j))*dnw(i-1,j)*Cp_w*(T0+T(i-1,j));
            Rt=Rt + twx1(i-1,j)*(dpth(i,j)-dpth(i-1,j))*dnw(i-1,j)*Cp_w*(T0+T(i-1,j));
        end
    end
    Rt=Rt - ttx(i-1,j)*(T(i,j)-T(i-1,j));
    
else
    if adv_opt==1
        if twx(i,j)*(pw(i+1,j)-pw(i,j)) > twx1(i,j)*(dpth(i+1,j)-dpth(i,j))
            Rt=Rt + twx(i,j)*(pw(i+1,j)-pw(i,j))*dnw(i+1,j)*Cp_w*(T0+T(i+1,j));
            Rt=Rt - twx1(i,j)*(dpth(i+1,j)-dpth(i,j))*dnw(i+1,j)*Cp_w*(T0+T(i+1,j));
        else
            Rt=Rt + twx(i,j)*(pw(i+1,j)-pw(i,j))*dnw(i,j)*Cp_w*(T0+T(i,j));
            Rt=Rt - twx1(i,j)*(dpth(i+1,j)-dpth(i,j))*dnw(i,j)*Cp_w*(T0+T(i,j));
        end
    end
    Rt=Rt + ttx(i,j)*(T(i+1,j)-T(i,j));
    
    if adv_opt==1
        if twx(i-1,j)*(pw(i,j)-pw(i-1,j)) > twx1(i-1,j)*(dpth(i,j)-dpth(i-1,j))
            Rt=Rt - twx(i-1,j)*(pw(i,j)-pw(i-1,j))*dnw(i,j)*Cp_w*(T0+T(i,j));
            Rt=Rt + twx1(i-1,j)*(dpth(i,j)-dpth(i-1,j))*dnw(i,j)*Cp_w*(T0+T(i,j));
        else
            Rt=Rt - twx(i-1,j)*(pw(i,j)-pw(i-1,j))*dnw(i-1,j)*Cp_w*(T0+T(i-1,j));
            Rt=Rt + twx1(i-1,j)*(dpth(i,j)-dpth(i-1,j))*dnw(i-1,j)*Cp_w*(T0+T(i-1,j));
        end
    end
    Rt=Rt - ttx(i-1,j)*(T(i,j)-T(i-1,j));
    
end


%% Boundary Condition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note that if flow into the grid, use "+" sign; if flow out, use "-"
%---------------------Top boundary condition-----------------------------
if i == 1
    if adv_opt==1
        twx_temp = dy(i,j)*dz(i,j)/(0.5*dx(i,j))*kx(i,j)*krw(i,j)/vsw(i,j);
        twx1_temp = twx_temp * dnw(i,j) * g;     
        if (pw(i,j) - dnw(i,j)*g*dpth(i,j)) > (pw_top(1,j) - dnw(i,j)*g*dpth_groundsurface(1,j)) % if hydraulic head at node i+1 is greater than at node i
            Rt=Rt - twx_temp*(pw(i,j)-pw_top(1,j))*dnw(i,j)*Cp_w*(T0+T(i,j));
            Rt=Rt + twx1_temp*(dpth(i,j)-dpth_groundsurface(1,j))*dnw(i,j)*Cp_w*(T0+T(i,j));
        else
            Rt=Rt - twx_temp*(pw(i,j)-pw_top(1,j))*dnw_top(1,j)*Cp_w*(T0+T_top(1,j));
            Rt=Rt + twx1_temp*(dpth(i,j)-dpth_groundsurface(1,j))*dnw_top(1,j)*Cp_w*(T0+T_top(1,j));
        end
    end
    ttx_temp = dy(i,j)*dz(i,j)/(0.5*dx(i,j))*lambda(i,j);
    Rt=Rt - ttx_temp*(T(i,j)-T_top(1,j));
end
%-------------------Bottom boundary condition -----------------------------
if i == COL
    Rt=Rt + dy(i,j)*dz(i,j)*qt(1,j);   % Flow into the domain
end


%% Accumulation Term
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sw_unfrozen = calc_sw_unfrozen(i,j);
sw_temp = (INDC2(i,j)==1) * 1.0 + (INDC2(i,j)>1) * sw_unfrozen;    
si_temp = (INDC2(i,j)==1) * 0.0 + (INDC2(i,j)>1) * (1-sw_unfrozen);

switch T_opt
    case 1
        Rt=-Rt + ...
            vb(i,j)*( (phi(i,j)*(sw_temp*dnw(i,j)*Cp_w + si_temp*dni*Cp_i) + (1-phi(i,j))*dns*Cp_s)*(T(i,j)-T_0(i,j))/dt - ...
            phi(i,j)*dni*Li*(si_temp-si_0(i,j))/dt);
    case 2
        Rt=-Rt + ...
            vb(i,j)*((phi(i,j)*(sw_temp*dnw(i,j)*Cp_w + si_temp*dni*Cp_i ) + (1-phi(i,j))*dns*Cp_s)*(T0+T(i,j))/dt - ...
            (phi_0(i,j)*(sw_0(i,j)*dnw_0(i,j)*Cp_w + si_0(i,j)*dni*Cp_i) + (1-phi_0(i,j))*dns*Cp_s)*(T0+T_0(i,j))/dt - ...
            phi(i,j)*dni*Li*(si_temp-si_0(i,j))/dt);
end


