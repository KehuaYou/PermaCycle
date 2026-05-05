function Rw=calc_Rw(i,j)
global pw  dpth
global twx twx1
global sw_0 dnw dnw_0
global cl cl_0
global phi phi_0 si_0 dni
global dt
global INDC2 
global vb
global pw_top  cl_top dnw_top dpth_groundsurface
global dy dz
global dx kx krw vsw g
global COL

%% Initialize Residual
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rw=0;

%% Flow and Transport
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if i==1
    if twx(i,j)*(pw(i+1,j)-pw(i,j)) > twx1(i,j)*(dpth(i+1,j)-dpth(i,j))
        Rw=Rw + twx(i,j)*(pw(i+1,j)-pw(i,j))*dnw(i+1,j)*(1-cl(i+1,j));
        Rw=Rw - twx1(i,j)*(dpth(i+1,j)-dpth(i,j))*dnw(i+1,j)*(1-cl(i+1,j));
    else
        Rw=Rw + twx(i,j)*(pw(i+1,j)-pw(i,j))*dnw(i,j)*(1-cl(i,j));
        Rw=Rw - twx1(i,j)*(dpth(i+1,j)-dpth(i,j))*dnw(i,j)*(1-cl(i,j));
    end

elseif i==COL
    if twx(i-1,j)*(pw(i,j)-pw(i-1,j)) > twx1(i-1,j)*(dpth(i,j)-dpth(i-1,j))
        Rw=Rw - twx(i-1,j)*(pw(i,j)-pw(i-1,j))*dnw(i,j)*(1-cl(i,j));
        Rw=Rw + twx1(i-1,j)*(dpth(i,j)-dpth(i-1,j))*dnw(i,j)*(1-cl(i,j));
    else
        Rw=Rw - twx(i-1,j)*(pw(i,j)-pw(i-1,j))*dnw(i-1,j)*(1-cl(i-1,j));
        Rw=Rw + twx1(i-1,j)*(dpth(i,j)-dpth(i-1,j))*dnw(i-1,j)*(1-cl(i-1,j));
    end

else
    if twx(i,j)*(pw(i+1,j)-pw(i,j)) > twx1(i,j)*(dpth(i+1,j)-dpth(i,j))
        Rw=Rw + twx(i,j)*(pw(i+1,j)-pw(i,j))*dnw(i+1,j)*(1-cl(i+1,j));
        Rw=Rw - twx1(i,j)*(dpth(i+1,j)-dpth(i,j))*dnw(i+1,j)*(1-cl(i+1,j));
    else
        Rw=Rw + twx(i,j)*(pw(i+1,j)-pw(i,j))*dnw(i,j)*(1-cl(i,j));
        Rw=Rw - twx1(i,j)*(dpth(i+1,j)-dpth(i,j))*dnw(i,j)*(1-cl(i,j));
    end
    if twx(i-1,j)*(pw(i,j)-pw(i-1,j)) > twx1(i-1,j)*(dpth(i,j)-dpth(i-1,j))
        Rw=Rw - twx(i-1,j)*(pw(i,j)-pw(i-1,j))*dnw(i,j)*(1-cl(i,j));
        Rw=Rw + twx1(i-1,j)*(dpth(i,j)-dpth(i-1,j))*dnw(i,j)*(1-cl(i,j));
    else
        Rw=Rw - twx(i-1,j)*(pw(i,j)-pw(i-1,j))*dnw(i-1,j)*(1-cl(i-1,j));
        Rw=Rw + twx1(i-1,j)*(dpth(i,j)-dpth(i-1,j))*dnw(i-1,j)*(1-cl(i-1,j));
    end

end


%% Boundary Condition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note that if flow into the grid, use "+" sign; if flow out, use "-"

%++++++++++++++++++++ Top Boundary Condition +++++++++++++++++++++++++
%++++++++++++++++++Fixed outside water pressure++++++++++++++++++++++++
if i == 1
    twx_temp = dy(i,j)*dz(i,j)/(0.5*dx(i,j))*kx(i,j)*krw(i,j)/vsw(i,j);
    twx1_temp = twx_temp * dnw(i,j) * g;
    if (pw(i,j) - dnw(i,j)*g*dpth(i,j)) > (pw_top(1,j) - dnw_top(1,j)*g*dpth_groundsurface(1,j)) % if hydraulic head at node i+1 is greater than at node i
        Rw=Rw - twx_temp*(pw(i,j)-pw_top(1,j))*dnw(i,j)*(1-cl(i,j));
        Rw=Rw + twx1_temp*(dpth(i,j)-dpth_groundsurface(1,j))*dnw(i,j)*(1-cl(i,j));
    else
        Rw=Rw - twx_temp*(pw(i,j)-pw_top(1,j))*dnw_top(1,j)*(1-cl_top(1,j));
        Rw=Rw + twx1_temp*(dpth(i,j)-dpth_groundsurface(1,j))*dnw_top(1,j)*(1-cl_top(1,j));
    end
end

%% Accumulation Term
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sw_unfrozen = calc_sw_unfrozen(i,j);
sw_temp = (INDC2(i,j)==1) * 1.0 + (INDC2(i,j)>1) * sw_unfrozen; 
si_temp = (INDC2(i,j)==1) * 0.0 + (INDC2(i,j)>1) * (1-sw_unfrozen);

Rw=-Rw+vb(i,j)*( (phi(i,j)*sw_temp*dnw(i,j)*(1-cl(i,j))+phi(i,j)*si_temp*dni)-...
    (phi_0(i,j)*sw_0(i,j)*dnw_0(i,j)*(1-cl_0(i,j))+phi_0(i,j)*si_0(i,j)*dni) )/dt;
