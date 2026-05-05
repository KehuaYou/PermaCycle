function Rc=calc_Rc(i,j)
global pw  dpth vb
global twx twx1 
global sw sw_0 dnw dnw_0 cl cl_0 
global phi phi_0 dt
global INDC2 
global Dsaltx 
global pw_top  cl_top dnw_top dpth_groundsurface
global dy dz 
global dx kx krw vsw g 
global  COL 

%% Initialize Residual
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rc=0;

%% Flow and Transport
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if i==1
    tmp=twx(i,j)*(pw(i+1,j)-pw(i,j)) - twx1(i,j)*(dpth(i+1,j)-dpth(i,j));
    if tmp > 0
        Rc=Rc + tmp*dnw(i+1,j)*cl(i+1,j);
    else
        Rc=Rc + tmp*dnw(i,j)*cl(i,j);
    end
    Rc=Rc + Dsaltx(i,j)*2*phi(i+1,j)*sw(i+1,j)*phi(i,j)*sw(i,j)/(phi(i+1,j)*sw(i+1,j)+phi(i,j)*sw(i,j))*(dnw(i+1,j)+dnw(i,j))/2*(cl(i+1,j)-cl(i,j));
    
elseif i==COL
    if twx(i-1,j)*(pw(i,j)-pw(i-1,j)) > twx1(i-1,j)*(dpth(i,j)-dpth(i-1,j))
        Rc=Rc - twx(i-1,j)*(pw(i,j)-pw(i-1,j))*dnw(i,j)*cl(i,j);
        Rc=Rc + twx1(i-1,j)*(dpth(i,j)-dpth(i-1,j))*dnw(i,j)*cl(i,j);
    else
        Rc=Rc - twx(i-1,j)*(pw(i,j)-pw(i-1,j))*dnw(i-1,j)*cl(i-1,j);
        Rc=Rc + twx1(i-1,j)*(dpth(i,j)-dpth(i-1,j))*dnw(i-1,j)*cl(i-1,j);
    end
    Rc=Rc - Dsaltx(i-1,j)*2*phi(i,j)*sw(i,j)*phi(i-1,j)*sw(i-1,j)/(phi(i,j)*sw(i,j)+phi(i-1,j)*sw(i-1,j))*(dnw(i,j)+dnw(i-1,j))/2*(cl(i,j)-cl(i-1,j));
    
else
    if twx(i,j)*(pw(i+1,j)-pw(i,j)) > twx1(i,j)*(dpth(i+1,j)-dpth(i,j))
        Rc=Rc + twx(i,j)*(pw(i+1,j)-pw(i,j))*dnw(i+1,j)*cl(i+1,j);
        Rc=Rc - twx1(i,j)*(dpth(i+1,j)-dpth(i,j))*dnw(i+1,j)*cl(i+1,j);
    else
        Rc=Rc + twx(i,j)*(pw(i+1,j)-pw(i,j))*dnw(i,j)*cl(i,j);
        Rc=Rc - twx1(i,j)*(dpth(i+1,j)-dpth(i,j))*dnw(i,j)*cl(i,j);
    end
    Rc=Rc + Dsaltx(i,j)*2*phi(i+1,j)*sw(i+1,j)*phi(i,j)*sw(i,j)/(phi(i+1,j)*sw(i+1,j)+phi(i,j)*sw(i,j))*(dnw(i+1,j)+dnw(i,j))/2*(cl(i+1,j)-cl(i,j));
    
    if twx(i-1,j)*(pw(i,j)-pw(i-1,j)) > twx1(i-1,j)*(dpth(i,j)-dpth(i-1,j))
        Rc=Rc - twx(i-1,j)*(pw(i,j)-pw(i-1,j))*dnw(i,j)*cl(i,j);
        Rc=Rc + twx1(i-1,j)*(dpth(i,j)-dpth(i-1,j))*dnw(i,j)*cl(i,j);
    else
        Rc=Rc - twx(i-1,j)*(pw(i,j)-pw(i-1,j))*dnw(i-1,j)*cl(i-1,j);
        Rc=Rc + twx1(i-1,j)*(dpth(i,j)-dpth(i-1,j))*dnw(i-1,j)*cl(i-1,j);
    end
    Rc=Rc - Dsaltx(i-1,j)*2*phi(i,j)*sw(i,j)*phi(i-1,j)*sw(i-1,j)/(phi(i,j)*sw(i,j)+phi(i-1,j)*sw(i-1,j))*(dnw(i,j)+dnw(i-1,j))/2*(cl(i,j)-cl(i-1,j));
    
end


%% Boundary Condition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note that if flow into the grid, use "-" sign; if flow out, use "+"

%++++++++++++++++++++ Top Boundary Condition +++++++++++++++++++++++++
%+++++++++++++++++++fixed outside water pressure++++++++++++++++++++++++++++++++
if i == 1
    twx_temp = dy(i,j)*dz(i,j)/(0.5*dx(i,j))*kx(i,j)*krw(i,j)/vsw(i,j);
    twx1_temp = twx_temp * dnw(i,j) * g;
    tmp=twx_temp*(pw(i,j)-pw_top(1,j)) - twx1_temp*(dpth(i,j)-dpth_groundsurface(1,j));
    
    if (pw(i,j) - dnw(i,j)*g*dpth(i,j)) > (pw_top(1,j) - dnw(i,j)*g*dpth_groundsurface(1,j)) % if hydraulic head at node i+1 is greater than at node i
        Rc=Rc - tmp*dnw(i,j)*cl(i,j);
    else
        Rc=Rc - tmp*dnw_top(1,j)*cl_top(1,j);
    end
    Dsalt_temp = Dsaltx(i,j) * phi(i,j) * sw(i,j) * dy(i,j)*dz(i,j)/(0.5*dx(i,j));
    Rc=Rc - Dsalt_temp*phi(i,j)*sw(i,j)*dnw(i,j)*(cl(i,j)-cl_top(1,j));
end


%% Accumulation Term
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sw_unfrozen = calc_sw_unfrozen(i,j);
sw_temp = (INDC2(i,j)==1) * 1.0 + (INDC2(i,j)>1) * sw_unfrozen;            
Rc=-Rc+vb(i,j)*(phi(i,j)*sw_temp*dnw(i,j)*cl(i,j)-phi_0(i,j)*sw_0(i,j)*dnw_0(i,j)*cl_0(i,j))/dt;




