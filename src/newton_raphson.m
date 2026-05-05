global COL ROW total_grid_num
global pw  si  sw  cl  T
global INDC2 
global residual
global pw_ini T_ini pw_scale T_scale cl_scale 
global rsidw0  rsidc0  rsidt0
global pw_dimensionless T_dimensionless cl_dimensionless 
global iteration back_track delta

%% Prepare Some Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%   Nondimensionalize Primary Variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pw_dimensionless=(pw-pw_ini)./pw_scale;
T_dimensionless=(T-T_ini)./T_scale;
cl_dimensionless=cl./cl_scale;

%------------------- Update Flow operators -------------------
for j=1:ROW
    for i=1:COL
        calc_trans(i,j);
    end
end

%% Residual and Jacobian Matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%   Initialize Residual and Jacobian Matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B=zeros(1,total_grid_num*3);
jac=spalloc(total_grid_num*3,total_grid_num*3,total_grid_num*3*(3*5));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%   Calculate Residual and Jacobian Matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
l=0; % Initialize position index
for j=1:ROW  
    for i=1:COL
        l=l+3;  % Position index in residual matrix
        top = l-2;
        bottom = l;

        rsidw0=calc_Rw(i,j);
        rsidc0=calc_Rc(i,j);
        rsidt0=calc_Rt(i,j);

        B(l-2)=-rsidw0;
        B(l-1)=-rsidc0;
        B(l)=-rsidt0;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%                   i-1, j (upper neighbor)                      %%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if i > 1
            left = l-5;
            right = l-3;
            jac(top:bottom, left:right) = calc_derivative(i,j,i-1,j);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%                              i, j (self)                      %%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        left = l-2;
        right = l;
        jac(top:bottom, left:right) = calc_derivative(i,j,i,j);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%           i+1, j   (lower neighbor)             %%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if i < COL
            left = l+1;
            right = l+3;
            jac(top:bottom, left:right) = calc_derivative(i,j,i+1,j);
        end

    end
end 

%% Variable Change
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   Calculate the change of primary variable %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delta = real(jac) \ real(B');
delta = real(delta);
residual = real(B');

f_damp = back_track^(iteration-1);
delta = delta * f_damp;
%% Virable Update
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   UPdate the primary variable %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
l=0;
for j=1:ROW
    for i=1:COL

        l=l+3;
        if INDC2(i,j) == 1 % L
            dpw=delta(l-2);
            dcl=delta(l-1);
            dT=delta(l);

            pw(i,j)=(pw_dimensionless(i,j)+dpw)*pw_scale+pw_ini(i,j);
            cl(i,j)=(cl_dimensionless(i,j)+dcl)*cl_scale;
            T(i,j)=(T_dimensionless(i,j)+dT)*T_scale+T_ini(i,j);

            sw(i,j)=1.0;
            si(i,j)=0;

        else    % L+I
            dpw=delta(l-2);
            dcl=delta(l-1);
            dT=delta(l);

            pw(i,j)=(pw_dimensionless(i,j)+dpw)*pw_scale+pw_ini(i,j);
            cl(i,j)=(cl_dimensionless(i,j)+dcl)*cl_scale;
            T(i,j)=(T_dimensionless(i,j)+dT)*T_scale+T_ini(i,j);

            sw(i,j) = calc_sw_unfrozen(i,j);
            si(i,j) = 1- sw(i,j);

        end

    end 
end 




