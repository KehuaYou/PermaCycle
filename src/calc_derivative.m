function jac_temp = calc_derivative(i_cal,j_cal,i,j)
global pw  cl  T 
global krw krg
global eps
global INDC2 
global lambda
global pw_ini T_ini pw_scale T_scale cl_scale
global pw_dimensionless T_dimensionless cl_dimensionless 
global rsidw0 rsidc0  rsidt0
global phi 

jac_temp = zeros(3,3);
ice_condition = INDC2(i,j)>1;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%          Change pw               %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pw(i,j)=(pw_dimensionless(i,j)+eps)*pw_scale+pw_ini(i,j);
calc_trans(i,j);

rsidw1= calc_Rw(i_cal,j_cal);
rsidc1=calc_Rc(i_cal,j_cal);
rsidt1=calc_Rt(i_cal,j_cal);

jac_temp(1,1)=(rsidw1-rsidw0)/eps;
jac_temp(2,1)=(rsidc1-rsidc0)/eps;
jac_temp(3,1)=(rsidt1-rsidt0)/eps;

pw(i,j)=pw_dimensionless(i,j)*pw_scale+pw_ini(i,j);
calc_trans(i,j);

%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%          Change cl                  %%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cl(i,j)=(cl_dimensionless(i,j)+eps)*cl_scale;

if ice_condition
    temp3= lambda(i,j);
    temp4 = krw(i,j);
    temp5 = krg(i,j);
end

if ice_condition
    [sw_temp, si_temp, sg_temp] = calc_saturations(i,j);
    lambda(i,j) = bulk_thermal_conductivity(phi(i,j), sw_temp, si_temp, sg_temp);
    [krw(i,j),krg(i,j)] = relative_perm(sw_temp,si_temp,sg_temp);
end

calc_trans(i,j);

rsidw1= calc_Rw(i_cal,j_cal);
rsidc1=calc_Rc(i_cal,j_cal);
rsidt1=calc_Rt(i_cal,j_cal);

jac_temp(1,2)=(rsidw1-rsidw0)/eps;
jac_temp(2,2)=(rsidc1-rsidc0)/eps;
jac_temp(3,2)=(rsidt1-rsidt0)/eps;

cl(i,j)=cl_dimensionless(i,j)*cl_scale;
if ice_condition
    lambda(i,j)=temp3;
    krw(i,j)=temp4;
    krg(i,j)=temp5;
end
calc_trans(i,j);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%          Change T                %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T(i,j)=(T_dimensionless(i,j)+eps)*T_scale+T_ini(i,j);

if ice_condition
    temp4= lambda(i,j);
    temp5 = krw(i,j);
    temp6 = krg(i,j);
end

if ice_condition
    [sw_temp, si_temp, sg_temp] = calc_saturations(i,j);
    lambda(i,j) = bulk_thermal_conductivity(phi(i,j), sw_temp, si_temp, sg_temp);
    [krw(i,j),krg(i,j)] = relative_perm(sw_temp,si_temp,sg_temp);
end

calc_trans(i,j);

rsidw1= calc_Rw(i_cal,j_cal);
rsidc1=calc_Rc(i_cal,j_cal);
rsidt1=calc_Rt(i_cal,j_cal);

jac_temp(1,3)=(rsidw1-rsidw0)/eps;
jac_temp(2,3)=(rsidc1-rsidc0)/eps;
jac_temp(3,3)=(rsidt1-rsidt0)/eps;

T(i,j)=T_dimensionless(i,j)*T_scale+T_ini(i,j);
if ice_condition
    lambda(i,j)=temp4;
    krw(i,j)=temp5;
    krg(i,j)=temp6;
end
calc_trans(i,j);

end