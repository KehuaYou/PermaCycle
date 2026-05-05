function [krw_calc, krg_calc] = relative_perm(sw_temp,si_temp,sg_temp)
global wn gn swr sgr

%======= reduction due to the presence of ice =======
if si_temp>0
    kh=(1-si_temp)^2;
else
    kh=1;
end

%======= reduction due to the mutliphase flow =======
sw_eff = sw_temp/(1-si_temp);
sg_eff = sg_temp/(1-si_temp);

sw_temp_new = (sw_eff-swr)/(1-swr);
sg_temp_new = (sg_eff-sgr)/(1-swr);

if sw_eff <= swr
    krw_calc = 0;
else
    krw_calc = kh*sw_temp_new^wn;
end

if sg_eff >= sgr
    krg_calc = kh*sg_temp_new^gn;
else
    krg_calc =0;
end
end