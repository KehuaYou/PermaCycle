function lambda_calc = bulk_thermal_conductivity(phi_temp, sw_temp, si_temp, sg_temp)
global lambda_s lambda_w lambda_g lambda_i
global conductivity_index

if conductivity_index == 0
    lambda_calc=(1-phi_temp)*lambda_s+phi_temp*(sw_temp*lambda_w+sg_temp*lambda_g+si_temp*lambda_i);
elseif conductivity_index == 1
    lambda_calc=1/((1-phi_temp)/lambda_s+phi_temp*sw_temp/lambda_w+phi_temp*sg_temp/lambda_g+phi_temp*si_temp/lambda_i);
else
    lambda_calc = lambda_s^(1-phi_temp) * lambda_w^(phi_temp*sw_temp) * lambda_i^(phi_temp*si_temp) * lambda_g^(phi_temp*sg_temp);
end

end