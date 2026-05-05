function result=calc_sw_unfrozen(i,j)
global cl T
global Ax Swr_freeze

cl_use = cl(i,j); 
T_use = T(i,j); 
T_freeze = -cl_use*(164.49*cl_use+49.462);
result = (T_use<=T_freeze) * (exp(Ax*(T_use-T_freeze))*(1-Swr_freeze) + Swr_freeze) +...
    (T_use>T_freeze) * 1.0 ;
end