function [sw_deriv, si_deriv, sg_deriv] = calc_saturations(i,j)
global INDC2

sw_unfrozen = calc_sw_unfrozen(i,j);

sw_deriv = (INDC2(i,j)==1) * 1.0 + (INDC2(i,j)>1) * sw_unfrozen;

si_deriv = (INDC2(i,j)==1) * 0.0 + (INDC2(i,j)>1) * (1-sw_unfrozen);

sg_deriv = 0;

end