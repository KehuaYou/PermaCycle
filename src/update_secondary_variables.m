global ROW 
global krw krg  
global lambda 
global sw sg si 
global phi 


for j=1:ROW
    for i=1:COL

        %-------- update bulk thermal conductivity ---------------
        lambda(i,j)= bulk_thermal_conductivity(phi(i,j), sw(i,j), si(i,j), sg(i,j));
        
        %-------- update relative permeability ---------------
        [krw(i,j),krg(i,j)] = relative_perm(sw(i,j),si(i,j),sg(i,j));       
        
    end
end


