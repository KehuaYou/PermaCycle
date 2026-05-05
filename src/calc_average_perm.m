global COL ROW
global kx 
global cnstx cnstTx 
global dx dy dz

%% Harmonic averaging of intrinsic permeability
for i=1:COL
    for j=1:ROW
        if i < COL
            cnstx(i,j)=2.0*dz(i+1,j)*dy(i+1,j)*dz(i,j)*dy(i,j)*kx(i+1,j)*kx(i,j)/ ...
                (dx(i,j)*dz(i+1,j)*dy(i+1,j)*kx(i+1,j) + dx(i+1,j)*dz(i,j)*dy(i,j)*kx(i,j)); %the average is about kx*dy*dz/dx
            cnstTx(i,j)=2.0*dz(i+1,j)*dy(i+1,j)*dz(i,j)*dy(i,j)/ ...
                (dx(i,j)*dz(i+1,j)*dy(i+1,j) + dx(i+1,j)*dz(i,j)*dy(i,j));
        else
            cnstx(i,j)=0.0;
            cnstTx(i,j)=0.0;
        end
    end
end