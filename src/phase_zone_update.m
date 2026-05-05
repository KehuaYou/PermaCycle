global COL ROW
global cl T pw
global INDC3 INDC1
global p_salinity T_salinity salinity

for j=1:ROW
    for i=1:COL

        if T(i,j) > -cl(i,j)*(164.49*cl(i,j)+49.462)
            INDC3(i,j)=1;
        else
            INDC3(i,j)=4;
        end

        if cl(i,j) <interpolation2(p_salinity,T_salinity,salinity,pw(i,j)/1e6,T(i,j))
            INDC1(i,j)=1;
        else
            INDC1(i,j)=0;
        end

    end
end