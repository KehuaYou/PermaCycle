global ROW COL
global sg si sw cm

for i=1:COL
    for j=1:ROW
        if sg(i,j)<0
            sg(i,j)=0;
        end

        if si(i,j)<0
            si(i,j)=0;
        end


        sw(i,j) = 1- sg(i,j) - si(i,j);


        if cm(i,j)<0
            cm(i,j)=0;
        end
    end
end