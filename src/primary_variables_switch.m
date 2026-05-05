global COL ROW
global INDC2 INDC3
global T si si_0 cl
global go_back

T_use = T;
cl_use = cl;

for i=1:COL
    for j=1:ROW
        %*************************** No phase change ***************
        if INDC2(i,j) == 2
            INDC2(i,j)=3;
        end
        if  INDC3(i,j)==3
            INDC3(i,j)=4;
        end

        %****************** Phase number increases****************
        if INDC2(i,j)==1 && INDC3(i,j) == 1
            if T_use(i,j)<=-cl_use(i,j)*(164.49*cl_use(i,j)+49.462)     % L changes to L+I
                INDC2(i,j)=2;
                INDC3(i,j)=3;
                % go_back=go_back+1;
                go_back=1;
                break
            end
        end

        %**************Phase number decreases*********************
        if si(i,j)<=0 && si_0(i,j)>0   % ice disappear
            INDC2(i,j)=1;
            INDC3(i,j)=1;
        end
    end

    if go_back==1
        break
    end
end