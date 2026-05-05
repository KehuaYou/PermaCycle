function y=calc_trans(i,j)
global COL
global cnstx 
global twx twx1 ttx 
global krw 
global vsw dnw
global lambda
global g
global pw dpth
global dx dy dz


if i < COL
    
    if pw(i+1,j)-(dnw(i+1,j)+dnw(i,j))/2*g*dpth(i+1,j) > pw(i,j)-(dnw(i+1,j)+dnw(i,j))/2*g*dpth(i,j)
        krup = krw(i+1,j);
    else
        krup = krw(i,j);
    end
    visav=(vsw(i+1,j)+vsw(i,j))/2.0;
    denav=(dnw(i+1,j)+dnw(i,j))/2.0;
    twx(i,j)=cnstx(i,j)*krup/visav;
    twx1(i,j)=twx(i,j)*denav*g;
    
    ttx(i,j) = 2 / ( dx(i,j)/(lambda(i,j)*dy(i,j)*dz(i,j)) +...
        dx(i+1,j)/(lambda(i+1,j)*dy(i+1,j)*dz(i+1,j)) );
    
end

if i > 1
    
    if pw(i-1,j) -(dnw(i-1,j)+dnw(i,j))/2*g*dpth(i-1,j) > pw(i,j)-(dnw(i-1,j)+dnw(i,j))/2*g*dpth(i,j)
        krup = krw(i-1,j);
    else
        krup = krw(i,j);
    end
    visav=(vsw(i-1,j)+vsw(i,j))/2.0;
    denav=(dnw(i-1,j)+dnw(i,j))/2.0;
    twx(i-1,j)=cnstx(i-1,j)*krup/visav;
    twx1(i-1,j)=twx(i-1,j)*denav*g;
    
    ttx(i-1,j)=2 / ( dx(i,j)/(lambda(i,j)*dy(i,j)*dz(i,j)) +...
        dx(i-1,j)/(lambda(i-1,j)*dy(i-1,j)*dz(i-1,j)) );
    
end

