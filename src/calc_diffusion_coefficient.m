global COL ROW
global kx krw vsw pw dnw g dpth pw_top dpth_groundsurface dx dy dz twx twx1
global cnstTx phi sw
global Dsaltx

dispersivity = 0.1;
D_molecular = 1e-9;

%% ------------- Vertical Velocity -----------------
% Positive Upward
vw_x = zeros(COL,ROW);
vw_x(1,:) = kx(1,:).*krw(1,:)./vsw(1,:).*((pw(1,:)-dnw(1,:).*g.*dpth(1,:)) -...
    (pw_top(1,:)-dnw(1,:).*g.*dpth_groundsurface))./(0.5*dx(1,:));

for i=2:COL
    vw_x(i,:) = (twx(i-1,:).*(pw(i,:)-pw(i-1,:)) -...
        twx1(i-1,:).*(dpth(i,:)-dpth(i-1,:)) )./(dy(i-1,:).*dz(i-1,:));
end

velocity = abs(vw_x);
D0 = D_molecular + dispersivity.*velocity;

%% ------------- Diffusion Operator -----------------
for i = 1:COL-1

    Dsaltx(i,1) = cnstTx(i,1)*(D0(i,1) * phi(i,1) * sw(i,1) +...
        D0(i+1,1) * phi(i+1,1) * sw(i+1,1)) / 2.0;
end

Dsaltx(COL,1) = Dsaltx(COL-1,1);

