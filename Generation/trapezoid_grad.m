function [pwa,pw,pwd,amp,g] = trapezoid_grad(area,smax,gmax)

disp('Trapezoid calculation')
disp(['Area = ',num2str(area)]);
disp(['Smax = ',num2str(smax)]);
disp(['Gmax = ',num2str(gmax)]);

risetime = rup_grd(gmax/smax);
target_area = abs(area); % + area_ramp_pre  + area_ramp_post;

actual_area = 0;
gtrap = 0;
pwa = 0;
pw  = 4;
pwd = 0;


% - Calculate Area
up    = ( 0:(risetime/4))*4*smax;
flat  = gmax;
down  = fliplr(up);
ramp_grad = [up(:); flat(:); down(:)];
ramp_area = sum(4*ramp_grad);


if ramp_area > target_area
    %Triangle
    pwa = rup_grd( sqrt(target_area/smax));
    pw  = 4;
    pwd = pwa;
    amp = smax*pwa;
else
    pwa = risetime;
    pw  = rup_grd( (target_area-ramp_area)/gmax);
    pwd = risetime;
    amp = gmax;
end

%Trapezoid
up    = ( 0:(pwa/4))*4*smax;
flat  = amp*ones(1,pw/4);
down  = fliplr(up);

g = [up(:); flat(:); down(:)];
g = g * area / sum(4*g);
    
amp = max(g);

function w =  rup_grd( w )

w = ceil( w / 4)*4;





