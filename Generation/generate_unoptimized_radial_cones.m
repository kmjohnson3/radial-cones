function [g_readout,g_rewind,k_readout,k_rewind,R,xdir,ydir,zdir] = generate_unoptimized_radial_cones( bw_readout , gmax, smax, ileaves,mat,fov,T, readout_time, cone_angle)

% Inputs (all scalars):
%   bw_readout = readout in kHz
%   gmax = maximum gradient stength in G/cm (5)
%   smax = maximum slew in G/cm/us
%   mat  = matrix size in pixels
%   fov  = fov in mm
%   ileaves = number of interleaves
%   T = update time of gradients in ms
%   readout_time = Length of readout in ms

% Outputs:
%   g_readout - [N x 3 ] Gradients of a base cone
%   k_readout - [N x 3 ] Kspace coordinates of base cone
%   g_rewind  - [N x 3 ] Gradients of a base cone + rewinder
%   k_rewind  - [N x 3 ] Kspace coordinates of base cone + rewinder
%   R - [3 x 3 x ileaves ] - Rotation matrix


% Precalc some constants
res = fov/mat; %mm
kmax = 5/res;
conv = fov/res/kmax/2;
bw = bw_readout;
frac = 0.5;

%---------------------------------------------
%   Generate Base K-Space Cone
%---------------------------------------------
if cone_angle==0
    % This is 3D radial
    error('This is Radial');
else
    [k,time,g,pre_area] = make_base_cone_flex(mat,fov,125  ,gmax,smax,readout_time,cone_angle,T,frac);
end
size(g)
g = [zeros(3,1) g']'; % Add aditional dwell point
g_readout = g;

%---------------------------------------------
%   Generate Rewiners 
%     This is not the fasteest way to rewind gradients since it limits slew
%     to gmax/sqrt(3) on each axis 
%---------------------------------------------

gmax3 = gmax/sqrt(3);
smax3 = smax/sqrt(3);
max_pts = 0;
readout_pts = size(k,1);
max_rewindtime = 0;
done = [0 0 0];
for pass = 1:2
    for dir =1:3
        
        if done(dir)==1
            continue;
        end
        
        gX = g(:,dir);
        gx_ramp_width= abs(gX(end)) / smax3;
        gx_ramp_area = 0.5*gx_ramp_width*gX(end);
        areax_end = -sum(gX*T)
        
        disp(['Area end = ',num2str(areax_end),' gx_ramp_end=',num2str(gx_ramp_area)]);
        
        if  sign( areax_end - gx_ramp_area) == sign(gx_ramp_area)
            if pass==2
                disp(['Dir=',num2str(dir),'Gradient end = ',num2str(gX(end)),' Area end = ',num2str(areax_end)]);
                disp('Gradient should not be fully rewound - calculating bridged gradient');
                
                w = 2*abs(areax_end/gX(end));
                if w > max_rewindtime
                    
                    % Calculate trap + square
                    w1 = 999999;
                    w2 = 999999;
                    w3 = 999999;
                    gend  = gX(end);
                    gmaxR =  gend + smax3*T*sign(gend);
                    
                    while ( (w2+w1+w3) > max_rewindtime ) && ( abs(gmaxR) < gmax3);
                        disp(['W2 is ',num2str(w2+w),' of ',num2str(max_rewindtime)]);
                        gmaxR = gmaxR + smax3*T*sign(gend);
                        
                        %Ramp down
                        w1  = abs(gmaxR / smax3);
                        pts = ceil( w1 / T);
                        gX_rampdn = linspace(gmaxR,0,pts+1);
                        
                        %Ramp up
                        w2  = abs(gmaxR-gX(end)) / smax3;
                        pts = ceil( w2 / T);
                        gX_rampup = linspace(gX(end),gmaxR,pts+1);
                        %gX_rampup = gX_rampup(1:end);
                        
                        %Flat top
                        area_target = T*( -sum(gX_rampup)-sum(gX)-sum(gX_rampdn));
                        w3 = area_target / gmaxR;
                        pts = ceil(w3/T);
                        gx_rewind_trap = ones(1,pts)*gmaxR;
                    end
                    
                    g_rewind = [gX_rampup'; gx_rewind_trap'; gX_rampdn'];
                    area_adjust = areax_end - sum(T*g_rewind);
                    %Scale just the amount stronger than gX(end)
                    ramp_down = linspace(g_rewind(1),0,numel(g_rewind))';
                    gdiff = g_rewind - ramp_down;
                    gdiff = gdiff * ( areax_end - sum(ramp_down*T))/sum(T*gdiff);
                    g_rewind = gdiff + ramp_down;
                    
                    %calculate constant an
                    gX = [gX; g_rewind];
                    
                    
                else
                    disp('Ramp down only')
                    % Gradient is a ramp down only
                    pts = 1+ceil( w / T);
                    gX_rewind = linspace(1,0,pts);
                    a_rewind  = areax_end / (T*sum(gX_rewind));
                    gX = [gX; a_rewind*gX_rewind'];
                end
                
                g2{dir} = gX;
                max_pts = max(numel(gX),max_pts);
                
            end
            
            disp(['Sum G = ',num2str(sum(T*gX))])
        elseif pass==1
            disp(['Dir=',num2str(dir),'Gradient end = ',num2str(g(end,dir)),' Area end = ',num2str(areax_end)]);
            disp(['Gradient Rewound + Trapezoid']);
            % Gradient maust be rewound first
            time_rewind = abs(g(end,dir)) / smax3;
            pts = ceil( time_rewind/ T);
            gX_rewind = linspace(g(end,dir),0,pts+1);
            gX_rewind = gX_rewind(2:end);
            gX = [gX; gX_rewind'];
            
            % Area
            area_trap =   -sum(gX*T);
            time_trap = sqrt( abs(area_trap) / smax3);
            pts = ceil( time_trap / T);
            gX_trap = [ linspace(0,1,pts) linspace(1,0,pts)];
            a_trap =  area_trap / (T*sum(gX_trap));
            if( abs( a_trap) > gmax3)
                pts_up = ceil(gmax3/smax3/T);
                pwa = pts_up*T;
                pw = abs(area_trap) - pwa*gmax3;
                pts = ceil(pw/T);
                gX_trap = [ linspace(0,1,pts_up) ones(1,pts) linspace(1,0,pts_up)];
                a_trap =  area_trap / (T*sum(gX_trap));
            end
            
            gX = [gX; a_trap*gX_trap'];
            
            max_rewindtime = max( 2*time_trap+time_rewind, max_rewindtime);
            
            g2{dir} = gX;
            max_pts = max(numel(gX),max_pts);
            done(dir)=1;
            
            disp(['Sum G = ',num2str(sum(T*gX))])
        end
        
    end
end

% Append the Rewinder
g3 = zeros(max_pts,3);
g3(1:numel(g2{1}),1) = g2{1};
g3(1:numel(g2{2}),2) = g2{2};
g3(1:numel(g2{3}),3) = g2{3};
g3 = [g3' zeros(3,1) ]'; % Add aditional dwell point
g = g3;
g_rewind = g; 

%Generate K-space from gradient
k = 4258*cumtrapz(g)*T/1000;

k_l = sum( sqrt(sum(diff(k(1:readout_pts,:)).^2,2))) / kmax;
disp(['Read Time = ',num2str(T*numel(g(:,1)))]);
disp(['Length Factor = ',num2str(k_l)]);
time = ( 1:length(g(:,1)) )* T;

% Generate K-Space / Gradient Plots
figure
subplot(221), plot(time,sqrt(g(:,1).^2 + g(:,2).^2 + g(:,3).^2)); title('gradient magnitude');
subplot(222), plot(time(1:end-1),1/T*sqrt(diff(g(:,1)).^2 + diff(g(:,2)).^2 + diff(g(:,3)).^2)); title('gradient slew');
subplot(223), plot(time,g(:,1),time,g(:,2),time,g(:,3)); title('Raw Gradients');
subplot(224), plot3(k(:,1),k(:,2),k(:,3),'-k'); title('Trajectory'); daspect([1 1 1]);

%Plot the Gradients
figure
plot(time,g(:,3),time,g(:,1),time,g(:,2),'LineWidth',2.5);
legend('G_r','G_t_x','G_t_y','Orientation','Horizontal');
legend('boxoff');
xlabel('Time (ms)','FontSize',16);
ylabel('Gradient (G/cm)','FontSize',16);
ylim([-gmax gmax]);
set(gca,'FontSize',16,'LineWidth',2.5);

% Generate and Plot First Moment
figure
m1 = [];
m1(:,1) = cumsum(time'.*g(:,1));
m1(:,2) = cumsum(time'.*g(:,2));
m1(:,3) = cumsum(time'.*g(:,3));
subplot(221), plot(time,g(:,1),time,g(:,2),time,g(:,3)); title('Raw Gradients');
subplot(222), plot(time,cumsum(g(:,1)),time,cumsum(g(:,2)),time,cumsum(g(:,3))); title('M0');
subplot(223), plot(time,m1(:,1),time,m1(:,2),time,m1(:,3)); title('M1');

%Recalc K-Space
k2 = [];
k2(:,1) = cumsum(4.258*T*g(:,1));
k2(:,2) = cumsum(4.258*T*g(:,2));
k2(:,3) = cumsum(4.258*T*g(:,3));
k = k2;
k_readout = k(1:size(g_readout,1),:);
k_rewind = k;

%---------------------------------------------
%   Put Cones on Radial Lines and rotate
%---------------------------------------------
pr = 0:(ileaves-1);
zdir = (2*pr-ileaves+1)/ileaves;
xdir = cos(  sqrt(ileaves*pi)*asin(zdir) ).*sqrt( 1-(zdir.^2));
ydir = sin(  sqrt(ileaves*pi)*asin(zdir) ).*sqrt( 1-(zdir.^2));

% This is very slow on some machines excluded from release but
% available:
%     http://www.mathworks.com/matlabcentral/fileexchange/37004-suite-of-functions-to-perform-uniform-sampling-of-a-sphere
% It does not have a significant effect on image quality for
% reasonable numbers of ileaves
%     [V,Tri,~,Ue]=ParticleSampleSphere('N',ileaves);
%     xdir = V(:,1);
%     ydir = V(:,2);
%     zdir = V(:,3);

kx = single(zeros(length(k(:,1)),ileaves));
ky = single(zeros(length(k(:,1)),ileaves));
kz = single(zeros(length(k(:,1)),ileaves));

omega_all = 2*pi*rand(ileaves,1);
k = k*fov/10;

for pos = 1:ileaves
    
    phi   = acos( zdir(pos)/ ( sqrt( xdir(pos)^2 + ydir(pos)^2 + zdir(pos)^2)));
    theta = pi/2+ atan2( ydir(pos),xdir(pos));
    
    %Unit vectors which describe the cone axis r'
    ux = xdir(pos);
    uy = ydir(pos);
    uz = zdir(pos);
    
    % Create a rotation matrix from the angles
    %about X
    Rphi = [ 1 0 0;
        0 cos(phi)  -sin(phi);
        0 sin(phi)  cos(phi)];
    Rtheta = [ cos(theta) -sin(theta) 0;
        sin(theta)   cos(theta) 0;
        0 0 1];
    
    omega = omega_all(pos);
    
    utu = [ ux^2  ux*uy ux*uz;
        ux*uy uy^2  uy*uz;
        ux*uz uy*uz  uz^2];
    usk = [ 0 -uz uy;
        uz  0 -ux;
        -uy ux 0];
    Raxis = utu + cos(omega)*(eye(3)-utu) + sin(omega)*usk;
    Rnet = Raxis*Rtheta*Rphi;
    
    R(:,:,pos) = Rnet;
    
end



