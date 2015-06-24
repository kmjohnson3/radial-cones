function [kA,time,gA,pre_area] = make_base_cone_flex(opxres,fov,maxbw,gmax,smax,readout_time,cone_angle,T,frac)

%%-------------------------------------------------------------
% Code to make a cone for Radial Cones Trajectory
%
% Inputs:
%   opxres = resolution
%   fov    = field of view in mm
%   maxbw  = max readout bandwidth (+/-)
%   smax   = max slew
%   readout_time = time for readout in ms
%   cone_angle = angle at edge of k-space in degrees
%   T = update time in s
%
%Outputs:
%   k = 3d k space trajectory
%   time = time points in s
%   g = required gradients

res = fov/opxres; %mm
kmax = 5/res;
conversion = fov/res/kmax/2;
gamma = 4258; %Hz/G
bwmax = 2e3*maxbw; %Bandwidth in Hz
gmax = min( gmax, (bwmax)/( gamma *fov/10)); %Max G/cm


disp('Cone Design Inputs: ');
disp(['  Xres = ',num2str(opxres)]);
disp(['  Fov  = ',num2str(fov)]);
disp(['  Maxbw= ',num2str(maxbw)]);
disp(['  Smax = ',num2str(smax)]);
disp(['  Readout Time = ',num2str(readout_time)]);
disp(['  Cone Angle = ',num2str(cone_angle)]);
disp(['  Kmax = ',num2str(kmax)]);

%% Make a Cone
costX = 100;
p = 0.0;
%while costX > 0.001

for pass = 1
   % ------------Calculate K-Space Trajectory------------------------------
    if pass ==1
        p = 0.1:0.3:5;
    else
        p = (-0.5:0.05:0.5)*(1/5^(pass-2)) + p_ideal;
    end
    
    t = linspace(1-2*frac,1,1000);
    for p_pos = 1:numel(p)
        
        fact = 0.0;
        for step = [1 0.5 0.1 0.05 0.01 0.005]
            rtime =0;
            while rtime < readout_time
                
                kmax_rad = kmax*cos(cone_angle/180*pi);
                fact = fact+step;
                k(:,3) = t*kmax_rad;
%                kTr = ( 2*(t-0.15).^2 )*kmax*sin(cone_angle/180*pi);
                kTr = ( abs(t).^2 )*kmax*sin(cone_angle/180*pi);

%                 kTr( (t-0.15) <0) =0;
%                 phi = (t-0.15).^p(p_pos);
                phi = sign(t).*abs(t).^p(p_pos);

%                 phi( (t-0.15) <0) =0;
                k(:,1) = cos( fact*phi*2*pi ).*kTr;
                k(:,2) = sin( fact*phi*2*pi ).*kTr;
                
                %plot3(k(:,1),k(:,2),k(:,3))
                %daspect([1 1 1])
                %drawnow
                
                if frac > 0.5
                    [C,time,gA,s,kA, phi, sta, stb] = minTimeGradient(k,[],[],gmax,smax,16e-3);;
                else
                    [C,time,gA,s,kA, phi, sta, stb] = minTimeGradient(k,0.0,[],gmax,smax,16e-3);;
                end
                rtime = max(time);
            end
            fact = fact - step;
        end
        a_l = sqrt(sum(diff(kA).^2,2));
        k_l = sum( a_l) / ( 2*frac*kmax);
               
        %------------Calculate Cost--------------------------------------------
        k_r = sqrt( kA(:,1).^2 + kA(:,2).^2 + kA(:,3).^2 );
        rr = linspace(min(k_r),max(k_r),1000)';
        
        %plot(k_r,kA(:,1),k_r,kA(:,2),k_r,kA(:,3));
        %pause
        
        kR(:,1)=interp1(k_r,kA(:,1),rr,'linear');
        kR(:,2)=interp1(k_r,kA(:,2),rr,'linear');
        kR(:,3)=interp1(k_r,kA(:,3),rr,'linear');
        arc_length = sqrt(sum(diff(kR).^2,2));
        
        
%         clf
%         subplot(211),plot(kR(:,1),'k'); hold on; plot(kR(:,2),'r'); plot(kR(:,3),'b');
%         subplot(212),plot(kA(:,1),'k'); hold on; plot(kA(:,2),'r'); plot(kA(:,3),'b');
%         drawnow
        
        A = [rr(1:end-1).^2];
        b = linsolve(A,arc_length);
        
        f = A*b;
        cost(p_pos) = sum( abs(f - arc_length).^2);
        disp(['E=',num2str(p(p_pos)),' Length Factor = ',num2str(k_l),' E = ',num2str(cost(p_pos))]);
        drawnow
    end
    [c,idx] = min(cost);
    p_ideal = p(idx);
    
    
 %   plot(p,cost);
    p = polyfit(p,cost,2);
    p_ideal = -p(2)/(2*p(1));
    
    
    %------------Recalculate with ideal------------------------------------
    fact = 0.0;
    for step = [1 0.5 0.1 0.05 0.01 0.005]
        rtime =0;
        while rtime < readout_time
                kmax_rad = kmax*cos(cone_angle/180*pi);
                fact = fact+step;
                kIn(:,3) = t*kmax_rad;
%                kTr = ( 2*(t-0.15).^2 )*kmax*sin(cone_angle/180*pi);
                kTr = ( abs(t).^2 )*kmax*sin(cone_angle/180*pi);

%                 kTr( (t-0.15) <0) =0;
%                 phi = (t-0.15).^p(p_pos);
                phi = sign(t).*abs(t).^p_ideal;

%                 phi( (t-0.15) <0) =0;
                kIn(:,1) = cos( fact*phi*2*pi ).*kTr;
                kIn(:,2) = sin( fact*phi*2*pi ).*kTr;
                if frac > 0.5
                    [C,time,gA,s,kA, phi, sta, stb] = minTimeGradient(kIn,[],[],gmax,smax,T);
                else
                    [C,time,gA,s,kA, phi, sta, stb] = minTimeGradient(kIn,0,[],gmax,smax,T);
                end
                rtime = max(time);
                disp(rtime)
        end
        fact = fact - step;
    end
    a_l = sqrt(sum(diff(kA).^2,2));
    k_l = sum( a_l) / kmax;
    
    
%     figure,plot3(kIn(:,1),kIn(:,2),kIn(:,3))
%     figure,plot3(kA(:,1),kA(:,2),kA(:,3))
%     pause
    
    %------------Calculate Cost--------------------------------------------
    k_r = sqrt( kA(:,1).^2 + kA(:,2).^2 + kA(:,3).^2 );
    rr = linspace(min(k_r),max(k_r),1000)';
    
    kR(:,1)=interp1(k_r,kA(:,1),rr);
    kR(:,2)=interp1(k_r,kA(:,2),rr);
    kR(:,3)=interp1(k_r,kA(:,3),rr);
    arc_length = sqrt(sum(diff(kR).^2,2));
    
    A = [rr(1:end-1).^2];
    b = linsolve(A,arc_length);
    
    f = A*b;
 %   cost(p_pos) = sum( abs(f - arc_length).^2);
    
%     %------------Plot          --------------------------------------------
%     plot(rr(1:end-1),arc_length,rr(1:end-1),f);
%     hold on
%     title(num2str(p))
%     drawnow
%     
end

%Need to figure out areas required to prewind and rewind
if frac == 0.5
    pre_area = [0 0 0];
elseif frac == 1.0
    [m,idx] = max(gA(:,3));
    pre_area = T*sum(gA(1:idx,:),1);
end
    





return
