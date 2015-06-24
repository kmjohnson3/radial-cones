clear 
clc

%Design Parameters
mat = 256; %matrix size
fov = 220; %mm
res = fov/mat; %mm
ileaves = 10000; % Number of interleaves
bw_readout = 125e3*2; % +/- 125 kHz
cone_angle = 30;
readout_time = 2.1;

% System
gmax = 4.9; % G/cm
smax = 15; %G/cm/us
T = 4e-3; % gradient sample time

% For export
BaseName = 'RadialCones';
dirname = pwd;

% Main function call
tic 
[g,gR,k,kR,R,xdir,ydir,zdir] = generate_unoptimized_radial_cones( bw_readout , gmax, smax, ileaves,mat,fov,T, readout_time, cone_angle);
toc
% Get aproximate density compensation for optimization
kr_line = sqrt(sum(k(:,3).^2,2));
gr_line = sqrt(sum(g(:,3).^2,2));
kw_line = gr_line.*( kr_line.^2);

%-------------------------------
%  Export Gradients and Rotation
%--------------------------------

% These would export gradients to a file for the scanner
%export_gradient( g(:,1),T*1e3,0.5,fullfile(dirname,[BaseName,'.gradx']),fov,mat,bwmax*1e-3,smax);
%export_gradient( g(:,2),T*1e3,0.5,fullfile(dirname,[BaseName,'.grady']),fov,mat,bwmax*1e-3,smax);
%export_gradient( g(:,3),T*1e3,0.5,fullfile(dirname,[BaseName,'.gradz']),fov,mat,bwmax*1e-3,smax);
%
%fid = fopen(fullfile(dirname,[BaseName,'.rotation']),'wb');
%for pos = 1:ileaves
%    fwrite(fid,R(:,:,pos),'float');
%end
%fclose(fid);


%-------------------------------
%  Export K-Space for optimization code
%--------------------------------

fid = fopen(fullfile(dirname,[BaseName,'.struct']),'w');
fwrite(fid,size(k,1),'int');
fwrite(fid,k(:,1),'float');
fwrite(fid,k(:,2),'float');
fwrite(fid,k(:,3),'float');
fwrite(fid,kw_line,'float');
fwrite(fid,ileaves,'int');
fwrite(fid,xdir(1:ileaves),'float');
fwrite(fid,ydir(1:ileaves),'float');
fwrite(fid,zdir(1:ileaves),'float');
fclose(fid);