%% Here to calculate the 3D directional variance of a 3D image
% This is the main program
% Here to load images and create a 3D stack
clear;
sz = 240;
sz2 = 240;
sz3 = 100; % modify 1, define sz, sz2, and sz3 according to the size of the 3D image

shgstack = zeros(sz,sz2,sz3);

for i = 1:sz3
    shgstack(:,:,i) = imread(['C:\Users\97916\Desktop\3D code\Supplementary Software\Example\',num2str(i),'.tif']); % modify 2 
    % 'shgstack' is the 3D image used for directional variance analysis.
    % Define the path of the folder where the images are saved
end

%% Here to create the mask selecting vessel regions

threshint = 5;  % modify 3 
% 'threshint' is the raw background intensity acquired by averaging intensity of several regions identified as background
maska3d = shgstack > threshint; % maska3d is a raw binary mask based on the raw background intensity

meanint = mean(shgstack(maska3d)); % calculate the mean intensity within the raw binary mask
threshinta = 0.45*meanint; % modify 4 
% 'threshinta' is an inmproved intensity threshold, and the factor '0.45' can be adjusted according to different sample
maska3da = shgstack > threshinta; % 'maska3da' is an improved binary mask

% Here to create the finalmask
finalmaska = maska3d.*maska3da; % 'finalmaska' is the binary mask acquired based on the above two masks

finalmask = zeros(sz,sz2,sz3);
for j = 1:sz3
    Maskf1 = finalmaska(:,:,j); 
    Maskf2 = miprmdebrisn(Maskf1,15); % 'miprmdebrisn.m' is a function which is able to remove very tiny structures
    finalmask(:,:,j) = Maskf2;
end

finalmask = logical(finalmask); % 'finalmask' is the final binary mask selecting the vessel-only regions

%% Here to calculate the voxel-wise 3D orientation
armdx = 5; % modify 5 
% the width of window in both 'x' and 'y' dimensions is both '2*armdx+1'.
% 'armdx' is chosen based on the diamater of vessels. Typically '2*armdx+1'
% is 2 to 3 times the diameter of vessel so as to provide optimal accuracy,
% as published in Biomed. Opt. Express 2015 (6), 2294-2310
armdz = 2; % modify 6
% '2*armdz+1' is the width of the window in 'z' dimension. We would expect
% the actual width of the window in 'x', 'y' dimensions is equal to the one
% in 'z' dimension. Therefore, 'armdz' is determined by the sampling 
% frequency in the 'xy' dimension and 'z' dimension  
filton = 1; % modify 7
% 'filton' is a binary variable that tells program whether to filter the
% data prior to analysis. 1: filtering; 0: no filtering
para = 0.58; % modify 8
% 'para' is the ratio of sampling frequency between 'xy' dimension and 'z'
% dimension
[aSDE,pSDER,bSDER,gSDER] = calcfibang3D(shgstack,armdx,armdz,filton,para); 
% 'calcfibang3D.m' is the function which calculates the voxel-wise 3D
% orientation of a 3D image. 'aSDE' is the calculated theta stack, 'pSDER'
% is the calculated phi stack, 'bSDER' is the calculated beta stack, and
% 'gSDER' is the calculated gamma stack

%% Here to calculate the voxel-wise 3D directional variance based on a localized window
bwhole = sqrt(1./(tan(2*bSDER*pi/180).^2)+1./(tan(2*gSDER*pi/180).^2));
Cwhole = (bwhole./sqrt(1+bwhole.^2)).*cos(2*aSDE*pi/180);
Swhole = (bwhole./sqrt(1+bwhole.^2)).*sin(2*aSDE*pi/180);
Zwhole = zeros(sz,sz2,sz3);
% 'bwhole', 'Cwhole', 'Swhole' and 'Zwhole' are all assisting variables
% used to acquire the directional variance results. The algorithm is 
% detailed in the 'Results' section of the paper 
for m = 1:sz
   for n = 1:sz2
     for o = 1:sz3
        if bSDER(m,n,o) <= 90   
            Zwhole(m,n,o) = 1/sqrt(1+bwhole(m,n,o)^2);
        else
            Zwhole(m,n,o) = -1/sqrt(1+bwhole(m,n,o)^2);
        end
     end
   end
end

Cwholemean = zeros(sz,sz2,sz3);
Swholemean = zeros(sz,sz2,sz3);
Zwholemean = zeros(sz,sz2,sz3);
ll = 0;

armdxx = 5; % modify 9
armdzz = 2; % modify 10
% 'armdxx' and 'armdzz' are defining the size of the window used to
% acquired the 3D directional variance (xy: 2*armdxx+1; z: 2*armdzz+1).
% Typically, these two parameters are defined the same as 'armdx' and
% 'armdz', so that the window used to calculate the voxel-wise orientation 
% and the one used to acquire localized directional variance have the same 
% size. However, they can be set to other values according to the purpuse 
% of analysis 

h = waitbar(0,'Please Wait');
ijk = 0;
nijk = (2*armdxx+1)*(2*armdxx+1)*(2*armdzz+1);
tic

for  i = -armdxx:armdxx
    for j = -armdxx:armdxx
        for k = -armdzz:armdzz
            ijk = ijk+1;
            waitbar(ijk/nijk)
            Cim = circshift(Cwhole,[i j k]);
            Sim = circshift(Swhole,[i j k]);
            Zim = circshift(Zwhole,[i j k]);
                        
            Cwholemean = Cwholemean+Cim;
            Swholemean = Swholemean+Sim;
            Zwholemean = Zwholemean+Zim;
            ll = ll+1;
        end
    end
end
    
toc
close(h)
Cwholemean = Cwholemean/ll;
Swholemean = Swholemean/ll;
Zwholemean = Zwholemean/ll;
Rwholemean = sqrt(Cwholemean.^2+Swholemean.^2+Zwholemean.^2);
Vmatr = 1-Rwholemean;
% 'Vmatr' is the voxel-wise 3D directional variance stack, with every voxel
% filled with a directional variance value showing the vessel organization
% within the localized window region surrounding it
Vmatr(isnan(Vmatr)) = 0;
Vvaluelocal = mean(Vmatr(finalmask));
% 'Vvaluelocal' is a mean value of the directional variance of this 3D 
% stack, with directional variance acquired by the localized window. Notice
% that only the regions identified as vessels are considered and contribute
% to the mean value
%% Here to calculate the 3D directional variance based on the entire 3D image
Cvalue = mean(Cwhole(finalmask));
Svalue = mean(Swhole(finalmask));
Zvalue = mean(Zwhole(finalmask));
Rvalue = sqrt(Cvalue^2+Svalue^2+Zvalue^2);
Vvalueentire = 1-Rvalue;
% 'Vvalueentire' is the 3D directional variance value acquired from the
% orientation information of all the vessels within the entire 3D stack
%% Here to do post-processing

% For post-processing, we prepare the 'pretty' images of orientation and
% directional variance. In these 'pretty' images, the raw intensity image 
% is used to provide the contrast of vessel features, and the orientation or
% directional maps are labeled by different colors to show the orientation 
% or directional variance information

% first, prepare the pretty theta orientation stack
shgstackre = shgstack;
shgstackre = shgstackre/max(max(max(shgstackre)));
prettytheta = zeros(sz,sz2,3,sz3);

uplim = 180;
botlim = 0;
bright = 0.99;
dark = 0.01;

for mm = 1:sz3
    shgima = shgstackre(:,:,mm);
    thetaima = aSDE(:,:,mm);
    thetaprettyima = prettyImage(thetaima,shgima,'none',jet(64),uplim,botlim,bright,dark);
    % 'prettyImage.m' is the function used to create 'pretty' images.
    % 'thetaima' is the theta orientation map, 'shgima' is the raw SHG
    % image. 'jet(64)' designates the color scheme. 'uplim' and 'botlim' 
    % are the upper and bottom limits of the orientation index. The range 
    % of both theta and phi is from 0 to 180. 'bright' and 'dark' are used 
    % to enhance the contrast of the image
    prettytheta(:,:,:,mm) = thetaprettyima;
end

% second, prepare the pretty phi orientation stack
prettyphi = zeros(sz,sz2,3,sz3);
for mm = 1:sz3
    shgima = shgstackre(:,:,mm);
    phiima = pSDER(:,:,mm);
    phiprettyima = prettyImage(phiima,shgima,'none',jet(64),uplim,botlim,bright,dark);
    % Here 'phiima' is the phi orientation map 
    prettyphi(:,:,:,mm) = phiprettyima;
end

% third, prepare the pretty directional variance stack, with directional
% variance acquired with localized window
uplim = 1;
botlim = 0;
prettyvar = zeros(sz,sz2,3,sz3);
for mm = 1:sz3
    shgima = shgstackre(:,:,mm);
    varima = Vmatr(:,:,mm);
    varprettyima = prettyImage(varima,shgima,'none',jet(64),uplim,botlim,bright,dark);
    % Here 'varima' is the directional variance map. For directional
    % variance, the range is from 0 to 1. Therefore, 'uplim' and 'botlim'
    % are modified accordingly
    prettyvar(:,:,:,mm) = varprettyima;
end