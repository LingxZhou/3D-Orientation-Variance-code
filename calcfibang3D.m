function [aSDE,pSDER,bSDER,gSDER]=calcfibang3D(shgim,armdx,armdz,filton,para)
% This function calculates the voxel-wise orientation. The output and input
% indices are explained in the main program
im = shgim;

im = double(im);
sz = size(im,1);
sz2 = size(im,2);
sz3 = size(im,3);

% Here starts the determination of theta
terimc = zeros(sz,sz2,3*sz3);
terimc(:,:,1:sz3) = im;
terimc(:,:,(sz3+1):(2*sz3)) = im;
terimc(:,:,(2*sz3+1):(3*sz3)) = im;

imc = zeros(sz,sz2,sz3);

for i = 1:sz3
imc(:,:,i) = mean(terimc(:,:,(sz3+i-armdz):(sz3+i+armdz)),3);
end

clear terimc

meanc = zeros(sz,sz2,sz3);
h = waitbar(0,'Please Wait')
ij = 0;
nij = sz3;
for i = 1:sz3
    ij = ij+1;
    waitbar(ij/nij)
    meanc(:,:,i) = calcfibangspeed(imc(:,:,i),armdx,filton);
end
close(h)
clear imc

% Here starts the determination of beta
terimj = zeros(3*sz,sz2,sz3);
terimj(1:sz,:,:) = im;
terimj((sz+1):(2*sz),:,:) = im;
terimj((2*sz+1):(3*sz),:,:) = im;

imj = zeros(sz,sz2,sz3);
for i = 1:sz
    imj(i,:,:) = mean(terimj((sz+i-armdz):(sz+i+armdz),:,:),1);
end
clear terimj
% To calculate the orientation in 2D manner
meanj = zeros(sz,sz2,sz3);
h = waitbar(0,'Please Wait')
ij = 0;
nij = sz;
for i = 1:sz
    ij = ij+1;
    waitbar(ij/nij)
    meanj(i,:,:) = calcfibangspeed(squeeze(imj(i,:,:)),armdx,filton);
end
close(h)
clear imj
meanj = (pi/2-meanj).*(meanj<=pi/2)+(3*pi/2-meanj).*(meanj>pi/2);

% Here starts the determination of gamma

terimg = zeros(sz,3*sz2,sz3);
terimg(:,1:sz2,:) = im;
terimg(:,(sz2+1):(2*sz2),:) = im;
terimg(:,(2*sz2+1):(3*sz2),:) = im;

img = zeros(sz,sz2,sz3);
for i = 1:sz
    img(:,i,:) = mean(terimg(:,(sz2+i-armdz):(sz2+i+armdz),:),2);
end
clear terimg
% To calculate the orientation in 2D manner
meang = zeros(sz,sz2,sz3);
h = waitbar(0,'Please Wait')
ij = 0;
nij = sz2;
for i = 1:sz2
    ij = ij+1;
    waitbar(ij/nij)
    meang(:,i,:) = calcfibangspeed(squeeze(img(:,i,:)),armdx,filton);
end
close(h)
clear img
meang = (pi/2+meang).*(meang<=pi/2)+(meang-(pi/2)).*(meang>pi/2);

% Here to acquire the orientation of phi
meanp = atan(sqrt(1./(tan(meanj).^2)+1./(tan(meang).^2)));
meanplim = meanp;
for i = 1:sz
    for j = 1:sz2
        for k = 1:sz3
            if meanj(i,j,k) <= pi/2
                meanp(i,j,k) = meanp(i,j,k);
            else
                meanp(i,j,k) = pi-meanp(i,j,k);
            end
        end
    end
end

% Here to acquire the real phi orientation regardless of Z resolution
meanpreal = (pi/2)-(atan(para*tan((pi/2)-meanplim)));
for i = 1:sz
    for j = 1:sz2
        for k = 1:sz3
            if meanj(i,j,k) <= pi/2
                meanpreal(i,j,k) = meanpreal(i,j,k);
            else
                meanpreal(i,j,k) = pi-meanpreal(i,j,k);
            end
        end
    end
end

% Here to calculate the real beta angle regardless of Z resolution
meanjreal = atan(para*tan(meanj));
meanjreal = meanjreal.*(meanjreal>=0)+(meanjreal+pi).*(meanjreal<0);

% Here to calculate the real gamma angle regardless of Z resolution
meangreal = atan(para*tan(meang));
meangreal = meangreal.*(meangreal>=0)+(meangreal+pi).*(meangreal<0);

% Here prepares the final output
aSDE = meanc*180/pi;
pSDER = meanpreal*180/pi;
bSDER = meanjreal*180/pi;
gSDER = meangreal*180/pi;



  
  



