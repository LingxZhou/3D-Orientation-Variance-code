function [metric2] = prettyImage(metric,avgi,filename,lut,uplim,botlim,bright,dark)
% Creates pretty images by modulating with total fluorescence. 
% OUTPUT IS NOT INTENDED FOR QUANTITATIVE USE. 
% 
% Output:
% metric2 = modulated metric map 
% avgi2 = top+bottom saturated intensity image.
% If 'filename.tif' already exists, filename is iterated, i.e. 'filename1'.
% If filename is specified as 'none', image is not saved.
% 
% Input:
% metric = metric map 
% avgi = modulating mask for metric, i.e. total SHG intensity image
% Last three inputs are optional:
% lut is the color-map/look-up-table for converting redox ratio to color,
% while uplim and botlim are the upper and bottom limits of the metric 
% for plotting.  Default values are lut = jet(64); 
%% Iterate filename if file already exists
k = 0;
base = filename;
while exist([filename,'.tif'],'file') == 2
  k = k+1;
  filename = [base,num2str(k)];
end   

%% Take total intensity and rescale across 1st to 99th percentile.
ImR=sort(reshape(nonzeros(avgi),1,[]),'descend'); 
maxFlr=ImR(round(dark*length(ImR))); %pick up 99th percentile 
minFlr=ImR(round(bright*length(ImR))); % and 1st values
avgi2=single((avgi-minFlr)/(maxFlr-minFlr)); %setting 1% to 0 and 99% to 1
avgi2=avgi2.*(avgi2<1) + (avgi2>=1); %saturate top 1%
avgi2=avgi2.*(avgi2>=0); %saturate bottom 1%

%% Set default color map and corresponding normalization and offset
if nargin < 6, botlim=0; end
if nargin < 5, uplim=1; end
if nargin < 4, lut = jet(64); end
bitScale = size(lut,1)-1;

%% Rescale metric map to bit scale
metric=(metric-botlim)/(uplim-botlim);
metric=metric.*(metric<1) + (metric>=1);
metric=round(bitScale*(metric.*(metric>=0)))+1;
    
%% Apply color map, convert to RGB, and modulate by total intensity
imageSize1=size(avgi,1);
imageSize2=size(avgi,2);
stackSize=size(avgi,3);
metric2 = single(ones(imageSize1,imageSize2,3,stackSize));
for istack=1:stackSize    
    metric2(:,:,:,istack) = ind2rgb(metric(:,:,istack),lut);
    for k = 1:3
        metric2(:,:,k,istack) = metric2(:,:,k,istack).*avgi2(:,:,istack);
    end %for k = 1:3
    
    % Save redox stack as multipage tiff
    if ~strcmpi(filename,'none')
        imwrite(metric2(:,:,:,istack),[filename,'.tif'],'Compression','none','WriteMode','append');
    end %if ~strcmpi(filename,'none')
  
end %for istack=1:stack

    if ~strcmpi(filename,'none')
        disp([filename,'.tif saved.']); 
    end %if ~strcmpi(filename,'none')

return