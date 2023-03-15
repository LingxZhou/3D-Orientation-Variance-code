function oimg = miprmdebrisn(img,excSize)

%   This function removes structures smaller than EXCSIZE 
%   in the binary image 'img'. The output 'oimg' is also a binary image

L = bwlabeln(img,18); % 18 is the neighborhood size
S = regionprops(L,'Area');
oimg = ismember(L,find([S.Area] >= excSize));
