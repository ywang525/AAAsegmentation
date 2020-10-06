


close all
clear all
clc

[filename, pathname] = uigetfile( ...
{'*.ima','All Image Files';...
'*.*','All Files' },...
'select images', ...
'MultiSelect', 'on');

for i1=1:length(filename)
    image = dicomread(strcat(pathname,filename{i1}));
    info(:,:,i1) = dicominfo(strcat(pathname,filename{i1}));
    IMG(:,:,i1) = image;
end

thickness=abs(info(1,1,1).ImagePositionPatient-info(1,1,2).ImagePositionPatient);
solutionz=thickness(3,1);
solutionx=info(1,1,1). PixelSpacing(1,1);
solutiony=info(1,1,1). PixelSpacing(2,1);

Ic = IMG-1024;
%Ic = IMG;

imshow(uint8(Ic(:,:,200)));
hold on
[x,y]=ginput(2);
Icut=Ic(y(1):y(2),x(1):x(2),:);
region=(find(Icut>120));
seg_inner=zeros(size(Icut));
seg_inner(region)=1;
plot_3Dto2D(seg_inner(:,:,151:250),10);

[sizex1, sizey1, sizez1]=size(Icut);
r=normrnd(50,10,sizex1,sizey1,sizez1);
Icut(region)=r(region);
plot_3Dto2D(Icut(:,:,151:250),10);

a=input('minimum #slice number: ');
b=input('maximum #slice number: ');
c=input('select #slice number: ');

seg_in=seg_inner(:,:,a-20:b+50);
imagD = double(Icut(:,:,a-20:b+50));
m = max(max(imagD(:))) ; n = min(min(imagD(:))) ; 
Img = uint8((imagD-n)*255/(m-n));

for i = 1:b-a+71
    PreImg = imadjust(Img(:,:,i),[0 0.7],[],1);
    I(:,:,i) = PreImg;
end
plot_3Dto2D(I,10);

im=I(:,:,c-a+20);
figure, imagesc(im), colormap gray, axis image, hold on

% - initialization
PointList = CreateContour(im);
m = roipoly(im,PointList(:,2),PointList(:,1));

% - plot initialization
imagesc(im), colormap gray, axis image, hold on,
hold on
plot(PointList(:,2),PointList(:,1),'g--','LineWidth',2);

% descending aorta segmentation
seg = levelset_ellipse(im, m, 0.8, 200);  %-- run segmentation
res=zeros(size(I));
res(:,:,c-a+20)=seg;
segn=seg;
for i=c-a+21:b-a+70
    SE = strel('disk',1);
    seg = imerode(seg,SE);
    m=seg;
    im=I(:,:,i);
    seg = levelset_ellipse(im, m, 0.8, 200);
    res(:,:,i)=seg;
end
seg=segn;
for i=c-a+19:-1:1
    SE = strel('disk',1);
    seg = imerode(seg,SE);
    m=seg;
    im=I(:,:,i);
    seg = levelset_ellipse(im, m, 0.8, 200);
    res(:,:,i)=seg;
end

V=size(find(res(:,:,21:b-a+21)),1)*solutionx*solutiony*solutionz*0.001;
vtkImageWriter(res,'changhai-ct-020.vtk',[solutionx, solutiony, solutionz]);
save changhai-ct-020.mat

     