


close all
clear all
clc

I = analyze75read('Changhai020bb-coronal.img');
I=permute(I,[2,3,1]);
Icut = I;
[sizex,sizey,sizez]=size(Icut);

imagD = double(Icut);
m = max(max(imagD(:))) ; n = min(min(imagD(:))) ; 
Img = uint8((imagD-n)*255/(m-n));
%Img = abs(255-Img);

for i = 1:sizez
    PreImg = imadjust(Img(:,:,i),[0 1],[],1);
    ImgA(:,:,i) = PreImg;
end
plot_3Dto2D(ImgA,5);
break
a=input('minimum #slice number: ');
b=input('maximum #slice number: ');
c=input('select #slice number: ');

im=ImgA(:,:,c);
figure, imagesc(im), colormap gray, axis image, hold on

% - initialization
PointList = CreateContour(im);
m = roipoly(im,PointList(:,2),PointList(:,1));

% - plot initialization
imagesc(im), colormap gray, axis image, hold on,
hold on
plot(PointList(:,2),PointList(:,1),'g--','LineWidth',2);

% descending aorta segmentation
seg = levelset_ellipse(im, m, 0.6, 300);  %-- run segmentation
res_inner(:,:,c)=seg;
segn=seg;
for i=c+1:b
    SE = strel('disk',1);
    seg = imerode(seg,SE);
    m=seg;
    im=ImgA(:,:,i);
    seg = levelset_ellipse(im, m, 0.6, 300);
    res_inner(:,:,i)=seg;
end
seg=segn;
for i=c-1:-1:a
    SE = strel('disk',1);
    seg = imerode(seg,SE);
    m=seg;
    im=ImgA(:,:,i);
    seg = levelset_ellipse(im, m, 0.6, 300);
    res_inner(:,:,i)=seg;
end

ImgB=Ireplacement(ImgA, res_inner, sizex, sizey,sizez);

im=ImgB(:,:,c);
figure, imagesc(im), colormap gray, axis image, hold on

% - initialization
PointList = CreateContour(im);
m = roipoly(im,PointList(:,2),PointList(:,1));

% - plot initialization
imagesc(im), colormap gray, axis image, hold on,
hold on
plot(PointList(:,2),PointList(:,1),'g--','LineWidth',2);

% - descending aorta segmentation
seg = levelset_ellipse(im, m, 0.9, 300);  %-- run segmentation
res_outer(:,:,c)=seg;
segn=seg;
for i=c+1:b
    SE = strel('disk',1);
    seg = imerode(seg,SE);
    m=seg;
    im=ImgB(:,:,i);
    seg = levelset_ellipse(im, m, 0.9, 300);
    res_outer(:,:,i)=seg;
end
seg=segn;
for i=c-1:-1:a
    SE = strel('disk',1);
    seg = imerode(seg,SE);
    m=seg;
    im=ImgB(:,:,i);
    seg = levelset_ellipse(im, m, 0.9, 300);
    res_outer(:,:,i)=seg;
end
