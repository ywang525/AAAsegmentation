function plot_3Dto2D(I, L, newfigoff, mi, ma)

% plot a 3D matrix as 2D
% input: matrix, number of images per row, window level low, high

%clear tp Itp
Itp=squeeze(I);
%Itp = squeeze(I(61:200,61:200,:,4));
%Itp = squeeze(I(:,:,end/2,:));

tp=reshape(Itp,[size(Itp,1)*size(Itp,2),size(Itp,3)]);
tp=reshape(tp,[size(Itp,1),size(Itp,2)*size(Itp,3)]);

if nargin==1 
   L=1
else 
   tp2=tp(:,1:end/L);
   for ii=2:L
       tp2=[tp2;tp(:,end/L*(ii-1)+1:end/L*ii)];
   end
   tp=tp2;
end

if imag(tp)~=0
 tp =abs(tp);
end

if nargin<3
figure
else
if newfigoff==0
figure
end
end

if nargin >3
	imshow(tp, [mi ma])
else
	imshow((tp),[])
end