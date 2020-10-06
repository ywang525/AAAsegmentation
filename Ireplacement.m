
function Inew=Ireplacement(I, resi,sizex1,sizey1,sizez1)

    a1=I;
    b1=resi; 
    b11=b1;
    SE=strel('disk',1);
    b12=imdilate(b11,SE);
    bregion=find(b12);
    bbbregion=find(b12-b11);
    [x,y]=size(bbbregion);
    total=sum(a1(bbbregion));
    ave=total/x;
    r=normrnd(ave-20,20,sizex1,sizey1,sizez1);
    a1(bregion)=r(bregion);
Inew=a1;

