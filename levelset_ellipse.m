function seg=levelset_ellipse(img, IM, init_mask, lambda, numiter)

% Check image type -> double precision
[dimy, dimx, c] = size(img);
if(isfloat(img)) % image is a double
    if(c==3) 
      img = rgb2gray(uint8(img)); 
    end
else           % image is a int
    if(c==3) 
      img = rgb2gray(img); 
    end
    img = double(img);
end

imagesc(img), colormap gray, axis image, hold on


alpha = .1; % -- curvature term, usefull ?
Refresh = 20; % -- refresh frequency
rad = 10;       % -- radius of the localization of the data attachment term
method = 1;     % -- 1-Chanvese, 2- Yezzi

seg = segment_one_ellipse(img,IM,init_mask,numiter,lambda,alpha,rad,method,Refresh);


