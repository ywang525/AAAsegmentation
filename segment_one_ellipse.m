
%------------------------------------------------------------------------
function seg = segment_one_ellipse(img,IM,init_mask,max_its,lambda,alpha,rad,method,displayUpdate)
  
    %-- Create a signed distance map (SDF) from mask
    phi = mask2phi(init_mask);
    
%     display = 1;
%     S.Xc=1;
%     S.Yc=1;
%     S.A=2;
%     S.B=2;
%     S.Phi=3;
%     S.P = [1;1;1;1;1;1];
%     varargout=S;
    
    %--main loop    
    for its=0:max_its   % Note: no automatic convergence test        
        
        % get the curve's narrow band 
%         idx = find( phi <= 0.8 & phi >= -0.8 );  % get the curve's narrow band 
%         [y x] = ind2sub(size(phi),idx);

        % shape constraint
        % [S,parms] = get_pose_force(phi(idx),y,x,P);
        
        % -- ellipse fitting
        [y,x,idx,parms] = Improved_LS_fit(phi);
        psi = sampson_dist( [x(:) y(:)] , parms.P );
        % -- build psi function
        S = - 2*(phi(idx) - psi);
        
        %-- Feature part based on gradient information
        F = get_feature_function(img,phi,idx,rad,method);
        
        %-- get forces from curvature penalty
        phiidx = phi(idx);
        curvature = get_curvature(phi,idx,x,y);  
        
        %-- gradient descent to minimize energy

        F = F./max(abs(F));
        S = S./max(abs(S));
        dphidt = F + lambda*S + alpha*curvature;
        
%         dphidt = F + alpha*curvature;
%         dphidt = S;
                       
        %-- maintain the CFL condition
        dt = .45/(max(dphidt)+eps);
        
        %-- evolve the curve
        phi(idx) = phiidx + dt.*dphidt;
        
        %-- Keep SDF smooth
        phi = sussman(phi,.5);
        
        %-- intermediate output
%         if((display>0)&&(mod(its,displayUpdate) == 0))
          if((mod(its,displayUpdate) == 0))
            showCurveAndPhi(IM,phi,its);
        end
    end

    %-- final output
%     if(display)
        showCurveAndPhi(IM,phi,its);
%     end

    %-- make mask from SDF
    seg = double(phi<=0); %-- Get mask from levelset

  
%---------------------------------------------------------------------
%---------------------------------------------------------------------
%-- AUXILIARY FUNCTIONS ----------------------------------------------
%---------------------------------------------------------------------
%---------------------------------------------------------------------
  
  
%-- Displays the image with curve superimposed
function showCurveAndPhi(IM, phi, i)
    
    imagesc(IM); colormap(gray); axis image; axis on; hold on;
    contour(phi, [0 0], 'g','LineWidth',1);
    %contour(phi, [0 0], 'k','LineWidth',2);
    if(~exist('i','var')) 
        hold off; drawnow;
    else
        hold off; title([num2str(i) ' Iterations']); drawnow;
    end
   
%-- converts a mask to a SDF
function phi = mask2phi(init_a)
    phi=bwdist(init_a)-bwdist(1-init_a)+im2double(init_a)-.5;
 
   
%-- Feature part based on gradient information
function [F] = get_feature_function(I,phi,idx,rad,method)

    
% %     upts = find(phi<=0);                 % interior points
% %     vpts = find(phi>0);                  % exterior points
% %     u = sum(I(upts))/(length(upts)+eps); % interior mean
% %     v = sum(I(vpts))/(length(vpts)+eps); % exterior mean
% %     
% %     F = (I(idx)-u).^2-(I(idx)-v).^2;         % force from image information


    %-- get the curve's narrow band
    [dimy dimx] = size(phi);
    [y x] = ind2sub(size(phi),idx);
    
    
    
    %-- get windows for localized statistics
    xneg = x-rad; xpos = x+rad;      %get subscripts for local regions
    yneg = y-rad; ypos = y+rad;
    xneg(xneg<1)=1; yneg(yneg<1)=1;  %check bounds
    xpos(xpos>dimx)=dimx; ypos(ypos>dimy)=dimy;

    %-- re-initialize u,v,Ain,Aout
    u=zeros(size(idx)); v=zeros(size(idx)); 
    Ain=zeros(size(idx)); Aout=zeros(size(idx)); 
    
    %-- compute local stats
    for i = 1:numel(idx)  % for every point in the narrow band
        img = I(yneg(i):ypos(i),xneg(i):xpos(i)); %sub image
        P = phi(yneg(i):ypos(i),xneg(i):xpos(i)); %sub phi
        upts = find(P<=0);            %local interior
        Ain(i) = length(upts)+eps;
        u(i) = sum(img(upts))/Ain(i);
        vpts = find(P>0);             %local exterior
        Aout(i) = length(vpts)+eps;
        v(i) = sum(img(vpts))/Aout(i);
    end   
 
    %-- get image-based forces
    switch method  %-choose which energy is localized
    case 1,                 %-- CHAN VESE
        F = -(u-v).*(2.*I(idx)-u-v);
    otherwise,              %-- YEZZI
        F = -((u-v).*((I(idx)-u)./Ain+(I(idx)-v)./Aout));
    end
    
  
%-- level set re-initialization by the sussman method
function D = sussman(D, dt)
    
    % forward/backward differences
    a = D - shiftR(D); % backward
    b = shiftL(D) - D; % forward
    c = D - shiftD(D); % backward
    d = shiftU(D) - D; % forward

    a_p = a;  a_n = a; % a+ and a-
    b_p = b;  b_n = b;
    c_p = c;  c_n = c;
    d_p = d;  d_n = d;

    a_p(a < 0) = 0;
    a_n(a > 0) = 0;
    b_p(b < 0) = 0;
    b_n(b > 0) = 0;
    c_p(c < 0) = 0;
    c_n(c > 0) = 0;
    d_p(d < 0) = 0;
    d_n(d > 0) = 0;

    dD = zeros(size(D));
    D_neg_ind = find(D < 0);
    D_pos_ind = find(D > 0);
    dD(D_pos_ind) = sqrt(max(a_p(D_pos_ind).^2, b_n(D_pos_ind).^2) ...
                       + max(c_p(D_pos_ind).^2, d_n(D_pos_ind).^2)) - 1;
    dD(D_neg_ind) = sqrt(max(a_n(D_neg_ind).^2, b_p(D_neg_ind).^2) ...
                       + max(c_n(D_neg_ind).^2, d_p(D_neg_ind).^2)) - 1;

    D = D - dt .* sussman_sign(D) .* dD;   
    %-- compute curvature along SDF
function curvature = get_curvature(phi,idx,x,y)
    [dimy, dimx] = size(phi);        

    %-- get subscripts of neighbors
    ym1 = y-1; xm1 = x-1; yp1 = y+1; xp1 = x+1;

    %-- bounds checking  
    ym1(ym1<1) = 1; xm1(xm1<1) = 1;              
    yp1(yp1>dimy)=dimy; xp1(xp1>dimx) = dimx;    

    %-- get indexes for 8 neighbors
    idup = sub2ind(size(phi),yp1,x);    
    iddn = sub2ind(size(phi),ym1,x);
    idlt = sub2ind(size(phi),y,xm1);
    idrt = sub2ind(size(phi),y,xp1);
    idul = sub2ind(size(phi),yp1,xm1);
    idur = sub2ind(size(phi),yp1,xp1);
    iddl = sub2ind(size(phi),ym1,xm1);
    iddr = sub2ind(size(phi),ym1,xp1);
    
    %-- get central derivatives of SDF at x,y
    phi_x  = -phi(idlt)+phi(idrt);
    phi_y  = -phi(iddn)+phi(idup);
    phi_xx = phi(idlt)-2*phi(idx)+phi(idrt);
    phi_yy = phi(iddn)-2*phi(idx)+phi(idup);
    phi_xy = -0.25*phi(iddl)-0.25*phi(idur)...
             +0.25*phi(iddr)+0.25*phi(idul);
    phi_x2 = phi_x.^2;
    phi_y2 = phi_y.^2;
    
    %-- compute curvature (Kappa)
    curvature = ((phi_x2.*phi_yy + phi_y2.*phi_xx - 2*phi_x.*phi_y.*phi_xy)./...
              (phi_x2 + phi_y2 +eps).^(3/2)).*(phi_x2 + phi_y2).^(1/2);        
  
%-- whole matrix derivatives
function shift = shiftD(M)
    shift = shiftR(M')';

function shift = shiftL(M)
    shift = [ M(:,2:size(M,2)) M(:,size(M,2)) ];

function shift = shiftR(M)
  shift = [ M(:,1) M(:,1:size(M,2)-1) ];

function shift = shiftU(M)
    shift = shiftL(M')';
  
function S = sussman_sign(D)
    S = D ./ sqrt(D.^2 + 1);    

  
   

  

