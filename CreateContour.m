%% Contour Drawing Program 

function [ctr] = CreateContour(img)

    but = 0;
    ctr = zeros(0,2);
%     imagesc(img); colormap(gray); axis image; hold on;
    text(6,6,'Left = Add point, Right = End ROI','FontSize',[8],'Color', 'r');        

    % % set(handles.Figure, 'String', 'Left = Add point, Middle = Quit');

    while ( but~= 3 )
        [x,y,but] = ginput(1);
        if ( but ~= 2 )
            ctr(end+1,1) = y;
            ctr(end,2) = x;
            imagesc(img); colormap gray; axis image;
            hold on; plot(ctr(:,2),ctr(:,1),'b','linewidth',2);                        
            tmp = ctr(1,:); tmp(end+1,:) = ctr(end,:);
            plot(tmp(:,2), tmp(:,1), 'b', 'linewidth', 2);            
            plot(ctr(:, 2), ctr(:,1), 'or', 'linewidth', 2); hold off;
            hold off;
            text(6,6,'Left = Add point, Middle = End ROI','FontSize',[8],'Color', 'r');
        end
    end

    ctr(end+1, :) = ctr(1,:);   
       
