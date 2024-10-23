function anim_bd_spectrum(run_id,filename)
%ANIM_BD_ORB Animate the progression of Floquet multipliers along a branch
% obtained via COCO
% Input:
%   run_id: string identifier of the continuation run
%   fs: functions for plotting, default: fs1(t,y)=t wrt fs2(t,y)=y
%   filename: if given save .avi file with this name
% Output:
%   video file with the name filename.avi

% Save as video file
if nargin>1
    myVideo = VideoWriter(filename,'Motion JPEG AVI');
    myVideo.Quality = 95; 
    myVideo.FrameRate = 10;
   open(myVideo)
end

% Plot each orbit one by one
labs = coco_bd_labs(run_id);
figure();
for lab = labs
    [orb, ~] = orb_get_data(run_id,lab);
    plot_spectrum(orb.mu)
    title(sprintf('Orbit: %i',lab));
    drawnow
    % Save frame to video file
    if nargin>1
        frame = getframe(gcf);
        writeVideo(myVideo, frame);
    else
        pause(0.1);
    end
end

if nargin>1
    close(myVideo);
end

end

