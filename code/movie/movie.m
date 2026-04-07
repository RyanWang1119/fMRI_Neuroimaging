INDIR = fullfile(pwd, 'data', 'task_residual');
OUTDIR = fullfile(INDIR, 'movie_output');

videoFileName = 'Relational_error_vs_match.mp4';
v = VideoWriter(videoFileName, 'MPEG-4'); 
v.FrameRate = 2; 
v.Quality = 95; 
open(v);

time_steps = 1:41;

for i = 1:length(time_steps)
    t = time_steps(i);
    

    visualize_spatial_patterns('Relational', 'error_vs_match', t);
    
    figHandle = gcf; 
    figure(figHandle); 
    set(figHandle, 'Color', 'w'); 
    set(figHandle, 'Renderer', 'painters'); 
    drawnow; 
    pause(0.1); 
    
    frame = getframe(figHandle); 
   
    if ~isempty(frame.cdata)
        writeVideo(v, frame);
    else
        warning('Frame at time %d was empty.', t);
    end
end

close(v);
fprintf('Video generated: %s\n', fullfile(pwd, videoFileName));