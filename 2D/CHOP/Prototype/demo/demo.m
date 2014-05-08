%> Name: runDemo
%>
%> Description: Run a demonstration of the capabilities of CHOP. The image
%> obtained from one of the cameras is processed and results are shown in a
%> frame-by-frame basis. If inference is sped up, this can be a real-time
%> demonstration.
%>
%> @param datasetName Name of the dataset.
%> 
%> @retval options Program options.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 01.02.2014
%> Ver 1.1 on 12.01.2014 Timing is added by Mete
function [ ] = demo( datasetName )
    %% Demo parameters.
    usedCam = 1;

    %% Obtain video from camera. This 
    % Set-up webcam video input
    try
      vid = videoinput('macvideo', usedCam);
    catch
      errordlg('No webcam available');
    end
    
    %% Initialize video params.
    % Set parameters for video # 1
    % Acquire only one frame each time
    set(vid,'FramesPerTrigger',1);
    % Go on forever until stopped
    set(vid,'TriggerRepeat',Inf);
    % Get a grayscale image
    set(vid,'ReturnedColorSpace','RGB');
    triggerconfig(vid, 'Manual');
    % Start video
    start(vid);
    
    %% Get images from both cameras and save them iteratively.
    handle = gui();
    while(true)
        trigger(vid);
        img=getdata(vid,1,'uint8');
    
        % Show image 1
        imshow(img);
        title('Image from left camera');
        
        %% Save the images in their corresponding positions.
        imwrite(img, [objectDir '/left/' num2str(objectID) '_r' num2str(imgAngle) '_l.png']);
        
        %% Step forward.
        ApplyRTBlocking(rt, ['!1m1f' num2str(stepsPerAngle) '\n']);
    end
    
    % Stop updating the figure.
    hold off;
    
    % Clean up everything and return!
    stop(vid);
    delete(vid);
    close(handle);



end

