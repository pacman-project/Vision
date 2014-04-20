%> Name: collectMultiviewData
%>
%> Description: Given the class name and object id, this function collects
%> images of the object from different viewpoints. We use a turntable and
%> two firewire cameras to obtain multi-view images.
%>
%> @param className Name of the class the object belongs to.
%> @param objectID Identifier of the object.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 15.04.2014
function [ ] = collectMultiviewData( className, objectID )
    % If a new object is trained, objectID can be skipped.
    if nargin < 2
        % Find a new object id.
        objectID = 1;
        objectDir = [pwd '/data/' className '/1'];
        while exist(objectDir, 'dir')
            objectID = objectID + 1;
            objectDir = [pwd '/data/' className '/' num2str(objectID)];
        end
    end

    %% Initialize parameters
    angleStepper = 5;
    stepsPerFullRotation = 7200;
    stepsPerAngle = (angleStepper/360) * stepsPerFullRotation;
%    persistent rt;
%    persistent vid;
%    persistent vid2;
    
    %% Check hardware
    try
        rt=serial('/dev/tty.usbserial');
        fopen(rt);
    catch
        errordlg('Rotary table not connected. Perhaps forgot to turn it on?');
    end
    
    % Open figure
    figure, hold on;
    set(gcf, 'Visible', 'off');

    % Set-up webcam video input
    try
      % Get images from both cams.
      vid = videoinput('macvideo', 3);
 %     vid2 = videoinput('macvideo', 3);
    catch
      errordlg('No webcam available');
    end
    
    %% Create folder structure.
    objectDir = [pwd '/data/' className '/' num2str(objectID)];
    if exist(objectDir, 'dir')
        rmdir(objectDir, 's');
    end
    mkdir(objectDir);
    mkdir([objectDir '/left']);
    mkdir([objectDir '/right']);
    
    %% Initialize video params.
    % Set parameters for video # 1
    % Acquire only one frame each time
    set(vid,'FramesPerTrigger',1);
    % Go on forever until stopped
    set(vid,'TriggerRepeat',Inf);
    % Get a grayscale image
    set(vid,'ReturnedColorSpace','RGB');
    triggerconfig(vid, 'Manual');
%     
%     % Set parameters for video # 2
%     % Acquire only one frame each time
%     set(vid2,'FramesPerTrigger',1);
%     % Go on forever until stopped
%     set(vid2,'TriggerRepeat',Inf);
%     % Get a grayscale image
%     set(vid2,'ReturnedColorSpace','RGB');
%     triggerconfig(vid2, 'Manual');

    %% Initialize rotary table.
    % Energize the motors.
    fprintf(rt, sprintf('!1we11\n'));
    % Remove limits on each side.
    pause(1);
    % Remove limits on each side.
    fprintf(rt, sprintf('!1wl0\n'));
    pause(1);
    fprintf(rt, sprintf('!1wr0\n'));
    pause(1);
    
    % Go home.
    ApplyRTBlocking(rt, '!1h1\n');
%     display('Now, please put the object onto the rotary table and press any key.');
%     f=figure;
%     set(f, 'Visible', 'off');
%     waitforbuttonpress;
%     close(f);
%     display('The turntable will now do a full rotation for you to verify the object is centered.'); 
%     display('Feel free to move the object as you wish to center it.');
%     pause(5);
%     ApplyRTBlocking(rt, ['!1m1f' num2str(stepsPerFullRotation) '\n']);
%     ApplyRTBlocking(rt, '!1h1\n');
    
    % Get confirmation, and we're ready to go.
%     display('Press any key to continue.');
%     f=figure;
%     set(f, 'Visible', 'off');
%     waitforbuttonpress;
%     close(f);
    
    % Start videos and timer object
    start(vid);
%   start(vid2);

    set(gcf, 'Visible', 'on');
    %% Get images from both cameras and save them iteratively.
    for imgAngle = 0:angleStepper:(360-angleStepper)
 %       start(vid);
        trigger(vid);
        img=getdata(vid,1,'uint8');
 %       stop(vid);
%         start(vid2);
%         trigger(vid2);
%         img2=getdata(vid2,1,'uint8');
%         stop(vid2);
    
        % Show image 1
%         subplot(1,2,1);
        imshow(img);
        title('Image from left camera');

%         % Show image 2
%         subplot(1,2,2);
%         imshow(img2);
%         title('Image from right camera');
        
        %% Save the images in their corresponding positions.
        imwrite(img, [objectDir '/left/' num2str(objectID) '_r' num2str(imgAngle) '_l.png']);
%         imwrite(img2, [objectDir '/right/' num2str(objectID) '_r' num2str(imgAngle) '_r.png']);
        
        %% Step forward.
        ApplyRTBlocking(rt, ['!1m1f' num2str(stepsPerAngle) '\n']);
    end
    
    % Go home and stay there.
    ApplyRTBlocking(rt, '!1h1\n');
    
    % Stop updating the figure.
    hold off;
    
    % Clean up everything and return!
    fclose(rt);
    delete(rt);
    stop(vid);
    delete(vid);
%    stop(vid2);
%     delete(vid2);
    close(gcf);
end

%% Function to apply rotary table commands and block until they are finished.
function [] = ApplyRTBlocking(rt, action)
    prevState = rt.BytesAvailable;
    fprintf(rt, sprintf(action));
    newState = prevState;
    while newState == prevState
       newState =  rt.ValuesReceived;
       pause(0.05);
    end
end