function [ ] = collectMultiviewData(  )
    % Define frame rate
    NumberFrameDisplayPerSecond=10;

    % Open figure
    hFigure=figure(1);

    % Set-up webcam video input
    try
       % For windows
       vid = videoinput('winvideo', 1);
    catch
       try
          % For macs.
          vid = videoinput('macvideo', 3);
       catch
          errordlg('No webcam available');
       end
    end

    % Set parameters for video
    % Acquire only one frame each time
    set(vid,'FramesPerTrigger',1);
    % Go on forever until stopped
    set(vid,'TriggerRepeat',Inf);
    % Get a grayscale image
    set(vid,'ReturnedColorSpace','RGB');
    triggerconfig(vid, 'Manual');

    % set up timer object
    TimerData=timer('TimerFcn', {@FrameRateDisplay,vid},'Period',1/NumberFrameDisplayPerSecond,'ExecutionMode','fixedRate','BusyMode','drop');

    % Start video and timer object
    start(vid);
    start(TimerData);

    % We go on until the figure is closed
    uiwait(hFigure);

    % Clean up everything
    stop(TimerData);
    delete(TimerData);
    stop(vid);
    delete(vid);
    % clear persistent variables
    clear functions;

    % This function is called by the timer to display one frame of the figure
    function FrameRateDisplay(obj, event,vid)
        persistent IM;
        trigger(vid);
        IM=getdata(vid,1,'uint8');

        % if first execution, we create the figure objects
        subplot(3,1,1);
        imshow(IM);
        title('Image from camera 1');

        % if first execution, we create the figure objects
        subplot(3,1,2);
        imshow(IM);
        title('Image from camera 2');

        % if first execution, we create the figure objects
        subplot(3,1,3);
        imshow(IM);
        title('Image from camera 3');
    end
end