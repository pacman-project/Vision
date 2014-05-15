function [] = ReviewDemo(datasetName)
    %% Program initializations
    global state;
    % Recognition-related stuff.
    global vocabulary, global redundantVocabulary, global modes, global highLevelModes;
    
    %% Set parameters for demonstration.
    usedCam = 1; % Which cam to use. If you are using the head, most likely 
    % this number is 2 or 3.
    sideNeutralZone = 1; % This amount of pixels from each side is not processed.
    imResizeFactor = 0.5; % The factor of downsampling.
    
    %% Prepare input video.
    % Stop video objects if they still exist.
    out = imaqfind;
    delete(out);
    
    state = 'Start'; % 'Start', 'Active', 'Exit', 'Pending'
    try
      vid = videoinput('macvideo', usedCam);
    catch
      errordlg('No webcam available');
    end
    
    % Acquire only one frame each time
    set(vid,'FramesPerTrigger',1);
    % Go on forever until stopped
    set(vid,'TriggerRepeat',Inf);
    % Get a grayscale image
    set(vid,'ReturnedColorSpace','RGB');
    triggerconfig(vid, 'Manual');
    start(vid);
   
    %% Read vocabulary.
    % Set program options, and create relevant folders.
    options = SetParameters(datasetName, false);
    createFolders(options);
    
    % Read vocabulary from the file.
    if exist([pwd '/output/' datasetName '/vb.mat'], 'file')
        load([pwd '/output/' datasetName '/vb.mat'], 'vocabulary', 'redundantVocabulary', 'modes', 'highLevelModes');
    else
        display('No vocabulary exists!');
    end
    % Read models for predictions from file
    if exist([pwd '/demo/Category_Pose/models/' datasetName '_models.mat'], 'file')
        load([pwd '/demo/Category_Pose/models/' datasetName '_models.mat']);
    else
        display('No prediction models exist!');
        models = [];
    end
    options.prediction.models = models;
    options.prediction.feature_params = feature_params;
    testFileName = [pwd '/input/' datasetName '/test/lastFrameTestDemo.png'];
    resultFileName = [pwd '/output/' datasetName '/reconstruction/test/lastFrameTestDemo_level4_leaf.png'];
    
    %% Start gui here.
    [hFig, hAxes] = createFigureAndAxes();
    insertButtons(hFig, hAxes, vid);

    %% Main loop.
    % Process input image, and show the results in GUI.
    oldImg = [];
    partImg = [];
    while true
        set(hFig, 'Interruptible', 'off');
        if ~strcmp(state,'Exit') && isvalid(vid)
            % Get image from the camera.
 %           trigger(vid);
 %           img=getdata(vid,1,'uint8');
            img = imread('4_r10_l_test.png');
            
            % Crop the image and resize it.
            img = img(sideNeutralZone:(end-sideNeutralZone), sideNeutralZone:(end-sideNeutralZone), :);
            img = imresize(img, imResizeFactor);
            
            % Process input image and get the result.
            if strcmp(state, 'Active')
                oldImg = img;
                partImg = processImage(img, datasetName, testFileName, options);
                state = 'Pending';
            elseif strcmp(state, 'Start')
                oldImg = img;
                partImg = img;
            end
            
            % Display input video frame on axis
            axes(hAxes.axis1);
            imshow(oldImg);
            
            % Show processed frame.
            axes(hAxes.axis2);
            imshow(partImg);
        elseif strcmp(state,'Exit')
            % Close the figure window
            stop(vid);
            % Close the video file
            delete(vid);
            close(hFig);
            break;
        end
        set(hFig, 'Interruptible', 'on');
    end
end

%% Process the image and get level 4 part detections.
function partImg = processImage(img, datasetName, testFileName, options)
    % Copy file to the input folder for test.
    inputImageName = 'lastFrameTestDemo';
    if ~exist([pwd '/input/' datasetName '/test'], 'dir')
        mkdir([pwd '/input/' datasetName '/test']);
    end
    imwrite(img, [pwd '/input/' datasetName '/test/' inputImageName '.png']);
    
    % Run test, and return the output image.
    singleTestImage(testFileName, options);
    partImg = img;
    
    [~, fileName, ~] = fileparts(testFileName);
    load([pwd '/output/' datasetName '/test/inference/' fileName '_test.mat'], 'exportArr');
    
    % Find the category/pose of the object.
    [detections, windowSizes] = FindDetections(exportArr, options.prediction.models, options.prediction.feature_params, [size(img,1), size(img,2)]);
    
    % Visualize detection frames.
    overlayImg = zeros(size(partImg,1), size(partImg,2));
    for detItr = 1:size(detections,1)
        halfSizes = round((windowSizes{detections(detItr, 3)} / 2) - 1);
        overlayImg((detections(detItr, 1)-halfSizes(1)):(detections(detItr, 1)+halfSizes(1)), ...
            [(detections(detItr, 2)-halfSizes(2)), (detections(detItr, 2)+halfSizes(2))]) = detections(detItr,4);
        overlayImg([(detections(detItr, 1)-halfSizes(1)),(detections(detItr, 1)+halfSizes(1))], ...
            (detections(detItr, 2)-halfSizes(2)):(detections(detItr, 2)+halfSizes(2))) = detections(detItr,4);
    end
    rgbOverlayImg = label2rgb(overlayImg, 'jet', 'k');
    
    % Overlay detections with rgb image.
    for bandItr = 1:size(img,3)
        bandImg = img(:,:,bandItr);
        bandOverlayImg = rgbOverlayImg(:,:,bandItr);
        bandImg(overlayImg>0) = bandOverlayImg(overlayImg>0);
        rgbOverlayImg(:,:,bandItr) = bandImg;
    end
    partImg = rgbOverlayImg;
end

%% Function to generate gui.
function [hFig, hAxes] = createFigureAndAxes()

    % Close figure opened by last run
    figTag = 'CVST_VideoOnAxis_9804532';
    close(findobj('tag',figTag));

    % Create new figure
    hFig = figure('numbertitle', 'off', ...
           'name', 'Video In Custom GUI', ...
           'menubar','none', ...
           'toolbar','none', ...
           'resize', 'on', ...
           'tag',figTag, ...
           'renderer','painters', ...
           'position',[680 678 480 240]);

    % Create axes and titles
    hAxes.axis1 = createPanelAxisTitle(hFig,[0.1 0.2 0.36 0.6],'Original Video');%[X Y W H]
    hAxes.axis2 = createPanelAxisTitle(hFig,[0.5 0.2 0.36 0.6],'Part detections');
end

%% Name the axes correctly.
function hAxis = createPanelAxisTitle(hFig, pos, axisTitle)

    % Create panel
    hPanel = uipanel('parent',hFig,'Position',pos,'Units','Normalized');

    % Create axis
    hAxis = axes('position',[0 0 1 1],'Parent',hPanel);
    set(hAxis,'xtick',[],'ytick',[],'xcolor',[1 1 1],'ycolor',[1 1 1]);

    % Set video title using uicontrol. uicontrol is used so that text
    % can be positioned in the context of the figure, not the axis.
    titlePos = [pos(1)+0.02 pos(2)+pos(3)+0.3 0.3 0.07];
    uicontrol('style','text',...
        'String', axisTitle,...
        'Units','Normalized',...
        'Parent',hFig,'Position', titlePos,...
        'BackgroundColor',get(hFig,'Color'));
end

%% Insert buttons to control the gui.
function insertButtons(hFig,hAxes,vid)

    set(hFig, 'Interruptible', 'off');

    % Play button with text Start/Pause/Continue
    uicontrol(hFig,'unit','pixel','style','pushbutton','string','Snapshot!',...
            'position',[10 10 75 25], 'tag','PBButton123','callback',...
            {@playCallback,vid,hAxes});

    % Exit button with text Exit
    uicontrol(hFig,'unit','pixel','style','pushbutton','string','Exit',...
            'position',[100 10 50 25],'callback', ...
            {@exitCallback,vid,hFig});
end

%% Get the frame, process it and show the result.
function playCallback(hObject,~,~,~)
   global state;
   try
        % Check the status of play button
        isTextStart = strcmp(get(hObject,'string'),'Snapshot!');
        if (isTextStart)
            set(hObject,'string','Continue');
            state = 'Active';
        else
            set(hObject,'string','Snapshot!');
            state = 'Start';
        end
   catch ME
       % Re-throw error message if it is not related to invalid handle
       if ~strcmp(ME.identifier, 'MATLAB:class:InvalidHandle')
           rethrow(ME);
       end
   end
end

%% Callback for exit button.
function exitCallback(~,~,~,~)
    global state;
    state = 'Exit';
end