function [] = ReviewDemo(datasetName)
    %% Program initializations
    global state;
    % Recognition-related stuff.
    global vocabulary, global redundantVocabulary, global modes, global highLevelModes, global vocabCounts;
    
    %% Set parameters for demonstration.
    usedCam = 1; % Which cam to use. If you are using the head, most likely 
    % this number is 2 or 3.
%    sideNeutralZone = 1; % This amount of pixels from each side is not processed.
%    imResizeFactor = 0.5; % The factor of downsampling.
    
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
        vocabCounts = cellfun(@(x) numel(x), vocabulary);
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
           trigger(vid);
           img=getdata(vid,1,'uint8');
   %         img = imread('/Users/rusi/Desktop/FreshGIT/Vision/2D/CHOP/Prototype/input/Swans/vocab/swimming.jpg');
 %           img = imread('/Users/rusi/Desktop/FreshGIT/Vision/2D/CHOP/Prototype/input/bham_45/test/aa.png');
      %       img = imread('/Users/rusi/Desktop/FreshGIT/Vision/2D/CHOP/Prototype/input/bham_45/wideTest/Mug/13/left/13_r0_l.png');
     %        img = imread('/Users/rusi/Dropbox/same_mug/same_mug_3.png');
            emptyResultImg = zeros(size(img,1), size(img,2));
            
            % Crop the image and resize it.
 %           img = img(sideNeutralZone:(end-sideNeutralZone), sideNeutralZone:(end-sideNeutralZone), :);
 %           img = imresize(img, imResizeFactor);
            
            % Process input image and get the result.
            if strcmp(state, 'Active')
                oldImg = img;
                [activationImg, detectionImg, partImg] = processImage(img, datasetName, testFileName, options);
                state = 'Pending';
            elseif strcmp(state, 'Start')
                oldImg = img;
                activationImg = emptyResultImg;
                detectionImg = emptyResultImg;
                partImg = emptyResultImg;
            end
            
            % Display input video frame on axis
            axes(hAxes.axis1);
            imshow(oldImg);
            
            % Show processed frame.
            axes(hAxes.axis2);
            imshow(activationImg);
            
            axes(hAxes.axis3);
            imshow(detectionImg);
            
            % Show processed frame.
            axes(hAxes.axis4);
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
function [activationImg, detectionImg, partImg] = processImage(img, datasetName, testFileName, options)
    global vocabCounts;
    % Copy file to the input folder for test.
    categoryStrs = {'Mug', 'Bowl', 'Tool', 'Pan'};
    inputImageName = 'lastFrameTestDemo';
    if ~exist([pwd '/input/' datasetName '/test'], 'dir')
        mkdir([pwd '/input/' datasetName '/test']);
    end
    imwrite(img, [pwd '/input/' datasetName '/test/' inputImageName '.png']);
    
    [~, fileName, ~] = fileparts(testFileName);
    delete([pwd '/output/' datasetName '/test/inference/' fileName '_test.mat']);
    
    % Run test, and return the output image, along with realizations.
    singleTestImage(testFileName, options);
    smoothedImg = imread([pwd '/output/' datasetName '/smoothed/' fileName '.png']);
    imageSize = [size(smoothedImg,1), size(smoothedImg,2)];
    orgImg = imread([pwd '/output/' datasetName '/original/' fileName '.png']);
    
    if exist([pwd '/output/' datasetName '/test/inference/' fileName '_test.mat'], 'file')
      load([pwd '/output/' datasetName '/test/inference/' fileName '_test.mat'], 'exportArr');
    else
      exportArr = [];
    end
    
    % Update feature parameters.
    options.prediction.feature_params.isTesting = 1;
    overlayImg = zeros(imageSize);
    
    % If no nodes exist, exit.
    if isempty(exportArr)
        partImg = overlayImg;
        detectionImg = partImg;
        return;
    end
    
    % Visualize the activated parts.
    activationImg = GetActivationImage(exportArr, vocabCounts);
    
    % Detect objects in the scene.
    [detections, ~, groupMask] = FindDetections(exportArr, options.prediction.models, options.prediction.feature_params, imageSize);
    
    if isempty(detections)
        partImg = overlayImg;
        detectionImg = partImg;
        return;
    end
    
    % Visualize detection frames.
    detMask = zeros(imageSize);
    objectVisImg = zeros(size(orgImg), 'uint8');
    for detItr = 1:size(detections,1)
        detMask(detections(detItr, 1), detections(detItr, 2)) = detections(detItr,4);
        
        % Visualize results.
        objectMaskStr = [options.currentFolder '/render/' categoryStrs{detections(detItr,4)} ...
            '/' categoryStrs{detections(detItr,4)} num2str(detections(detItr,5)*30) '.png'];
        objectMask = imread(objectMaskStr);
        binaryObjectMask = rgb2gray(objectMask)>0;
        objectMaskSize = [size(objectMask,1), size(objectMask,2)];
        lowBounds = detections(detItr, 1:2) - round([size(objectMask,1), size(objectMask,2)]/2);
        for bandItr = 1:size(orgImg,3)
            overlapImg = objectVisImg(lowBounds(1):(lowBounds(1)+objectMaskSize(1)-1), lowBounds(2):(lowBounds(2)+objectMaskSize(2)-1), bandItr);
            objectMaskBand = objectMask(:,:,bandItr);
            overlapImg(binaryObjectMask) = objectMaskBand(binaryObjectMask);
            objectVisImg(lowBounds(1):(lowBounds(1)+objectMaskSize(1)-1), lowBounds(2):(lowBounds(2)+objectMaskSize(2)-1), bandItr) = overlapImg;
        end
    end
    
    % TODO: Show objectVisImg in gui.
    detectionImg = objectVisImg;
    
    % Get highest-level clean reconstructed image and multiply it with
    % group image. Each object's realizations will thus have a different
    % colour.
    maxLevelId = min(1, (max(exportArr(:,4))-2));
    cleanRealizationMask = double(imread([options.outputFolder '/reconstruction/test/' fileName '_level' num2str(maxLevelId) 'clean.png']));
    cleanRealizationMask = cleanRealizationMask / max(max(cleanRealizationMask));
    groupImg = double(label2rgb(bwlabel(groupMask), 'jet', 'k', 'shuffle'));
    for bandItr = 1:size(groupImg,3)
        groupImg(:,:,bandItr) = groupImg(:,:,bandItr) .* cleanRealizationMask;
    end
    groupImg = uint8(round(groupImg));
    
    % Combine group image with smoothed image.
    cleanRealizationMask = uint8(cleanRealizationMask>0);
    for bandItr = 1:size(groupImg,3)
        groupImgBand = groupImg(:,:,bandItr);
        groupImg(:,:,bandItr) = groupImgBand .* cleanRealizationMask + smoothedImg .* uint8(~cleanRealizationMask);
    end
    partImg = groupImg;
end

%% Function to generate gui.
function [hFig, hAxes] = createFigureAndAxes()

    % Close figure opened by last run
    figTag = 'CVST_VideoOnAxis_9804532';
    close(findobj('tag',figTag));

    % Create new figure
    hFig = figure('numbertitle', 'off', ...
           'name', 'CHOP in Object Detection / Pose Estimation', ...
           'menubar','none', ...
           'toolbar','none', ...
           'resize', 'on', ...
           'tag',figTag, ...
           'renderer','painters', ...
           'position',[200 200 730 600]);

    % Create axes and titles
    hAxes.axis1 = createPanelAxisTitle(hFig,[0.041 0.55 0.438 0.40],'Original Image');%[X Y W H]
    hAxes.axis2 = createPanelAxisTitle(hFig,[0.041 0.1 0.438 0.40],'Activated Compositions');%[X Y W H]
    hAxes.axis3 = createPanelAxisTitle(hFig,[0.52 0.55 0.438 0.40],'Detected Objects');%[X Y W H]
    hAxes.axis4 = createPanelAxisTitle(hFig,[0.52 0.1 0.438 0.40],'Object Decomposition');%[X Y W H]
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
     titlePos = [pos(1)+0.06, pos(2)+pos(3)-0.038, 0.3, 0.033];
     uicontrol('style','text',...
         'String', axisTitle,...
         'Units','Normalized',...
         'Parent',hFig,'Position', titlePos,...
         'BackgroundColor',get(hFig,'Color'));
end

%% Insert buttons to control the gui.
function insertButtons(hFig,hAxes,vid)

    set(hFig, 'Interruptible', 'off');

    % Play button with text Snapshot!/Continue.
    uicontrol(hFig,'unit','pixel','style','pushbutton','string','Snapshot!',...
            'position',[210 15 150 30], 'tag','PBButton123','callback',...
            {@playCallback,vid,hAxes});

    % Exit button with text Exit
    uicontrol(hFig,'unit','pixel','style','pushbutton','string','Exit',...
            'position',[370 15 150 30],'callback', ...
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