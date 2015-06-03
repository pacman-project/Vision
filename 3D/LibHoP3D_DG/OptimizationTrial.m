commonRoot = 'D:/'; 
root = [commonRoot, 'LibHoP3D/'];
addPathsYalmip(root);

fileListPrecomputed = false; 
depthPathDefault = [];
el_path = 'D:\Input Data\VladislavSTD\Vladislav_STD/depth_layer4_after_inh1';
lenDPW = length(el_path);

outRoot = 'D:\Input Data\VladislavSTD\Vladislav_STD/depth_layer4_after_inh2';
subset_len = 1;
is_subset = false;
[list_El, lenF] = extractFileList(fileListPrecomputed, el_path, depthPathDefault, false, subset_len);


parfor i = 1:lenF

    marks = imread(list_El{i});
%     imtool(marks, [0, 50]);
    


    maxGroupSize = 400;
    layerID = 2;
    alpha = 3.0;

    partSize = [15,15];


    [r,c] = size(marks);
    % so far we proceed line by line

    [rows, cols] = find(marks>0);

    [rows, inds] = sort(rows, 'ascend');
    cols = cols(inds);




    % format of the datastructure we create on each step
    %  rows
    %  cols
    %  elements


    rowSel = [];
    colSel = [];
    outImage = zeros(r,c);
    curRow = rows(1);
    numGroup = 0;


    while curRow < rows(end)

        % collect a bin to be optimized   
        numElsInBin = 0;
        numGroup = numGroup + 1;
        prevNumEls = 0;

        if numGroup ~= 1 % attach also 4 previous rows
            indsPG = find(group(1, :) >= curRow - 6);
            group = group(:, indsPG);
            numElsInBin = length(indsPG);
        else
            group = [];
        end

        while numElsInBin < maxGroupSize && curRow < rows(end)

            indsC = find(rows == curRow);
            lenIC = length(indsC);
            if lenIC > 0
                adder = zeros(3, lenIC);
                adder(1, :) = curRow;
                adder(2, :) = cols(indsC);
                inds2 = sub2ind(size(marks), adder(1, :), adder(2, :));
                adder(3, :) = marks(inds2);
                group = [group, adder];
                numElsInBin = numElsInBin + lenIC;
            end

            curRow = curRow + 1;
        end

        curRow  % display a value

        if numElsInBin > 0 % if something is in the bin => optimize
            % at the moment use coverage and overlap constraints only

            standCov = partSize(1)* partSize(2); % standard coverage
            coverages = standCov * ones(1, numElsInBin);
            C = diag(coverages);

            % compute overlap matrix (numElsInBin, numElsInBin)
            listRect = zeros(numElsInBin, 4);
            listRect(:,1) = group(2,:)'; % x-coordinate
            listRect(:,2) = group(1,:)'; % y-coordinate
            listRect(:,3) = partSize(1); % size in x direction
            listRect(:,4) = partSize(2); % size in y direction

            overlaps = rectint(listRect, listRect);
            O = overlaps - C;
            C = C - O/2;

%             apply optimization
            x = binvar(numElsInBin, 1);

            tic

            obj = x'*(O-alpha*C)*x; 
            ops = sdpsettings('solver','bnb', 'bnb.solver', 'fmincon');
            if numGroup == 1
                sol = optimize(integer(x),obj, ops);
            else
                sol = optimize([x(1:prevNumEls) == 1, integer(x)],obj, ops);
            end

            if sol.problem == 0
             % Extract and display value
             xx = value(x);
            else
             display('Hmm, something went wrong!');
             sol.info
             yalmiperror(sol.problem)
            end

            ii = find(xx==1);
            rowSel = [rowSel; group(1,ii)'];
            colSel = [colSel; group(2,ii)'];

            group = group(:, ii);

        end


    end

    % now perform inhibition for these group

    inds3 = sub2ind(size(marks), rowSel, colSel);
    outImage(inds3) = marks(inds3);
    
    marks2 = uint8(outImage);
    curStr = list_El{i};
%       ll = strfind(curStr, '/');
%       fileName = curStr(ll:end);

    fileName = curStr(lenDPW+1:end);
    outFile = [outRoot, fileName];

    ll = strfind(outFile, '/');
    lll = strfind(outFile, '\');
    ll = [ll, lll];
    
    ll = max(ll); % last position
    folderName = outFile(1:ll);
    b = exist(folderName,'dir');

    if b == 0
        mkdir(folderName);
    end
    imwrite(marks2, outFile, 'png');
    
    % imtool(outImage, [1, 50]);

end












