% performs inference of the first layer parts

function performInference2DG(list_input, lenF, inputDataType, receptiveField, inPath, zScale, filtOptions, is_overwrite, outRoot)

    options = GetOptions(filtOptions);
    options.elementRadius = receptiveField;
    
    indsSS = randperm(lenF);
    lenDPW = length(inPath);

    parfor i = 1:lenF
        
        curStr = list_input{indsSS(i)};
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
        
        % chech if file exists
        b = exist(outFile, 'file');

        if ~b || is_overwrite
        
            if inputDataType == 1  % depth images
%                 if ~ strcmp(list_input{indsSS(i)}, 'D:/Input Data\VladislavSTD\Vladislav_STD\depthEdgesCorners/99_1_1_2_3_3.png')
%                     continue;
%                 end

                I = imread(list_input{indsSS(i)});
                I = double(I(:,:,1));
                I = I * zScale;
                mask = [];
                [I, ~, ~, mask, is_successfull] = preliminaryProcessing(I, mask, filtOptions);
                if ~ is_successfull
                    continue;
                end

%                 % visualize this image in 3D
%                 surf(I, 'FaceColor',[0.3, 0.3 0.3], 'FaceAlpha', 0.4, 'EdgeColor', 'none', 'FaceLighting', 'phong');
%                 camlight left
%                 axis equal;
%                 hold on

                [marksOut, Normals] = performInference1DG(I, mask, options);
%                 imtool(marksOut, [0, max(max(marksOut))]);

                marksOut = uint8(marksOut);
                normFile = [outFile(1:end - 4), '_N.mat'];
                
                parSave(normFile, Normals); 
                imwrite(marksOut, outFile, 'png');

            elseif inputDataType == 2  % meshes

                [V, F, N] = meshRead(list_input{indsSS(i)});

                % compute all differential parameters in every point
                [Umin,Umax,Cmin,Cmax,Cmean,Cgauss,faceNs] = compute_curvature(V,F,options);

                %chech if this patch can be considered as a planar one
                a = 2;
            end
        end
        
        i
    end

end


function parSave(normFile, Normals)
    save(normFile, 'Normals'); 
end

