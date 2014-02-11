%> Name: discoverSubs
%>
%> Description: This function runs SUBDUE (executable) with the input object 
%> graph level and parameters defined in options. 
%>
%> @param graphFileName The input graph's path.
%> @param resultFileName The result file's path.
%> @param options Program options.
%> @param currentFolder Path to the workspace folder.
%> @param preDefinedFileName The file that includes pre-defined subs.
%>
%> @retval []
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 05.02.2014
function [vocabLevel, graphLevel] = runSubdueExec(graphFileName, resultFileName, options, currentFolder, preDefinedFileName)
    % Set SUBDUE options
    subdueOptions = ' ';
    if options.subdue.diverse
       subdueOptions = [subdueOptions '-diverse '];
    end
    if options.subdue.valuebased
       subdueOptions = [subdueOptions '-valuebased '];
    end
    if options.subdue.overlap
       subdueOptions = [subdueOptions '-overlap '];
    end
    if options.subdue.threshold > 0.0001
       subdueOptions = [subdueOptions '-threshold ' num2str(options.subdue.threshold) ' '];
    end
    if ~isempty(preDefinedFileName)
       subdueOptions = [subdueOptions '-ps ' preDefinedFileName ' -discovery 1 '];
    end

    subdueOptions = [subdueOptions '-nsubs ' num2str(options.subdue.nsubs) ...
                        ' -minsize ' num2str(options.subdue.minSize) ...
                        ' -maxsize ' num2str(options.subdue.maxSize) ...
                        ' -beam ' num2str(options.subdue.beam) ...
                        ' -fileoutput 3 ' ...
                        ' -out ' resultFileName ' ' graphFileName];

    % Form the command and run subdue with specified options
    if ismac
        command = [currentFolder '/miners/subdueMac' subdueOptions];
    elseif isunix
        command = [currentFolder '/miners/subdueUnix' subdueOptions];
    elseif ispc
        command = [currentFolder '/miners/subdueWin' subdueOptions];

        % If Windows, replace all / in command with \.
        command = strrep(command, '/', options.subdue.winSep);
    end

    if isempty(preDefinedFileName)
   %     command
        [~, ~] = system(command);
    else
        [~, ~] = system(command);
    end

    [vocabLevel, graphLevel] = parseResultFile(resultFileName, options);
end