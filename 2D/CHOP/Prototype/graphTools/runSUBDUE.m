%> Name: runSUBDUE
%>
%> Description: This function runs SUBDUE with the input image and
%> parameters defined in options. The result file is saved under ./
%>
%> @param graphFileName The input graph's path.
%> @param resultFileName The result file's path.
%> @param options Program options.
%> @param currentFolder Path to the workspace folder.
%>
%> @retval []
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 28.11.2013
function [ ] = runSUBDUE( graphFileName, resultFileName, options, currentFolder)

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

    subdueOptions = [subdueOptions ' -nsubs ' num2str(options.subdue.nsubs) ...
                        ' -minsize ' num2str(options.subdue.minSize) ...
                        ' -maxsize ' num2str(options.subdue.maxSize) ...
                        ' -beam ' num2str(options.subdue.beam) ...
                        ' -fileoutput 3 ' ...
                        ' -out ' resultFileName ' ' graphFileName];

    % Form the command and run subdue with specified options
    command = [currentFolder '/subdue' subdueOptions];
    system(command);
end
