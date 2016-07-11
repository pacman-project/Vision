% This script saves the workspace into relevant folders.
% Start with the easy way.
% savedFlag = true;
% try
%     save([pwd '/Workspace.mat'], '-v7');
%     copyfile([pwd '/Workspace.mat'], [pwd '/workspaces/' options.datasetName '/level' num2str(levelItr) '.mat']);
% catch
%     savedFlag = false;
% end

% if savedFlag
%    return;
% end

% Secondly, we try saving each variable individually.
savedFlag = false; % TODO REMOVE.
varList = whos;
if exist([pwd '/workspaces/' options.datasetName '/level' num2str(levelItr)], 'dir')
   rmdir([pwd '/workspaces/' options.datasetName '/level' num2str(levelItr)], 's');
end
mkdir([pwd '/workspaces/' options.datasetName '/level' num2str(levelItr)]);
save([pwd '/Workspace.mat'], 'varList', 'savedFlag','datasetName','levelItr');
copyfile([pwd '/Workspace.mat'], [pwd '/workspaces/' options.datasetName '/level' num2str(levelItr) '.mat']);
for varItr = 1:numel(varList)
   varName = varList(varItr).name;
   if ~strcmp(varName, 'ans')
       try
           save([pwd '/workspaces/' options.datasetName '/level' num2str(levelItr) '/' varName '.mat'], varName, '-v7');
       catch
           display(['Saving variable ' varName ' in latest format. This may take a while.']);
           save([pwd '/workspaces/' options.datasetName '/level' num2str(levelItr) '/' varName '.mat'], varName, '-v7.3');
       end
   end
end