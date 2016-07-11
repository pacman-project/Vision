% This script loads the workspace from relevant folders.
load([pwd '/Workspace.mat']);
if savedFlag
   return;
end

% We try loading each variable individually.
for varItr = 1:numel(varList)
   varName = varList(varItr).name;
   if ~strcmp(varName, 'ans')
       try
           load([pwd '/workspaces/' datasetName '/level' num2str(levelItr) '/' varName '.mat']);
       end
   end
end