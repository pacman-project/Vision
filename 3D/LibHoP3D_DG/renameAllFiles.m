inFolder = 'D:\Warehouse_Clean_Renamed_Vladislav_Rusen_Final\';
outFolder =  'D:\OrderedModels\';

is_subset = false;
subsetPercent = 1.0;

list_in = ReadImageNames(inFolder);
lenF = length(list_in);

% removedModels = [22 44 52 53 54 65 67 69 72 73 55 48 64 82 307 83 94 116 130 153 179 198 212 234 251 238 300 282 278 309 366 374 383 367 334 320 331 326 336 321 328 405 393 433 421 428 410 430 417];
categoryNames = cell(lenF,1);
for i = 1:lenF
%    if ~ismember(i, removedModels)
        fileString = list_in{i};
        idx = strfind(fileString, '\');
        startIdx = idx(end-1);
        endIdx = idx(end);
        categoryName = fileString(startIdx+1:endIdx-1);
%         if ~exist([outFolder '\' categoryName], 'dir')
%             mkdir([outFolder '\' categoryName]);
%         end
        categoryNames{i} = categoryName;
        fileNameNew = [[outFolder '\'], 'D_', num2str(i) ,'.obj'];
        copyfile(list_in{i}, fileNameNew);
 %   end
end

