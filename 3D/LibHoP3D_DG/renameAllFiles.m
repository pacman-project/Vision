inFolder = 'D:\MissingWarehouse\MissingWarehouse\';
outFolder =  'D:\Warehouse_Clean_Renamed_Vladislav_1\';

is_subset = false;
subsetPercent = 1.0;

list_in = ReadImageNames(inFolder);
lenF = length(list_in);

for i = 1:lenF
    
    fileNameNew = [outFolder, 'D_', num2str(i) ,'.obj'];
    copyfile(list_in{i}, fileNameNew);
end

