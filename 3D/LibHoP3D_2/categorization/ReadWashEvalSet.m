% this file reads an evaluation set for the Washington dataset

s={};

path = 'D:\Input Data\Washington\Wash-rgbd-datasetCleaned_Reserv1';
subfolders = list_of_subfolders(path);

testSet = zeros(51, 10);
fid = fopen('WashingtonTestSet.txt');
tline = fgetl(fid);
i = 1;
j = 1;
while ischar(tline)
   ll = strfind(tline, '_');
   if ~isempty(ll)
       substr = tline(ll(end)+1:end);
       testSet(i, j) = str2num(substr);
       i = i+1;
   end
   if i == 52
       j = j + 1;
       i = 1;
   end
   tline = fgetl(fid);
end

a = 2;