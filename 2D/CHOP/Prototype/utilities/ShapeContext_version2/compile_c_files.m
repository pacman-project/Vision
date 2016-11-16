% This script will compile all the C files of the registration methods
cd('functions');
files=dir('*.c');
for i=1:length(files)
    filename=[files(i).name];
    if(length(filename)>10)
        disp(['compiling : ' filename]);
        mex(filename,'-v');
    end
end
cd('..');
