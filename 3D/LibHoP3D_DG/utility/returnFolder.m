% this function retirns a folder name from a file name

function [ folderName ] = returnFolder( outFile )

    ll = strfind(outFile, '/');
    lll = strfind(outFile, '\');
    ll = [ll, lll];
    ll = max(ll); % last position
    folderName = outFile(1:ll);
    
end

