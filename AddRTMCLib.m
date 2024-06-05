function AddRTMCLib(baseFolder)
if nargin == 0
    baseFolder = fileparts( which( mfilename ) );
end

addpath(baseFolder);

addingfolders = dir( baseFolder );
for folderid = 3 : size( addingfolders, 1 )
    if addingfolders(folderid).isdir
        if addingfolders(folderid).name(1) == '.' || addingfolders(folderid).name(1) == '+' || addingfolders(folderid).name(1) == '@'
            %.gitignore and .git folder
        else
            AddRTMCLib([ baseFolder, filesep, addingfolders(folderid).name])
        end
    end
end
savepath