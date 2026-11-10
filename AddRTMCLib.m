function AddRTMCLib(baseFolder, saveIt)
% AddRTMCLib  Add the library and all its subfolders to the MATLAB path.
%   AddRTMCLib                 repository of this file, path saved (savepath)
%   AddRTMCLib(folder)         given folder, path saved
%   AddRTMCLib(folder, false)  path not saved (used by scripts and tests)
if nargin < 1 || isempty(baseFolder)
    baseFolder = fileparts( which( mfilename ) );
end
if nargin < 2
    saveIt = true;
end

addpath(baseFolder);

% Packages (+TMatrix, +Comparison, lamb/+lamb) and the class folder
% (@MaterialClass) are found through their parent folder, so they are
% skipped below; adding lamb/ makes the Lamb f-k package visible.
addingfolders = dir( baseFolder );
for folderid = 3 : size( addingfolders, 1 )
    if addingfolders(folderid).isdir
        if addingfolders(folderid).name(1) == '.' || addingfolders(folderid).name(1) == '+' || addingfolders(folderid).name(1) == '@'
            %.gitignore and .git folder
        else
            AddRTMCLib([ baseFolder, filesep, addingfolders(folderid).name], false)
        end
    end
end
if saveIt
    savepath
end
