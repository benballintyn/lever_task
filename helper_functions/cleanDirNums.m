function [] = cleanDirNums(basedir)
DELETED_COUNT=0;
MOVED_COUNT=0;
try
    run_num = load([basedir '/run_num.mat']); run_num=run_num.run_num;
catch e
    error('Could not load run_num.mat')
end

dirContents = dir(basedir);
dirContents = dirContents(~ismember({dirContents.name},{'..','.'}));
dirContents = dirContents([dirContents.isdir]);
for i=1:length(dirContents)
    folderNums(i) = str2num(dirContents(i).name);
end
[~,folderOrder] = sort(folderNums);
dirContents = dirContents(folderOrder);

for i=1:length(dirContents)
    curdir = [dirContents(i).folder '/' dirContents(i).name];
    curfiles = dir(curdir);
    curfiles = curfiles(~ismember({curfiles.name},{'..','.'}));
    if (isempty(curfiles))
        rmdir(curdir)
        DELETED_COUNT = DELETED_COUNT + 1;
    end
end

clear folderNums
dirContents = dir(basedir);
dirContents = dirContents(~ismember({dirContents.name},{'..','.'}));
dirContents = dirContents([dirContents.isdir]);
for i=1:length(dirContents)
    folderNums(i) = str2num(dirContents(i).name);
end
[~,folderOrder] = sort(folderNums);
dirContents = dirContents(folderOrder);
folderNums = folderNums(folderOrder);
if (folderNums(1) ~= 1)
    warning('First folder in the directory is not 1')
end
for i=2:length(dirContents)
    if (folderNums(i) < folderNums(i-1))
        error('Folder numbers are not sorted.')
    elseif (folderNums(i) < 1)
        error('Folder numbers should be positive integers')
    elseif (folderNums(i) > (folderNums(i-1) + 1))
        warning(['Folder ' num2str(dirContents(i-1).name) ' followed by folder ' num2str(dirContents(i).name)])
        if (str2num(dirContents(i-1).name) ~= folderNums(i-1))
            error('Previous directory name and folderNum do not match')
        elseif (str2num(dirContents(i).name) ~= folderNums(i))
            error('Current directory name and folderNum do not match')
        end
        newFolderNum = str2num(dirContents(i-1).name) + 1;
        curdir = [dirContents(i).folder '/' dirContents(i).name];
        newdir = [dirContents(i).folder '/' num2str(newFolderNum)];
        movefile(curdir,newdir)
        folderNums(i) = newFolderNum;
        dirContents(i).name = num2str(newFolderNum);
        MOVED_COUNT = MOVED_COUNT + 1;
    end    
end
old_run_num = run_num;
run_num = str2num(dirContents(end).name);
save([basedir '/run_num.mat'],'run_num','-mat')

disp(sprintf('%u directories deleted',DELETED_COUNT))
disp(sprintf('%u directories moved',MOVED_COUNT))
if (old_run_num == run_num)
    disp('run_num unchanged')
else
    disp(sprintf('run_num changed from %1$u to %2$u',old_run_num,run_num))
end
end

