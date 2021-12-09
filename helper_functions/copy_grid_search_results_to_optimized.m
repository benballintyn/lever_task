function [] = copy_grid_search_results_to_optimized(datadir)
OPTIMIZED_DIR = [datadir '/optimized/'];
GRID_SEARCH_DIR = [datadir '/grid_search/'];
if (~exist(OPTIMIZED_DIR,'dir'))
    error(['optimized directory does not exist inside ' datadir])
elseif (~exist(GRID_SEARCH_DIR,'dir'))
    warning('No grid search results to transfer')
    return;
else
    directories = dir(OPTIMIZED_DIR);
    directories = directories(~ismember({directories.name},{'.','..'}));
    directories = directories([directories.isdir]);
    for i=1:length(directories)
        nums(i) = str2num(directories(i).name);
    end
    [~,sortedOrder] = sort(nums);
    directories = directories(sortedOrder);
    
    run_num = load([OPTIMIZED_DIR '/run_num.mat']); run_num=run_num.run_num;
    
    if (run_num ~= max(nums))
        error('run_num and maximum folder number do not agree. Try cleanDirNums.')
    end
    
    gridDirs = dir(GRID_SEARCH_DIR);
    gridDirs = gridDirs(~ismember({gridDirs.name},{'.','..'}));
    gridDirs = gridDirs([gridDirs.isdir]);
    for i=1:length(gridDirs)
        nums(i) = str2num(gridDirs(i).name);
    end
    [~,sortedOrder] = sort(nums);
    gridDirs = gridDirs(sortedOrder);
    
    for i=1:length(gridDirs)
        copyfile([GRID_SEARCH_DIR gridDirs(i).name],[OPTIMIZED_DIR num2str(run_num + i)])
        disp(['copied grid_search/' gridDirs(i).name ' to optimized/' num2str(run_num + i)])
    end
    
    run_num = run_num + length(gridDirs);
    save([OPTIMIZED_DIR '/run_num.mat'],'run_num','-mat')
end

