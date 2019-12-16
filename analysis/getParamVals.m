function [paramVals,paramNames,scores,subscores] = getParamVals(basedir,modelDir,subdir)
datadir = [basedir '/' modelDir '/' subdir];
files=dir(datadir);
count=0;
for i=1:length(files)
    if (~(strcmp(files(i).name,'.') || strcmp(files(i).name,'..')) && files(i).isdir)
        count=count+1;
        if (exist([datadir '/' num2str(count) '/score.mat'],'file'))
            p = load([datadir '/' num2str(count) '/params.mat']); p=p.params;
            paramVals(count,:) = p;
            score = load([datadir '/' num2str(count) '/score.mat']); score=score.score;
            scores(count) = score;
            if (exist([datadir '/' num2str(count) '/subscores.mat'],'file'))
                ss = load([datadir '/' num2str(count) '/subscores.mat']); ss=ss.subscores;
                subscores(count,:) =  ss;
            else
                subscores = [];
            end
        else
            count=count-1;
            continue;
        end
    end
end

info = load([datadir '/' modelDir '.mat']); info=info.(modelDir);
paramNames = info.agentParams;
paramNames = [paramNames info.envParams];
end

