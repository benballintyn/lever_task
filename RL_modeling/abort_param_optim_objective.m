function [] = abort_param_optim_objective(x,abortProcessType,allData)

sessionTypeCount = zeros(1,4);
for i=1:length(allData)
    for j=1:length(allData{1})
        sessType = [num2str(allData{i}{j}.Reward_L / allData{i}{j}.Reward_S) 'xFR' num2str(allData{i}{j}.NumPressRequired_S)];
        switch sessType
            case '2xFR6'
                sessInd = 1;
                FRind = 1;
            case '2xFR12'
                sessInd = 2;
                FRind = 2;
            case '5xFR6'
                sessInd = 3;
                FRind = 1;
            case '5xFR12'
                sessInd = 4;
                FRind = 2;
        end
    end
end

end

