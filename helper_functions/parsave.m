function [] = parsave(fname,savevar,varname)
switch varname
    case 'subscores'
        subscores = savevar;
        save(fname,'subscores','-mat')
    case 'run_num'
        run_num = savevar;
        save(fname,'run_num','-mat')
    case 'nS'
        nS = savevar;
        save(fname,'nS','-mat')
    case 'nL'
        nL = savevar;
        save(fname,'nL','-mat')
end

