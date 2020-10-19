function [vR,vS,choices,switchPoint] = theoreticalPredictSwitch(Lr,Sr,Nl,Ns,itiTime,uDecayConstant)
startNl = Nl;
nT = 200;
choices = zeros(4,nT);
vR = zeros(4,nT);
vS = zeros(4,nT);
noSwitch = true;
for i=1:nT
    vR(1,i) = Lr/Nl;
    vS(1,i) = Sr/Ns;
    if (vR(1,i) >= vS(1,i))
        Nl = Nl + 1;
        choices(1,i) = 2;
    else
        if (noSwitch)
            switchPoint(1) = i;
            noSwitch = false;
        end
        choices(1,i) = 1;
    end
end

noSwitch = true;
Nl = startNl;
%itiTime = 5;
for i=1:nT
    vR(2,i) = Lr/(itiTime + Nl);
    vS(2,i) = Sr/(itiTime + Ns);
    if (vR(2,i) >= vS(2,i))
        Nl = Nl + 1;
        choices(2,i) = 2;
    else
        if (noSwitch)
            switchPoint(2) = i;
            noSwitch = false;
        end
        choices(2,i) = 1;
    end
end

noSwitch = true;
Nl = startNl;
%uDecayConstant = 2000;
u = @(x) exp(-x/uDecayConstant);
totalR = 0;
for i=1:nT
    vR(3,i) = (u(totalR)*Lr)/Nl;
    vS(3,i) = (u(totalR)*Sr)/Ns;
    if (vR(3,i) >= vS(3,i))
        Nl = Nl + 1;
        totalR = totalR + Lr;
        choices(3,i) = 2;
    else
        if (noSwitch)
            switchPoint(3) = i;
            noSwitch = false;
        end
        totalR = totalR + Sr;
        choices(3,i) = 1;
    end
end

noSwitch = true;
Nl = startNl;
totalR = 0;
for i=1:nT
    vR(4,i) = (u(totalR)*Lr)/(itiTime + Nl);
    vS(4,i) = (u(totalR)*Sr)/(itiTime + Ns);
    if (vR(4,i) >= vS(4,i))
        Nl = Nl + 1;
        totalR = totalR + Lr;
        choices(4,i) = 2;
    else
        if (noSwitch)
            switchPoint(4) = i;
            noSwitch = false;
        end
        totalR = totalR + Sr;
        choices(4,i) = 1;
    end
end
end

