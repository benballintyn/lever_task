function [hitBound,boundcross] = logisticAbortProcess(Qcur,Qalt,req,logisticParams,ansFunc,useFullANS,ansSigma,noAbortANS)
hitBound = 0;
boundcross = nan;
f = @(Qcur,Qalt,temp,offset) 1 - 1/(1 + exp(-(Qcur - Qalt + offset)/temp));
pAbort = f(Qcur,Qalt,logisticParams.temp,logisticParams.offset);
r = rand;
if (r <= pAbort)
    hitBound = 1;
    boundcross = 1;
    return
else
    Qs(1) = Qcur;
    for i=2:req
        if (useFullANS & ~noAbortANS)
            Qs(i) = Qs(i-1) + Qs(i-1)/(req - ansFunc(i,ansSigma) + 1);
        else
            Qs(i) = Qs(i-1) + Qs(i-1)/(req - i + 1);
        end
        pAbort = f(Qs(i),Qalt,logisticParams.temp,logisticParams.offset);
        r = rand;
        if (r <= pAbort)
            hitBound = 1;
            boundcross = i;
            return;
        end
    end
end
end

