function [trajectory,hitBound,boundcross,Qs] = drift_valueUpdateProcess(Qcur,Qalt,driftRateCoeff,req,noiseAmplitude,startVal,boundaryVal,nsteps)
noise = noiseAmplitude*randn(1,nsteps);
trajectory = nan(1,nsteps);
trajectory(1) = startVal;
hitBound = 0;
boundcross = nan;
Qs(1) = Qcur;
for i=2:nsteps
    Qs(i) = Qs(i-1) + Qs(i-1)/(req - i + 1);
    driftRate = -driftRateCoeff*(Qs(i) - Qalt);
    trajectory(i) = trajectory(i-1) + driftRate + noise(i);
    if (trajectory(i) > boundaryVal)
        hitBound = 1;
        boundcross = i;
        break;
    end
end
end

