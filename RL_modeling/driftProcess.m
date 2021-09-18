function [trajectory,hitBound,boundcross] = driftProcess(driftRate,noiseAmplitude,startVal,boundaryVal,nsteps)
noise = noiseAmplitude*randn(1,nsteps);
trajectory = nan(1,nsteps);
trajectory(1) = startVal;
hitBound = 0;
boundcross = nan;
for i=2:nsteps
    trajectory(i) = trajectory(i-1) + driftRate + noise(i);
    if (trajectory(i) > boundaryVal)
        hitBound = 1;
        boundcross = i;
        break;
    end
end
end

