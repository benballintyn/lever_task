function [trajectory,hitBound,boundcross] = driftProcess(driftRate,c,nsteps)
noise = c*randn(1,nsteps);
trajectory = nan(1,nsteps);
trajectory(1) = 0;
hitBound = 0;
boundcross = nan;
for i=2:nsteps
    trajectory(i) = trajectory(i-1) + driftRate + noise(i);
    if (trajectory(i) > 1)
        hitBound = 1;
        boundcross = i;
        break;
    end
end
end

