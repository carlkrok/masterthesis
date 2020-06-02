function PlotEulerAngles( quaternions, timeVec, refTime, xRef, yRef, zRef )

numSteps = length(quaternions(:,1));

xAxisRot = zeros(numSteps,1);
yAxisRot = zeros(numSteps,1);
zAxisRot = zeros(numSteps,1);

for stepIter = 1:numSteps
    thisQuaternion = quaternions(stepIter,:)';
    thisRotMat = QuaternionToRotMat(thisQuaternion);
    [ xAxisRot(stepIter), yAxisRot(stepIter), zAxisRot(stepIter) ] = ...
        RotMatToEuler( thisRotMat );
end

figure
hold on
grid on
x = timeVec;
tiledlayout(3,1)
ax1 = nexttile;
plot(x, xAxisRot.*(180/pi),refTime, xRef.*(180/pi))
title('Euler angle X-Axis orientation')
ylabel('[deg]')
xlabel('[s]')
ax2 = nexttile;
plot(x, yAxisRot.*(180/pi),refTime, yRef.*(180/pi))
title('Euler angle Y-Axis orientation')
ylabel('[deg]')
xlabel('[s]')
ax3 = nexttile;
plot(x, zAxisRot.*(180/pi),refTime, zRef.*(180/pi))
title('Euler angle Z-Axis orientation')
ylabel('[deg]')
xlabel('[s]')
grid(ax1,'on')
grid(ax2,'on')
grid(ax3,'on')
hold off

end

