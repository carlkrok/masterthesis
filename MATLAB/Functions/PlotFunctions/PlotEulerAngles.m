function PlotEulerAngles( quaternions )

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
x = 0:numSteps-1;
tiledlayout(3,1)
ax1 = nexttile;
plot(x, xAxisRot.*(180/pi))
title('Euler angle X-Axis orientation')
ax2 = nexttile;
plot(x, yAxisRot.*(180/pi))
title('Euler angle Y-Axis orientation')
ax3 = nexttile;
plot(x, zAxisRot.*(180/pi))
title('Euler angle Z-Axis orientation')
grid(ax1,'on')
grid(ax2,'on')
grid(ax3,'on')
hold off

end

