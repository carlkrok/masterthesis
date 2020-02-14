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
nexttile
plot(x, xAxisRot.*(180/pi))
title('X-Axis')
nexttile
plot(x, yAxisRot.*(180/pi))
title('Y-Axis')
nexttile
plot(x, zAxisRot.*(180/pi))
title('Z-Axis')
hold off

end

