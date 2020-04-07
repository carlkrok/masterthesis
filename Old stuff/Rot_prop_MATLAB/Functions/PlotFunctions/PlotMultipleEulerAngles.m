function PlotMultipleEulerAngles( eph1, eph2, eph3 )

numSteps = length(eph1(:,1));

xAxisRot_1 = zeros(numSteps,1);
yAxisRot_1 = zeros(numSteps,1);
zAxisRot_1 = zeros(numSteps,1);

for stepIter = 1:numSteps
    thisQuaternion = eph1(stepIter,8:11)';
    thisRotMat = QuaternionToRotMat(thisQuaternion);
    [ xAxisRot_1(stepIter), yAxisRot_1(stepIter), zAxisRot_1(stepIter) ] = ...
        RotMatToEuler( thisRotMat );
end

xAxisRot_2 = zeros(numSteps,1);
yAxisRot_2 = zeros(numSteps,1);
zAxisRot_2 = zeros(numSteps,1);

for stepIter = 1:numSteps
    thisQuaternion = eph2(stepIter,8:11)';
    thisRotMat = QuaternionToRotMat(thisQuaternion);
    [ xAxisRot_2(stepIter), yAxisRot_2(stepIter), zAxisRot_2(stepIter) ] = ...
        RotMatToEuler( thisRotMat );
end

xAxisRot_3 = zeros(numSteps,1);
yAxisRot_3 = zeros(numSteps,1);
zAxisRot_3 = zeros(numSteps,1);

for stepIter = 1:numSteps
    thisQuaternion = eph3(stepIter,8:11)';
    thisRotMat = QuaternionToRotMat(thisQuaternion);
    [ xAxisRot_3(stepIter), yAxisRot_3(stepIter), zAxisRot_3(stepIter) ] = ...
        RotMatToEuler( thisRotMat );
end

figure
hold on
grid on
x = 0:numSteps-1;
tiledlayout(3,1)
nexttile
plot(x', xAxisRot_1.*(180/pi),x', xAxisRot_2.*(180/pi),x', xAxisRot_3.*(180/pi))
title('Euler angle X-Axis orientation')
legend('Free motion', 'Perturbed motion', 'Actuated motion')
nexttile
plot(x', yAxisRot_1.*(180/pi),x', yAxisRot_2.*(180/pi),x', yAxisRot_3.*(180/pi))
title('Euler angle Y-Axis orientation')
legend('Free motion', 'Perturbed motion', 'Actuated motion')
nexttile
plot(x', zAxisRot_1.*(180/pi),x', zAxisRot_2.*(180/pi),x', zAxisRot_3.*(180/pi))
title('Euler angle Z-Axis orientation')
legend('Free motion', 'Perturbed motion', 'Actuated motion')
hold off


% xAxisRot = zeros(numSteps,1);
% yAxisRot = zeros(numSteps,1);
% zAxisRot = zeros(numSteps,1);
% 
% figure
% hold on
% grid on
% x = 0:numSteps-1;
% tiledlayout(3,1)
% nexttile
% for ephIter = 1:3
%     thisEph = allEphs(ephIter, :, :);
%     for stepIter = 1:numSteps
%         thisQuaternion = [thisEph(1,stepIter,8), thisEph(1,stepIter,9), ...
%             thisEph(1,stepIter,10), thisEph(1,stepIter,11)]';
%         thisRotMat = QuaternionToRotMat(thisQuaternion);
%         [ xAxisRot(stepIter), yAxisRot(stepIter), zAxisRot(stepIter) ] = ...
%             RotMatToEuler( thisRotMat );
%     end
%     plot(x, xAxisRot.*(180/pi))
% end
% title('Euler angle X-Axis orientation')
% nexttile
% for ephIter = 1:3
%     thisEph = allEphs(ephIter, :, :);
%     for stepIter = 1:numSteps
%         thisQuaternion = [thisEph(1,stepIter,8), thisEph(1,stepIter,9), ...
%             thisEph(1,stepIter,10), thisEph(1,stepIter,11)]';
%         thisRotMat = QuaternionToRotMat(thisQuaternion);
%         [ xAxisRot(stepIter), yAxisRot(stepIter), zAxisRot(stepIter) ] = ...
%             RotMatToEuler( thisRotMat );
%     end
%     plot(x, yAxisRot.*(180/pi))
% end
% title('Euler angle Y-Axis orientation')
% nexttile
% for ephIter = 1:3
%     thisEph = allEphs(ephIter, :, :);
%     for stepIter = 1:numSteps
%         thisQuaternion = [thisEph(1,stepIter,8), thisEph(1,stepIter,9), ...
%             thisEph(1,stepIter,10), thisEph(1,stepIter,11)]';
%         thisRotMat = QuaternionToRotMat(thisQuaternion);
%         [ xAxisRot(stepIter), yAxisRot(stepIter), zAxisRot(stepIter) ] = ...
%             RotMatToEuler( thisRotMat );
%     end
%     plot(x, zAxisRot.*(180/pi))
% end
% title('Euler angle Z-Axis orientation')
% hold off

end

