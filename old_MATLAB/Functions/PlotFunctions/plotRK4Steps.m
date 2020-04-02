function plotRK4Steps( prevTimeStep, thisTimeStepLength, stepTimeIter, ...
    Y, k1, k2, k3, k4, Yfinal )

%     figure(1)
%     hold on
%     grid on
%     plot([prevTimeStep, prevTimeStep + thisTimeStepLength], [norm(Y(stepTimeIter-1,11:13)), norm(Y(stepTimeIter-1,11:13) + k1(11:13))])
%     plot([prevTimeStep, prevTimeStep + thisTimeStepLength], [norm(Y(stepTimeIter-1,11:13)), norm(Y(stepTimeIter-1,11:13) + k2(11:13))])
%     plot([prevTimeStep, prevTimeStep + thisTimeStepLength], [norm(Y(stepTimeIter-1,11:13)), norm(Y(stepTimeIter-1,11:13) + k3(11:13))])
%     plot([prevTimeStep, prevTimeStep + thisTimeStepLength], [norm(Y(stepTimeIter-1,11:13)), norm(Y(stepTimeIter-1,11:13) + k4(11:13))])
%     plot([prevTimeStep, prevTimeStep + thisTimeStepLength], [norm(Y(stepTimeIter-1,11:13)), norm(Yfinal(11:13))])
%     legend("k1", "k2", "k3", "k4","Yfinal")
%     hold off

figure(1)
hold on
grid on
title("RK4 Steps in \omega^{body}_{sat}(1)")
plot([prevTimeStep, prevTimeStep + thisTimeStepLength], [Y(stepTimeIter-1,11), Y(stepTimeIter-1,11) + k1(11)])
plot([prevTimeStep, prevTimeStep + thisTimeStepLength], [Y(stepTimeIter-1,11), Y(stepTimeIter-1,11) + k2(11)])
plot([prevTimeStep, prevTimeStep + thisTimeStepLength], [Y(stepTimeIter-1,11), Y(stepTimeIter-1,11) + k3(11)])
plot([prevTimeStep, prevTimeStep + thisTimeStepLength], [Y(stepTimeIter-1,11), Y(stepTimeIter-1,11) + k4(11)])
plot([prevTimeStep, prevTimeStep + thisTimeStepLength], [Y(stepTimeIter-1,11), Y(stepTimeIter-1,11) + Yfinal(11)])
legend("k1", "k2", "k3", "k4","Yfinal")
hold off

figure(2)
hold on
grid on
title("RK4 Steps in \omega^{body}_{sat}(2)")
plot([prevTimeStep, prevTimeStep + thisTimeStepLength], [Y(stepTimeIter-1,12), Y(stepTimeIter-1,12) + k1(12)])
plot([prevTimeStep, prevTimeStep + thisTimeStepLength], [Y(stepTimeIter-1,12), Y(stepTimeIter-1,12) + k2(12)])
plot([prevTimeStep, prevTimeStep + thisTimeStepLength], [Y(stepTimeIter-1,12), Y(stepTimeIter-1,12) + k3(12)])
plot([prevTimeStep, prevTimeStep + thisTimeStepLength], [Y(stepTimeIter-1,12), Y(stepTimeIter-1,12) + k4(12)])
plot([prevTimeStep, prevTimeStep + thisTimeStepLength], [Y(stepTimeIter-1,12), Y(stepTimeIter-1,12) + Yfinal(12)])
legend("k1", "k2", "k3", "k4","Yfinal")
hold off

figure(3)
hold on
grid on
title("RK4 Steps in \omega^{body}_{sat}(3)")
plot([prevTimeStep, prevTimeStep + thisTimeStepLength], [Y(stepTimeIter-1,13), Y(stepTimeIter-1,13) + k1(13)])
plot([prevTimeStep, prevTimeStep + thisTimeStepLength], [Y(stepTimeIter-1,13), Y(stepTimeIter-1,13) + k2(13)])
plot([prevTimeStep, prevTimeStep + thisTimeStepLength], [Y(stepTimeIter-1,13), Y(stepTimeIter-1,13) + k3(13)])
plot([prevTimeStep, prevTimeStep + thisTimeStepLength], [Y(stepTimeIter-1,13), Y(stepTimeIter-1,13) + k4(13)])
plot([prevTimeStep, prevTimeStep + thisTimeStepLength], [Y(stepTimeIter-1,13), Y(stepTimeIter-1,13) + Yfinal(13)])
legend("k1", "k2", "k3", "k4","Yfinal")
hold off
    
end

