function PlotQuaternionError( qRef, quaternions)

qError = zeros(length(quaternions), 4);
for quatIter = 1:length(quaternions)
    qError(quatIter,:) = QuaternionProduct( QuaternionInverse( qRef ), quaternions(quatIter,:)' );
end

t = 1:length(quaternions);

figure
hold on
tiledlayout(4,1)
nexttile
grid on
plot(t, qError(:,1))
title("Quat error (1)")
nexttile
grid on
plot(t, qError(:,2))
title("Quat error (2)")
nexttile
grid on
plot(t, qError(:,3))
title("Quat error (3)")
nexttile
grid on
plot(t, qError(:,4))
title("Quat error (4)")
hold off


end

