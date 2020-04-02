function PlotQuaternionError( qRef, quaternions)

qError = zeros(length(quaternions), 4);
for quatIter = 1:length(quaternions)
    qError(quatIter,:) = QuaternionProduct( QuaternionInverse( qRef ), quaternions(quatIter,:)' );
end

t = 1:length(quaternions);

figure
hold on
tiledlayout(4,1)
ax1 = nexttile;
grid on
plot(t, qError(:,1))
title("Quat error (1)")
ax2 = nexttile;
grid on
plot(t, qError(:,2))
title("Quat error (2)")
ax3 = nexttile;
grid on
plot(t, qError(:,3))
title("Quat error (3)")
ax4 = nexttile;
grid on
plot(t, qError(:,4))
title("Quat error (4)")
grid(ax1,'on')
grid(ax2,'on')
grid(ax3,'on')
grid(ax4,'on')
hold off


end

