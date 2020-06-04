function PlotQuaternionError( qRefVec, quaternions, timeVec)

qError = zeros(length(quaternions), 4);
for quatIter = 1:length(quaternions(:,1))
    qError(quatIter,:) = QuaternionProduct( QuaternionInverse( qRefVec(:, quatIter) ), quaternions(quatIter,:)' );
end

t = timeVec;

figure
hold on
tiledlayout(4,1)
ax1 = nexttile;
grid on
plot(t, qError(:,1))
xlabel('[s]')
title("Quaternion Error (1)")
ax2 = nexttile;
grid on
plot(t, qError(:,2))
xlabel('[s]')
title("Quaternion Error (2)")
ax3 = nexttile;
grid on
plot(t, qError(:,3))
xlabel('[s]')
title("Quaternion Error (3)")
ax4 = nexttile;
grid on
plot(t, qError(:,4))
xlabel('[s]')
title("Quaternion Error (4)")
grid(ax1,'on')
grid(ax2,'on')
grid(ax3,'on')
grid(ax4,'on')
hold off


end

