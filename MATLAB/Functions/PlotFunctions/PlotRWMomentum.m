function PlotRWMomentum( rw_momentum )

numSteps = length(rw_momentum(:,1));

figure
hold on
grid on
x = 0:numSteps-1;
tiledlayout(4,1)
nexttile
plot(x, rw_momentum(:,1))
title('RW 1 momentum')
nexttile
plot(x, rw_momentum(:,2))
title('RW 2 momentum')
nexttile
plot(x, rw_momentum(:,3))
title('RW 3 momentum')
nexttile
plot(x, rw_momentum(:,4))
title('RW 4')
hold off

end

