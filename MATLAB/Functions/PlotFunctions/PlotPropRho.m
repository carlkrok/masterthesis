function PlotPropRho( rhoVec, timeVec )


figure
hold on
title('Propellant Densities')
ylabel('kg/m^3')
time = timeVec;
for i = 1:length(rhoVec(1,:))
    plot(time, rhoVec(:,i))
end
legend('$T_1$','$T_2$','$T_3$','$T_4$','$T_5$','$T_6$')
hold off

end

