function PlotPropMass( rhoVec, timeVec )

figure
hold on
title('Propellant Consumption')
ylabel('[ \mu g ]')
xlabel('[ s ]')
time = timeVec;
for propIter = 1:length(rhoVec(1,:))
    plot(time, (ones(size(rhoVec(:,propIter))).*0.01 - ...
        rhoVec(:,propIter) .* 0.01^3)*1e6)
end
legend('T_1','T_2','T_3','T_4','T_5','T_6')
hold off

end

