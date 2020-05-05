function PlotPropMass( rhoVec, timeVec )

figure
hold on
title('Propellant Consumption')
ylabel('[ \mu g ]')
xlabel('[ s ]')
time = timeVec;

plot(time, (ones(size(rhoVec(:,1))).*0.01 - rhoVec(:,1) .* 0.01^3)*1e6)
plot(time, (ones(size(rhoVec(:,2))).*0.01 - rhoVec(:,2) .* 0.01^3)*1e6, '--' )
plot(time, (ones(size(rhoVec(:,3))).*0.01 - rhoVec(:,3) .* 0.01^3)*1e6, '--' )
plot(time, (ones(size(rhoVec(:,4))).*0.01 - rhoVec(:,4) .* 0.01^3)*1e6, '--' )
plot(time, (ones(size(rhoVec(:,5))).*0.01 - rhoVec(:,5) .* 0.01^3)*1e6, '--')
plot(time, (ones(size(rhoVec(:,6))).*0.01 - rhoVec(:,6) .* 0.01^3)*1e6, '--' )


legend('T_1','T_2','T_3','T_4','T_5','T_6')
hold off

end

