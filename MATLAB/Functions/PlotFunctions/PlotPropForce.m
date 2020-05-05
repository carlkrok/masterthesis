function PlotPropForce( forceVec )


figure
hold on
grid on
title('Propulsion Force Applied')
ylabel('[N]')
xlabel('[s]')
time = forceVec(:,1);
for i = 2:length(forceVec(1,:))
    plot(time, forceVec(:,i))
end
legend('T_1','T_2','T_3','T_4','T_5','T_6')
hold off

end
