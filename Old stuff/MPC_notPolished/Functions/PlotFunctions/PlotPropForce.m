function PlotPropForce( forceVec )


figure
hold on
grid on
title('Propulsion Force Applied')
ylabel('N')
time = 1:length(forceVec(:,1));
for i = 1:length(forceVec(1,:))
    plot(time, forceVec(:,i))
end
hold off

end
