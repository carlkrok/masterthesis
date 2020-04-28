function PlotPropRho( rhoVec )


figure
hold on
title('Propulsion module densities')
ylabel('kg/m^3')
time = 1:length(rhoVec(:,1));
for i = 1:length(rhoVec(1,:))
    plot(time, rhoVec(:,i))
end
hold off

end

