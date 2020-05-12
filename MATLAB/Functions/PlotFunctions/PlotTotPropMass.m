function PlotTotPropMass( rhoVec, timeVec )

global satelliteConfiguration

if satelliteConfiguration == 1
    propMass = 0.01;
    propVol = 0.01^3;
elseif satelliteConfiguration == 2
    propMass = 0.1;
    propVol = 0.1*0.1*0.01;
end

figure
hold on
title('Propellant Consumption')
ylabel('[ \mu g ]')
xlabel('[ s ]')
time = timeVec;
totMass = zeros(1,length(rhoVec(:,1)));
for propIter = 1:length(rhoVec(:,1))
    totMass(1,propIter) = totMass(1,propIter) + ((ones(size(rhoVec(propIter,:))).*propMass - ...
        rhoVec(propIter,:) .* propVol)*1e6) * ones(length(rhoVec(1,:)),1);
end
plot(time, totMass)
legend('Tot')
hold off

end

