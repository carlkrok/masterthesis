function PlotOrbitRadiusDeviation( t, xyz, isCircular )

r0 = norm( [xyz(1,1), xyz(1,2), xyz(1,3)] );

numSteps = length(xyz(:,1));

radiusDeviations = zeros(numSteps, 1);
for stepIter = 1:numSteps
    radiusDeviations(stepIter) = norm( [xyz(stepIter,1), ...
        xyz(stepIter,2), xyz(stepIter,3)] );
    if isCircular
        radiusDeviations(stepIter) = radiusDeviations(stepIter) - r0;
    end 
end

figure
hold on
grid on
xlabel('Time [s]')
ylabel('Deviation [m]')
plot(t, radiusDeviations)

end

