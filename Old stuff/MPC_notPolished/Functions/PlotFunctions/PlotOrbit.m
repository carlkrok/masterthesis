function PlotOrbit( xyz, includeEarth )

figure
hold on
axis equal
grid on
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')
plot3( xyz(:,1), xyz(:,2), xyz(:,3) )
hold off

end

