function [ refFrame ] = getReferenceFrame( target, v )
    
zAxisRef = target / norm(target);
yAxisRef = cross( target, v ) / norm(cross( target, v ));
xAxisRef = cross( yAxisRef, zAxisRef ) / norm(cross( yAxisRef, zAxisRef ));

refFrame = [xAxisRef'; yAxisRef'; zAxisRef'];

if cross( target, v ) == 0
    warning("Reference Rotation Matrix not obtainable");
    refFrame = eye(3,3);

end

