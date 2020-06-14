

timestep = 1;
firstRefTime = 30;
secondRefTime = 90;
velSettlingTolerance = 0.99999;


rmsTotal = ErrorVelRMS( velVec, refVel, 1 )

firstSettleTime = VelSettlingTime( qError, velSettlingTolerance, firstRefTime, firstRefTime + 39, timestep )

secondSettleTime = VelSettlingTime( qError, velSettlingTolerance, secondRefTime, secondRefTime + 45, timestep )

averageSettlingTime = (firstSettleTime + secondSettleTime)/2

settledRms = ErrorVelRMS( velVec, refVel, (secondRefTime + secondSettleTime)/timestep )


function rms = ErrorVelRMS( velVec, refVel, startingIter )

rms = sqrt(mean((velVec(startingIter:end) - refVel .* ones(length(velVec(startingIter:end)),1)).^2));

end


function settlingTime = VelSettlingTime( velVec, velSettlingTolerance, refVel, refTime, endTime, timestep )

isSettled = false;
settlingTime = -1;

for indexIter = 0:1:((endTime-refTime)/timestep)
    
    if ( abs(velVec(refTime/timestep + indexIter) - refVel) <= velSettlingTolerance ) && isSettled == false

        settlingTime = indexIter*timestep;
        isSettled = true;

    elseif ( abs(velVec(refTime/timestep + indexIter) - refVel) > velSettlingTolerance ) && isSettled == true
        
        settlingTime = -1;
        isSettled = false;
        
    end

end

end

