

timestep = 0.1;
firstRefTime = 30;
secondRefTime = 180;
velSettlingTolerance = 0.00005;

omegaRef = [simConfig.firstReferenceOmega(1) * ones((simConfig.secondReferenceOmegaTime/timestep), 1); ...
    simConfig.secondReferenceOmega(1) * ones((simConfig.thirdReferenceOmegaTime-simConfig.secondReferenceOmegaTime)/timestep, 1); ...
    simConfig.thirdReferenceOmega(1) * ones((duration-simConfig.thirdReferenceOmegaTime)/timestep+1, 1)];


rmsX = ErrorVelRMS( omega(:,1), omegaRef, 1 )
rmsY = ErrorVelRMS( omega(:,2), 0.*omegaRef, 1 )
rmsZ = ErrorVelRMS( omega(:,3), 0.*omegaRef, 1 )

rmsTotal = rmsX + rmsY + rmsZ

firstSettleTime = VelSettlingTime( omega(:,1), omegaRef, velSettlingTolerance, firstRefTime, firstRefTime + 100, timestep )

secondSettleTime = VelSettlingTime( omega(:,1), omegaRef, velSettlingTolerance, secondRefTime, secondRefTime + 50, timestep )

averageSettlingTime = (firstSettleTime + secondSettleTime)/2

settledRmsX = ErrorVelRMS( omega(:,1), omegaRef, (secondRefTime + secondSettleTime)/timestep )
settledRmsY = ErrorVelRMS( omega(:,2), omegaRef, (secondRefTime + secondSettleTime)/timestep )
settledRmsZ = ErrorVelRMS( omega(:,3), omegaRef, (secondRefTime + secondSettleTime)/timestep )
% settledRmsX = ErrorVelRMS( omega(:,1), omegaRef, (240)/timestep )
% settledRmsY = ErrorVelRMS( omega(:,2), omegaRef, (240)/timestep )
% settledRmsZ = ErrorVelRMS( omega(:,3), omegaRef, (240)/timestep )

settledRmsTotal = settledRmsX + settledRmsY + settledRmsZ


function rms = ErrorVelRMS( velVec, refVel, startingIter )

rms = sqrt(mean((velVec(startingIter:end) - refVel(startingIter:end)).^2));

end


function settlingTime = VelSettlingTime( velVec, refVec, velSettlingTolerance, refTime, endTime, timestep )

isSettled = false;
settlingTime = -1;

for indexIter = 0:1:((endTime-refTime)/timestep)
    
    if ( abs(velVec(refTime/timestep + indexIter) - refVec(refTime/timestep + indexIter)) <= velSettlingTolerance ) && isSettled == false

        settlingTime = indexIter*timestep;
        isSettled = true;

    elseif ( abs(velVec(refTime/timestep + indexIter) - refVec(refTime/timestep + indexIter)) > velSettlingTolerance ) && isSettled == true
        
        settlingTime = -1;
        isSettled = false;
        
    end

end

end

