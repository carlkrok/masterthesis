

timestep = 1;
firstRefTime = 300;
secondRefTime = 900;
quatSettlingTolerance = 0.99999;

quatRef = [simConfig.firstReferenceQuaternion * ones(1,(simConfig.secondReferenceQuaternionTime/timestep_controller)+1), ...
    simConfig.secondReferenceQuaternion * ones(1,(simConfig.thirdReferenceQuaternionTime-simConfig.secondReferenceQuaternionTime)/timestep_controller), ...
    simConfig.thirdReferenceQuaternion * ones(1,(duration-simConfig.thirdReferenceQuaternionTime)/timestep_controller)];


qError = zeros(length(quaternions), 4);
for quatIter = 1:length(quaternions(:,1))
    qError(quatIter,:) = QuaternionProduct( QuaternionInverse( quatRef(:, quatIter) ), quaternions(quatIter,:)' );
end



rmsTotal = ErrorQuatRMS( qError, 1 )

firstSettleTime = QuatSettlingTime( qError, quatSettlingTolerance, firstRefTime, firstRefTime + 300, timestep )

secondSettleTime = QuatSettlingTime( qError, quatSettlingTolerance, secondRefTime, secondRefTime + 450, timestep )

averageSettlingTime = (firstSettleTime + secondSettleTime)/2

settledRms = ErrorQuatRMS( qError, (secondRefTime + secondSettleTime)/timestep )



function rms = ErrorQuatRMS( errorQuat, startingIter )

rms = sqrt(mean((errorQuat(startingIter:end,1) - ones(length(errorQuat(startingIter:end,1)),1)).^2));

end


function settlingTime = QuatSettlingTime( errorQuat, quatSettlingTolerance, refTime, endTime, timestep )

isSettled = false;
settlingTime = -1;

for indexIter = 0:1:((endTime-refTime)/timestep)
    
    if ( errorQuat(refTime/timestep + indexIter, 1) >= quatSettlingTolerance ) && isSettled == false

        settlingTime = indexIter*timestep;
        isSettled = true;

    elseif ( errorQuat(refTime/timestep + indexIter, 1) < quatSettlingTolerance ) && isSettled == true
        
        settlingTime = -1;
        isSettled = false;
        
    end

end

end

