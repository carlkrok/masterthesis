

timestep = 0.1;
firstRefTime = 30;
secondRefTime = 90;
quatSettlingTolerance = 0.99999;


qError = zeros(length(quaternions), 4);
for quatIter = 1:length(quaternions(:,1))
    qError(quatIter,:) = QuaternionProduct( QuaternionInverse( quatRef(:, quatIter) ), quaternions(quatIter,:)' );
end



rmsTotal = ErrorQuatRMS( qError, 1 )

firstSettleTime = QuatSettlingTime( qError, quatSettlingTolerance, firstRefTime, firstRefTime + 39, timestep )

secondSettleTime = QuatSettlingTime( qError, quatSettlingTolerance, secondRefTime, secondRefTime + 45, timestep )

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

