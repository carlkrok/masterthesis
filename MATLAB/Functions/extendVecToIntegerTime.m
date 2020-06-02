function intVec = extendVecToIntegerTime(contVec)

intVec = [contVec(1); contVec(1)];

for contIter = 2:length(contVec)

    intVec = [intVec; contVec(contIter); contVec(contIter)];
    
end

end