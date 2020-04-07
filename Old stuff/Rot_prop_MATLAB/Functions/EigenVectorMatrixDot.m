function rightEigenMat_dot = EigenVectorMatrixDot( rightEigenMat, leftEigenMat, eigenVec, AMatDot )

rightEigenMat_dot = zeros(3,3);

for kIter = 1:3

    for jIter = 1:3
        
        if kIter == jIter
            
            c = zeros(3,1);
            
            for iIter = 1:3
                
                if iIter ~= kIter
                
                    c = c - c_kj(rightEigenMat, leftEigenMat, eigenVec, AMatDot, ...
                    kIter, iIter) * rightEigenMat_dot(:, kIter);
                
                end
                
            end
            
        else
            
            c = c_kj(rightEigenMat, leftEigenMat, eigenVec, AMatDot, ...
                    kIter, jIter);
        end
        
        rightEigenMat_dot(:, kIter) = rightEigenMat_dot(:, kIter) ...
            + ( c .* rightEigenMat_dot(:, jIter) );
        
    end

end

end

function c_kj = c_kj( rightEigenMat, leftEigenMat, eigenVec, AMatDot, ...
        k_index, j_index)
    
    c_kj = ( leftEigenMat(:, j_index)' * AMatDot * ...
            rightEigenMat(:, k_index) ) / ...
            ( ( eigenVec(k_index) - eigenVec(j_index) ) * ...
            leftEigenMat(:, j_index)' * rightEigenMat(:, j_index) );

end