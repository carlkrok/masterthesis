function [ qInv ] = QuaternionInverse( q )

qLength = norm( q );

qInv = [q(1)/qLength; ...
        -q(2)/qLength; ...
        -q(3)/qLength; ...
        -q(4)/qLength];

end

