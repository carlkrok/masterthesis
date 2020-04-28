function [ qError ] = QuaternionError( qRef, qCurr )
% qError describes the rotation from the current quaternion (qCurr) to the
% desired quaternion(qRef).

qError = QuaternionProduct( QuaternionInverse( qRef ), qCurr );

end

