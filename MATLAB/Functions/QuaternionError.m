function [ qError ] = QuaternionError( qRef, qCurr )

qError = QuaternionProduct( QuaternionInverse( qRef ), qCurr );

end

