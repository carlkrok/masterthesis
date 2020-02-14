function [ qSum ] = QuaternionSummation( q1, q2 )

etaSum = q1(1)*q2(1)-q2(2:4)'*q1(2:4);
epsilonSum = q2(1)*q1(2:4) + q1(1)*q2(2:4)+cross(q1(2:4), q2(2:4));

qSum = [etaSum; epsilonSum];

end

