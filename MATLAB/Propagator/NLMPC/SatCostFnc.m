function J = SatCostFnc(Y,U,e,data, t, mjd)

global simConfig

qRef_ECI = simConfig.referenceQuaternion;

J = 0;

Orient_cost = 1e3;

MTQ_cost = 1;
RW_cost = 1;
PROP_cost = 1;

for ephIter = 1:length(Y(:,1))
    
    qError_ECI = QuaternionError(qRef_ECI, Y(ephIter,7:10));
    
    J = J + Orient_cost * qError_ECI(2:4)' * qError_ECI(2:4);
    
end

% for controlIter = 1:length(U(:,1))
%     
%     J = J + U(controlIter,1:3) * ((MTQ_cost .* eye(3,3)) ...
%         * U(controlIter,1:3)');
%     
%     J = J + U(controlIter,4:7) * ((RW_cost .* eye(4,4)) ...
%         * U(controlIter,4:7)');
%     
%     J = J + U(controlIter,8:13) * ((PROP_cost .* eye(6,6)) ...
%         * U(controlIter,8:13)');
%     
% end

end