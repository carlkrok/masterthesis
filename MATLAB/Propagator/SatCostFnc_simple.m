function J = SatCostFnc_simple(Y,U,e,data,I_mat_body, com_struct, tot_mass, ...
    b_earth_eci)

%global simConfig

%qRef_ECI = simConfig.referenceQuaternion;

J = 0;

Orient_cost = 1e6;

MTQ_cost = 0.01;
RW_cost = 0.1;
PROP_cost = 10000;

for ephIter = 1:length(Y(:,1))
    
    %qError_ECI = QuaternionError(qRef_ECI, Y(ephIter,7:10));
    
    J = J + Orient_cost * ( norm( Y(ephIter,8:10) ) ); ...qError_ECI(2:4)' * qError_ECI(2:4);
    
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