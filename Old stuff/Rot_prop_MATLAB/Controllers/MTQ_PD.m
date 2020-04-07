function [m, t_ref_sat, t_can_sat] = MTQ_PD( q_error, omega_error, ...
    q_error_integrated, b_earth, ...
    sat_I_mat, K_p, K_d, K_i, maxDipoleMoment, t )

global plotData

t_ref_sat =  - K_p * q_error(2:4) ...
     - K_d * omega_error; ...
     %- K_i * q_error_integrated(2:4);
 
 
 
 b_unit = b_earth / norm(b_earth);
 %t_ref_unit = t_ref_sat / norm( t_ref_sat );
 
 angleDiff = acos(dot(b_unit, [0;0;1]));
 angleAxis = cross( b_unit, [0;0;1] );
 angleAxis = angleAxis / norm(angleAxis);
 
 RotMat_bToTemp = AngleAxisToRotMat( angleDiff, angleAxis );
 
 t_ref_tempFrame = RotMat_bToTemp * t_ref_sat;
 
 t_can_tempframe = t_ref_tempFrame .* [1;1;0];
 
 t_can_sat = RotMat_bToTemp' * t_can_tempframe;
 
m = cross( b_earth, t_can_sat) / norm( b_earth )^2;


% t_mtq_lb = -maxDipoleMoment .* ones(3,1);
% t_mtq_ub = maxDipoleMoment .* ones(3,1);

% m = lsqlin(S_crossProdMat(-b_earth),t_can_sat,[],[],[],[], ...
%     t_mtq_lb, ... % Lower bound
%     t_mtq_ub); % Upper bound

% if t>1e3
%     figure(1)
%     hold on
%     axis equal
%     grid on
%     plot3([0 b_earth(1)/norm(b_earth)], [0 b_earth(2)/norm(b_earth)], [0 b_earth(3)/norm(b_earth)], 'r')
%     plot3([0 t_ref_sat(1)/norm(t_ref_sat)], [0 t_ref_sat(2)/norm(t_ref_sat)], [0 t_ref_sat(3)/norm(t_ref_sat)], 'g')
%     plot3([0 t_can_sat(1)/norm(t_can_sat)], [0 t_can_sat(2)/norm(t_can_sat)], [0 t_can_sat(3)/norm(t_can_sat)], 'b')
%     hold off
%     disp('here')
% end

scaling_factor = 1;
for mIter = 1:length(m)
    if abs(m(mIter)) > maxDipoleMoment 
        this_scaling_factor = abs(m(mIter)) / maxDipoleMoment;
        if this_scaling_factor > scaling_factor
            scaling_factor = this_scaling_factor;
        end
    end
end
% 
% m = m./scaling_factor;

% for mIter = 1:length(m)
%     if m(mIter) > maxDipoleMoment 
%         m(mIter) = maxDipoleMoment;
%     elseif m(mIter) < -maxDipoleMoment 
%         m(mIter) = -maxDipoleMoment;
%     end
% end



end

