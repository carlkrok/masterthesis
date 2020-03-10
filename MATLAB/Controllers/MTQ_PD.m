function [m, t_ref_sat] = MTQ_PD( q_error, omega_error, b_earth, ...
    sat_I_mat, K_p, K_d, maxDipoleMoment, sat_omega, rw_A_mat, rw_momentum )

t_ref_sat =  - K_p * sat_I_mat * q_error(2:4) ...
     - K_d * sat_I_mat * omega_error ...
     + cross( sat_omega, (sat_I_mat * sat_omega + rw_A_mat * rw_momentum ));

 
m = cross( b_earth, t_ref_sat) / norm( b_earth )^2;

% t_mtq_lb = -maxDipoleMoment .* ones(3,1);
% t_mtq_ub = maxDipoleMoment .* ones(3,1);
% 
% m = lsqlin(S_crossProdMat(-b_earth),t_ref_sat,[],[],[],[], ...
%     t_mtq_lb, ... % Lower bound
%     t_mtq_ub); % Upper bound



% scaling_factor = 1;
% for mIter = 1:length(m)
%     if abs(m(mIter)) > maxDipoleMoment 
%         this_scaling_factor = abs(m(mIter)) / maxDipoleMoment;
%         if this_scaling_factor > scaling_factor
%             scaling_factor = this_scaling_factor;
%         end
%     end
% end
% 
% m = m./scaling_factor;

for mIter = 1:length(m)
    if m(mIter) > maxDipoleMoment 
        m(mIter) = maxDipoleMoment;
    elseif m(mIter) < -maxDipoleMoment 
        m(mIter) = -maxDipoleMoment;
    end
end

end

