function m = MTQ_PD( q_error, omega_error, b_earth, ...
    sat_I_mat, K_p, K_d, maxDipoleMoment, sat_omega, rw_A_mat, rw_momentum )

%t_ref_sat =  - K_p * q_error(2:4) + K_d * omega_error;

t_ref_sat =  - K_p .* (sat_I_mat * q_error(2:4)) ...
     - K_d .* (sat_I_mat * omega_error); ...
     %+ cross( sat_omega, (sat_I_mat * sat_omega + rw_A_mat * rw_momentum ));

m = cross( b_earth, t_ref_sat) / norm( b_earth )^2;

scaling_factor = 1;
for mIter = 1:length(m)
    if abs(m(mIter)) > maxDipoleMoment 
        this_scaling_factor = abs(m(mIter)) / maxDipoleMoment;
        if this_scaling_factor > scaling_factor
            scaling_factor = this_scaling_factor;
        end
    end
end

m = m./scaling_factor;

% res = cross( m, b_earth);
% figure(1)
% hold on
% axis equal
% grid on
% xlabel("X");
% ylabel("Y");
% zlabel("Z");
% plot3([0, b_earth(1)], [0, b_earth(2)], [0, b_earth(3)],'r')
% plot3([0, t_ref_sat(1)], [0, t_ref_sat(2)], [0, t_ref_sat(3)],'g')
% plot3([0, m(1)], [0, m(2)], [0, m(3)],'b')
% plot3([0, res(1)], [0, res(2)], [0, res(3)],'c')
% % plot3([0, b_earth(1)/norm(b_earth)], [0, b_earth(2)/norm(b_earth)], [0, b_earth(3)/norm(b_earth)],'r')
% % plot3([0, t_ref_sat(1)/norm(t_ref_sat)], [0, t_ref_sat(2)/norm(t_ref_sat)], [0, t_ref_sat(3)]/norm(t_ref_sat),'g')
% % plot3([0, m(1)/norm(m)], [0, m(2)/norm(m)], [0, m(3)/norm(m)],'b')
% % plot3([0, res(1)/norm(res)], [0, res(2)/norm(res)], [0, res(3)/norm(res)],'c')
% legend("Earth magnetic field", "t_ref_sat", "m", "t_result_mtq")
% hold off

end

