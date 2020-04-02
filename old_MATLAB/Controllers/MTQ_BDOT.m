function m = MTQ_BDOT( b_earth, time, maxDipoleMoment )

persistent prev_b_earth
persistent prev_time

if isempty( prev_b_earth )
    prev_b_earth = b_earth;
end
if isempty( prev_time )
    prev_time = time-1;
end

b_earth_dot = ( b_earth - prev_b_earth ) / ( time - prev_time );

if any(isnan(b_earth_dot), 'all')
    warning("NaN found in b_earth_dot")
    b_earth_dot = zeros(3,1);
end

k = 1;

m = - k * maxDipoleMoment .* sign( b_earth_dot );

prev_b_earth = b_earth;
prev_time = time;

end

% figure(1)
% hold on
% axis equal
% grid on
% xlabel("X");
% ylabel("Y");
% zlabel("Z");
% plot3([0, b_earth_body(1)], [0, b_earth_body(2)], [0, b_earth_body(3)],'r')
% plot3([0, m_mtq_body(1)], [0, m_mtq_body(2)], [0, m_mtq_body(3)],'b')
% plot3([0, t_mtq_body(1)], [0, t_mtq_body(2)], [0, t_mtq_body(3)],'c')
% legend("Earth magnetic field", "m", "t_result_mtq")
% hold off

