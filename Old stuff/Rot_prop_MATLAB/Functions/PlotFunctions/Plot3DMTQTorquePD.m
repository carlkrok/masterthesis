function Plot3DMTQTorquePD( mtq_t, mtq_t_cmd, mtq_m, mtq_b )

figure
hold on
axis equal
grid on
xlabel("X");
ylabel("Y");
zlabel("Z");
for tIter = 1:length(mtq_t(:,4))
    plot3([0, mtq_b(tIter,1)], [0, mtq_b(tIter,2)], [0, mtq_b(tIter,3)],'r')
    %plot3([0, mtq_m(tIter,1)], [0, mtq_m(tIter,2)], [0, mtq_m(tIter,3)],'g')
    plot3([0, mtq_t(tIter,1)], [0, mtq_t(tIter,2)], [0, mtq_t(tIter,3)],'b')
    plot3([0, mtq_t_cmd(tIter,1)], [0, mtq_t_cmd(tIter,2)], [0, mtq_t_cmd(tIter,3)],'c')
    %legend("Earth magnetic field", "m", "t_result_mtq", "t_ref_sat")
end
legend("Earth magnetic field", "Resulting \tau", "Reference \tau")
hold off

end

