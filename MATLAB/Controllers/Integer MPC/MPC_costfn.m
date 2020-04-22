function res = MPC_costfn( U, Y_lb_tot, Y_ub_tot, Y_lb_dotVec_tot, ...
    Y_ub_dotVec_tot, B_chi, chi_0, H, W)

global satData

U(7:13) = satData.propulsion.maxThrust .* U(7:13);

Y = chi_0 + B_chi*U';

res = Y'*H*Y; % + U*W*U'; % .Fval

% res.Ineq = 0;
% if any(Y.*Y_ub_dotVec_tot - Y_ub_tot)
%     res.Ineq = 1;
% elseif any(-1 .* Y.*Y_lb_dotVec_tot - Y_lb_tot)
%     res.Ineq = 1;
% end

end

