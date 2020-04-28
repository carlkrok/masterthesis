function [f_cmd, t_ref_sat_body] = PROP_PD( q_error, omega_error, ...
        K_p, K_d, com_struct, RotMat_structToBody, thrusters, maxThrust, t )

f_cmd = zeros(length(thrusters),1);
    
t_ref_sat_body =  - K_p * q_error(2:4) ...
     - K_d * omega_error; 
 
%t_ref_sat_struct = RotMat_structToBody' * t_ref_sat_body;
 
A_mat = zeros(3,length(f_cmd));

 for propIter = 1:length(f_cmd)
    %this_t_dir = cross( (thrusters(propIter).structArm - com_struct), ...
    %    -thrusters(propIter).structureThrustDir);
    %this_t_dir = this_t_dir / norm(this_t_dir);
    
    this_t_dir = cross( RotMat_structToBody * ...
        (thrusters(propIter).structArm - com_struct), ...
        RotMat_structToBody * (-thrusters(propIter).structureThrustDir));
    
    A_mat(:,propIter) = this_t_dir;
    
%     this_t_can = ( t_ref_sat_struct ) .* this_t_dir
%     if any( this_t_can > 0 )
%         f_cmd(propIter) = norm(this_t_can) / ...
%         norm(thrusters(propIter).structArm - com_struct);
%         t_ref_sat = t_ref_sat - this_t_can;
%     end
 end
 
A_MPinv_mat = A_mat'/(A_mat*A_mat');
 
 f_cmd = A_MPinv_mat * t_ref_sat_body;

%  if t > 1300
%     disp('here')
%  end

for propIter = 1:length(f_cmd)
    if f_cmd(propIter) < 0
        f_cmd(propIter) = 0;
    elseif f_cmd(propIter) > maxThrust
        f_cmd(propIter) = maxThrust;
    end
end
 
%f_cmd = lsqlin(A_mat, t_ref_sat_struct,[],[],[],[],zeros(6,1), ...
%    []); %maxThrust.*ones(6,1));

end

