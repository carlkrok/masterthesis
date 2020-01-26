function [ moment_body ] = Forces_RW( rw_w, rw_rotMat, rw_I_mat )

moment_body = zeros(3,1);

rwIter = 1;
for thisRW_w = rw_w
    moment_body = moment_body + rw_rotMat(:,:,rwIter) * (rw_I_mat * thisRW_w);
    rwIter = rwIter + 1;
end

