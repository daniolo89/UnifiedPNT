function k = cal_Kv(sinEl, cosAz, sinAz, cosEl, lambda)
    % Calculate the wavenumber vector 
    k = -(2*pi/lambda) * [ sinEl .* cosAz, sinEl .* sinAz, cosEl];
end
