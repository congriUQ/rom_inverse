function [out, d_out] = constant(x, mask)
%Dummy feature function for constant bias. Returns 1 as feature function
%value and 0 as gradient
%Mask is just a dummy parameter every feature function gets

out = 1;
d_out = 0*x(:);
end

