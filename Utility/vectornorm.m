function vecnorm =  vectornorm(A,B,dim)
vecnorm = sqrt(sum((A-B).^2,dim));
end
