function mHom = tmat2hom(m)
%TMAT2HOM converts the input transformation matrix to homogeneous form.

s = size(m);
if ~all(s == s(1))
    error('The input matrix must be a square matrix.');
end

mHom = eye(s(1) + 1);
mHom(1:s(1),1:s(2)) = m;