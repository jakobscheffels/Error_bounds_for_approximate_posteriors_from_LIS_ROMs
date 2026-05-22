function B = BMatrix(J,psi,eta)
    B = zeros(3,8);

    for i=1:4
        grad = partialShapeFunction(i,psi,eta);
        grad = J\grad;
        B(1,2*i-1)=grad(1);
        B(2,2*i) = grad(2);
        B(3,2*i-1:2*i) = [grad(2) grad(1)];
    end
end
