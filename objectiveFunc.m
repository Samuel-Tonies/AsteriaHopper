function J = objectiveFunc(DV,N)
    nNodes = N+1;
    F = DV(4*nNodes+1:5*nNodes);
    tf = DV(end);
    dt = tf / N;
    J = 0.5 * dt * sum(F(1:end-1).^2 + F(2:end).^2);
end