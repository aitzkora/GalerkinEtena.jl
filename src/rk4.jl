function rk4(op, u₀::Array{Float64,2}, tInit::Float64, tFinal::Float64, dt::Float64)
    # Low storage Runge-Kutta coefficients
    rk4a = [ 0.0 
            -567301805773.0/1357537059087.0
            -2404267990393.0/2016746695238.0
            -3550918686646.0/2091501179385.0
            -1275806237668.0/842570457699.0];
    rk4b = [ 1432997174477.0/9575080441755.0
            5161836677717.0/13612068292357.0
            1720146321549.0/2090206949498.0
            3134564353537.0/4481467310338.0
            2277821191437.0/14882151754819.0];
    rk4c = [             0.0  
            1432997174477.0/9575080441755.0
            2526269341429.0/6820363962896.0
            2006345519317.0/3224310063776.0
            2802321613138.0/2924317926251.0];
    t = tInit
    u = copy(u₀)
    resu = zeros(size(u))
    Nsteps = Int64(ceil((tFinal - tInit) / dt))
    dt = (tFinal-tInit) /Nsteps
    for tstep = 1:Nsteps
        for INTRK=1:5
            tloc = t + rk4c[INTRK] .* dt
            rhsu = op(u, tloc)
            resu = rk4a[INTRK] .* resu .+ dt .* rhsu
            u = u + rk4b[INTRK] .* resu
        end
        t += dt
    end
    return u
end

