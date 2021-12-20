using ControlSystems

# Example from the paper
L = tf(25, [1,10,10,10])
dm = diskmargin(L, 0)
@test dm.ω0 ≈ 1.94   atol=0.02
@test dm.γmin ≈ 0.63    atol=0.02
@test dm.γmax ≈ 1.59    atol=0.02
@test dm.α ≈ 0.46   atol=0.02
show(dm)
plot(dm)
nyquistplot(L)
plot!(dm, nyquist=true)
plot!(Disk(dm), nyquist=true)




## Frequency-dependent margin
w = exp10.(LinRange(-2, 2, 500))
dms = diskmargin(L, 0, w)
plot(w, dms)

##
s = tf("s")
L = 6.25*(s + 3)*(s + 5) / (s*(s + 1)^2 *(s^2 + 0.18s + 100))

## αmax > 2
dm = diskmargin(L, 0, 200)
@test dm.γmax < dm.γmin
@test dm.γmin ≈ 0 atol = 1e-6
@test dm.ϕm ≈ 90 atol = 1e-4


w = exp10.(LinRange(-1, 2, 500))
dms = diskmargin(L, 0, w)
plot(w, dms) 


## MIMO loop at a time
a = 10
P = [
        tf([1,-a^2], [1, 0, a^2]) tf([a, a], [1, 0, a^2])
        -tf([a, a], [1, 0, a^2]) tf([1,-a^2], [1, 0, a^2])
    ]
P = minreal(ss(P))
K = ss(1.0I(2))
Li = K*P
Lo = P*K

@test tf(minreal(ControlSystems.broken_feedback(Li, 1))) ≈ tf(1, [1, 0])
@test tf(minreal(ControlSystems.broken_feedback(Lo, 1))) ≈ tf(1, [1, 0])

dm = diskmargin(Li)[1]
@test dm.α ≈ 2
@test dm.γmax > 1e10
@test dm.ϕm ≈ 90

dmm = diskmargin(P, K)
dmm.input[1].α == dm.α


# using IntervalArithmetic
# δ(a=1) = Complex(-a..a, -a..a)
# Δ(n, a) = diagm([δ(a) for _ in 1:n])
# M = [0 1; -0.1 -0.1]
# D = Δ(2, 1)
# 0 ∈ det(I-M*D)


L3 = let
    tempA = [1.0 0.0 9.84 0.0 0.0; 0.0 1.0 0.01 2.14634e-6 0.0; 0.0 0.0 1.0 0.0 0.0; 0.0 0.0 0.0 1.0 -1.73983959887; 0.0 0.0 0.0 0.0 0.56597684805]
    tempB = [0.0 4.901416e-5 0.00019883877999999998; 0.0 0.0 0.0; 0.0 0.0 4.0414389999999996e-5; 0.0 -0.02004649371 0.0; 0.0 -0.00490141631 0.0]
    tempC = [0.0 -0.83516488404 0.0 0.0 0.0; 186.74725411661 0.0 0.0 0.0 0.0; -7.44299057498 0.0 7035.08410814126 0.0 0.0]
    tempD = [0.0 0.0 0.0; 34875.36444283988 0.0 0.0; 48304.01940122544 0.0 0.0]
    ss(tempA, tempB, tempC, tempD, 0.01)
end

w = 2π .* exp10.(LinRange(-2, 2, 300))

dm = diskmargin(L3, 0, 4.05)
@test dm[1].α ≈ 0.794418036911981 rtol=1e-3

dm = diskmargin(L3, 0)

dm = diskmargin(L3, ss(I(3), L3.Ts), 0, w)
plot(w, dm)


