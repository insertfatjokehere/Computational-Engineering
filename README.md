# Computational-Engineering

## FDTD_First_Example.m

This was the first assignment in MATLAB. The code initially produces an Electro-Magnetic Wave, then hits a permeable surface.
To reflect the wave back from the surface, you change

```matlab
    %Absorbing boundary condition at right-hand side
    u3(NX) = u2(NX-1) + a2*(u3(NX-1)-u2(NX));
```

into this:

```matlab
    %Absorbing boundary condition at right-hand side
    u3(NX) = 0; % u2(NX-1) + a2*(u3(NX-1)-u2(NX));
```
## FDTD_1D.m

This was the second assignment. This works on the same principle as the first assignment, except you change the Stability Criterian. To manipulate the stability of the waves, you change `dt` to be within `0` and `1`.

```matlab
dx=lambda/20.0; %space increment of 1-D lattice
dt=dx/cc; %time step
```
