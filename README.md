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
