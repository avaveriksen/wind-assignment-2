function F = myfunctions(x)
    VS1 = x(1);
    VS2 = x(2);
   
    F(1) = VS1 - Z21*(P1/VS1 - (VS1 - (P1*Z11)/VS1)/Z31) - (P1*Z11)/VS1 - VB;
    F(2) = VS2 - Z22*(P2/VS2 - (VS2 - (P2*Z12)/VS2)/Z32) - (P2*Z12)/VS2 - VB;
    F(3) = (VPOC + Z5*(P1/VS1 + P2/VS2 - VA1/Z31 - VA2/Z32))/(Z5/Z41 + 1) - VB;
end