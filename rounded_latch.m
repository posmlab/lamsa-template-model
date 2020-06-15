function latch = rounded_latch(R, m_L, v_0, coeff_fric)
latch.max_width=R;
latch.mass=m_L;
latch.v_0=v_0;
yL=@(x) latch.max_width*(1-sqrt(1-x^2/latch.max_width^2));
syms x;
yL_prime = diff(yL(x));
yL_doubleprime = diff(yL(x),2);
latch.y_L = {yL, matlabFunction(yL_prime), matlabFunction(yL_doubleprime)};
latch.coeff_fric = coeff_fric;
end 