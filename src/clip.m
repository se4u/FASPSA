function m = clip(m, val)
m(m > val) = val;
m(m < -val) = -val;