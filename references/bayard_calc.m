function error = bayard_calc(gyro, st, dt);

  
  
q1 = gyro.random_walk^2;

q2 = gyro.bias_stability^2/3600;


l = sqrt(q1 + 2 * sqrt(st.r * q2));

p11 = sqrt(st.r)*l;

p12 = sqrt(st.r*q2);

p22 = sqrt(q2)*l;

p = q2/3*dt^3 + p22*dt^2 + (2*p12 + q1)*dt + p11 + st.b^2;

error = sqrt(p);
