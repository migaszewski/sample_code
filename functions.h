/* ------------- metoda Taylora ------------ */
void taylor3d ( double x0[], double *t, double h, struct params_str *params )
{
 int N = (*params).npl, i, j, kk, l, n = (*params).order;
 int i0, j0, ij0, ji0;

 for (i=0; i<N; i++) {
     i0 = n*i;
     (*params).X[i0] = x0[6*i];
     (*params).Y[i0] = x0[6*i+1];
     (*params).Z[i0] = x0[6*i+2];
     (*params).U[i0] = x0[6*i+3];
     (*params).V[i0] = x0[6*i+4];
     (*params).W[i0] = x0[6*i+5];
     (*params).R[i0] = sqrt(x0[6*i]*x0[6*i] + x0[6*i+1]*x0[6*i+1] + x0[6*i+2]*x0[6*i+2]);
     (*params).S[i0] = pow((*params).R[i0], -3.);
     for (j=i+1; j<N; j++) {
         ij0 = n*(N*i + j);
         ji0 = n*(N*j + i);
         (*params).D[ij0] = sqrt(pow(x0[6*i]-x0[6*j], 2.) + pow(x0[6*i+1]-x0[6*j+1], 2.) + pow(x0[6*i+2]-x0[6*j+2], 2.));
         (*params).P[ij0] = pow((*params).D[ij0], -3.);
         (*params).D[ji0] = (*params).D[ij0];
         (*params).P[ji0] = (*params).P[ij0];
     }
 }
 
 for (kk=0; kk<n-1; kk++) {

     for (i=0; i<N; i++) {

         i0 = n*i;
         
         (*params).X[i0+kk+1] = (*params).U[i0+kk]/(kk+1.);
         (*params).Y[i0+kk+1] = (*params).V[i0+kk]/(kk+1.);
         (*params).Z[i0+kk+1] = (*params).W[i0+kk]/(kk+1.);

         (*params).U[i0+kk+1] = 0.0;
         (*params).V[i0+kk+1] = 0.0;
         (*params).W[i0+kk+1] = 0.0;
         
         for (l=0; l<=kk; l++) {

             (*params).U[i0+kk+1] -= (*params).mu[i]*(*params).S[i0+l]*(*params).X[i0+kk-l];
             (*params).V[i0+kk+1] -= (*params).mu[i]*(*params).S[i0+l]*(*params).Y[i0+kk-l];
             (*params).W[i0+kk+1] -= (*params).mu[i]*(*params).S[i0+l]*(*params).Z[i0+kk-l];
             
             for (j=0; j<N; j++) {
                 if (j != i) {
                     ij0 = n*(N*i + j);
                     j0  = n*j;
                     
                     (*params).U[i0+kk+1] -= k2*(*params).m[j]*((*params).P[ij0+l]*((*params).X[i0+kk-l] - (*params).X[j0+kk-l]) + (*params).S[j0+l]*(*params).X[j0+kk-l]);
                     
                     (*params).V[i0+kk+1] -= k2*(*params).m[j]*((*params).P[ij0+l]*((*params).Y[i0+kk-l] - (*params).Y[j0+kk-l]) + (*params).S[j0+l]*(*params).Y[j0+kk-l]);
                     
                     (*params).W[i0+kk+1] -= k2*(*params).m[j]*((*params).P[ij0+l]*((*params).Z[i0+kk-l] - (*params).Z[j0+kk-l]) + (*params).S[j0+l]*(*params).Z[j0+kk-l]);
                 }
             }
             
         }

         (*params).U[i0+kk+1] /= kk+1.;
         (*params).V[i0+kk+1] /= kk+1.;
         (*params).W[i0+kk+1] /= kk+1.;
         
         (*params).R[i0+kk+1] = 0.0;
         for (l=0; l<=kk; l++) (*params).R[i0+kk+1] += (*params).X[i0+l]*(*params).U[i0+kk-l] + (*params).Y[i0+l]*(*params).V[i0+kk-l] + (*params).Z[i0+l]*(*params).W[i0+kk-l];             
         for (l=1; l<=kk; l++) (*params).R[i0+kk+1] -= (kk-l+1)*(*params).R[i0+l]*(*params).R[i0+kk-l+1];             
         (*params).R[i0+kk+1] /= (kk+1.)*(*params).R[i0];
         
         (*params).S[i0+kk+1] = 0.0;
         for (l=0; l<=kk; l++) (*params).S[i0+kk+1] -= 3*(kk-l+1)*(*params).S[i0+l]*(*params).R[i0+kk-l+1];
         for (l=1; l<=kk; l++) (*params).S[i0+kk+1] -= (kk-l+1)*(*params).R[i0+l]*(*params).S[i0+kk-l+1];
         (*params).S[i0+kk+1] /= (kk+1.)*(*params).R[i0];
     
         for (j=i+1; j<N; j++) {
         
             ij0 = n*(N*i + j);
             ji0 = n*(N*j + i);
             j0  = n*j;
                          
             (*params).D[ij0+kk+1] = 0.0;
             for (l=0; l<=kk; l++) {
                 (*params).D[ij0+kk+1] += ((*params).X[i0+l] - (*params).X[j0+l])*((*params).U[i0+kk-l] - (*params).U[j0+kk-l]);
                 (*params).D[ij0+kk+1] += ((*params).Y[i0+l] - (*params).Y[j0+l])*((*params).V[i0+kk-l] - (*params).V[j0+kk-l]);
                 (*params).D[ij0+kk+1] += ((*params).Z[i0+l] - (*params).Z[j0+l])*((*params).W[i0+kk-l] - (*params).W[j0+kk-l]);
             }
             for (l=1; l<=kk; l++) (*params).D[ij0+kk+1] -= (kk-l+1)*(*params).D[ij0+l]*(*params).D[ij0+kk-l+1];
             (*params).D[ij0+kk+1] /= (kk+1.)*(*params).D[ij0];
             (*params).D[ji0+kk+1] = (*params).D[ij0+kk+1];
             
             (*params).P[ij0+kk+1] = 0.0;
             for (l=0; l<=kk; l++) (*params).P[ij0+kk+1] -= 3*(kk-l+1)*(*params).P[ij0+l]*(*params).D[ij0+kk-l+1];
             for (l=1; l<=kk; l++) (*params).P[ij0+kk+1] -= (kk-l+1)*(*params).D[ij0+l]*(*params).P[ij0+kk-l+1];
             (*params).P[ij0+kk+1] /= (kk+1.)*(*params).D[ij0];
             (*params).P[ji0+kk+1] = (*params).P[ij0+kk+1];
         
         }
     
     }
     
 }
 
 for (i=0; i<N; i++) {

     i0 = n*i;
     
     x0[6*i]   = (*params).X[i0+n-1];
     x0[6*i+1] = (*params).Y[i0+n-1];
     x0[6*i+2] = (*params).Z[i0+n-1];
     x0[6*i+3] = (*params).U[i0+n-1];
     x0[6*i+4] = (*params).V[i0+n-1];
     x0[6*i+5] = (*params).W[i0+n-1];

     for (kk=n-2; kk>=0; kk--) {
         x0[6*i]   = (*params).X[i0+kk] + h*x0[6*i];
         x0[6*i+1] = (*params).Y[i0+kk] + h*x0[6*i+1];
         x0[6*i+2] = (*params).Z[i0+kk] + h*x0[6*i+2];
         x0[6*i+3] = (*params).U[i0+kk] + h*x0[6*i+3];
         x0[6*i+4] = (*params).V[i0+kk] + h*x0[6*i+4];
         x0[6*i+5] = (*params).W[i0+kk] + h*x0[6*i+5];
     }
 }
 
 *t += h;

}

double kepler_equation (double e, double M)
{
 double E, E0, f0, f1, f2, f3, d1, d2, d3, delta, nu;
 int i;

  E = M;
  
  delta=1.;
 
  i=0;
  while ((delta>1e-12) && (i<100)) {
     ++i;
     E0=E;
     f0=E-e*sin(E)-M;
     f1=1-e*cos(E);
     f2=e*sin(E);
     f3=e*cos(E);
     d1=-f0/f1;
     d2=-f0/(f1+d1*f2/2);
     d3=-f0/(f1+d2*f2/2+d2*d2*f3/6);
     E=E+d3;
     delta=fabs(E0-E);
  }
 
  nu = 2.*atan2(sqrt(1. + e)*sin(E/2.), sqrt(1. - e)*cos(E/2.));

 return nu;

}

void kepl2cart ( double kepl[], double cart[], struct params_str params )
{
 int N = params.npl;
 double a[N], e[N], I[N], Om[N], om[N], f[N];
 double P11, P12, P21, P22, P31, P32;
 double r, h, rd, X, Y, Xd, Yd;
 double x[N], y[N], z[N];
 double vx[N], vy[N], vz[N];
 double xd[N], yd[N], zd[N];
 double Vx, Vy, Vz;
 int i, ii;
 
 for (i=0; i<N; i++) {
    ii = 6*i;
    a[i]  = kepl[ii];
    e[i]  = kepl[ii+1];
    I[i]  = kepl[ii+2];
    Om[i] = kepl[ii+3];
    om[i] = kepl[ii+4];
    f[i]  = kepl[ii+5];  // Mmean
    f[i]  = kepler_equation(e[i], f[i]);
 }
 
 for (i=0; i<N; i++) {
 
    ii = 6*i;

    P11 = cos(Om[i])*cos(om[i]) - sin(Om[i])*cos(I[i])*sin(om[i]);
    P12 = -cos(Om[i])*sin(om[i]) - sin(Om[i])*cos(I[i])*cos(om[i]);
    P21 = sin(Om[i])*cos(om[i]) + cos(Om[i])*cos(I[i])*sin(om[i]);
    P22 = -sin(Om[i])*sin(om[i]) + cos(Om[i])*cos(I[i])*cos(om[i]);
    P31 = sin(I[i])*sin(om[i]);
    P32 = sin(I[i])*cos(om[i]);
    
    r = a[i]*(1.0 - e[i]*e[i]) / (1 + e[i]*cos(f[i]));
    h = sqrt(params.mu[i]*a[i]*(1.0 - e[i]*e[i]));
    rd = h*e[i]*sin(f[i])/(a[i]*(1.0 - e[i]*e[i]));
    
    X = r*cos(f[i]);
    Y = r*sin(f[i]);
    Xd = rd*cos(f[i]) - (h/r)*sin(f[i]);
    Yd = rd*sin(f[i]) + (h/r)*cos(f[i]);
    
    x[i] = P11*X + P12*Y;
    y[i] = P21*X + P22*Y;
    z[i] = P31*X + P32*Y;
    
    vx[i] = P11*Xd + P12*Yd;
    vy[i] = P21*Xd + P22*Yd;
    vz[i] = P31*Xd + P32*Yd;
    
    cart[ii]   =  x[i];
    cart[ii+1] =  y[i];
    cart[ii+2] =  z[i];
    cart[ii+3] = vx[i];
    cart[ii+4] = vy[i];
    cart[ii+5] = vz[i];
    
 }

}

void cart2kepl_a ( double cart[], double kepl[], struct params_str params )
{
 int N = params.npl;
 double x[N], y[N], z[N];
 double xd[N], yd[N], zd[N];
 double vx[N], vy[N], vz[N];
 double h[N], r[N], v2[N];
 double Vx, Vy, Vz, M;
 double a[N], e[N], I[N], Om[N], om[N], f[N];
 double sinOm, cosOm, sinomf, cosomf, sinf, cosf, rd, omf;
 double hx[N], hy[N], hz[N];
 int i, ii;
 double rv[N], E[N];
 double cosE, sinE, cosw, sinw, x0, y0, Mmean[N];
 
 for (i=0; i<N; i++) {
    ii = 6*i;
    x[i]  = cart[ii];
    y[i]  = cart[ii+1];
    z[i]  = cart[ii+2];
    xd[i] = cart[ii+3];
    yd[i] = cart[ii+4];
    zd[i] = cart[ii+5];
 }

 for (i=0; i<N; i++) {
    vx[i] = xd[i];
    vy[i] = yd[i];
    vz[i] = zd[i];
 }
 
 for (i=0; i<N; i++) {
    hx[i] = y[i]*vz[i] - z[i]*vy[i];
    hy[i] = z[i]*vx[i] - x[i]*vz[i];
    hz[i] = x[i]*vy[i] - y[i]*vx[i];
    h[i] = sqrt(hx[i]*hx[i] + hy[i]*hy[i] + hz[i]*hz[i]);
    r[i] = sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i]);
    v2[i] = vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
 }
 
 for (i=0; i<N; i++) {
    a[i] = 1./(2./r[i] - v2[i]/params.mu[i]);
    e[i] = sqrt(1. - h[i]*h[i]/(params.mu[i]*a[i]));
    if ((z[i] == 0.0) && (zd[i] == 0.0)) {
        I[i] = 0.0;
        Om[i] = 0.0;
        rd = (x[i]*vx[i] + y[i]*vy[i] + z[i]*vz[i])/r[i];
        sinf = a[i]*(1. - e[i]*e[i])*rd/(h[i]*e[i]);
        cosf = (a[i]*(1. - e[i]*e[i])/r[i] - 1.)/e[i];
        f[i] = atan2(sinf, cosf); // prawdziwa
        sinomf = y[i]/r[i];
        cosomf = x[i]/r[i];
        omf = atan2(sinomf, cosomf);
        om[i] = omf - f[i];
        om[i] = atan2(sin(om[i]),cos(om[i]));
        f[i] = 2.*atan2(sqrt(1 - e[i])*sin(f[i]/2.), sqrt(1 + e[i])*cos(f[i]/2.));  // mimosrodowa
        f[i] = f[i] - e[i]*sin(f[i]);  // srednia
        f[i] = atan2(sin(f[i]), cos(f[i]));
    }
    else {
        I[i] = acos(hz[i]/h[i]);
        sinOm = hx[i]/(h[i]*sin(I[i]));
        cosOm = -hy[i]/(h[i]*sin(I[i]));
        Om[i] = atan2(sinOm, cosOm);
        sinomf = z[i]/(r[i]*sin(I[i]));
        cosomf = (x[i]/r[i] + sinOm*sinomf*cos(I[i]))/cosOm;
        omf = atan2(sinomf, cosomf);
        rd = (x[i]*vx[i] + y[i]*vy[i] + z[i]*vz[i])/r[i];
        sinf = a[i]*(1 - e[i]*e[i])*rd/(h[i]*e[i]);
        cosf = (a[i]*(1 - e[i]*e[i])/r[i] - 1.)/e[i];
        f[i] = atan2(sinf, cosf); // prawdziwa
        om[i] = omf - f[i];
        om[i] = atan2(sin(om[i]),cos(om[i]));
        f[i] = 2.*atan2(sqrt(1 - e[i])*sin(f[i]/2.), sqrt(1 + e[i])*cos(f[i]/2.));  // mimosrodowa
        f[i] = f[i] - e[i]*sin(f[i]);  // srednia
        f[i] = atan2(sin(f[i]), cos(f[i]));
    }
 }
 
 for (i=0; i<N; i++) {
    ii =6*i;
    kepl[ii]   = a[i];
    kepl[ii+1] = e[i];
    kepl[ii+2] = I[i];
    kepl[ii+3] = Om[i];
    kepl[ii+4] = om[i];
    kepl[ii+5] = f[i];
 }

}

void cart2kepl_P ( double cart[], double kepl[], struct params_str params )
{
 int N = params.npl;
 double x[N], y[N], z[N];
 double xd[N], yd[N], zd[N];
 double vx[N], vy[N], vz[N];
 double h[N], r[N], v2[N];
 double Vx, Vy, Vz, M;
 double a[N], e[N], I[N], Om[N], om[N], f[N];
 double sinOm, cosOm, sinomf, cosomf, sinf, cosf, rd, omf;
 double hx[N], hy[N], hz[N];
 int i, ii;
 double rv[N], E[N];
 double cosE, sinE, cosw, sinw, x0, y0, Mmean[N];
 
 for (i=0; i<N; i++) {
    ii = 6*i;
    x[i]  = cart[ii];
    y[i]  = cart[ii+1];
    z[i]  = cart[ii+2];
    xd[i] = cart[ii+3];
    yd[i] = cart[ii+4];
    zd[i] = cart[ii+5];
 }

 Vx = Vy = Vz = 0.0;
 M = params.m0;
 for (i=0; i<N; i++) {
    Vx += params.m[i]*xd[i];
    Vy += params.m[i]*yd[i];
    Vz += params.m[i]*zd[i];
    M  += params.m[i];
 }
  
 for (i=0; i<N; i++) {
    vx[i] = (params.m[i]/params.beta[i])*(xd[i] - Vx/M);
    vy[i] = (params.m[i]/params.beta[i])*(yd[i] - Vy/M);
    vz[i] = (params.m[i]/params.beta[i])*(zd[i] - Vz/M);
 }
 
 for (i=0; i<N; i++) {
    hx[i] = y[i]*vz[i] - z[i]*vy[i];
    hy[i] = z[i]*vx[i] - x[i]*vz[i];
    hz[i] = x[i]*vy[i] - y[i]*vx[i];
    h[i] = sqrt(hx[i]*hx[i] + hy[i]*hy[i] + hz[i]*hz[i]);
    r[i] = sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i]);
    v2[i] = vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
 }
 
 for (i=0; i<N; i++) {
    a[i] = 1./(2./r[i] - v2[i]/params.mu[i]);
    e[i] = sqrt(1. - h[i]*h[i]/(params.mu[i]*a[i]));
    if ((z[i] == 0.0) && (zd[i] == 0.0)) {
        I[i] = 0.0;
        Om[i] = 0.0;
        rd = (x[i]*vx[i] + y[i]*vy[i] + z[i]*vz[i])/r[i];
        sinf = a[i]*(1. - e[i]*e[i])*rd/(h[i]*e[i]);
        cosf = (a[i]*(1. - e[i]*e[i])/r[i] - 1.)/e[i];
        f[i] = atan2(sinf, cosf); // prawdziwa
        sinomf = y[i]/r[i];
        cosomf = x[i]/r[i];
        omf = atan2(sinomf, cosomf);
        om[i] = omf - f[i];
        om[i] = atan2(sin(om[i]),cos(om[i]));
        f[i] = 2.*atan2(sqrt(1 - e[i])*sin(f[i]/2.), sqrt(1 + e[i])*cos(f[i]/2.));  // mimosrodowa
        f[i] = f[i] - e[i]*sin(f[i]);  // srednia
        f[i] = atan2(sin(f[i]), cos(f[i]));
    }
    else {
        I[i] = acos(hz[i]/h[i]);
        sinOm = hx[i]/(h[i]*sin(I[i]));
        cosOm = -hy[i]/(h[i]*sin(I[i]));
        Om[i] = atan2(sinOm, cosOm);
        sinomf = z[i]/(r[i]*sin(I[i]));
        cosomf = (x[i]/r[i] + sinOm*sinomf*cos(I[i]))/cosOm;
        omf = atan2(sinomf, cosomf);
        rd = (x[i]*vx[i] + y[i]*vy[i] + z[i]*vz[i])/r[i];
        sinf = a[i]*(1 - e[i]*e[i])*rd/(h[i]*e[i]);
        cosf = (a[i]*(1 - e[i]*e[i])/r[i] - 1.)/e[i];
        f[i] = atan2(sinf, cosf); // prawdziwa
        om[i] = omf - f[i];
        om[i] = atan2(sin(om[i]),cos(om[i]));
        f[i] = 2.*atan2(sqrt(1 - e[i])*sin(f[i]/2.), sqrt(1 + e[i])*cos(f[i]/2.));  // mimosrodowa
        f[i] = f[i] - e[i]*sin(f[i]);  // srednia
        f[i] = atan2(sin(f[i]), cos(f[i]));
    }
 }
 
 for (i=0; i<N; i++) {
    ii =6*i;
    kepl[ii]   = a[i];
    kepl[ii+1] = e[i];
    kepl[ii+2] = I[i];
    kepl[ii+3] = Om[i];
    kepl[ii+4] = om[i];
    kepl[ii+5] = f[i];
 }

}

double energy(double cart[], struct params_str params )
{
 int N = params.npl;
 double V[3], H, v2, r, d[3], del, M;
 int i, j, ii, jj;
 double Vx, Vy, Vz;
 
 
 M = params.m0;
 for (i=0; i<N; i++)  M += params.m[i];
 
 V[0] = 0.0;
 V[1] = 0.0;
 V[2] = 0.0;
 for (i=0; i<N; i++) {
     ii = 6*i;
     V[0] += params.m[i]*cart[ii+3];
     V[1] += params.m[i]*cart[ii+4];
     V[2] += params.m[i]*cart[ii+5];
 }
 V[0] = -V[0]/M;
 V[1] = -V[1]/M;
 V[2] = -V[2]/M;
 
 H = 0.5*params.m0*(V[0]*V[0] + V[1]*V[1] + V[2]*V[2]);
 
 for (i=0; i<N; i++) {
     ii = 6*i;
     Vx = cart[ii+3] + V[0];
     Vy = cart[ii+4] + V[1];
     Vz = cart[ii+5] + V[2];
     v2 = Vx*Vx + Vy*Vy + Vz*Vz;
     r  = sqrt(pow(cart[ii],2.) + pow(cart[ii+1],2.) + pow(cart[ii+2],2.));
     H += 0.5*params.m[i]*v2 - k2*params.m0*params.m[i]/r;
     for (j=i+1; j<N; j++) {
	 jj = 6*j;
         d[0] = cart[ii] - cart[jj];
	 d[1] = cart[ii+1] - cart[jj+1];
 	 d[2] = cart[ii+2] - cart[jj+2];
	 del = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
	 del = sqrt(del);
	 H -= k2*params.m[i]*params.m[j]/del;
     }
 }

 return H;

}
