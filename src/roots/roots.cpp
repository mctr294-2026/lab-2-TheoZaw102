#include "roots.hpp"
#include <cmath>
#include <algorithm>

const double TOL = 1e-6;
const int Max_Iter = 1e6;

static bool bracketed(double fa, double fb) {
    if (fa == 0.0 || fb == 0.0) return true;
    return (fa>0) != (fb>0);
}

bool bisection(std::function<double(double)> f,
               double a, double b,
               double *root) {
    if (!root) return false;
    if (a > b) std::swap(a,b);

    double fa = f(a);
    double fb = f(b);

    if (std::abs(fa) <= TOL) {*root=a; return true;}
    if (std::abs(fb) <= TOL) {*root=b; return true;}

    if (!bracketed(fa, fb)) return false;

    for (int i = 0; i < Max_Iter; i++) {
        double mid = (a + b) / 2.0;
        double fm = f(mid);

        if (std::abs(fm) <= TOL) {
            *root = mid;
            return true;
        }
        if (std::abs(b-a) <= TOL) {
            *root = mid;
            return true;
        }

        if ((fm > 0)==(fa > 0)) {
            a = mid;
            fa = fm;
        } else {
            b= mid;
            fb = fm;
        }
    }
    return false;
}

bool regula_falsi(std::function<double(double)> f,
                  double a, double b,
                  double *root) {
    if (!root) return false;
    if (a>b) std::swap(a,b);

    double fa = f(a);
    double fb = f(b);

    if (std::abs(fa) <= TOL) {*root=a; return true;}
    if (std::abs(fb) <= TOL) {*root=b; return true;}

    if (!bracketed(fa, fb)) return false;

    for (int i = 0; i < Max_Iter; i++) {
        double botfrac = (fb - fa);
        if (std::abs(botfrac) < TOL) return false;

        double c = a - (fa*(b-a)/botfrac);
        double fc = f(c);

        if (std::abs(fc) <= TOL) {
            *root = c;
            return true;
        }
        if ((fc>0) == (fa>0)) {
            a = c;
            fa = fc;
        } else {
            b = c;
            fb = fc;
        }
        if (std::abs(b-a) <= TOL) {
            *root =c;
            return true;
        } 
    }
    return false;
    
}

bool newton_raphson(std::function<double(double)> f,
                    std::function<double(double)> g,
                    double a, double b, double c,
                    double *root) {
    if (!root) return false;
    if (a>b) std::swap(a,b);

    double x = c;
    if (x<a) x = a;
    if (x>b) x = b;

    for (int i = 0; i < Max_Iter; i++) {
        double fx = f(x);
        if (std::abs(fx) <= TOL) {
            *root = x;
            return true;
        }

        double gx = g(x);
        if (std::abs(gx) < TOL) return false;

        double x2 = x - fx / gx;

        if (x2 < a || x2 > b) return false;

        if (std::abs(x2 - x) <= TOL) {
            *root = x2;
            return true;
        }
        x = x2;
    }
    return false;
}

bool secant(std::function<double(double)> f,
            double a, double b, double c,
            double *root) {
    if (!root) return false;
    if (a>b) std::swap(a,b);

    double x0 = c;
    if (x0 < a) x0 = a;
    if (x0 > b) x0 = b;

    double x1 = b;
    if (x1 == x0) x1 = a;

    double f0 = f(x0);
    double f1 = f(x1);

    for (int i = 0; i < Max_Iter; i++) {
        if (std::abs(f1) <= TOL) {
            *root = x1;
            return true;
        }
        double botfrac = (f1 - f0);
        if (std::abs(botfrac) < TOL) return false;

        double x2 = x1 - (f1*(x1-x0)/botfrac);
        if (x2 < a || x2 > b) return false;

        if (std::abs(x2 - x1) <= TOL) {
            *root = x2;
            return true;
        }
        x0 = x1;
        f0 = f1;
        x1 = x2;
        f1 = f(x1);
    }
    return false;
}

