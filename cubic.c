#include <stdio.h>
#include <math.h>

typedef struct complex {
    double real, imag;
    } COMPLEX;

COMPLEX Complex(),      /* returns a complex made from the real args */
   ctimes(),                    /* product -- dyadic */
   cplus();                     /* sum -- dyadic */

#define Pi 3.14159265358979323846

/*
 * cubic -- apply closed-form solution to an arbitrary cubic
 *          in real coefficients
 * based on the exposition of the algorithm in the CRC Standard
 * Mathematical Tables, 27'th edition, CRC Press, Boca Raton FL,
 * 1984, on page 9, and Chrystal: Textbook of Algebra
 *
 * syntax: cubic c3 c2 c1 c0
 * Where the eqn is c3 x^3 + c2 x^2 + c1 x + c0 = 0
 * ci are arbitrary constants in format suitable
 * for conversion using "%f" format in scanf().
 */
/*
 * Written by Bennett Todd
 * Date: 10/15/85 -- initial release
 * This code is released to the public domain; do with it as you wish.
 * If you use it, credit would be appreciated.
 * Please send any bug fixes, enhancements, or general comments to
 *      Bennett Todd
 *      Duke Computation Center
 *      Durham, NC 27706-7756
 *      +1 919 684 3695
 *      UUCP: ...{decvax,seismo,philabs,ihnp4,akgua}!mcnc!ecsvax!duccpc!bet
 *      BITNET: dbtodd@tucc
 *
 * Changed by: Alan Wendt, U of AZ CS, arizona!wendt.
 */
main(argc, argv)
    int argc;
    char **argv;
    {
    /* Sorry about the variable names; they are chosen to match CRC */
    double o, p, q, r, a, b, t, l, m;
    COMPLEX root_1, root_2, root_3, omega, omegasquared;

    omega = Complex(-0.5, sqrt(3.0) / 2.0);     /* third roots of unity */
    omegasquared = Complex(-0.5, -sqrt(3.0) / 2.0);

    /* (attempt to) parse command line */
    if (argc != 5) {
        fprintf(stderr, "syntax: cubic c3 c2 c1 c0\n");
        fprintf(stderr, "Indicate missing terms with zeroes.\n");
        exit(1);
        }
    if (sscanf(argv[1], "%lf", &o) != 1) {
        fprintf(stderr, "cubic: cannot parse %s\n", argv[1]);
        exit(1);
        }
    if (sscanf(argv[2], "%lf", &p) != 1) {
        fprintf(stderr, "cubic: cannot parse %s\n", argv[2]);
        exit(1);
        }
    if (sscanf(argv[3], "%lf", &q) != 1) {
        fprintf(stderr, "cubic: cannot parse %s\n", argv[3]);
        exit(1);
        }
    if (sscanf(argv[4], "%lf", &r) != 1) {
        fprintf(stderr, "cubic: cannot parse %s\n", argv[4]);
        exit(1);
        }

    /* sanity check */
    if (o == 0) {
        fprintf(stderr, "cubic: sorry buddy, that's a quadratic.\n");
        exit(1);
        }

    p /= o;
    q /= o;
    r /= o;

    /* Change of variable from "x^3 + px^2 + qx + r = 0" to
     * "x'^3 + ax' + b = 0" by propagating the substitution
     * "x = x' - p/3".
     */
    a = q - p * p / 3.0;
    b = (2.0 * p * p * p - 9.0 * p * q + 27.0 * r) / 27.0;

    t = b * b / 4.0 + a * a * a / 27.0;

    if (t < 0.0) {
        /*  as b^2 / 4 + a^3/27 < 0, a < 0, -a and rho > 0 */
        double rho, theta;
        rho = pow(-a, 1.5) / sqrt(27.0);
        theta = acos(b / rho / 2.0);
        rho = -2.0 * pow(rho, 1.0 / 3.0);
        root_1 = Complex(cos(theta / 3.0) * rho, 0.0);
        root_2 = Complex(cos((2.0 * Pi + theta) / 3.0) * rho, 0.0);
        root_3 = Complex(cos((4.0 * Pi + theta) / 3.0) * rho, 0.0);
        }

    else {
        t = sqrt(t);

        l = b / 2 + t;
        m = b / 2 - t;

        /* pow doesn't work for negative bases */
        if (l >= 0) l = - pow(l, 1.0 / 3.0);
        else l = pow(-l, 1.0 / 3.0);

        if (m >= 0) m = -pow(m, 1.0 / 3.0);
        else m = pow(-m, 1.0 / 3.0);

        root_1 = Complex(l + m, 0.0);

        root_2 = cplus(
                ctimes(omega, Complex(l, 0.0)),
                ctimes(omegasquared, Complex(m, 0.0))
                );

        root_3 = cplus(
                ctimes(omegasquared, Complex(l, 0.0)),
                ctimes(omega, Complex(m, 0.0))
                );
        }

    p /= 3.0;

    printf("%.16g", root_1.real - p);
    if (root_1.imag != 0) printf("%s%.16g i",
        root_1.imag < 0 ? " - " : " + ",
        fabs(root_1.imag));

    printf("\n%.16g", root_2.real - p);
    if (root_2.imag != 0) printf("%s%.16g i",
        root_2.imag < 0 ? " - " : " + ",
        fabs(root_2.imag));

    printf("\n%.16g", root_3.real - p);
    if (root_3.imag != 0) printf("%s%.16g i",
        root_3.imag < 0 ? " - " : " + ",
        fabs(root_3.imag));

    printf("\n");
    }

/* complex math routines.  */
COMPLEX Complex(x, y)
    double x, y;
    {
    COMPLEX result;

    result.real = x;
    result.imag = y;
    return result;
    }

COMPLEX ctimes(a, b)
    COMPLEX a, b;
    {
    COMPLEX result;

    result.real = a.real*b.real - a.imag*b.imag;
    result.imag = a.real*b.imag + a.imag*b.real;
    return result;
    }

COMPLEX cplus(a, b)
    COMPLEX a, b;
    {
    COMPLEX result;
    result.real = a.real + b.real;
    result.imag = a.imag + b.imag;
    return result;
    } 
