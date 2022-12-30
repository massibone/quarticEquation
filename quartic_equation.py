from cmath import *

import numpy as np

P = np.poly1d([1, 2, 3, 4], True)


def ferrari(a, b, c, d, e):
    "resolution of P=ax^4+bx^3+cx^2+dx+e=0"
    "CN all coeffs real."
    "First shift : x= z-b/4/a  =>  P=z^4+pz^2+qz+r"
    z0 = b/4/a
    a2, b2, c2, d2 = a*a, b*b, c*c, d*d
    p = -3*b2/(8*a2)+c/a
    q = b*b2/8/a/a2 - 1/2*b*c/a2 + d/a
    r = -3/256*b2*b2/a2/a2 + c*b2/a2/a/16-b*d/a2/4+e/a
    "Second find X so P2=AX^3+BX^2+C^X+D=0"
    A = 8
    B = -4*p
    C = -8*r
    D = 4*r*p-q*q
    y0, y1, y2 = cardan(A, B, C, D)
    if abs(y1.imag) < abs(y0.imag):
        y0 = y1
    if abs(y2.imag) < abs(y0.imag):
        y0 = y2
    a0 = (-p+2*y0.real)**.5
    if a0 == 0:
        b0 = y0**2-r
    else:
        b0 = -q/2/a0
    r0, r1 = roots2(1, a0, y0+b0)
    r2, r3 = roots2(1, -a0, y0-b0)
    return (r0-z0, r1-z0, r2-z0, r3-z0)


J = exp(2j*pi/3)
Jc = 1/J


def cardan(a, b, c, d):
    u = np.empty(2, np.complex128)
    z0 = b/3/a
    a2, b2 = a*a, b*b
    p = -b2/3/a2 + c/a
    q = (b/27*(2*b2/a2-9*c/a)+d)/a
    D = -4*p*p*p-27*q*q
    r = sqrt(-D/27+0j)
    u = ((-q-r)/2)**0.33333333333333333333333
    v = ((-q+r)/2)**0.33333333333333333333333
    w = u*v
    w0 = abs(w+p/3)
    w1 = abs(w*J+p/3)
    w2 = abs(w*Jc+p/3)
    if w0 < w1:
        if w2 < w0:
            v *= Jc
    elif w2 < w1:
        v *= Jc
    else:
        v *= J
    return u+v-z0, u*J+v*Jc-z0, u*Jc+v*J-z0


def roots2(a, b, c):
    bp = b/2
    delta = bp*bp-a*c
    u1 = (-bp-delta**.5)/a
    u2 = -u1-b/a
    return u1, u2


print(" Quartic equation   ax4+bx3+cx2+dx+e=0")
print(ferrari(1, -7, 5, 31, -30))
