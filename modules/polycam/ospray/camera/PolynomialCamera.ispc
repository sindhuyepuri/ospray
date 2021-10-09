// Copyright 2009-2019 Intel Corporation
// SPDX-License-Identifier: Apache-2.0

#include "PolynomialCamera.ih"

void PolynomialCamera_initRay(const Camera *uniform _self,
    varying Ray &ray,
    const varying CameraSample &sample)
{
  const PolynomialCamera *uniform self =
      (const PolynomialCamera *uniform)_self;

  const vec2f screen = Camera_subRegion(_self, sample.screen);

  const vec3f dir = self->dir;
  const vec3f org =
      self->pos_00 + screen.x * self->pos_du + screen.y * self->pos_dv;

  const float time = Camera_shutterTime(_self, sample.time);
  setRay(ray, org, dir, self->super.nearClip, inf, time);
}

// TODO: PT4fData (coef, e0, e1, e2, e3, e4)
// TODO: Poly4f (trunc_deg, PolyTerm terms, consolidated?)
//     TODO: PolyTerm (scalar factor, power of each term)
// TODO: Transform4f (Poly4f equations, trunc_deg)
// TODO: get_system
//      TODO: two_plane_5x
//          TODO: make_TPT4f
//          TODO: refract_spherical_5
//          TODO: propagate_5
// TODO: Optical Material

PT4fData* initPoly4(float coeff, uint32_t e0, uint32_t e1, uint32_t e2, uint32_t e3) {
    PT4fData newData;
    newData.coeff = coeff;
    newData.e0 = e0;
    newData.e1 = e1;
    newData.e2 = e2;
    newData.e3 = e3;
}

void propagate_5(PT4fData **polynomials, int *num_terms, float l) {
    
    num_terms[0] = 7;

    polynomials[0][0] = initPoly4(1, 
        1, 0, 0, 0);
    polynomials[0][1] = initPoly4(l, 
        0, 0, 1, 0);
    polynomials[0][2] = initPoly4(l/2., 
        0, 0, 3, 0);
    polynomials[0][3] = initPoly4(l/2., 
        0, 0, 1, 1);
    polynomials[0][4] = initPoly4(3*l/8.,
        0, 0, 1, 4);
    polynomials[0][5] = initPoly4(3*l/4.,
        0, 0, 3, 2);
    polynomials[0][6] = initPoly4(3*l/8.,
        0, 0, 5, 0);

    num_terms[1] = 7;

    polynomials[1][0] = initPoly4(1, 
        0, 1, 0, 0);
    polynomials[1][1] = initPoly4(l, 
        0, 0, 0, 1);
    polynomials[1][2] = initPoly4(l/2.,
        0, 0, 0, 3);
    polynomials[1][3] = initPoly4(l/2.,
        0, 0, 2, 1);
    polynomials[1][4] = initPoly4(3*l/8.,
        0, 0, 4, 1);
    polynomials[1][5] = initPoly4(3*l/4.,
        0, 0, 2, 3);
    polynomials[1][6] = initPoly4(3*l/8.,
        0, 0, 0, 5);

    num_terms[2] = 1;

    polynomials[2][0] = initPoly4(1,
        0, 0, 1, 0);

    num_terms[3] = 1;

    polynomials[3][9] = initPoly4(1,
        0, 0, 0, 1);
}

Transform4f refract_spherical_5(float R, float n1, float n2, int trunc = 5) {
    int num_terms[4];
    PT4fData *element[4];
    for (int i = 0; i < 4; ++i) {
        element[i] = new PT4fData[1024];
    }

    PT4fData **polynomials = element;

    float R2 = R*R;
    float R3 = R2*R;
    float R4 = R2*R2;
    float R5 = R4*R;
    float n12 = n1*n1;
    float n13 = n12*n1;
    float n14 = n13*n1;
    float n22 = n2*n2;
    float n23 = n22*n2;
    float n24 = n23*n2;
    // Maple output     // cf. Spherical5.mw
    
    // How does x' depend?

    num_terms[0] = 0;

    // Linear (paraxial) term:
    float x1_xCoeff    = 1;
    polynomials[0][num_terms[0]++] = initPoly4( x1_xCoeff,
                            1,0,0,0);  // x coefficient

    // Degree-3 terms:
    float x1_x3Coeff   = 0.5*(n2-n1)/(R2*n2);
    float x1_xy2Coeff  = x1_x3Coeff;
    float x1_x2dxCoeff = 0.5*(n2-n1)/(R*n2);
    float x1_y2dxCoeff = x1_x2dxCoeff;
    polynomials[0][num_terms[0]++] = initPoly4( x1_x3Coeff,
                            3,0,0,0);  // x^3 coefficient
    polynomials[0][num_terms[0]++] = initPoly4( x1_xy2Coeff,
                            1,2,0,0);  // x * y^2 coefficient
    polynomials[0][num_terms[0]++] = initPoly4( x1_x2dxCoeff,
                            2,0,1,0);  // x^2 * dx coefficient
    polynomials[0][num_terms[0]++] = initPoly4( x1_y2dxCoeff,
                            0,2,1,0);  // y^2 * dx coefficient

    // Degree-5 terms:
    polynomials[0][num_terms[0]++] = initPoly4(-0.25*(-n23+n13)/(n23*R),
                        0,2,1,2);
    polynomials[0][num_terms[0]++] = initPoly4(-0.25*(-n23+n13)/(n23*R),
                        0,2,3,0);
    polynomials[0][num_terms[0]++] = initPoly4(-0.5*(-n23+n1*n22-n12*n2+n13)/(n23*R2),
                        0,3,1,1);
    polynomials[0][num_terms[0]++] = initPoly4(-0.125*(-3*n23+5*n1*n22-4*n12*n2+2*n13)/(n23*R3),
                        0,4,1,0);
    polynomials[0][num_terms[0]++] = initPoly4(-0.25*n1*(-n22+n12)/(n23*R2),
                        1,2,0,2);
    polynomials[0][num_terms[0]++] = initPoly4(-0.25*(-2*n23+n1*n22-2*n12*n2+3*n13)/(n23*R2),
                        1,2,2,0);
    polynomials[0][num_terms[0]++] = initPoly4(-0.5*(-n2+n1)*(n12+n22)/(n23*R3),
                        1,3,0,1);
    polynomials[0][num_terms[0]++] = initPoly4(-0.125*(-n2+n1)*(-2*n1*n2+3*n22+2*n12)/(n23*R4),
                        1,4,0,0);
    polynomials[0][num_terms[0]++] = initPoly4(-0.25*(-n23+n13)/(n23*R),
                        2,0,1,2);
    polynomials[0][num_terms[0]++] = initPoly4(-0.25*(-n23+n13)/(n23*R),
                        2,0,3,0);
    polynomials[0][num_terms[0]++] = initPoly4(-0.5*(-n23+n1*n22-n12*n2+n13)/(n23*R2),
                        2,1,1,1);
    polynomials[0][num_terms[0]++] = initPoly4(-0.25*(-5*n23-6*n12*n2+4*n13+7*n1*n22)/(n23*R3),
                        2,2,1,0);
    polynomials[0][num_terms[0]++] = initPoly4(-0.25*n1*(-n22+n12)/(n23*R2),
                        3,0,0,2);
    polynomials[0][num_terms[0]++] = initPoly4(-0.25*(-2*n23+n1*n22-2*n12*n2+3*n13)/(n23*R2),
                        3,0,2,0);
    polynomials[0][num_terms[0]++] = initPoly4(-0.5*(-n2+n1)*(n12+n22)/(n23*R3),
                        3,1,0,1);
    polynomials[0][num_terms[0]++] = initPoly4(-0.25*(-n2+n1)*(-2*n1*n2+3*n22+2*n12)/(n23*R4),
                        3,2,0,0);
    polynomials[0][num_terms[0]++] = initPoly4(-0.125*(-7*n23+9*n1*n22-8*n12*n2+6*n13)/(n23*R3),
                        4,0,1,0);
    polynomials[0][num_terms[0]++] = initPoly4(-0.125*(-n2+n1)*(-2*n1*n2+3*n22+2*n12)/(n23*R4),
                        5,0,0,0);

    
    // How does y' depend?
    // Axial symmetry - just copy stuff from x'
    
    num_terms[1] = 5;

    // Linear term:
    float y1_yCoeff    = x1_xCoeff;
    polynomials[1][0] = initPoly4( y1_yCoeff,
                    0,1,0,0);  // y coefficient

    // Degree-3 terms:
    float y1_y3Coeff   = x1_x3Coeff;
    float y1_x2yCoeff  = x1_x3Coeff;
    float y1_y2dyCoeff = x1_x2dxCoeff;
    float y1_x2dyCoeff = x1_x2dxCoeff;
    polynomials[1][1] = initPoly4( y1_y3Coeff,
                    0,3,0,0);  // y^3 coefficient
    polynomials[1][2] = initPoly4( y1_x2yCoeff,
                    2,1,0,0);  // x^2 * y coefficient
    polynomials[1][3] = initPoly4( y1_y2dyCoeff,
                    0,2,0,1);  // y^2 * dy coefficient
    polynomials[1][4] = initPoly4( y1_x2dyCoeff,
                    2,0,0,1);  // x^2 * dy coefficient

    // Degree-5 Terms:
    polynomials[1][num_terms[1]++] = initPoly4(-0.25*(-n23+n13)/(n23*R),
                        2,0,2,1);
    polynomials[1][num_terms[1]++] = initPoly4(-0.25*(-n23+n13)/(n23*R),
                        2,0,0,3);
    polynomials[1][num_terms[1]++] = initPoly4(-0.5*(-n23+n1*n22-n12*n2+n13)/(n23*R2),
                        3,0,1,1);
    polynomials[1][num_terms[1]++] = initPoly4(-0.125*(-3*n23+5*n1*n22-4*n12*n2+2*n13)/(n23*R3),
                        4,0,0,1);
    polynomials[1][num_terms[1]++] = initPoly4(-0.25*n1*(-n22+n12)/(n23*R2),
                        2,1,2,0);
    polynomials[1][num_terms[1]++] = initPoly4(-0.25*(-2*n23+n1*n22-2*n12*n2+3*n13)/(n23*R2),
                        2,1,0,2);
    polynomials[1][num_terms[1]++] = initPoly4(-0.5*(-n2+n1)*(n12+n22)/(n23*R3),
                        3,1,1,0);
    polynomials[1][num_terms[1]++] = initPoly4(-0.125*(-n2+n1)*(-2*n1*n2+3*n22+2*n12)/(n23*R4),
                        4,1,0,0);
    polynomials[1][num_terms[1]++] = initPoly4(-0.25*(-n23+n13)/(n23*R),
                        0,2,2,1);
    polynomials[1][num_terms[1]++] = initPoly4(-0.25*(-n23+n13)/(n23*R),
                        0,2,0,3);
    polynomials[1][num_terms[1]++] = initPoly4(-0.5*(-n23+n1*n22-n12*n2+n13)/(n23*R2),
                        1,2,1,1);
    polynomials[1][num_terms[1]++] = initPoly4(-0.25*(-5*n23-6*n12*n2+4*n13+7*n1*n22)/(n23*R3),
                        2,2,0,1);
    polynomials[1][num_terms[1]++] = initPoly4(-0.25*n1*(-n22+n12)/(n23*R2),
                        0,3,2,0);
    polynomials[1][num_terms[1]++] = initPoly4(-0.25*(-2*n23+n1*n22-2*n12*n2+3*n13)/(n23*R2),
                        0,3,0,2);
    polynomials[1][num_terms[1]++] = initPoly4(-0.5*(-n2+n1)*(n12+n22)/(n23*R3),
                        1,3,1,0);
    polynomials[1][num_terms[1]++] = initPoly4(-0.25*(-n2+n1)*(-2*n1*n2+3*n22+2*n12)/(n23*R4),
                        2,3,0,0);
    polynomials[1][num_terms[1]++] = initPoly4(-0.125*(-7*n23+9*n1*n22-8*n12*n2+6*n13)/(n23*R3),
                        0,4,0,1);
    polynomials[1][num_terms[1]++] = initPoly4(-0.125*(-n2+n1)*(-2*n1*n2+3*n22+2*n12)/(n23*R4),
                        0,5,0,0);
    
    // How does dx' depend?
    
    num_terms[2] = 9;

    // Linear terms (paraxial):
    float dx1_dxCoeff   = n1/n2;
    float dx1_xCoeff    = (n1-n2)/(R*n2);
    polynomials[2][0] = initPoly4( dx1_dxCoeff,
                    0,0,1,0);  // dx coefficient
    polynomials[2][1] = initPoly4( dx1_xCoeff,
                    1,0,0,0);  // x coefficient

    // Degree-3 terms:
    float dx1_x3Coeff   = 0.5*n1*(n1-n2)/(n22*R3);
    float dx1_x2dxCoeff = 0.5*(n1-n2)*(2*n1+n2)/(n22*R2);
    float dx1_xy2Coeff  = dx1_x3Coeff;
    float dx1_xydyCoeff = n1*(n1-n2)/(n22*R2);
    float dx1_y2dxCoeff = 0.5*(n1-n2)/(R2*n2);
    float dx1_xdx2Coeff = 0.5*n1*(n1-n2)/(R*n22);
    float dx1_xdy2Coeff = dx1_xdx2Coeff;
    polynomials[2][2] = initPoly4( dx1_x3Coeff,
                    3,0,0,0);  // x^3 coefficient
    polynomials[2][3] = initPoly4( dx1_x2dxCoeff,
                    2,0,1,0);  // x^2 * dx coefficient
    polynomials[2][4] = initPoly4( dx1_xy2Coeff,
                    1,2,0,0);  // x * y^2 coefficient
    polynomials[2][5] = initPoly4( dx1_xydyCoeff,
                    1,1,0,1);  // x * y * dy coefficient
    polynomials[2][6] = initPoly4( dx1_y2dxCoeff,
                    0,2,1,0);  // y^2 * dx coefficient
    polynomials[2][7] = initPoly4( dx1_xdx2Coeff,
                    1,0,2,0);  // x * dx^2 coefficient
    polynomials[2][8] = initPoly4( dx1_xdy2Coeff,
                    1,0,0,2);  // x * dy^2 coefficient


    // Degree-5 Terms:

    polynomials[2][num_terms[2]++] = initPoly4(0.25*(-n22+n12)/(n22*R2),
                        0,2,1,2);
    polynomials[2][num_terms[2]++] = initPoly4(0.25*(-n22+n12)/(n22*R2),
                        0,2,3,0);
    polynomials[2][num_terms[2]++] = initPoly4(0.5*(-n22+n12)/(n22*R3),
                        0,3,1,1);
    polynomials[2][num_terms[2]++] = initPoly4(0.125*(-n2+n1)*(2*n1+n2)/(n22*R4),
                        0,4,1,0);
    polynomials[2][num_terms[2]++] = initPoly4(0.125*n1*(-n23+n13)/(R*n24),
                        1,0,0,4);
    polynomials[2][num_terms[2]++] = initPoly4(0.25*n1*(-n23+n13)/(R*n24),
                        1,0,2,2);
    polynomials[2][num_terms[2]++] = initPoly4(0.125*n1*(-n23+n13)/(R*n24),
                        1,0,4,0);
    polynomials[2][num_terms[2]++] = initPoly4(0.5*n12*(-n22+n12)/(R2*n24),
                        1,1,0,3);
    polynomials[2][num_terms[2]++] = initPoly4(0.5*n12*(-n22+n12)/(R2*n24),
                        1,1,2,1);
    polynomials[2][num_terms[2]++] = initPoly4(0.25*n1*(-n23+3*n13-2*n1*n22)/(R3*n24),
                        1,2,0,2);
    polynomials[2][num_terms[2]++] = initPoly4(0.25*(-n1*n23+2*n12*n22+n14-2*n24)/(R3*n24),
                        1,2,2,0);
    polynomials[2][num_terms[2]++] = initPoly4(0.5*n1*(-n23+n13)/(R4*n24),
                        1,3,0,1);
    polynomials[2][num_terms[2]++] = initPoly4(0.125*n1*(-n23+n13)/(R5*n24),
                        1,4,0,0);
    polynomials[2][num_terms[2]++] = initPoly4(0.25*(-n12*n22+2*n14-n24)/(R2*n24),
                        2,0,1,2);
    polynomials[2][num_terms[2]++] = initPoly4(0.25*(-n12*n22+2*n14-n24)/(R2*n24),
                        2,0,3,0);
    polynomials[2][num_terms[2]++] = initPoly4(0.5*(-n12*n22+2*n14-n24)/(R3*n24),
                        2,1,1,1);
    polynomials[2][num_terms[2]++] = initPoly4(0.25*(-3*n1*n23+2*n12*n22+2*n14-n24)/(R4*n24),
                        2,2,1,0);
    polynomials[2][num_terms[2]++] = initPoly4(0.25*n1*(-n23+n13)/(R3*n24),
                        3,0,0,2);
    polynomials[2][num_terms[2]++] = initPoly4(0.25*(-n1*n23+3*n14-2*n24)/(R3*n24),
                        3,0,2,0);
    polynomials[2][num_terms[2]++] = initPoly4(0.5*n1*(-n23+n13)/(R4*n24),
                        3,1,0,1);
    polynomials[2][num_terms[2]++] = initPoly4(0.25*n1*(-n23+n13)/(R5*n24),
                        3,2,0,0);
    polynomials[2][num_terms[2]++] = initPoly4(0.125*(-5*n1*n23+4*n14+2*n12*n22-n24)/(R4*n24),
                        4,0,1,0);
    polynomials[2][num_terms[2]++] = initPoly4(0.125*n1*(-n23+n13)/(R5*n24),
                        5,0,0,0);


    // How does dy' depend?


    num_terms[3] = 9;

    // Linear terms (paraxial):
    float dy1_dyCoeff   = n1/n2;
    float dy1_yCoeff    = (n1-n2)/(R*n2);
    polynomials[3][0] = initPoly4( dy1_dyCoeff,
                    0,0,0,1);  // dy coefficient
    polynomials[3][1] = initPoly4( dy1_yCoeff,
                    0,1,0,0);  // y coefficient

    // Degree-3 terms:
    float dy1_y3Coeff   = 0.5*n1*(n1-n2)/(n22*R3);
    float dy1_y2dyCoeff = 0.5*(n1-n2)*(2*n1+n2)/(n22*R2);
    float dy1_x2yCoeff  = dy1_y3Coeff;
    float dy1_xydxCoeff = n1*(n1-n2)/(n22*R2);
    float dy1_x2dyCoeff = 0.5*(n1-n2)/(R2*n2);
    float dy1_ydy2Coeff = 0.5*n1*(n1-n2)/(R*n22);
    float dy1_ydx2Coeff = dy1_ydy2Coeff;
    polynomials[3][2] = initPoly4( dy1_y3Coeff,
                    0,3,0,0);  // y^3 coefficient
    polynomials[3][3] = initPoly4( dy1_y2dyCoeff,
                    0,2,0,1);  // y^2 * dy coefficient
    polynomials[3][4] = initPoly4( dy1_x2yCoeff,
                    2,1,0,0);  // x^2 * y coefficient
    polynomials[3][5] = initPoly4( dy1_xydxCoeff,
                    1,1,1,0);  // x * y * dx coefficient
    polynomials[3][6] = initPoly4( dy1_x2dyCoeff,
                    2,0,0,1);  // x^2 * dy coefficient
    polynomials[3][7] = initPoly4( dy1_ydy2Coeff,
                    0,1,0,2);  // y * dy^2 coefficient
    polynomials[3][8] = initPoly4( dy1_ydx2Coeff,
                    0,1,2,0);  // y * dx^2 coefficient

    // Degree-5 Terms:

    polynomials[3][num_terms[3]++] = initPoly4(0.25*(-n22+n12)/(n22*R2),
                        2,0,2,1);
    polynomials[3][num_terms[3]++] = initPoly4(0.25*(-n22+n12)/(n22*R2),
                        2,0,0,3);
    polynomials[3][num_terms[3]++] = initPoly4(0.5*(-n22+n12)/(n22*R3),
                        3,0,1,1);
    polynomials[3][num_terms[3]++] = initPoly4(0.125*(-n2+n1)*(2*n1+n2)/(n22*R4),
                        4,0,0,1);
    polynomials[3][num_terms[3]++] = initPoly4(0.125*n1*(-n23+n13)/(R*n24),
                        0,1,4,0);
    polynomials[3][num_terms[3]++] = initPoly4(0.25*n1*(-n23+n13)/(R*n24),
                        0,1,2,2);
    polynomials[3][num_terms[3]++] = initPoly4(0.125*n1*(-n23+n13)/(R*n24),
                        0,1,0,4);
    polynomials[3][num_terms[3]++] = initPoly4(0.5*n12*(-n22+n12)/(R2*n24),
                        1,1,3,0);
    polynomials[3][num_terms[3]++] = initPoly4(0.5*n12*(-n22+n12)/(R2*n24),
                        1,1,1,2);
    polynomials[3][num_terms[3]++] = initPoly4(0.25*n1*(-n23+3*n13-2*n1*n22)/(R3*n24),
                        2,1,2,0);
    polynomials[3][num_terms[3]++] = initPoly4(0.25*(-n1*n23+2*n12*n22+n14-2*n24)/(R3*n24),
                        2,1,0,2);
    polynomials[3][num_terms[3]++] = initPoly4(0.5*n1*(-n23+n13)/(R4*n24),
                        3,1,1,0);
    polynomials[3][num_terms[3]++] = initPoly4(0.125*n1*(-n23+n13)/(R5*n24),
                        4,1,0,0);
    polynomials[3][num_terms[3]++] = initPoly4(0.25*(-n12*n22+2*n14-n24)/(R2*n24),
                        0,2,2,1);
    polynomials[3][num_terms[3]++] = initPoly4(0.25*(-n12*n22+2*n14-n24)/(R2*n24),
                        0,2,0,3);
    polynomials[3][num_terms[3]++] = initPoly4(0.5*(-n12*n22+2*n14-n24)/(R3*n24),
                        1,2,1,1);
    polynomials[3][num_terms[3]++] = initPoly4(0.25*(-3*n1*n23+2*n12*n22+2*n14-n24)/(R4*n24),
                        2,2,0,1);
    polynomials[3][num_terms[3]++] = initPoly4(0.25*n1*(-n23+n13)/(R3*n24),
                        0,3,2,0);
    polynomials[3][num_terms[3]++] = initPoly4(0.25*(-n1*n23+3*n14-2*n24)/(R3*n24),
                        0,3,0,2);
    polynomials[3][num_terms[3]++] = initPoly4(0.5*n1*(-n23+n13)/(R4*n24),
                        1,3,1,0);
    polynomials[3][num_terms[3]++] = initPoly4(0.25*n1*(-n23+n13)/(R5*n24),
                        2,3,0,0);
    polynomials[3][num_terms[3]++] = initPoly4(0.125*(-5*n1*n23+4*n14+2*n12*n22-n24)/(R4*n24),
                        0,4,0,1);
    polynomials[3][num_terms[3]++] = initPoly4(0.125*n1*(-n23+n13)/(R5*n24),
                        0,5,0,0);
    Transform4f result; // = make_TPT4f(element, num_terms, trunc);
    result.trunc_deg = trunc;
    for (int i = 0; i < 4; i++) {
        g
    }
    for (int i = 0; i < 4; ++i) {
        delete[] element[i];
    }
    return result;
}

/*! create a new ispc-side version of a orthographiccamera - with given
    pointer to the c-class version - and return both class pointer and
    pointer to internal data back via the supplied reference
    parameters */
export void *uniform PolynomialCamera_create(void *uniform cppE)
{
  uniform PolynomialCamera *uniform cam =
      uniform new uniform PolynomialCamera;
  cam->super.cppEquivalent = cppE;
  cam->super.initRay = PolynomialCamera_initRay;
  cam->super.doesDOF = false;
  return cam;
}

export void PolynomialCamera_set(void *uniform _self,
    const uniform vec3f &dir,
    const uniform vec3f &pos_00,
    const uniform vec3f &pos_du,
    const uniform vec3f &pos_dv)
{
  uniform PolynomialCamera *uniform self =
      (uniform PolynomialCamera * uniform) _self;
  self->dir = dir;
  self->pos_00 = pos_00;
  self->pos_du = pos_du;
  self->pos_dv = pos_dv;
}
