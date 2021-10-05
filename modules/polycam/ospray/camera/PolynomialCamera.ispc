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

PT4fData* initPoly4(float coeff, uint32_t e0, uint32_t e1, uint32_t e2, uint32_t e3, uint32_t e4) {
    PT4fData newData;
    newData.coeff = coeff;
    newData.e0 = e0;
    newData.e1 = e1;
    newData.e2 = e2;
    newData.e3 = e3;
    newData.e4 = e4;
}

void propagate_5(PT4fData **polynomials, int *num_terms, float l) {
    
    num_terms[0] = 7;

    polynomials[0][0] = initPoly4(1, 
        1, 0, 0, 0);
    polynomials[0][1] = initPoly4(l, 
        0, 0, 1, 0);
    polynomials[0][2] = initPoly4(l/2.f, 
        0, 0, 3, 0);
    polynomials[0][3] = initPoly4(l/2.f, 
        0, 0, 1, 1);
    polynomials[0][4] = initPoly4(3*l/8.f,
        0, 0, 1, 4);
    polynomials[0][5] = initPoly4(3*l/4.f,
        0, 0, 3, 2);
    polynomials[0][6] = initPoly4(3*l/8.f,
        0, 0, 5, 0);

    num_terms[1] = 7;

    polynomials[1][0] = initPoly4(1, 
        0, 1, 0, 0);
    polynomials[1][1] = initPoly4(l, 
        0, 0, 0, 1);
    polynomials[1][2] = initPoly4(l/2.f,
        0, 0, 0, 3);
    polynomials[1][3] = initPoly4(l/2.f,
        0, 0, 2, 1);
    polynomials[1][4] = initPoly4(3*l/8.f,
        0, 0, 4, 1);
    polynomials[1][5] = initPoly4(3*l/4.f,
        0, 0, 2, 3);
    polynomials[1][6] = initPoly4(3*l/8.f,
        0, 0, 0, 5);

    num_terms[2] = 1;

    polynomials[2][0] = initPoly4(1,
        0, 0, 1, 0);

    num_terms[3] = 1;

    polynomials[3][9] = initPoly4(1,
        0, 0, 0, 1);
}

Transform4f refract_spherical_5(float R, float n1, float n2, int trunc = 5) {
    
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
