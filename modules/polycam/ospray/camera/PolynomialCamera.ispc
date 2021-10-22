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

PT4fData pt4fdataConstructor(float c, float exp0, float exp1, float exp2, float exp3) {
    PT4fData newStruct;
    newStruct.coef = c;
    newStruct.e0 = exp0;
    newStruct.e1 = exp1;
    newStruct.e2 = exp2;
    newStruct.e3 = exp3;
    return newStruct;
}

float** two_plane_5(float dz) {
    // int num_terms[4];
    // float* polynomials[4];
    // // for (int i = 0; i < 4; i++) {
    //     polynomials[i] = new float[5][1024];
    // }
    // num_terms[0] = 1;
    // polynomials[0][0] = 1;
    // polynomials[0][1] = 0;
    // polynomials[0][2] = 0;
    // polynomials[0][3] = 1;
    // polynomials[0][4] = 0;
    // num_terms[1] = 1;
    // polynomials[1][0] = 1;
    // polynomials[1][1] = 0;
    // polynomials[1][2] = 0;
    // polynomials[1][3] = 0;
    // polynomials[1][4] = 1;

    // num_terms[2] = 0;
    // polynomials[2][num_terms[2]++] =
    int num_terms[4];
    PT4fData *polynomials[4];
    for (int i = 0; i < 4; ++i) {
        polynomials[i] = new PT4fData[1024];
    }
    
    num_terms[0] = 1;
    polynomials[0][0] = pt4fdataConstructor(1, 0, 0, 1, 0);
    num_terms[1] = 1;
    polynomials[1][0] = pt4fdataConstructor(1, 0, 0, 0, 1);

    num_terms[2] = 0;
    polynomials[2][num_terms[2]++] = pt4fdataConstructor(1/dz)
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
