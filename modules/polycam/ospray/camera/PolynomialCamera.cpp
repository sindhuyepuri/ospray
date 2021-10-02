// Copyright 2009-2020 Intel Corporation
// SPDX-License-Identifier: Apache-2.0

#include "PolynomialCamera.h"
#include "camera/PolynomialCamera_ispc.h"


#include "TruncPoly/TruncPolySystem.hh"

#include "OpticalElements/OpticalMaterial.hh"
#include "OpticalElements/Spherical5.hh"
#include "OpticalElements/Cylindrical5.hh"
#include "OpticalElements/Propagation5.hh"
#include "OpticalElements/TwoPlane5.hh"

#include "OpticalElements/FindFocus.hh"

#include <iostream>

namespace ospray {

PolynomialCamera::PolynomialCamera()
{
  ispcEquivalent = ispc::PolynomialCamera_create(this);
}

std::string PolynomialCamera::toString() const
{
  return "ospray::PolynomialCamera";
}

void PolynomialCamera::commit()
{
  Camera::commit();

  // ------------------------------------------------------------------
  // first, "parse" the additional expected parameters
  // ------------------------------------------------------------------
  height = getParam<float>("height", 1.f); // imgPlane_size_y
  aspect = getParam<float>("aspect", 1.f);

  // ------------------------------------------------------------------
  // now, update the local precomputed values
  // ------------------------------------------------------------------
  dir = normalize(dir);
  vec3f pos_du = normalize(cross(dir, up));
  vec3f pos_dv = cross(pos_du, dir);

  pos_du *= height * aspect; // imgPlane_size_x
  pos_dv *= height;

  vec3f pos_00 = pos - 0.5f * pos_du - 0.5f * pos_dv;

  ispc::PolynomialCamera_set(getIE(),
      (const ispc::vec3f &)dir,
      (const ispc::vec3f &)pos_00,
      (const ispc::vec3f &)pos_du,
      (const ispc::vec3f &)pos_dv);
}

extern "C" OSPError OSPRAY_DLLEXPORT ospray_module_init_polycam(
    int16_t versionMajor, int16_t versionMinor, int16_t)
{
    auto status = moduleVersionCheck(versionMajor, versionMinor);
    if (status == OSP_NO_ERROR) {
        ospray::Camera::registerType<PolynomialCamera>("polycam");
    }
    return status;
}

} // namespace ospray

