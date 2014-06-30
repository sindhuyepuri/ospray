#undef NDEBUG

// ospray
#include "raycast.h"
#include "ospray/camera/perspectivecamera.h"
// embree
#include "common/sys/sync/atomic.h"
// ispc exports
#include "raycast_ispc.h"

namespace ospray {

  template<void *(*ISPC_CREATE)(void*,void*,void*)>
  void RayCastRenderer<ISPC_CREATE>::commit()
  {
    if (getIE()) ispc::RayCastRenderer_destroy(getIE());

    model = (Model *)getParamObject("world",NULL);
    model = (Model *)getParamObject("model",model);
    camera = (Camera *)getParamObject("camera",NULL);

    ispcEquivalent = ISPC_CREATE(this->getIE(),model->getIE(),camera->getIE()); 
    // switch(SHADE_MODE) {
    // case RC_EYELIGHT:
    //   ispcEquivalent = ispc::RayCastRenderer_create_eyeLight(this->getIE(),model->getIE(),camera->getIE()); 
    //   break;
    // case RC_EYELIGHT_GEOMID:
    //   ispcEquivalent = ispc::RayCastRenderer_create_eyeLight_geomID(this->getIE(),model->getIE(),camera->getIE()); 
    //   break;
    // default:
    //   throw std::runtime_error("shade mode not implemented...");
    // };
  }

#if 0
  typedef RayCastRenderer<RC_EYELIGHT> RayCastRenderer_EyeLight;
  typedef RayCastRenderer<RC_PRIMID>   RayCastRenderer_PrimID;
  typedef RayCastRenderer<RC_GEOMID>   RayCastRenderer_GeomID;
  typedef RayCastRenderer<RC_INSTID>   RayCastRenderer_InstID;
  typedef RayCastRenderer<RC_EYELIGHT_PRIMID>   RayCastRenderer_PrimID;
  typedef RayCastRenderer<RC_EYELIGHT_GEOMID>   RayCastRenderer_GeomID;
  typedef RayCastRenderer<RC_EYELIGHT_INSTID>   RayCastRenderer_InstID;
  typedef RayCastRenderer<RC_GNORMAL>  RayCastRenderer_Ng;
  typedef RayCastRenderer<RC_TESTSHADOW>  RayCastRenderer_Shadow;

  OSP_REGISTER_RENDERER(RayCastRenderer_EyeLight,raycast);
  OSP_REGISTER_RENDERER(RayCastRenderer_EyeLight,raycast_eyelight);
  OSP_REGISTER_RENDERER(RayCastRenderer_PrimID,  raycast_eyelight_primID);
  OSP_REGISTER_RENDERER(RayCastRenderer_GeomID,  raycast_eyelight_geomID);
  OSP_REGISTER_RENDERER(RayCastRenderer_InstID,  raycast_eyelight_instID);
  OSP_REGISTER_RENDERER(RayCastRenderer_PrimID,  raycast_primID);
  OSP_REGISTER_RENDERER(RayCastRenderer_GeomID,  raycast_geomID);
  OSP_REGISTER_RENDERER(RayCastRenderer_InstID,  raycast_instID);
  OSP_REGISTER_RENDERER(RayCastRenderer_Ng,      raycast_Ng);
  OSP_REGISTER_RENDERER(RayCastRenderer_Shadow,  raycast_shadow);
#endif

  typedef RayCastRenderer<ispc::RayCastRenderer_create_eyeLight> RayCastRenderer_eyeLight;
  OSP_REGISTER_RENDERER(RayCastRenderer_eyeLight,raycast);
  OSP_REGISTER_RENDERER(RayCastRenderer_eyeLight,eyeLight);

  typedef RayCastRenderer<ispc::RayCastRenderer_create_eyeLight_primID>
  RayCastRenderer_eyeLight_primID;
  OSP_REGISTER_RENDERER(RayCastRenderer_eyeLight,eyeLight_primID);
  OSP_REGISTER_RENDERER(RayCastRenderer_eyeLight,primID);

  typedef RayCastRenderer<ispc::RayCastRenderer_create_eyeLight_instID>
  RayCastRenderer_eyeLight_instID;
  OSP_REGISTER_RENDERER(RayCastRenderer_eyeLight,eyeLight_instID);
  OSP_REGISTER_RENDERER(RayCastRenderer_eyeLight,instID);

  typedef RayCastRenderer<ispc::RayCastRenderer_create_eyeLight_geomID>
  RayCastRenderer_eyeLight_geomID;
  OSP_REGISTER_RENDERER(RayCastRenderer_eyeLight,eyeLight_geomID);
  OSP_REGISTER_RENDERER(RayCastRenderer_eyeLight,geomID);

};

