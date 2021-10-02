
#pragma once

#include "fb/ImageOp.h"
#include "ospray_module_polyfilter_export.h"
#include "rkcommon/tasking/parallel_for.h"

namespace ospray {

struct OSPRAY_MODULE_POLYFILTER_EXPORT PolyFilterFrameOp : public FrameOp
{
    PolyFilterFrameOp();

    ~PolyFilterFrameOp() override;

    std::unique_ptr<LiveImageOp> attach(FrameBufferView &fbView) override;

    std::string toString() const override;
        
}; 

}
