
#pragma once

#include "fb/ImageOp.h"
#include "ospray_module_polyparallel_export.h"
#include "rkcommon/tasking/parallel_for.h"

namespace ospray {

struct OSPRAY_MODULE_POLYPARALLEL_EXPORT PolyParallelFrameOp : public FrameOp
{
    PolyParallelFrameOp();

    ~PolyParallelFrameOp() override;

    std::unique_ptr<LiveImageOp> attach(FrameBufferView &fbView) override;

    std::string toString() const override;
        
}; 

}
