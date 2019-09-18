// ======================================================================== //
// Copyright 2009-2019 Intel Corporation                                    //
//                                                                          //
// Licensed under the Apache License, Version 2.0 (the "License");          //
// you may not use this file except in compliance with the License.         //
// You may obtain a copy of the License at                                  //
//                                                                          //
//     http://www.apache.org/licenses/LICENSE-2.0                           //
//                                                                          //
// Unless required by applicable law or agreed to in writing, software      //
// distributed under the License is distributed on an "AS IS" BASIS,        //
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. //
// See the License for the specific language governing permissions and      //
// limitations under the License.                                           //
// ======================================================================== //

#include "TileOperation.h"
#include "DistributedFrameBuffer.h"
#include "DistributedFrameBuffer_ispc.h"

namespace ospray {

using DFB = DistributedFrameBuffer;

void TileOperation::attach(DistributedFrameBuffer *dfb) {}

TileDesc::TileDesc(const vec2i &begin, size_t tileID, size_t ownerID)
    : begin(begin), tileID(tileID), ownerID(ownerID)
{}

LiveTileOperation::LiveTileOperation(
    DFB *dfb, const vec2i &begin, size_t tileID, size_t ownerID)
    : TileDesc(begin, tileID, ownerID), dfb(dfb)
{}

void LiveTileOperation::accumulate(const ospray::Tile &tile)
{
  // accumulate, compute the finished normalized colors, and compute the
  // error, but don't write the finished colors into the color buffer yet
  // because pixel ops might modify them later Note: also used for FB_NONE
  error = DFB_accumulateTile((const ispc::VaryingTile *)&tile,
      (ispc::VaryingTile *)&finished,
      (ispc::VaryingTile *)&accum,
      (ispc::VaryingTile *)&variance,
      dfb->hasAccumBuffer,
      dfb->hasVarianceBuffer);

  if (dfb->hasNormalBuffer || dfb->hasAlbedoBuffer) {
    ispc::DFB_accumulateAuxTile((const ispc::VaryingTile *)&tile,
        (ispc::Tile *)&finished,
        (ispc::VaryingTile *)&accum);
  }
}

void LiveTileOperation::tileIsFinished()
{
  dfb->tileIsFinished(this);
}

} // namespace ospray
