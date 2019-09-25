// ======================================================================== //
// Copyright 2017-2019 Intel Corporation                                    //
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

#include "ospray_test_fixture.h"

extern OSPRayEnvironment *ospEnv;

namespace OSPRayTestScenes {

  // Helper functions /////////////////////////////////////////////////////////

  // creates a torus
  // volumetric data: stores data of torus
  // returns created ospvolume of torus
  static OSPVolume CreateTorus(std::vector<float> &volumetricData,
                               const int size)
  {
    volumetricData.resize(size * size * size, 0);

    const float r    = 30;
    const float R    = 80;
    const int size_2 = size / 2;

    for (int x = 0; x < size; ++x) {
      for (int y = 0; y < size; ++y) {
        for (int z = 0; z < size; ++z) {
          const float X = x - size_2;
          const float Y = y - size_2;
          const float Z = z - size_2;

          const float d = (R - std::sqrt(X * X + Y * Y));
          volumetricData[size * size * x + size * y + z] =
              r * r - d * d - Z * Z;
        }
      }
    }

    OSPVolume torus    = ospNewVolume("shared_structured_volume");
    OSPData voxelsData = ospNewData(size * size * size,
                                    OSP_FLOAT,
                                    volumetricData.data(),
                                    OSP_DATA_SHARED_BUFFER);
    ospSetObject(torus, "voxelData", voxelsData);
    ospRelease(voxelsData);
    ospSetVec3i(torus, "dimensions", size, size, size);
    ospSetInt(torus, "voxelType", OSP_FLOAT);
    ospSetVec2f(torus, "voxelRange", -10000.f, 10000.f);
    ospSetVec3f(torus, "gridOrigin", -0.5f, -0.5f, -0.5f);
    ospSetVec3f(torus, "gridSpacing", 1.f / size, 1.f / size, 1.f / size);
    return torus;
  }

  /////////////////////////////////////////////////////////////////////////////

  Base::Base()
  {
    const ::testing::TestCase *const testCase =
        ::testing::UnitTest::GetInstance()->current_test_case();
    const ::testing::TestInfo *const testInfo =
        ::testing::UnitTest::GetInstance()->current_test_info();
    imgSize = ospEnv->GetImgSize();

    {
      std::string testCaseName = testCase->name();
      std::string testInfoName = testInfo->name();
      size_t pos               = testCaseName.find('/');
      if (pos == std::string::npos) {
        testName = testCaseName + "_" + testInfoName;
      } else {
        std::string instantiationName = testCaseName.substr(0, pos);
        std::string className         = testCaseName.substr(pos + 1);
        testName = className + "_" + instantiationName + "_" + testInfoName;
      }
      for (char &byte : testName)
        if (byte == '/')
          byte = '_';
    }

    rendererType    = "scivis";
    frames          = 1;
    samplesPerPixel = 16;
    SetImageTool();
  }

  void Base::SetImageTool()
  {
    try {
      imageTool.reset(
          new OSPImageTools(imgSize, GetTestName(), frameBufferFormat));
    } catch (std::bad_alloc &) {
      FAIL() << "Failed to create image tool.\n";
    }
  }

  Base::~Base()
  {
    for (auto l : lightsList)
      ospRelease(l);

    for (auto i : instances)
      ospRelease(i);

    ospRelease(framebuffer);
    ospRelease(renderer);
    ospRelease(camera);
    ospRelease(world);
  }

  void Base::SetUp()
  {
    ASSERT_NO_FATAL_FAILURE(CreateEmptyScene());
  }

  void Base::AddLight(OSPLight light)
  {
    lightsList.push_back(light);
  }

  void Base::AddModel(OSPGeometricModel model, affine3f xfm)
  {
    ospCommit(model);
    OSPData data = ospNewData(1, OSP_GEOMETRIC_MODEL, &model);

    OSPGroup group = ospNewGroup();
    ospSetObject(group, "geometry", data);
    ospCommit(group);

    OSPInstance instance = ospNewInstance(group);
    ospSetParam(instance, "xfm", OSP_AFFINE3F, &xfm);
    ospCommit(instance);
    ospRelease(group);

    AddInstance(instance);

    ospRelease(model);
    ospRelease(data);
  }

  void Base::AddModel(OSPVolumetricModel model, affine3f xfm)
  {
    ospCommit(model);
    OSPData data = ospNewData(1, OSP_VOLUMETRIC_MODEL, &model);

    OSPGroup group = ospNewGroup();
    ospSetObject(group, "volume", data);
    ospCommit(group);

    OSPInstance instance = ospNewInstance(group);
    ospSetParam(instance, "xfm", OSP_AFFINE3F, &xfm);
    ospCommit(instance);
    ospRelease(group);

    AddInstance(instance);

    ospRelease(model);
    ospRelease(data);
  }

  void Base::AddInstance(OSPInstance instance)
  {
    ospCommit(instance);
    instances.push_back(instance);
  }

  void Base::PerformRenderTest()
  {
    SetLights();

    if (!instances.empty()) {
      OSPData insts =
          ospNewData(instances.size(), OSP_INSTANCE, instances.data());
      ospSetObject(world, "instance", insts);
      ospRelease(insts);
    }

    ospCommit(camera);
    ospCommit(world);
    ospCommit(renderer);

    ospResetAccumulation(framebuffer);

    RenderFrame();
    uint32_t *framebuffer_data =
        (uint32_t *)ospMapFrameBuffer(framebuffer, OSP_FB_COLOR);

    if (ospEnv->GetDumpImg()) {
      EXPECT_EQ(imageTool->saveTestImage(framebuffer_data), OsprayStatus::Ok);
    } else {
      EXPECT_EQ(imageTool->compareImgWithBaseline(framebuffer_data),
                OsprayStatus::Ok);
    }

    ospUnmapFrameBuffer(framebuffer_data, framebuffer);
  }

  void Base::CreateEmptyScene()
  {
    SetCamera();
    SetWorld();
    SetRenderer();
    SetFramebuffer();
  }

  void Base::SetCamera()
  {
    float cam_pos[]  = {0.f, 0.f, 0.f};
    float cam_up[]   = {0.f, 1.f, 0.f};
    float cam_view[] = {0.f, 0.f, 1.f};

    camera = ospNewCamera("perspective");

    ospSetFloat(camera, "aspect", imgSize.x / (float)imgSize.y);
    ospSetParam(camera, "position", OSP_VEC3F, cam_pos);
    ospSetParam(camera, "direction", OSP_VEC3F, cam_view);
    ospSetParam(camera, "up", OSP_VEC3F, cam_up);
  }

  void Base::SetWorld()
  {
    world = ospNewWorld();
  }

  void Base::SetLights()
  {
    if (!lightsList.empty()) {
      lights = ospNewData(lightsList.size(), OSP_LIGHT, lightsList.data());
      ospSetObject(world, "light", lights);
      ospRelease(lights);
    }
  }

  void Base::SetRenderer()
  {
    renderer = ospNewRenderer(rendererType.c_str());
    ospSetInt(renderer, "aoSamples", 0);
    ospSetFloat(renderer, "bgColor", 1.0f);
    ospSetObject(renderer, "model", world);
    ospSetObject(renderer, "camera", camera);
    ospSetInt(renderer, "spp", samplesPerPixel);
  }

  void Base::SetFramebuffer()
  {
    framebuffer = ospNewFrameBuffer(imgSize.x,
                                    imgSize.y,
                                    frameBufferFormat,
                                    OSP_FB_COLOR | OSP_FB_ACCUM | OSP_FB_DEPTH);
  }

  OSPMaterial Base::CreateMaterial(std::string type)
  {
    OSPMaterial material = ospNewMaterial(rendererType.c_str(), type.data());
    EXPECT_TRUE(material);
    if (type == "Glass") {
      ospSetFloat(material, "etaInside", 1.5);
      ospSetVec3f(material, "attenuationColor", 0.f, 1.f, 1.f);
      ospSetFloat(material, "attenuationDistance", 5.0);
    } else if (type == "Luminous") {
      ospSetFloat(material, "intensity", 3.0);
    }

    ospCommit(material);

    return material;
  }

  void Base::RenderFrame()
  {
    for (int frame = 0; frame < frames; ++frame)
      ospRenderFrameBlocking(framebuffer, renderer, camera, world);
  }

  SingleObject::SingleObject()
  {
    auto params  = GetParam();
    rendererType = std::get<0>(params);
    materialType = std::get<1>(params);
  }

  void SingleObject::SetUp()
  {
    Base::SetUp();

    SetMaterial();

    OSPLight distant = ospNewLight("distant");
    ASSERT_TRUE(distant) << "Failed to create lights";
    ospSetFloat(distant, "intensity", 1.0f);
    ospSetVec3f(distant, "direction", 1.0f, 1.0f, 1.0f);
    ospSetFloat(distant, "angularDiameter", 1.0f);
    ospCommit(distant);
    AddLight(distant);

    OSPLight ambient = ospNewLight("ambient");
    ASSERT_TRUE(ambient) << "Failed to create lights";
    ospSetFloat(ambient, "intensity", 0.2f);
    ospCommit(ambient);
    AddLight(ambient);
  }

  void SingleObject::SetMaterial()
  {
    material = CreateMaterial(materialType);
  }

  SpherePrecision::SpherePrecision()
  {
    auto params  = GetParam();
    radius       = std::get<0>(params);
    dist         = radius * std::get<1>(params);
    move_cam     = std::get<2>(params);
    rendererType = std::get<3>(params);
  }

  void SpherePrecision::SetUp()
  {
    Base::SetUp();

    float fov  = 160.0f * std::min(std::tan(radius / std::abs(dist)), 1.0f);
    float cent = move_cam ? 0.0f : dist + radius;

    ospSetVec3f(camera, "position", 0.f, 0.f, move_cam ? -dist - radius : 0.0f);
    ospSetVec3f(camera, "direction", 0.f, 0.f, 1.f);
    ospSetVec3f(camera, "up", 0.f, 1.f, 0.f);
    ospSetFloat(camera, "fovy", fov);

    ospSetInt(renderer, "spp", 16);
    ospSetVec4f(renderer, "bgColor", 0.2f, 0.2f, 0.4f, 1.0f);
    // scivis params
    ospSetInt(renderer, "aoSamples", 16);
    ospSetFloat(renderer, "aoIntensity", 1.f);
    // pathtracer params
    ospSetInt(renderer, "maxDepth", 2);

    OSPGeometry sphere      = ospNewGeometry("spheres");
    OSPGeometry inst_sphere = ospNewGeometry("spheres");

    float sph_center_r[]{-0.5f * radius,
        -0.3f * radius,
        cent,
        radius,
        0.8f * radius,
        -0.3f * radius,
        cent,
        0.9f * radius,
        0.8f * radius,
        1.5f * radius,
        cent,
        0.9f * radius,
        0.0f,
        -10001.3f * radius,
        cent,
        10000.f * radius,
        0.f,
        0.f,
        0.f,
        90.f * radius};
    OSPData positionSrc = ospNewSharedData(sph_center_r, OSP_VEC3F, 4, 16);
    OSPData positionData = ospNewData(OSP_VEC3F, 4);
    ospCopyData(positionSrc, positionData);
    ospRelease(positionSrc);
    OSPData radiusSrc = ospNewSharedData(sph_center_r + 3, OSP_FLOAT, 4, 16);
    OSPData radiusData = ospNewData(OSP_FLOAT, 4);
    ospCopyData(radiusSrc, radiusData);
    ospRelease(radiusSrc);
    ospCommit(positionData);
    ospCommit(radiusData);
    ospSetObject(sphere, "sphere.position", positionData);
    ospSetObject(sphere, "sphere.radius", radiusData);
    ospCommit(sphere);
    ospRelease(positionData);
    ospRelease(radiusData);

    positionSrc = ospNewSharedData(sph_center_r + 16, OSP_VEC3F, 1);
    positionData = ospNewData(OSP_VEC3F, 1);
    ospCopyData(positionSrc, positionData);
    ospRelease(positionSrc);
    radiusSrc = ospNewSharedData(sph_center_r + 19, OSP_FLOAT, 1);
    radiusData = ospNewData(OSP_FLOAT, 1);
    ospCopyData(radiusSrc, radiusData);
    ospRelease(radiusSrc);
    ospSetObject(inst_sphere, "sphere.position", positionData);
    ospSetObject(inst_sphere, "sphere.radius", radiusData);
    ospCommit(inst_sphere);
    ospRelease(positionData);
    ospRelease(radiusData);

    OSPGeometricModel model1 = ospNewGeometricModel(sphere);
    ospRelease(sphere);
    OSPGeometricModel model2 = ospNewGeometricModel(inst_sphere);
    ospRelease(inst_sphere);

    OSPMaterial sphereMaterial =
        ospNewMaterial(rendererType.c_str(), "default");
    ospSetFloat(sphereMaterial, "d", 1.0f);
    ospCommit(sphereMaterial);
    ospSetObject(model1, "material", sphereMaterial);
    ospSetObject(model2, "material", sphereMaterial);
    ospRelease(sphereMaterial);

    affine3f xfm(vec3f(0.01, 0, 0),
                 vec3f(0, 0.01, 0),
                 vec3f(0, 0, 0.01),
                 vec3f(-0.5f * radius, 1.6f * radius, cent));

    AddModel(model1);
    AddModel(model2, xfm);

    OSPLight distant = ospNewLight("distant");
    ASSERT_TRUE(distant) << "Failed to create lights";
    ospSetFloat(distant, "intensity", 3.0f);
    ospSetVec3f(distant, "direction", 0.3f, -4.0f, 0.8f);
    ospSetVec3f(distant, "color", 1.0f, 0.5f, 0.5f);
    ospSetFloat(distant, "angularDiameter", 1.0f);
    ospCommit(distant);
    AddLight(distant);

    OSPLight ambient = ospNewLight("ambient");
    ASSERT_TRUE(ambient) << "Failed to create lights";
    ospSetFloat(ambient, "intensity", 0.1f);
    ospCommit(ambient);
    AddLight(ambient);
  }

  Box::Box()
  {
    rendererType = "pathtracer";

    auto params = GetParam();

    cuboidMaterialType = std::get<0>(params);
    sphereMaterialType = std::get<1>(params);
  }

  void Box::SetUp()
  {
    Base::SetUp();

    SetMaterials();

    std::vector<OSPGeometricModel> models;
    OSPGroup group = ospNewGroup();

    float wallsVertices[] = {// left wall
                             1.f,
                             -1.f,
                             4.f,
                             1.f,
                             1.f,
                             4.f,
                             1.f,
                             1.f,
                             0.f,
                             1.f,
                             -1.f,
                             0.f,
                             // right wall
                             -1.f,
                             -1.f,
                             4.f,
                             -1.f,
                             1.f,
                             4.f,
                             -1.f,
                             1.f,
                             0.f,
                             -1.f,
                             -1.f,
                             0.f,
                             // bottom, back and top walls
                             -1.f,
                             -1.f,
                             0.f,
                             1.f,
                             -1.f,
                             0.f,
                             1.f,
                             -1.f,
                             4.f,
                             -1.f,
                             -1.f,
                             4.f,
                             -1.f,
                             1.f,
                             4.f,
                             1.f,
                             1.f,
                             4.f,
                             1.f,
                             1.f,
                             0.f,
                             -1.f,
                             1.f,
                             0.f};
    float wallsColors[]   = {
        1.f, 0.f, 0.f, 1.f, 1.f, 0.f, 0.f, 1.f, 1.f, 0.f, 0.f, 1.f, 1.f,
        0.f, 0.f, 1.f, 0.f, 1.f, 0.f, 1.f, 0.f, 1.f, 0.f, 1.f, 0.f, 1.f,
        0.f, 1.f, 0.f, 1.f, 0.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f,
        1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f,
        1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f};
    int32_t wallsIndices[] = {0,  1,  2,  2,  3,  0,  4,  5,  6,  6,
                              7,  4,  8,  9,  10, 10, 11, 8,  10, 11,
                              12, 12, 13, 10, 12, 13, 14, 14, 15, 12};
    OSPGeometry wallsMesh  = ospNewGeometry("triangles");
    OSPData data           = ospNewData(16, OSP_VEC3F, wallsVertices);
    ospCommit(data);
    ospSetObject(wallsMesh, "vertex.position", data);
    ospRelease(data);
    data = ospNewData(16, OSP_VEC4F, wallsColors);
    ospCommit(data);
    ospSetObject(wallsMesh, "vertex.color", data);
    ospRelease(data);
    data = ospNewData(10, OSP_VEC3UI, wallsIndices);
    ospCommit(data);
    ospSetObject(wallsMesh, "index", data);
    ospRelease(data);
    ospCommit(wallsMesh);

    OSPGeometricModel model = ospNewGeometricModel(wallsMesh);
    ospRelease(wallsMesh);
    OSPMaterial wallsMaterial = GetMaterial("OBJMaterial");
    ospSetObject(model, "material", wallsMaterial);
    ospRelease(wallsMaterial);
    models.push_back(model);

    float lightVertices[]   = {-0.3f,
                             0.999f,
                             2.7f,
                             0.3f,
                             0.999f,
                             2.7f,
                             0.3f,
                             0.999f,
                             3.3f,
                             -0.3f,
                             0.999f,
                             3.3f};
    int32_t lightIndices[]  = {0, 1, 2, 2, 3, 0};
    OSPGeometry lightSquare = ospNewGeometry("triangles");
    data                    = ospNewData(4, OSP_VEC3F, lightVertices);
    ospCommit(data);
    ospSetObject(lightSquare, "vertex.position", data);
    ospRelease(data);
    data = ospNewData(2, OSP_VEC3UI, lightIndices);
    ospSetObject(lightSquare, "index", data);
    ospRelease(data);
    OSPMaterial lightMaterial =
        ospNewMaterial(rendererType.c_str(), "Luminous");
    ospSetFloat(lightMaterial, "intensity", 20.f);
    ospSetVec3f(lightMaterial, "color", 1.f, 0.7f, 0.3f);
    ospCommit(lightMaterial);
    ospCommit(lightSquare);

    model = ospNewGeometricModel(lightSquare);
    ospRelease(lightSquare);
    ospSetObject(model, "material", lightMaterial);
    ospRelease(lightMaterial);
    models.push_back(model);

    float cuboidVertices[] = {
        0.7f, -0.999f, 3.3f,    0.3f, -0.999f, 3.6f, 0.0f, -0.999f,
        3.2f, 0.4f,    -0.999f, 2.9f, 0.4f,    0.2,  2.9f, 0.0f,
        0.2,  3.2f,    0.3f,    0.2,  3.6f,    0.7f, 0.2,  3.3f,
    };
    int32_t cuboidIndices[] = {0, 1, 2, 2, 3, 1, 4, 5, 6, 6, 7, 4,
                               0, 1, 6, 6, 7, 0, 1, 2, 5, 5, 6, 1,
                               2, 3, 4, 4, 5, 2, 3, 0, 7, 7, 4, 3};

    OSPGeometry cuboid = ospNewGeometry("triangles");
    data               = ospNewData(8, OSP_VEC3F, cuboidVertices);

    ospCommit(data);
    ospSetObject(cuboid, "vertex.position", data);
    ospRelease(data);
    data = ospNewData(12, OSP_VEC3UI, cuboidIndices);
    ospCommit(data);
    ospSetObject(cuboid, "index", data);
    ospRelease(data);
    ospCommit(cuboid);

    model = ospNewGeometricModel(cuboid);
    ospRelease(cuboid);
    ospSetObject(model, "material", cuboidMaterial);
    ospRelease(cuboidMaterial);
    models.push_back(model);

    float sphereVertex[] = {-0.3f, -0.55f, 2.5f};
    OSPGeometry sphere   = ospNewGeometry("spheres");
    data = ospNewData(1, OSP_VEC3F, sphereVertex);
    ospCommit(data);
    ospSetObject(sphere, "sphere.position", data);
    ospRelease(data);
    ospSetFloat(sphere, "radius", 0.45f);
    ospCommit(sphere);

    model = ospNewGeometricModel(sphere);
    ospRelease(sphere);
    ospSetObject(model, "material", sphereMaterial);
    ospRelease(sphereMaterial);
    models.push_back(model);

    for (auto &m : models)
      ospCommit(m);

    OSPData models_data =
        ospNewData(models.size(), OSP_GEOMETRIC_MODEL, models.data());
    ospSetObject(group, "geometry", models_data);
    ospCommit(group);
    ospRelease(models_data);

    OSPInstance instance = ospNewInstance(group);
    ospCommit(instance);
    ospRelease(group);

    AddInstance(instance);

    // NOTE(jda) - still need to set the world on the renderer for geom lights
    ospSetObject(renderer, "world", world);
  }

  void Box::SetMaterials()
  {
    cuboidMaterial = GetMaterial(cuboidMaterialType);
    sphereMaterial = GetMaterial(sphereMaterialType);
  }

  OSPMaterial Box::GetMaterial(std::string type)
  {
    OSPMaterial newMaterial = ospNewMaterial(rendererType.c_str(), type.data());
    // defaults for "OBJMaterial"
    if (type == "Glass") {
      ospSetFloat(newMaterial, "eta", 1.4);
    } else if (type == "Luminous") {
      ospSetFloat(newMaterial, "intensity", 0.7f);
    }
    ospCommit(newMaterial);

    return newMaterial;
  }

  Sierpinski::Sierpinski()
  {
    auto params      = GetParam();
    renderIsosurface = std::get<1>(params);
    level            = std::get<2>(params);

    rendererType = std::get<0>(params);
  }

  void Sierpinski::SetUp()
  {
    Base::SetUp();

    float cam_pos[]  = {-0.5f, -1.f, 0.2f};
    float cam_up[]   = {0.f, 0.f, -1.f};
    float cam_view[] = {0.5f, 1.f, 0.f};
    ospSetParam(camera, "position", OSP_VEC3F, cam_pos);
    ospSetParam(camera, "direction", OSP_VEC3F, cam_view);
    ospSetParam(camera, "up", OSP_VEC3F, cam_up);

    int size = 1 << level;

    volumetricData.resize(size * size * size, 0);

    auto getIndex = [size](int layer, int x, int y) {
      return size * size * layer + size * x + y;
    };

    volumetricData[getIndex(0, size - 1, size - 1)] = 255;
    for (int layer = 1; layer < size; ++layer) {
      for (int x = size - layer - 1; x < size; ++x) {
        for (int y = size - layer - 1; y < size; ++y) {
          int counter = 0;
          if (volumetricData[getIndex(layer - 1, x, y)]) {
            counter++;
          }
          if (volumetricData[getIndex(layer - 1, (x + 1) % size, y)]) {
            counter++;
          }
          if (volumetricData[getIndex(layer - 1, x, (y + 1) % size)]) {
            counter++;
          }

          volumetricData[getIndex(layer, x, y)] = (counter == 1) ? 255 : 0;
        }
      }
    }

    OSPVolume pyramid = ospNewVolume("shared_structured_volume");
    OSPData voxelsData =
        ospNewData(size * size * size, OSP_UCHAR, volumetricData.data());
    ospSetObject(pyramid, "voxelData", voxelsData);
    ospRelease(voxelsData);
    ospSetVec3i(pyramid, "dimensions", size, size, size);
    ospSetInt(pyramid, "voxelType", OSP_UCHAR);
    ospSetVec2f(pyramid, "voxelRange", 0, 255);
    ospSetVec3f(pyramid, "gridOrigin", -0.5f, -0.5f, -0.5f);
    ospSetVec3f(pyramid, "gridSpacing", 1.f / size, 1.f / size, 1.f / size);
    ospCommit(pyramid);

    OSPVolumetricModel model = ospNewVolumetricModel(pyramid);
    ospSetFloat(model, "samplingRate", 1.f);

    OSPTransferFunction transferFun =
        ospNewTransferFunction("piecewise_linear");
    ospSetVec2f(transferFun, "valueRange", 0, 255);
    float colors[]      = {1.0f, 1.0f, 1.0f, 0.0f, 0.0f, 0.0f};
    float opacites[]    = {0.f, 1.0f};
    OSPData tfColorData = ospNewData(2, OSP_VEC3F, colors);
    ospSetObject(transferFun, "color", tfColorData);
    ospRelease(tfColorData);
    OSPData tfOpacityData = ospNewData(2, OSP_FLOAT, opacites);
    ospSetObject(transferFun, "opacity", tfOpacityData);
    ospRelease(tfOpacityData);
    ospCommit(transferFun);
    ospSetObject(model, "transferFunction", transferFun);
    ospCommit(model);

    if (renderIsosurface) {
      OSPGeometry isosurface = ospNewGeometry("isosurfaces");
      ospSetObject(isosurface, "volume", model);
      float isovalue        = 0.f;
      OSPData isovaluesData = ospNewData(1, OSP_FLOAT, &isovalue);
      ospSetObject(isosurface, "isovalue", isovaluesData);
      ospRelease(isovaluesData);

      ospCommit(isosurface);

      OSPGeometricModel isoInstance = ospNewGeometricModel(isosurface);
      ospRelease(isosurface);
      AddModel(isoInstance);
    } else {
      AddModel(model);
    }

    ospRelease(pyramid);
    ospRelease(transferFun);

    OSPLight ambient = ospNewLight("ambient");
    ASSERT_TRUE(ambient) << "Failed to create lights";
    ospSetFloat(ambient, "intensity", 0.5f);
    ospCommit(ambient);
    AddLight(ambient);
  }

  Torus::Torus()
  {
    rendererType = GetParam();
  }

  void Torus::SetUp()
  {
    Base::SetUp();

    float cam_pos[]  = {-0.7f, -1.4f, 0.f};
    float cam_up[]   = {0.f, 0.f, -1.f};
    float cam_view[] = {0.5f, 1.f, 0.f};
    ospSetParam(camera, "position", OSP_VEC3F, cam_pos);
    ospSetParam(camera, "direction", OSP_VEC3F, cam_view);
    ospSetParam(camera, "up", OSP_VEC3F, cam_up);

    OSPVolume torus = CreateTorus(volumetricData, 256);
    ospCommit(torus);

    OSPVolumetricModel volumetricModel = ospNewVolumetricModel(torus);
    ospRelease(torus);

    OSPTransferFunction transferFun =
        ospNewTransferFunction("piecewise_linear");
    ospSetVec2f(transferFun, "valueRange", -10000.f, 10000.f);
    float colors[]      = {1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f};
    float opacites[]    = {1.0f, 1.0f};
    OSPData tfColorData = ospNewData(2, OSP_VEC3F, colors);
    ospSetObject(transferFun, "color", tfColorData);
    ospRelease(tfColorData);
    OSPData tfOpacityData = ospNewData(2, OSP_FLOAT, opacites);
    ospSetObject(transferFun, "opacity", tfOpacityData);
    ospRelease(tfOpacityData);
    ospCommit(transferFun);
    ospSetObject(volumetricModel, "transferFunction", transferFun);
    ospRelease(transferFun);

    ospCommit(volumetricModel);

    OSPGeometry isosurface = ospNewGeometry("isosurfaces");
    ospSetObject(isosurface, "volume", volumetricModel);
    float isovalues[2]    = {-7000.f, 0.f};
    OSPData isovaluesData = ospNewData(2, OSP_FLOAT, isovalues);
    ospSetObject(isosurface, "isovalue", isovaluesData);
    ospRelease(isovaluesData);
    ospCommit(isosurface);

    OSPGeometricModel model = ospNewGeometricModel(isosurface);
    AddModel(model);

    OSPLight ambient = ospNewLight("ambient");
    ASSERT_TRUE(ambient) << "Failed to create lights";
    ospSetFloat(ambient, "intensity", 0.5f);
    ospCommit(ambient);
    AddLight(ambient);
  }

  SlicedCube::SlicedCube() {}

  void SlicedCube::SetUp()
  {
    Base::SetUp();

    float cam_pos[]  = {-0.7f, -1.4f, 0.f};
    float cam_up[]   = {0.f, 0.f, -1.f};
    float cam_view[] = {0.5f, 1.f, 0.f};
    ospSetParam(camera, "position", OSP_VEC3F, cam_pos);
    ospSetParam(camera, "direction", OSP_VEC3F, cam_view);
    ospSetParam(camera, "up", OSP_VEC3F, cam_up);

    int size = 100;

    volumetricData.resize(size * size * size, 0);

    for (int x = 0; x < size; ++x) {
      for (int y = 0; y < size; ++y) {
        for (int z = 0; z < size; ++z) {
          float X = (x / (float)size) - 0.5;
          float Y = (y / (float)size) - 0.5;
          float Z = (z / (float)size) - 0.5;

          volumetricData[size * size * x + size * y + z] =
              sin(30 * X + 3 * sin(30 * Y + 3 * sin(30 * Z)));
        }
      }
    }

    OSPVolume blob = ospNewVolume("shared_structured_volume");
    ASSERT_TRUE(blob);
    OSPData voxelsData = ospNewData(size * size * size,
                                    OSP_FLOAT,
                                    volumetricData.data(),
                                    OSP_DATA_SHARED_BUFFER);
    ASSERT_TRUE(voxelsData);
    ospSetObject(blob, "voxelData", voxelsData);
    ospSetVec3i(blob, "dimensions", size, size, size);
    ospSetInt(blob, "voxelType", OSP_FLOAT);
    ospSetVec2f(blob, "voxelRange", 0.f, 3.f);
    ospSetVec3f(blob, "gridOrigin", -0.5f, -0.5f, -0.5f);
    ospSetVec3f(blob, "gridSpacing", 1.f / size, 1.f / size, 1.f / size);
    ospCommit(blob);

    OSPVolumetricModel volumetricModel = ospNewVolumetricModel(blob);

    OSPTransferFunction transferFun =
        ospNewTransferFunction("piecewise_linear");
    ASSERT_TRUE(transferFun);
    ospSetVec2f(transferFun, "valueRange", -1.f, 1.f);
    float colors[]      = {0.85f, 0.85f, 1.0f, 0.1f, 0.0f, 0.0f};
    float opacites[]    = {0.0f, 1.0f};
    OSPData tfColorData = ospNewData(2, OSP_VEC3F, colors);
    ASSERT_TRUE(tfColorData);
    ospSetObject(transferFun, "color", tfColorData);
    OSPData tfOpacityData = ospNewData(2, OSP_FLOAT, opacites);
    ASSERT_TRUE(tfOpacityData);
    ospSetObject(transferFun, "opacity", tfOpacityData);
    ospCommit(transferFun);

    ospSetObject(volumetricModel, "transferFunction", transferFun);
    ospRelease(transferFun);
    ospCommit(volumetricModel);

    OSPGeometry slice = ospNewGeometry("slices");
    ASSERT_TRUE(slice);
    float planes[] = {
        1.f, 1.f, 1.f, 0.5f, 1.f, 1.f, 1.f, 0.f, 1.f, 1.f, 1.f, -0.5f};
    OSPData planesData = ospNewData(3, OSP_VEC4F, planes);
    ASSERT_TRUE(planesData);
    ospSetObject(slice, "plane", planesData);
    ospSetObject(slice, "volume", volumetricModel);
    ospCommit(slice);
    ospRelease(volumetricModel);

    OSPGeometricModel model = ospNewGeometricModel(slice);
    AddModel(model);

    OSPLight ambient = ospNewLight("ambient");
    ASSERT_TRUE(ambient);
    ospSetFloat(ambient, "intensity", 0.5f);
    ospCommit(ambient);
    AddLight(ambient);
  }

  MTLMirrors::MTLMirrors()
  {
    auto param = GetParam();
    Kd         = std::get<0>(param);
    Ks         = std::get<1>(param);
    Ns         = std::get<2>(param);
    d          = std::get<3>(param);
    Tf         = std::get<4>(param);

    rendererType = "pathtracer";
  }

  void MTLMirrors::SetUp()
  {
    Base::SetUp();

    ospSetVec3f(renderer, "bgColor", 0.5f, 0.5f, 0.45f);

    OSPData data;

    float mirrorsVertices[] = {
        2.f, 1.f, 5.f, 2.f, -1.f, 5.f, 0.f,  1.f, 2.f, 0.f,  -1.f, 2.f,
        0.f, 1.f, 5.f, 0.f, -1.f, 5.f, -2.f, 1.f, 2.f, -2.f, -1.f, 2.f,
    };
    int32_t mirrorsIndices[] = {0, 1, 2, 1, 2, 3, 4, 5, 6, 5, 6, 7};
    OSPGeometry mirrors      = ospNewGeometry("triangles");
    ASSERT_TRUE(mirrors);
    data = ospNewData(8, OSP_VEC3F, mirrorsVertices);
    ASSERT_TRUE(data);
    ospCommit(data);
    ospSetObject(mirrors, "vertex.position", data);
    data = ospNewData(4, OSP_VEC3UI, mirrorsIndices);
    ASSERT_TRUE(data);
    ospSetObject(mirrors, "index", data);
    OSPMaterial mirrorsMaterial =
        ospNewMaterial(rendererType.c_str(), "OBJMaterial");
    ASSERT_TRUE(mirrorsMaterial);
    ospSetVec3f(mirrorsMaterial, "Kd", Kd.x, Kd.y, Kd.z);
    ospSetVec3f(mirrorsMaterial, "Ks", Ks.x, Ks.y, Ks.z);
    ospSetFloat(mirrorsMaterial, "Ns", Ns);
    ospSetFloat(mirrorsMaterial, "d", d);
    ospSetVec3f(mirrorsMaterial, "Tf", Tf.x, Tf.y, Tf.z);
    ospCommit(mirrorsMaterial);
    ospCommit(mirrors);

    OSPGeometricModel model = ospNewGeometricModel(mirrors);
    ospSetObject(model, "material", mirrorsMaterial);
    AddModel(model);

    float sphereCenters[] = {1.f, 0.f, 7.f};
    OSPGeometry light     = ospNewGeometry("spheres");
    ASSERT_TRUE(light);
    ospSetFloat(light, "radius", 1.f);
    data = ospNewData(1, OSP_VEC3F, sphereCenters);
    ASSERT_TRUE(data);
    ospSetObject(light, "sphere.position", data);
    ospCommit(light);

    model = ospNewGeometricModel(light);
    ospSetObject(
        model, "material", ospNewMaterial(rendererType.c_str(), "Luminous"));
    AddModel(model);

    OSPLight ambient = ospNewLight("ambient");
    ASSERT_TRUE(ambient);
    ospSetFloat(ambient, "intensity", 0.01f);
    ospCommit(ambient);
    AddLight(ambient);
  }

  Pipes::Pipes()
  {
    auto params  = GetParam();
    rendererType = std::get<0>(params);
    materialType = std::get<1>(params);
    radius       = std::get<2>(params);

    frames = 10;
  }

  void Pipes::SetUp()
  {
    Base::SetUp();

    float cam_pos[]  = {-7.f, 2.f, 0.7f};
    float cam_view[] = {7.f, -2.f, -0.7f};
    float cam_up[]   = {0.f, 0.f, 1.f};
    ospSetParam(camera, "position", OSP_VEC3F, cam_pos);
    ospSetParam(camera, "direction", OSP_VEC3F, cam_view);
    ospSetParam(camera, "up", OSP_VEC3F, cam_up);

    float vertex[] = {
      -2.f, 2.f,  -2.f,  2.f,  2.f,  -2.f,  2.f,  -2.f, -2.f,
      -2.f, -2.f, -2.f,  -1.f, -1.f, -1.f,  1.f,  -1.f,
      -1.f,  1.f,  1.f,  -1.f,  -1.f, 1.f,  -1.f,  -1.f,
      1.f,  1.f,  1.f,  1.f,  1.f,  1.f,  -1.f, 1.f,
      -1.f, -1.f, 1.f,  -2.f, -2.f, 2.f,  2.f,  -2.f, 2.f,
      2.f,  2.f,  2.f,  -2.f, 2.f,  2.f
    };

    int index[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14};
    OSPGeometry streamlines = ospNewGeometry("streamlines");
    ASSERT_TRUE(streamlines);
    OSPData data = ospNewData(16, OSP_VEC3F, vertex);
    ASSERT_TRUE(data);
    ospSetObject(streamlines, "vertex.position", data);
    ospRelease(data);
    data = ospNewData(15, OSP_UINT, index);
    ASSERT_TRUE(data);
    ospSetObject(streamlines, "index", data);
    ospRelease(data);
    ospSetFloat(streamlines, "radius", radius);
    ospCommit(streamlines);

    OSPGeometricModel model = ospNewGeometricModel(streamlines);
    ospRelease(streamlines);
    OSPMaterial material = CreateMaterial(materialType);
    ospSetObject(model, "material", material);
    ospRelease(material);
    AddModel(model);

    OSPLight ambient = ospNewLight("ambient");
    ASSERT_TRUE(ambient);
    ospSetFloat(ambient, "intensity", 0.5f);
    ospCommit(ambient);
    AddLight(ambient);
  }

  TextureVolume::TextureVolume()
  {
    rendererType = GetParam();
  }

  void TextureVolume::SetUp()
  {
    Base::SetUp();

    float cam_pos[]  = {-0.7f, -1.4f, 0.f};
    float cam_up[]   = {0.f, 0.f, -1.f};
    float cam_view[] = {0.5f, 1.f, 0.f};
    ospSetParam(camera, "position", OSP_VEC3F, cam_pos);
    ospSetParam(camera, "direction", OSP_VEC3F, cam_view);
    ospSetParam(camera, "up", OSP_VEC3F, cam_up);

    OSPVolume torus = CreateTorus(volumetricData, 256);
    ospCommit(torus);

    OSPVolumetricModel volumetricModel = ospNewVolumetricModel(torus);
    ospRelease(torus);

    OSPTransferFunction transferFun =
        ospNewTransferFunction("piecewise_linear");
    ospSetVec2f(transferFun, "valueRange", -10000.f, 10000.f);
    float colors[]      = {1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f};
    float opacites[]    = {1.0f, 1.0f};
    OSPData tfColorData = ospNewData(2, OSP_VEC3F, colors);
    ospSetObject(transferFun, "color", tfColorData);
    ospRelease(tfColorData);
    OSPData tfOpacityData = ospNewData(2, OSP_FLOAT, opacites);
    ospSetObject(transferFun, "opacity", tfOpacityData);
    ospRelease(tfOpacityData);
    ospCommit(transferFun);
    ospSetObject(volumetricModel, "transferFunction", transferFun);
    ospRelease(transferFun);
    ospCommit(volumetricModel);

    OSPMaterial sphereMaterial =
        ospNewMaterial(rendererType.c_str(), "default");
    OSPTexture tex = ospNewTexture("volume");
    ospSetObject(tex, "volume", volumetricModel);
    ospRelease(volumetricModel);
    ospCommit(tex);
    ospSetObject(sphereMaterial, "map_Kd", tex);
    ospRelease(tex);
    ospCommit(sphereMaterial);

    float sphereVertex[] = {0.f, 0.f, 0.f};
    OSPGeometry sphere   = ospNewGeometry("spheres");
    auto data = ospNewData(1, OSP_VEC3F, sphereVertex);
    ospCommit(data);
    ospSetObject(sphere, "sphere.position", data);
    ospRelease(data);
    ospSetFloat(sphere, "radius", 0.51f);
    ospCommit(sphere);

    OSPGeometricModel model = ospNewGeometricModel(sphere);
    ospRelease(sphere);
    ospSetObject(model, "material", sphereMaterial);
    ospRelease(sphereMaterial);
    AddModel(model);

    OSPLight ambient = ospNewLight("ambient");
    ASSERT_TRUE(ambient) << "Failed to create lights";
    ospSetFloat(ambient, "intensity", 0.5f);
    ospCommit(ambient);
    AddLight(ambient);
  }

  DepthCompositeVolume::DepthCompositeVolume()
  {
    auto params  = GetParam();
    rendererType = std::get<0>(params);
    bgColor      = std::get<1>(params);
  }

  void DepthCompositeVolume::SetUp()
  {
    Base::SetUp();

    const float cam_pos[]  = {0.f, 0.f, 1.0f};
    const float cam_up[]   = {0.f, 1.f, 0.f};
    const float cam_view[] = {0.0f, 0.f, -1.f};
    ospSetParam(camera, "position", OSP_VEC3F, cam_pos);
    ospSetParam(camera, "direction", OSP_VEC3F, cam_view);
    ospSetParam(camera, "up", OSP_VEC3F, cam_up);

    OSPVolume torus = CreateTorus(volumetricData, 256);
    ospCommit(torus);

    OSPVolumetricModel model = ospNewVolumetricModel(torus);
    ospRelease(torus);

    OSPTransferFunction transferFun =
        ospNewTransferFunction("piecewise_linear");
    ospSetVec2f(transferFun, "valueRange", -10000.f, 10000.f);
    const float colors[]   = {1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f};
    const float opacites[] = {0.05f, 0.1f};
    OSPData tfColorData    = ospNewData(2, OSP_VEC3F, colors);
    ospSetObject(transferFun, "color", tfColorData);
    ospRelease(tfColorData);
    OSPData tfOpacityData = ospNewData(2, OSP_FLOAT, opacites);
    ospSetObject(transferFun, "opacity", tfOpacityData);
    ospRelease(tfOpacityData);
    ospCommit(transferFun);
    ospSetObject(model, "transferFunction", transferFun);
    ospRelease(transferFun);

    AddModel(model);

    OSPTexture depthTex = ospNewTexture("texture2d");
    static std::vector<float> data;
    data.resize(imgSize.x * imgSize.y);
    for (int y = 0; y < imgSize.y; y++) {
      for (int x = 0; x < imgSize.x; x++) {
        const size_t index = imgSize.x * y + x;
        if (x < imgSize.x / 3) {
          data[index] = 999.f;
        } else if (x < (imgSize.x * 2) / 3) {
          data[index] = 0.75f;
        } else {
          data[index] = 0.00001f;
        }
      }
    }
    auto ospData =
        ospNewSharedData(data.data(), OSP_FLOAT, imgSize.x, 0, imgSize.y);
    ospCommit(ospData);
    ospSetInt(depthTex, "format", OSP_TEXTURE_R32F);
    ospSetInt(depthTex, "filter", OSP_TEXTURE_FILTER_NEAREST);
    ospSetObject(depthTex, "data", ospData);
    ospCommit(depthTex);
    ospRelease(ospData);

    ospSetObject(renderer, "maxDepthTexture", depthTex);
    ospSetParam(renderer, "bgColor", OSP_VEC4F, &bgColor.x);
    ospRelease(depthTex);

    OSPLight ambient = ospNewLight("ambient");
    ASSERT_TRUE(ambient) << "Failed to create lights";
    ospSetFloat(ambient, "intensity", 0.5f);
    ospCommit(ambient);
    AddLight(ambient);
  }

  Subdivision::Subdivision()
  {
    auto params  = GetParam();
    rendererType = std::get<0>(params);
    materialType = std::get<1>(params);
    frames       = 10;
  }

  void Subdivision::SetUp()
  {
    Base::SetUp();

    float cube_vertices[8][3] = {{-1.0f, -1.0f, -1.0f},
                                 {1.0f, -1.0f, -1.0f},
                                 {1.0f, -1.0f, 1.0f},
                                 {-1.0f, -1.0f, 1.0f},
                                 {-1.0f, 1.0f, -1.0f},
                                 {1.0f, 1.0f, -1.0f},
                                 {1.0f, 1.0f, 1.0f},
                                 {-1.0f, 1.0f, 1.0f}};

    float cube_colors[8][4] = {{0.0f, 0.0f, 0.0f, 1.f},
                               {1.0f, 0.0f, 0.0f, 1.f},
                               {1.0f, 0.0f, 1.0f, 1.f},
                               {0.0f, 0.0f, 1.0f, 1.f},
                               {0.0f, 1.0f, 0.0f, 1.f},
                               {1.0f, 1.0f, 0.0f, 1.f},
                               {1.0f, 1.0f, 1.0f, 1.f},
                               {0.0f, 1.0f, 1.0f, 1.f}};

    unsigned int cube_indices[24] = {
        0, 4, 5, 1, 1, 5, 6, 2, 2, 6, 7, 3, 0, 3, 7, 4, 4, 7, 6, 5, 0, 1, 2, 3,
    };

    unsigned int cube_faces[6] = {4, 4, 4, 4, 4, 4};

    float cube_vertex_crease_weights[8] = {2, 2, 2, 2, 2, 2, 2, 2};

    unsigned int cube_vertex_crease_indices[8] = {0, 1, 2, 3, 4, 5, 6, 7};

    float cube_edge_crease_weights[12] = {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2};

    unsigned int cube_edge_crease_indices[24] = {
        0, 1, 1, 2, 2, 3, 3, 0, 4, 5, 5, 6, 6, 7, 7, 4, 0, 4, 1, 5, 2, 6, 3, 7,
    };

    int numIndices = 24;
    int numFaces   = 6;

    auto subd     = ospNewGeometry("subdivision");
    auto vertices = ospNewData(8, OSP_VEC3F, cube_vertices);
    ospSetObject(subd, "vertex.position", vertices);
    ospRelease(vertices);
    auto indices = ospNewData(numIndices, OSP_UINT, cube_indices);
    ospSetObject(subd, "index", indices);
    auto faces = ospNewData(numFaces, OSP_UINT, cube_faces);
    ospSetObject(subd, "face", faces);
    ospRelease(faces);
    auto edge_crease_indices =
        ospNewData(12, OSP_VEC2UI, cube_edge_crease_indices);
    ospSetObject(subd, "edgeCrease.index", edge_crease_indices);
    ospRelease(edge_crease_indices);
    auto edge_crease_weights =
        ospNewData(12, OSP_FLOAT, cube_edge_crease_weights);
    ospSetObject(subd, "edgeCrease.weight", edge_crease_weights);
    ospRelease(edge_crease_weights);
    auto vertex_crease_indices =
        ospNewData(8, OSP_UINT, cube_vertex_crease_indices);
    ospSetObject(subd, "vertexCrease.index", vertex_crease_indices);
    ospRelease(vertex_crease_indices);
    auto vertex_crease_weights =
        ospNewData(8, OSP_FLOAT, cube_vertex_crease_weights);
    ospSetObject(subd, "vertexCrease.weight", vertex_crease_weights);
    ospRelease(vertex_crease_weights);
    auto colors = ospNewData(8, OSP_VEC4F, cube_colors);
    ospSetObject(subd, "vertex.color", colors);
    ospRelease(colors);
    ospSetFloat(subd, "level", 128.0f);

    ospCommit(subd);

    OSPGeometricModel model = ospNewGeometricModel(subd);
    ospRelease(subd);
    OSPMaterial material = CreateMaterial(materialType);
    ospSetObject(model, "material", material);
    ospRelease(material);
    AddModel(model);

    float cam_pos[]  = {-1.5f, 2.f, 1.7f};
    float cam_view[] = {1.5f, -2.f, -1.7f};
    float cam_up[]   = {0.f, 1.f, 0.f};
    ospSetParam(camera, "position", OSP_VEC3F, cam_pos);
    ospSetParam(camera, "direction", OSP_VEC3F, cam_view);
    ospSetParam(camera, "up", OSP_VEC3F, cam_up);

    OSPLight directional = ospNewLight("directional");
    ASSERT_TRUE(directional);
    ospSetFloat(directional, "intensity", 0.5f);
    ospSetVec3f(directional, "direction", -.2f, -.3f, -.4f);
    ospCommit(directional);
    AddLight(directional);
  }
  
  /* heavily based on Perlin's Java reference implementation of
  * the improved perlin noise paper from Siggraph 2002 from here
  * https://mrl.nyu.edu/~perlin/noise/
  **/
  class PerlinNoise 
  {
    struct PerlinNoiseData
    {
      PerlinNoiseData()
      {
        const int permutation[256] = { 
          151,160,137,91,90,15,131,13,201,95,96,53,194,233,7,225,140,36,103,30,69,142,8,99,
          37,240,21,10,23,190,6,148,247,120,234,75,0,26,197,62,94,252,219,203,117,35,11,32,
          57,177,33,88,237,149,56,87,174,20,125,136,171,168,68,175,74,165,71,134,139,48,27,
          166,77,146,158,231,83,111,229,122,60,211,133,230,220,105,92,41,55,46,245, 40,244,
          102,143,54,65,25,63,161,1,216,80,73,209,76,132,187,208,89,18,169,200,196,135,130,
          116,188,159,86,164,100,109,198,173,186,3,64,52,217,226,250,124,123,5,202, 38,147,
          118,126,255,82,85,212,207,206,59,227,47,16,58,17,182,189, 28, 42,223,183,170,213,
          119,248,152,2,44,154,163,70,221,153,101,155,167,43,172,9,129,22,39,253,19,98,108,
          110,79,113,224,232,178,185,112,104,218,246,97,228,251, 34,242,193,238,210,144,12,
          191,179,162,241,81,51,145,235,249,14,239,107, 49,192,214, 31,181,199,106,157,184,
          84,204,176,115,121,50,45,127,4,150,254,138,236,205,93,222,114,67, 29, 24, 72,243,
          141,128,195,78,66,215,61,156,180
        };
        for (int i=0; i < 256 ; i++) p[256+i] = p[i] = permutation[i];
      }
      inline int operator[](size_t idx) const { return p[idx]; }
      int p[512];
    };

    static PerlinNoiseData p;
    static inline float smooth(float t) { return t * t * t * (t * (t * 6.f - 15.f) + 10.f); }
    static inline float lerp  (float t, float a, float b) { return a + t * (b - a); }
    static inline float grad(int hash, float x, float y, float z) {
      const int h   = hash & 15;      
      const float u = h < 8 ? x : y;
      const float v = h < 4 ? y : h == 12 || h == 14 ? x : z;
      return ((h & 1) == 0 ? u : -u) + ((h & 2) == 0 ? v : -v);
    }

  public:
    static float noise(vec3f q, float frequency = 8.f)
    {
      float x = q.x * frequency;
      float y = q.y * frequency;
      float z = q.z * frequency;
      const int X = (int)floor(x) & 255;
      const int Y = (int)floor(y) & 255; 
      const int Z = (int)floor(z) & 255;
      x -= floor(x);
      y -= floor(y);
      z -= floor(z);
      const float u = smooth(x);
      const float v = smooth(y);
      const float w = smooth(z);
      const int  A = p[X] + Y;
      const int  B = p[X + 1] + Y;
      const int AA = p[A] + Z;
      const int BA = p[B] + Z;
      const int BB = p[B + 1] + Z; 
      const int AB = p[A + 1] + Z;

      return lerp(w, lerp(v, lerp(u, grad(p[AA], x, y, z), 
                                    grad(p[BA], x - 1, y, z)),
                            lerp(u, grad(p[AB], x, y - 1, z),  
                                    grad(p[BB], x - 1, y - 1, z))),
                    lerp(v, lerp(u, grad(p[AA + 1], x, y, z - 1),  
                                    grad(p[BA + 1], x - 1, y, z - 1)),
                            lerp(u, grad(p[AB + 1], x, y - 1, z - 1),
                                    grad(p[BB + 1], x - 1, y - 1, z - 1))));
    }
  };
  PerlinNoise::PerlinNoiseData PerlinNoise::p;

  HeterogeneousVolume::HeterogeneousVolume()
  {
    rendererType = "pathtracer";

    auto params = GetParam();

    albedo             = std::get<0>(params);
    anisotropy         = std::get<1>(params);
    densityScale       = std::get<2>(params);
    ambientColor       = std::get<3>(params);
    enableDistantLight = std::get<4>(params);
    enableGeometry     = std::get<5>(params);
    constantVolume     = std::get<6>(params);
    samplesPerPixel    = std::get<7>(params);

    imgSize = imgSize / 2.f;
    Base::SetImageTool();
  }

  void HeterogeneousVolume::SetUp()
  {
    Base::SetUp();

    if (!constantVolume)
      densityScale *= 10.f;
    
    const float theta = 1.f/2.f*M_PI - M_PI/6;
    const float phi   = 3.f/2.f*M_PI - M_PI/6;
    const float r = 5.f;
    vec3f p = r * vec3f(
      sin(theta)*cos(phi), 
      cos(theta), 
      sin(theta)*sin(phi));
    vec3f v = vec3f(0.f, -0.2f, 0.f) - p;

    ospSetVec3fv(camera, "position", &p.x);
    ospSetVec3fv(camera, "direction", &v.x);
    ospSetVec3f(camera, "up", 0.f, 1.f, 0.f);
    
    // create volume
    vec3l dims;
    if (constantVolume)
      dims = vec3l(8, 8, 8);
    else
      dims = vec3l(64, 64, 64);
    
    const float spacing = 3.f/(reduce_max(dims)-1);
    OSPVolume volume = ospNewVolume("vkl_structured_volume");

    auto turbulence = [](const vec3f& p, float base_freqency, int octaves)
    {
      float value = 0.f;
      float scale = 1.f;
      for (int o = 0; o < octaves; ++o)
      {
        value += PerlinNoise::noise(scale*p, base_freqency)/scale;
        scale *= 2.f;
      }
      return value;
    };

    // generate volume data
    auto numVoxels = dims.product();
    std::vector<float> voxels(numVoxels, 0);
    for (int z = 0; z < dims.z; ++z)
    for (int y = 0; y < dims.y; ++y)
    for (int x = 0; x < dims.x; ++x)
    {
      if (constantVolume)
        voxels[dims.x * dims.y * z + dims.x * y + x] = 1.0f;
      else
      {
        vec3f p = vec3f(x+0.5f, y+0.5f, z+0.5f)/dims;
        vec3f X = 2.f * p - vec3f(1.f); 
        if (length((1.4f + 0.4 * turbulence(p, 12.f, 12)) * X) < 1.f) 
          voxels[dims.x * dims.y * z + dims.x * y + x] = 0.5f + 0.5f * PerlinNoise::noise(p, 12);
      }
    }
    voxels[0] = 0.0f;

    OSPData voxelData = ospNewData(numVoxels, OSP_FLOAT, voxels.data());
    ospSetObject(volume, "voxelData", voxelData);
    ospRelease(voxelData);

    ospSetInt(volume, "voxelType", OSP_FLOAT);
    ospSetVec3i(volume, "dimensions", dims.x, dims.y, dims.z);
    ospSetVec3f(volume, "gridOrigin", -1.5f, -1.5f, -1.5f);
    ospSetVec3f(volume, "gridSpacing", spacing, spacing, spacing);
    ospCommit(volume);

    // create transfer function
    OSPTransferFunction tfn = ospNewTransferFunction("piecewise_linear");
    ospSetVec2f(tfn, "valueRange", 0.0f, 1.0f);
    OSPData tfColorData = ospNewData(1, OSP_VEC3F, &albedo);
    ospSetData(tfn, "color", tfColorData);
    ospRelease(tfColorData);
    std::vector<float> opacities = { 0.f, 1.f };
    OSPData tfOpacityData = ospNewData(opacities.size(), OSP_FLOAT, opacities.data());
    ospSetData(tfn, "opacity", tfOpacityData);
    ospRelease(tfOpacityData);
    ospCommit(tfn);
    OSPVolumetricModel volumetricModel = ospNewVolumetricModel(volume);
    ospSetFloat(volumetricModel,  "densityScale",     densityScale);
    ospSetFloat(volumetricModel,  "anisotropy",       anisotropy);
    ospSetObject(volumetricModel, "transferFunction", tfn);
    ospCommit(volumetricModel);
    ospRelease(tfn);
    ospRelease(volume);
    ospCommit(volumetricModel);

    OSPData volumeData = ospNewData(1, OSP_VOLUMETRIC_MODEL, &volumetricModel);
    OSPGroup volumeGroup = ospNewGroup();
    ospSetData(volumeGroup, "volume", volumeData);
    ospCommit(volumeGroup);
    ospRelease(volumeData);
 
    OSPInstance volumeInstance = ospNewInstance(volumeGroup);
    if (!constantVolume) {
      AffineSpace3f xfm = AffineSpace3f::scale(vec3f(1.25f));
      ospSetAffine3fv(volumeInstance, "xfm", &xfm.l.vx.x);
    }
    ospCommit(volumeInstance);
    ospRelease(volumeGroup);
    AddInstance(volumeInstance);

    if (enableGeometry)
    {
      std::vector<OSPGeometricModel> models;
      std::vector<vec3f> planeVertices = { 
        vec3f(-8.f, -2.5f, -8.f),
        vec3f(+8.f, -2.5f, -8.f),
        vec3f(+8.f, -2.5f, +8.f),
        vec3f(-8.f, -2.5f, +8.f) };
      std::vector<vec4i> planeIndices = { vec4i(0, 1, 2, 3) };

      OSPData positionData = ospNewData(planeVertices.size(), OSP_VEC3F,  planeVertices.data());
      OSPData indexData    = ospNewData(planeIndices.size(),  OSP_VEC4UI, planeIndices.data());
      OSPGeometry mesh  = ospNewGeometry("quads");
      ospSetData(mesh, "vertex.position", positionData);
      ospSetData(mesh, "index",           indexData);
      ospCommit(mesh);

      OSPGeometricModel model = ospNewGeometricModel(mesh);
      ospRelease(mesh);
      OSPMaterial material = ospNewMaterial("pathtracer", "OBJMaterial");
      ospCommit(material);
      ospSetObject(model, "material", material);
      ospRelease(material);
      if (enableGeometry)
        models.push_back(model);

      for (auto &m : models)
        ospCommit(m);

      OSPData modelsData = ospNewData(models.size(), OSP_GEOMETRIC_MODEL, models.data());
      OSPGroup modelsGroup = ospNewGroup();
      ospSetData(modelsGroup, "geometry", modelsData);
      ospCommit(modelsGroup);
      ospRelease(modelsData);
      OSPInstance modelsInstance = ospNewInstance(modelsGroup);
      ospCommit(modelsInstance);
      ospRelease(modelsGroup);
      AddInstance(modelsInstance);
    }


    //vec4f backgroundColor(0.03f, 0.07f, 0.23f);
    vec4f backgroundColor(ambientColor, 1.f);

    std::vector<OSPLight> lightHandles;
    OSPLight ambientLight = ospNewLight("ambient");
    ospSetFloat(ambientLight, "intensity", 1.f);
    ospSetVec3fv(ambientLight, "color", &ambientColor.x);
    ospCommit(ambientLight);
    lightHandles.push_back(ambientLight);
    
    OSPLight distantLight = ospNewLight("distant");
    ospSetFloat(distantLight, "intensity", 2.6f);
    ospSetVec3f(distantLight, "color", 1.0f, 0.96f, 0.88f);
    ospSetFloat(distantLight, "angularDiameter", 1.f);
    ospSetVec3f(distantLight, "direction", -0.5826f, -0.7660f, -0.2717f);
    ospCommit(distantLight);
    if (enableDistantLight)
      lightHandles.push_back(distantLight);

    OSPData lights = ospNewData(lightHandles.size(), OSP_LIGHT, lightHandles.data(), 0);
    ospCommit(lights);
    ospSetData(renderer, "light", lights);
    ospRelease(lights);
    
    OSPData texelData = ospNewData(1, OSP_VEC3F, &ambientColor.x);

    OSPTexture backplateTexture = ospNewTexture("texture2d");
    ospSetVec2i(backplateTexture, "size", 1, 1);
    ospSetInt(backplateTexture, "format", OSP_TEXTURE_RGB32F);
    ospSetInt(backplateTexture, "filter", OSP_TEXTURE_FILTER_NEAREST);
    ospSetData(backplateTexture, "data", texelData);
    ospCommit(backplateTexture);

    ospRelease(texelData);

    ospSetObject(renderer, "backplate", backplateTexture);
    ospRelease(backplateTexture);

    // NOTE(jda) - still need to set the world on the renderer for geom lights
    ospSetObject(renderer, "world", world);
    ospSetInt(renderer, "maxDepth", std::max(20, samplesPerPixel));
    ospSetInt(renderer, "spp", samplesPerPixel);
  }
}  // namespace OSPRayTestScenes
