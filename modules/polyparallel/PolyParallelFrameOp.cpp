/* Polynomial Optics
 * (C) 2012, Matthias Hullin <hullin@cs.ubc.ca>
 * (University of British Columbia)
 * (C) 2012, Johannes Hanika <hanatos@gmail.com>
 * (Weta Digital)
 * http://www.cs.ubc.ca/labs/imager/tr/2012/PolynomialOptics/
 */

#include "TruncPoly/TruncPolySystem.hh"

#include "OpticalElements/Cylindrical5.hh"
#include "OpticalElements/FindFocus.hh"
#include "OpticalElements/OpticalMaterial.hh"
#include "OpticalElements/Propagation5.hh"
#include "OpticalElements/Spherical5.hh"
#include "OpticalElements/TwoPlane5.hh"

#include "include/spectrum.h"

#include "PolyParallelFrameOp.h"

namespace ospray {

struct OSPRAY_MODULE_POLYPARALLEL_EXPORT LivePolyParallelFrameOp
    : public LiveFrameOp
{
  LivePolyParallelFrameOp(FrameBufferView &_fbView) : LiveFrameOp(_fbView) {}
  void process(const Camera *) override
  {
    float *color = static_cast<float *>(fbView.colorBuffer);
    const int arr_size = fbView.fbDims.x * fbView.fbDims.y * 4;
    float *arr = new float[arr_size];
    execute_filter(32.0, fbView.fbDims.x, fbView.fbDims.y, color, arr);
  }

  void execute_filter(float f_stop, int fb_width, int fb_height, float *color, float *new_color)
  {
    // Initial distance set to 50 meters
    float d = 5000;
    int degree = 3;
    float sample_mul = 1000;
    // f-stop set to 2.0
    float r_entrance = 39 / f_stop;
    int num_lambdas = 12;

    const float sensor_width = 36;

    const float lambda_from = 440;
    const float lambda_to = 660;

    int sensor_xres = fb_width;
    int sensor_yres = fb_height;
    const float sensor_scaling = sensor_xres / sensor_width;

    // Frame set to 299 for animated defocus
    int frame = 299;
    float r_pupil = r_entrance;

    // Focus on lambda = 550nm (Green)
    Transform4f system = get_system(550, degree, d);
    float d3 = find_focus_X(system);
    float magnification = get_magnification_X(system >> propagate_5(d3));
    // Incorporates animated defocus
    Transform4f prop =
        propagate_5(d3 - ((frame >= 100) ? (0.02 * (frame - 100)) : 0), degree);
    system = system >> prop;

    // Precomputes spectrum for each lambda
    float rgb[3 * num_lambdas];
    for (int ll = 0; ll < num_lambdas; ++ll) {
      float lambda = lambda_from
          + (lambda_to - lambda_from) * (ll / (float)(num_lambdas - 1));
      if (num_lambdas == 1)
        lambda = 550;
      spectrum_p_to_rgb(lambda, 1, rgb + 3 * ll);
    }

    // Optical system at 2 spectral locations (Cyan and Orange)
    Transform4d system_spectral_center = get_system(500) >> prop;
    Transform4d system_spectral_right = get_system(600, degree) >> prop;
    // Obtains system that maps (xyworld + xyaperture) -> (ray)
    System54f system_spectral =
        system_spectral_center.lerp_with(system_spectral_right, 550, 600);

    // Shortcut to calculating Lambertian (Trig Rules)
    system_spectral[2] = (system_spectral[2] * system_spectral[2]
        + system_spectral[3] * system_spectral[3]);
    system_spectral[2] %= 2;
    System53d system_lambert_cos2 = system_spectral.drop_equation(3);

    float pixel_size = sensor_width / (float)fb_width / magnification;

    for (int ll = 0; ll < num_lambdas; ++ll) {
      float lambda = lambda_from
          + (lambda_to - lambda_from) * (ll / (float)(num_lambdas - 1));
      if (num_lambdas == 1)
        lambda = 550;
      cout << "[" << lambda << "nm]" << flush;

      // Bake lambda dependency
      System43f system_lambda =
          system_lambert_cos2.bake_input_variable(4, lambda);
      system_lambda %= degree;

      
      tasking::parallel_for(fb_height * fb_width, [&](int pixelID) {
        int i = pixelID % fb_width;
        int j = pixelID / fb_width;

        const float y_sensor =
            ((j - fb_height / 2) / (float)fb_width) * sensor_width;
        const float y_world = y_sensor / magnification;

        System33f system_y = system_lambda.bake_input_variable(1, y_world);
        const float x_sensor = (i / (float)fb_width - 0.5) * sensor_width;
        const float x_world = x_sensor / magnification;
        const float rgbin[3] = {(float)color[(j * fb_width + i) * 4 + 0],
            (float)color[(j * fb_width + i) * 4 + 1],
            (float)color[(j * fb_width + i) * 4 + 2]};

        float L_in = spectrum_rgb_to_p(lambda, rgbin);

        int num_samples = max(1, (int)(L_in * sample_mul));
        float sample_weight = L_in / num_samples;

          // Sampling aperture creates the fringing effect
        for (int sample = 0; sample < num_samples; ++sample) {
          float x_ap, y_ap;
          do {
            x_ap = (rand() / (float)RAND_MAX - 0.5) * 2 * r_pupil;
            y_ap = (rand() / (float)RAND_MAX - 0.5) * 2 * r_pupil;
          } while (x_ap * x_ap + y_ap * y_ap > r_pupil * r_pupil);

          float in[5], out[4];

          in[0] = x_world + pixel_size * (rand() / (float)RAND_MAX - 0.5);
          in[1] = x_ap;
          in[2] = y_ap;

          system_y.evaluate(in, out);

          out[0] = out[0] * sensor_scaling + sensor_xres / 2;
          out[1] = out[1] * sensor_scaling + sensor_yres / 2;

          float lambert = sqrt(1 - out[2]);
          if (lambert != lambert)
            lambert = 0; // NaN check

          set_linear_atXY(lambert * sample_weight * rgb[0 + 3 * ll], out[0], out[1], 0, 0,
              color, new_color, fb_height, fb_width);            
          set_linear_atXY(lambert * sample_weight * rgb[1 + 3 * ll], out[0], out[1], 0, 1,
              color, new_color, fb_height, fb_width);
          set_linear_atXY(lambert * sample_weight * rgb[2 + 3 * ll], out[0], out[1], 0, 2,
              color, new_color, fb_height, fb_width);
        }
      });
    }

    tasking::parallel_for(fb_width * fb_height, [&](int pixelID ){
      int i = pixelID % fb_width;
      int j = pixelID / fb_width;
      float max_value = max(new_color[(j * fb_width + i) * 4 + 0],
            max(new_color[(j * fb_width + i) * 4 + 1],
                new_color[(j * fb_width + i) * 4 + 2]));
      for (int c = 0; c < 3; c++) {
        new_color[(j * fb_width + i) * 4 + c] = 
          max(new_color[(j * fb_width + i) * 4 + c], 0.02f * max_value);
        color[(j * fb_width + i) * 4 + c] = new_color[(j * fb_width + i) * 4 + c];
      }
    });
  }

  void set_linear_atXY(float value,
      float fx,
      float fy,
      const int z,
      int c,
      float *color,
      float *new_color,
      int height,
      int width)
  {
    int x = (int)fx - (fx >= 0 ? 0 : 1), nx = x + 1,
        y = (int)fy - (fy >= 0 ? 0 : 1), ny = y + 1;
    float dx = fx - x, dy = fy - y;
    if (y >= 0 && y < height) {
      if (x >= 0 && x < width) {
        float w1 = (1 - dx) * (1 - dy);
        new_color[(y * width + x) * 4 + c] =
            (float)(w1 * value + new_color[(y * width + x) * 4 + c]);
      }
      if (nx >= 0 && nx < width) {
        float w1 = dx * (1 - dy);
        new_color[(y * width + nx) * 4 + c] =
            (float)(w1 * value + new_color[(y * width + nx) * 4 + c]);
      }
    }
    if (ny >= 0 && ny < height) {
      if (x >= 0 && x < width) {
        float w1 = (1 - dx) * dy;
        new_color[(ny * width + x) * 4 + c] =
            (float)(w1 * value + new_color[(ny * width + x) * 4 + c]);
      }
      if (nx >= 0 && nx < width) {
        float w1 = dx * dy;
        new_color[(ny * width + nx) * 4 + c] =
            (float)(w1 * value + new_color[(ny * width + nx) * 4 + c]);
      }
    }
  }

  Transform4f get_system_1(float lambda, int degree = 3, float distance = 50000) {
    OpticalMaterial glass1("N-LASF45", false);
    OpticalMaterial glass2("N-BK7", false);
    
    float d0 = distance;
    const float R1 = 18.47;
    const float d1 = 6.40;
    const float R2 = 31.25;
    const float d2 = 13.00;
    const float R3 = -31.25;
    const float d3 = 6.40;
    const float R4 = -18.47;

    return two_plane_5(d0, degree)
        >> refract_spherical_5(R1, 1.f, glass1.get_index(lambda), degree)
        >> propagate_5(d1, degree) >> refract_spherical_5(
          R2, glass1.get_index(lambda), glass2.get_index(lambda), degree)
        >> propagate_5(d2, degree) >> refract_spherical_5(
          R3, glass2.get_index(lambda), glass1.get_index(lambda), degree)
        >> propagate_5(d3, degree)
        >> refract_spherical_5(R4, glass1.get_index(lambda), 1.f, degree);
  
  }

  Transform4f get_system(float lambda, int degree = 3, float distance = 50000)
  {
    // Let's simulate Edmund Optics achromat #NT32-921:
    /* Clear Aperture CA (mm) 	39.00
      Eff. Focal Length EFL (mm) 	120.00
      Back Focal Length BFL (mm) 	111.00
      Center Thickness CT 1 (mm) 	9.60
      Center Thickness CT 2 (mm) 	4.20
      Radius R1 (mm) 	65.22
      Radius R2 (mm) 	-62.03
      Radius R3 (mm) 	-1240.67
      Substrate 	N-SSK8/N-SF10
    */

    OpticalMaterial glass1("N-SSK8", false);
    OpticalMaterial glass2("N-SF10", false);

    // Also try: const float d0 = 5000; // Scene is 5m away
    float d0 = distance; // Scene is 5km away
    const float R1 = 65.22;
    const float d1 = 9.60;
    const float R2 = -62.03;
    const float d2 = 4.20;
    const float R3 = -1240.67;

    return two_plane_5(d0, degree)
        >> refract_spherical_5(R1, 1.f, glass1.get_index(lambda), degree)
        >> propagate_5(d1, degree) >> refract_spherical_5(
            R2, glass1.get_index(lambda), glass2.get_index(lambda), degree)
        >> propagate_5(d2, degree)
        >> refract_spherical_5(R3, glass2.get_index(lambda), 1.f, degree);
  }
};

PolyParallelFrameOp::PolyParallelFrameOp() {}

PolyParallelFrameOp::~PolyParallelFrameOp() {}

std::unique_ptr<LiveImageOp> PolyParallelFrameOp::attach(FrameBufferView &fbView)
{
  return rkcommon::make_unique<LivePolyParallelFrameOp>(fbView);
}

std::string PolyParallelFrameOp::toString() const
{
  return "ospray::PolyParallelFrameOp";
}

} // namespace ospray

extern "C" OSPError OSPRAY_DLLEXPORT ospray_module_init_polyparallel(
    int16_t versionMajor, int16_t versionMinor, int16_t /*versionPatch*/)
{
  auto status = ospray::moduleVersionCheck(versionMajor, versionMinor);

  if (status == OSP_NO_ERROR)
    ospray::ImageOp::registerType<ospray::PolyParallelFrameOp>("polyparallel");

  return status;
}
