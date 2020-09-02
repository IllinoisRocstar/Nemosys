#include <gtest.h>
#include <math.h>

#include "Omega_h_adapt.hpp"
#include "Omega_h_array_ops.hpp"
#include "Omega_h_build.hpp"
#include "Omega_h_metric.hpp"

#include "oshGeoMesh.H"

#include <nvml.h>

TEST(oshGeoMeshCUDA, MemoryUsage) {
  pid_t pid = getpid();
  nvmlInit_v2();
  unsigned int deviceCount;
  nvmlDeviceGetCount_v2(&deviceCount);
  nvmlDevice_t device;
  unsigned long long oldUsedGpuMemory = 0;
  bool found_process = false;
  unsigned int infoCount = 0;
  std::vector<nvmlProcessInfo_t> infos;
  for (int i = 0; i < deviceCount; ++i) {
    nvmlDeviceGetHandleByIndex_v2(0, &device);
    nvmlDeviceGetComputeRunningProcesses(device, &infoCount, nullptr);
    infoCount *= 2;
    infos.resize(infoCount);
    auto err =
        nvmlDeviceGetComputeRunningProcesses(device, &infoCount, infos.data());
    for (int j = 0; j < infoCount; ++j) {
      if (infos[j].pid == pid) {
        // Only care about error for the device this process is running on.
        EXPECT_EQ(err, NVML_SUCCESS);
        oldUsedGpuMemory = infos[j].usedGpuMemory;
        found_process = true;
        break;
      }
    }
    if (found_process) {
      break;
    }
  }
  std::cout << "Current GPU memory used (MiB): "
            << oldUsedGpuMemory / std::pow(1024, 2) << std::endl;
  ASSERT_GT(oldUsedGpuMemory, 0);

  auto lib = NEM::MSH::OmegaHInterface::GetLibrary();
  auto oshMesh =
      Omega_h::build_box(lib->world(), OMEGA_H_SIMPLEX, 1.0, 1.0, 1.0, 5, 5, 5);
  auto metrics = Omega_h::get_implied_isos(&oshMesh);
  auto scalar = Omega_h::get_metric_scalar_for_nelems(&oshMesh, metrics, 100);
  metrics = Omega_h::multiply_each_by(metrics, scalar);
  oshMesh.add_tag(Omega_h::VERT, "metric", 1, metrics);
  auto opts = Omega_h::AdaptOpts(&oshMesh);
  Omega_h::adapt(&oshMesh, opts);
  oshMesh.remove_tag(Omega_h::VERT, "metric");
  NEM::MSH::oshGeoMesh ogm{&oshMesh};

  unsigned long long newUsedGpuMemory = 0;
  infoCount = 0;
  nvmlDeviceGetComputeRunningProcesses(device, &infoCount, nullptr);
  infoCount *= 2;
  infos.resize(infoCount);
  nvmlDeviceGetComputeRunningProcesses(device, &infoCount, infos.data());
  for (int i = 0; i < infoCount; ++i) {
    if (infos[i].pid == pid) {
      newUsedGpuMemory = infos[i].usedGpuMemory;
      break;
    }
  }
  std::cout << "GPU memory used after Omega_h test (MiB): "
            << newUsedGpuMemory / std::pow(1024, 2) << std::endl;
  ASSERT_GT(newUsedGpuMemory, oldUsedGpuMemory);
}

int main(int argc, char *argv[]) {
  ::testing::InitGoogleTest(&argc, argv);

  int res = RUN_ALL_TESTS();

  return res;
}
