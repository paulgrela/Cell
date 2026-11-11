
#ifndef CELL_ENGINE_THREAD_PINNING_TO_CPU_CORE_H
#define CELL_ENGINE_THREAD_PINNING_TO_CPU_CORE_H

#include <map>
#include <vector>
#include <fstream>
#include <filesystem>

namespace CellEngineThreadPinning
{
    inline int ReadIntFromFile(const std::string& Path)
    {
        std::ifstream File(Path);
        int Value = 0;
        if (File)
            File >> Value;
        return Value;
    }

    inline int GetNumaNodeOfCPU(const int Cpu)
    {
        std::error_code Error;
        for (const auto& Entry : std::filesystem::directory_iterator("/sys/devices/system/cpu/cpu" + std::to_string(Cpu), Error))
        {
            const std::string Name = Entry.path().filename().string();
            if (Name.starts_with("node") && Name.size() > 4 && Name.find_first_not_of("0123456789", 4) == std::string::npos)
                return std::stoi(Name.substr(4));
        }
        return 0;
    }

    inline std::vector<int> GetCPUsForThreads()
    {
        cpu_set_t AllowedCpus;
        CPU_ZERO(&AllowedCpus);
        sched_getaffinity(0, sizeof(AllowedCpus), &AllowedCpus);

        std::map<std::tuple<int, int, int>, std::vector<int>> CPUsOfPhysicalCore;
        for (int Cpu = 0; Cpu < CPU_SETSIZE; Cpu++)
            if (CPU_ISSET(Cpu, &AllowedCpus))
            {
                const std::string Topology = "/sys/devices/system/cpu/cpu" + std::to_string(Cpu) + "/topology/";
                CPUsOfPhysicalCore[{ GetNumaNodeOfCPU(Cpu), ReadIntFromFile(Topology + "physical_package_id"), ReadIntFromFile(Topology + "core_id") }].push_back(Cpu);
            }

        std::vector<int> CPUs;
        for (std::size_t HardwareThreadOfCore = 0; ; HardwareThreadOfCore++)
        {
            bool AnyAdded = false;
            for (const auto& LogicalCpus : CPUsOfPhysicalCore | views::values)
                if (HardwareThreadOfCore < LogicalCpus.size())
                {
                    CPUs.push_back(LogicalCpus[HardwareThreadOfCore]);
                    AnyAdded = true;
                }
            if (AnyAdded == false)
                break;
        }
        return CPUs;
    }

    inline bool BindCurrentThreadToCPU(const int Cpu)
    {
        cpu_set_t CpuSet;
        CPU_ZERO(&CpuSet);
        CPU_SET(Cpu, &CpuSet);
        return pthread_setaffinity_np(pthread_self(), sizeof(CpuSet), &CpuSet) == 0;
    }
}

#endif