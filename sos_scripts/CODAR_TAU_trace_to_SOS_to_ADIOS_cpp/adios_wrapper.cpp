/*
 * C++ program to extract TAU trace data from SOS listeners, and 
 * write it out to an ADIOS2 data file/stream.
 */

#include "adios_wrapper.hpp"
#include "adios2.h"
#include <unordered_set>
#include "sos_wrapper.hpp"

namespace extractor {

void adios::initialize() {
    PRINTSTACK
    /** ADIOS class factory of IO class objects, DebugON is recommended */
    ad = adios2::ADIOS(true);
    /*** IO class object: settings and factory of Settings: Variables,
     * Parameters, Transports, and Execution: Engines */
    bpIO = ad.DeclareIO("TAU trace data from SOS");
    // if not defined by user, we can change the default settings
    // BPFile is the default engine
    bpIO.SetEngine(config["adios"]["adios_method"].get<std::string>());
    bpIO.SetParameters({{"num_threads", "1"}});

    // ISO-POSIX file output is the default transport (called "File")
    // Passing parameters to the transport
    bpIO.AddTransport("File", {{"Library", "POSIX"}});
}

void adios::define_variables(void) {
    PRINTSTACK
    const std::size_t Nx = 1;
    const adios2::Dims shape{static_cast<size_t>(Nx)};
    const adios2::Dims start{static_cast<size_t>(Nx)};
    const adios2::Dims count{Nx};
    program_count = bpIO.DefineVariable<int>("program_count", shape, start, count, adios2::ConstantDims);
    comm_rank_count = bpIO.DefineVariable<int>("comm_rank_count", shape, start, count, adios2::ConstantDims);
    thread_count = bpIO.DefineVariable<int>("thread_count", shape, start, count, adios2::ConstantDims);
    event_type_count = bpIO.DefineVariable<int>("event_type_count", shape, start, count, adios2::ConstantDims);
    timer_count = bpIO.DefineVariable<int>("timer_count", shape, start, count, adios2::ConstantDims);
    timer_event_count = bpIO.DefineVariable<int>("timer_event_count", shape, start, count, adios2::ConstantDims);
    counter_count = bpIO.DefineVariable<int>("counter_count", shape, start, count, adios2::ConstantDims);
    counter_event_count = bpIO.DefineVariable<int>("counter_event_count", shape, start, count, adios2::ConstantDims);
    comm_count = bpIO.DefineVariable<int>("comm_count", shape, start, count, adios2::ConstantDims);

    event_timestamps = bpIO.DefineVariable<long>("event_timestamps", {1, 6}, {0, 0}, {1, 6});
    counter_values = bpIO.DefineVariable<long>("counter_values", {1, 6}, {0, 0}, {1, 6});
    comm_timestamps = bpIO.DefineVariable<long>("comm_timestamps", {1, 8}, {0, 0}, {1, 8});
}

void adios::open() {
    PRINTSTACK
    if (!opened) {
        std::stringstream ss;
        ss << config["adios"]["outputdir"].get<std::string>();
        ss << "/";
        ss << config["adios"]["filename"].get<std::string>();
        printf("Writing %s\n", ss.str().c_str());
        bpWriter = bpIO.Open(ss.str(), adios2::Mode::Write);
        opened = true;
    }
}

void adios::close() {
    PRINTSTACK
    if (opened) {
        bpWriter.Close();
        opened = false;
    }
};

void adios::define_attribute(std::string name, std::string value) {
    PRINTSTACK
    static std::unordered_set<std::string> seen;
    if (seen.count(name) == 0) {
        seen.insert(name);
        bpIO.DefineAttribute<std::string>(name, value);
    }
}

void adios::write_variables(sos& my_sos,
    int num_timer_values,
    int num_counter_values,
    int num_comm_values,
    std::vector<long>& timer_values_array,
    std::vector<long>& counter_values_array,
    std::vector<long>& comm_values_array) 
{
    int programs = my_sos.get_prog_count();
    int comm_ranks = my_sos.get_comm_rank_count();
    int threads = my_sos.get_thread_count();
    int event_types = my_sos.get_event_type_count();
    int timers = my_sos.get_timer_count();
    int counters = my_sos.get_counter_count();

    bpWriter.BeginStep();

    bpWriter.Put(program_count, &programs);
    bpWriter.Put(comm_rank_count, &comm_ranks);
    bpWriter.Put(thread_count, &threads);
    bpWriter.Put(event_type_count, &event_types);
    bpWriter.Put(timer_count, &timers);
    bpWriter.Put(timer_event_count, &num_timer_values);
    bpWriter.Put(counter_count, &counters);
    bpWriter.Put(counter_event_count, &num_counter_values);
    bpWriter.Put(comm_count, &num_comm_values);

    event_timestamps.SetShape({(size_t) num_timer_values});
    bpWriter.Put(event_timestamps, timer_values_array.data());

    counter_values.SetShape({(size_t) num_counter_values});
    bpWriter.Put(counter_values, counter_values_array.data());

    comm_timestamps.SetShape({(size_t) num_comm_values});
    bpWriter.Put(comm_timestamps, comm_values_array.data());

    bpWriter.EndStep();
}

}; // end namespace extractor
