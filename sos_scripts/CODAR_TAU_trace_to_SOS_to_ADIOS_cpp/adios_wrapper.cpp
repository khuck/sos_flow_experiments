/*
 * C++ program to extract TAU trace data from SOS listeners, and 
 * write it out to an ADIOS2 data file/stream.
 */

#include "adios_wrapper.hpp"
#include "adios2.h"
#include <unordered_set>
#include <vector>
#include "sos_wrapper.hpp"

namespace extractor {

void adios::initialize() {
    /** ADIOS class factory of IO class objects, DebugON is recommended */
    ad = adios2::ADIOS(true);
    /*** IO class object: settings and factory of Settings: Variables,
     * Parameters, Transports, and Execution: Engines */
    bpIO = ad.DeclareIO("TAU trace data from SOS");
    // if not defined by user, we can change the default settings
    // BPFile is the default engine
    bpIO.SetEngine(config["adios"]["adios_method"].get<std::string>());
    bpIO.SetParameters({{"num_threads", "2"}});

    // ISO-POSIX file output is the default transport (called "File")
    // Passing parameters to the transport
    bpIO.AddTransport("File", {{"Library", "POSIX"}});
}

void adios::define_variables(void) {
    const std::size_t Nx = 1;
    const adios2::Dims shape{static_cast<size_t>(Nx)};
    const adios2::Dims start{static_cast<size_t>(Nx)};
    const adios2::Dims count{Nx};
    program_count = bpIO.DefineVariable<int>("program_count");
    comm_rank_count = bpIO.DefineVariable<int>("comm_rank_count");
    thread_count = bpIO.DefineVariable<int>("thread_count");
    event_type_count = bpIO.DefineVariable<int>("event_type_count");
    timer_count = bpIO.DefineVariable<int>("timer_count");
    timer_event_count = bpIO.DefineVariable<size_t>("timer_event_count");
    counter_count = bpIO.DefineVariable<int>("counter_count");
    counter_event_count = bpIO.DefineVariable<size_t>("counter_event_count");
    comm_count = bpIO.DefineVariable<size_t>("comm_count");

    event_timestamps = bpIO.DefineVariable<unsigned long>("event_timestamps", {1, 6}, {0, 0}, {1, 6});
    counter_values = bpIO.DefineVariable<unsigned long>("counter_values", {1, 6}, {0, 0}, {1, 6});
    comm_timestamps = bpIO.DefineVariable<unsigned long>("comm_timestamps", {1, 8}, {0, 0}, {1, 8});
}

void adios::open() {
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
    if (opened) {
        bpWriter.Close();
        opened = false;
    }
};

void adios::define_attribute(std::string name, std::string value) {
    static std::unordered_set<std::string> seen;
    if (seen.count(name) == 0) {
        seen.insert(name);
        bpIO.DefineAttribute<std::string>(name, value);
    }
}

void adios::write_variables(sos& my_sos,
    size_t num_timer_values,
    size_t num_counter_values,
    size_t num_comm_values,
    std::vector<unsigned long>& timer_values_array,
    std::vector<unsigned long>& counter_values_array,
    std::vector<unsigned long>& comm_values_array) 
{
    TAU_SCOPED_TIMER_FUNC();
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

    if (num_timer_values > 0) {
        event_timestamps.SetShape({num_timer_values, 6});
        const adios2::Dims timer_start{0, 0};
        const adios2::Dims timer_count{static_cast<size_t>(num_timer_values), 6};
        const adios2::Box<adios2::Dims> timer_selection{timer_start, timer_count};
        event_timestamps.SetSelection(timer_selection);
        bpWriter.Put(event_timestamps, timer_values_array.data());
    }

    if (num_counter_values > 0) {
        counter_values.SetShape({num_counter_values, 6});
        const adios2::Dims counter_start{0, 0};
        const adios2::Dims counter_count{static_cast<size_t>(num_counter_values), 6};
        const adios2::Box<adios2::Dims> counter_selection{counter_start, counter_count};
        counter_values.SetSelection(counter_selection);
        bpWriter.Put(counter_values, counter_values_array.data());
    }

    if (num_comm_values > 0) {
        comm_timestamps.SetShape({num_comm_values, 8});
        const adios2::Dims comm_start{0, 0};
        const adios2::Dims comm_count{static_cast<size_t>(num_comm_values), 8};
        const adios2::Box<adios2::Dims> comm_selection{comm_start, comm_count};
        comm_timestamps.SetSelection(comm_selection);
        bpWriter.Put(comm_timestamps, comm_values_array.data());
    }

    bpWriter.EndStep();
}

}; // end namespace extractor
