/*
 * C++ program to extract TAU trace data from SOS listeners, and 
 * write it out to an ADIOS2 data file/stream.
 */

#include "adios_wrapper.hpp"
#include "adios2.h"

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

void adios::declare_variables(void) {
/*
        g = ad.declare_group("TAU_metrics", "", ad.FLAG.YES)
        ad.define_var(g, "program_count", "", ad.DATATYPE.unsigned_integer, "", "", "")
        ad.define_var(g, "comm_rank_count", "", ad.DATATYPE.unsigned_integer, "", "", "")
        ad.define_var(g, "thread_count", "", ad.DATATYPE.unsigned_integer, "", "", "")
        ad.define_var(g, "event_type_count", "", ad.DATATYPE.unsigned_integer, "", "", "")
        ad.define_var(g, "timer_count", "", ad.DATATYPE.unsigned_integer, "", "", "")
        ad.define_var(g, "timer_event_count", "", ad.DATATYPE.unsigned_integer, "", "", "")
        ad.define_var(g, "event_timestamps", "", ad.DATATYPE.unsigned_long, "timer_event_count,6", "timer_event_count,6", "0,0")
        ad.define_var(g, "counter_count", "", ad.DATATYPE.unsigned_integer, "", "", "")
        ad.define_var(g, "counter_event_count", "", ad.DATATYPE.unsigned_integer, "", "", "")
        ad.define_var(g, "counter_values", "", ad.DATATYPE.unsigned_long, "counter_event_count,6", "counter_event_count,6", "0,0")
        ad.define_var(g, "comm_count", "", ad.DATATYPE.unsigned_integer, "", "", "")
        ad.define_var(g, "comm_timestamps", "", ad.DATATYPE.unsigned_long, "comm_count,8", "comm_count,8", "0,0")
*/
    const std::size_t Nx = 1;
    const adios2::Dims shape{static_cast<size_t>(Nx)};
    const adios2::Dims start{static_cast<size_t>(Nx)};
    const adios2::Dims count{Nx};
    bpIO.DefineVariable<int>("program_count", shape, start, count, adios2::ConstantDims);
    bpIO.DefineVariable<int>("comm_rank_count", shape, start, count, adios2::ConstantDims);
    bpIO.DefineVariable<int>("thread_count", shape, start, count, adios2::ConstantDims);
    bpIO.DefineVariable<int>("event_type_count", shape, start, count, adios2::ConstantDims);
    bpIO.DefineVariable<int>("timer_count", shape, start, count, adios2::ConstantDims);
    bpIO.DefineVariable<int>("timer_event_count", shape, start, count, adios2::ConstantDims);
    bpIO.DefineVariable<int>("counter_count", shape, start, count, adios2::ConstantDims);
    bpIO.DefineVariable<int>("counter_event_count", shape, start, count, adios2::ConstantDims);
    bpIO.DefineVariable<int>("comm_count", shape, start, count, adios2::ConstantDims);

/*
    Nx = 6;
    bpIO.DefineVariable<int>("event_timestamps", shape, start, count, adios2::ConstantDims);
    bpIO.DefineVariable<int>("counter_values", shape, start, count, adios2::ConstantDims);
    Nx = 8;
    bpIO.DefineVariable<int>("comm_timestamps", shape, start, count, adios2::ConstantDims);
    */
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

void adios::define_attribute(std::string name, std::string value, int nid, int tid) {
    std::stringstream ss;
    ss << "TAU:" << nid << ":" << tid << ":MetaData:" << name;
    bpIO.DefineAttribute<std::string>(ss.str(), value.c_str());
}

void adios::define_attribute(std::string name, std::string value) {
    bpIO.DefineAttribute<std::string>(name, value);
}

}; // end namespace extractor
