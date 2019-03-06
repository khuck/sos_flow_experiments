/*
 * C++ program to extract TAU trace data from SOS listeners, and 
 * write it out to an ADIOS2 data file/stream.
 */

#include "utils.hpp"
#include "adios2.h"

namespace extractor {

/* Class containing ADIOS archive info */
class adios {
    private:
        bool opened;
        json config;
        adios2::ADIOS ad;
        adios2::IO bpIO;
        adios2::Engine bpWriter;
    public:
        adios(json& _config) : 
            opened(false),
            config(_config)
        {
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
            open();
        };
        ~adios() {
            PRINTSTACK
            close();
        };
		void open() {
			std::stringstream ss;
			ss << config["adios"]["outputdir"].get<std::string>();
			ss << "/";
			ss << config["adios"]["filename"].get<std::string>();
			printf("Writing %s\n", ss.str().c_str());
			bpWriter = bpIO.Open(ss.str(), adios2::Mode::Write);
			opened = true;
		}
        void close() {
            PRINTSTACK
    		if (opened) {
        		bpWriter.Close();
        		opened = false;
    		}
        };
};

}; // end namespace extractor
