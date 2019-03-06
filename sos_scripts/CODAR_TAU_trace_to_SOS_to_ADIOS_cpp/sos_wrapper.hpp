/*
 * C++ program to extract TAU trace data from SOS listeners, and 
 * write it out to an ADIOS2 data file/stream.
 */

#pragma once

#include "utils.hpp"
#include "sos.h"
#include <nlohmann/json.hpp>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string.h>

using json = nlohmann::json;

namespace extractor {

/* Class containing SOS connection info */
class sos {
    private:
        bool connected;
        SOS_runtime * my_sos;
        json config;
    public:
        sos(json& _config) : 
            connected(false), my_sos(nullptr), config(_config) {
            PRINTSTACK
            connect();
        };
        ~sos() {
            PRINTSTACK
            disconnect();
        };
        void connect() {
            PRINTSTACK
            SOS_init(&(my_sos), SOS_ROLE_ANALYTICS, SOS_RECEIVES_NO_FEEDBACK, nullptr);
            if (my_sos == nullptr) {
                fprintf(stderr, "ERROR: Unable to contact the SOS daemon. Terminating..\n");
                fflush(stderr);
                exit (EXIT_FAILURE);
            }
            connected = true;
            srandom(my_sos->my_guid);
            if (config["sosd"]["SOS_CMD_PORT"] != nullptr) {
                strcpy(my_sos->daemon->remote_port,
                    config["sosd"]["SOS_CMD_PORT"].get<std::string>().c_str());
            }
            test_connection();
        }
		void test_connection() {
			SOS_buffer *request;
			SOS_buffer *reply;
			SOS_buffer_init_sized_locking(my_sos, &request,
					SOS_DEFAULT_BUFFER_MAX, false);
			SOS_buffer_init_sized_locking(my_sos, &reply,
					SOS_DEFAULT_BUFFER_MAX, false);

			SOS_msg_header header;
			header.msg_size = -1;
			header.msg_type = SOS_MSG_TYPE_PROBE;
			header.msg_from = my_sos->my_guid;
			header.ref_guid = 0;

			int offset = 0;
			SOS_msg_zip(request, header, 0, &offset);
			header.msg_size = offset;
			offset = 0;
			SOS_msg_zip(request, header, 0, &offset);
			SOS_buffer_wipe(reply);
			SOS_send_to_daemon(request, reply);
			offset = 0;
			SOS_msg_unzip(reply, &header, 0, &offset);

			uint64_t queue_depth_local       = 0;
			uint64_t queue_depth_cloud       = 0;
			uint64_t queue_depth_db_tasks    = 0;
			uint64_t queue_depth_db_snaps    = 0;

			SOS_buffer_unpack(reply, &offset, "gggg",
					&queue_depth_local,
					&queue_depth_cloud,
					&queue_depth_db_tasks,
					&queue_depth_db_snaps);
			SOS_buffer_destroy(request);
			SOS_buffer_destroy(reply);

		}
		void disconnect() {
            PRINTSTACK
            if (connected && my_sos != nullptr) {
                SOS_finalize(my_sos);
            }
            connected = false;
            my_sos = nullptr;
        };
};

}; // end namespace extractor
