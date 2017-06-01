#include "globals.h"
#include <sstream>
#include <vector>
#include <iterator>
#include <unistd.h>

void initialize(int * argc, char *** argv) {
    static bool initialized = false;
    if (!initialized) {
        _runtime = NULL;
        //printf("init() trying to connect...\n");
        SOS_init(argc, argv, &_runtime, SOS_ROLE_CLIENT, SOS_RECEIVES_NO_FEEDBACK, NULL);
        if(_runtime == NULL) {
            printf("%d Unable to connect to SOS daemon. Determining whether to spawn...\n", _commrank);
            fork_exec_sosd();
        }
        int repeat = 10;
        while(_runtime == NULL) {
            sleep(2);
            _runtime = NULL;
            //printf("init() trying to connect...\n");
            SOS_init(argc, argv, &_runtime, SOS_ROLE_CLIENT, SOS_RECEIVES_NO_FEEDBACK, NULL);
            if (_runtime != NULL) {
                printf("%d Connected to SOS daemon. Continuing...\n", _commrank);
                break;
            } else if (--repeat < 0) { 
                printf("%d Unable to connect to SOS daemon. Failing...\n", _commrank);
                return;
            }
        }
        initialized = true;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    make_pub();
}

void make_pub() {
        char pub_name[SOS_DEFAULT_STRING_LEN] = {0};
        char app_version[SOS_DEFAULT_STRING_LEN] = {0};

        //printf("[make_pub]: Creating new pub...\n");

        _runtime->config.comm_rank = _commrank;
        _runtime->config.comm_size = _commsize;

        sprintf(pub_name, "call and response test");
        sprintf(app_version, "v0.alpha");
        SOS_pub_create(_runtime, &_sos_pub, pub_name, SOS_NATURE_DEFAULT);

        strcpy(_sos_pub->prog_ver, app_version);
        _sos_pub->meta.channel       = 1;
        _sos_pub->meta.layer         = SOS_LAYER_LIB;
        // sos_pub->meta.pri_hint      = SOS_PRI_IMMEDIATE;
        // sos_pub->meta.scope_hint    = SOS_SCOPE_SELF;
        // sos_pub->meta.retain_hint   = SOS_RETAIN_SESSION;

        //printf("[make_pub]:   ... done.  (pub->guid == %ld)\n", _sos_pub->guid);
        //printf("[make_pub]: Announcing the pub...\n");
        SOS_announce(_sos_pub);
}

void finalize(void) {
    if (_runtime == NULL) { return; }
    static bool finalized = false;
    //printf("%s\n", __func__); fflush(stdout);
    if (finalized) return;
    // shutdown the daemon, if necessary
    if (_shutdown_daemon) {
        // to make sure we have a clean shutdown, allow the queues to drain.
        sleep(2);
        send_shutdown_message();
        // shouldn't be necessary, but sometimes the shutdown message is ignored?
        //fork_exec_sosd_shutdown();
    }
    std::cout << "Finalize SOS" << std::endl;
    SOS_finalize(_runtime);
    finalized = true;
}

void sample_value(std::string name, double value) {
  //std::cout << "Sending data : " << name << ": " << value << std::endl;
  SOS_pack(_sos_pub, name.c_str(), SOS_VAL_TYPE_DOUBLE, &value);
}

void flush_it(void) {
  SOS_publish(_sos_pub);
}
