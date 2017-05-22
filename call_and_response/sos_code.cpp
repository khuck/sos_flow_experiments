#include "globals.h"
#include <sstream>
#include <vector>
#include <iterator>
#include <unistd.h>

void do_fork(std::string forkCommand) {
    std::istringstream iss(forkCommand);
    std::vector<std::string> tokens;
    copy(std::istream_iterator<std::string>(iss),
         std::istream_iterator<std::string>(),
         std::back_inserter(tokens));
    const char **args = (const char **)calloc(tokens.size()+1, sizeof(char*));
    for (int i = 0; i < tokens.size() ; i++) {
        args[i] = tokens[i].c_str();
    }
    int rc = execvp(args[0],const_cast<char* const*>(args));
    if (rc < 0) {
        perror("\nError in execvp");
    }
    // exit the daemon spawn!
    std::cout << "Daemon exited!" << std::endl;
    _exit(0);
}

void fork_exec_sosd_shutdown(void) {
    // first, figure out who should fork a daemon on this node
    int i;
    if (_commrank == _daemon_rank) {
        int pid = vfork();
        if (pid == 0) {
            char* forkCommand;
            forkCommand = getenv ("SOS_FORK_SHUTDOWN");
            if (forkCommand) {
                std::cout << "Rank " << _commrank << " stopping SOS daemon(s): " << forkCommand << std::endl;
                std::string foo(forkCommand);
                do_fork(foo);
            } else {
                std::cout << "Please set the SOS_FORK_SHUTDOWN environment variable to stop SOS in the background." << std::endl;
            }
        }
    }
    //
    // wait until it is running
    //
    //wait(2);
}

void send_shutdown_message(void) {
    int i;
    SOS_buffer     *buffer;
    SOS_msg_header  header;
    int offset;
    if (_commrank == _daemon_rank) {
    	setenv("SOS_SHUTDOWN", "1", 1);

        SOS_buffer_init(_runtime, &buffer);

        header.msg_size = -1;
        header.msg_type = SOS_MSG_TYPE_SHUTDOWN;
        header.msg_from = _runtime->my_guid;
        header.pub_guid = 0;

        offset = 0;
        const char * format = "iigg";
        SOS_buffer_pack(buffer, &offset, const_cast<char*>(format),
                header.msg_size,
                header.msg_type,
                header.msg_from,
                header.pub_guid);

        header.msg_size = offset;
        offset = 0;
        const char * format2 = "i";
        SOS_buffer_pack(buffer, &offset, const_cast<char*>(format2), header.msg_size);

        std::cout << "Sending SOS_MSG_TYPE_SHUTDOWN ..." << std::endl;

        SOS_send_to_daemon(buffer, buffer);

        SOS_buffer_destroy(buffer);
    }
}

void fork_exec_sosd(void) {
    // first, figure out who should fork a daemon on this node
    int i;
    // get my hostname
    const int hostlength = 128;
    char hostname[hostlength] = {0};
    gethostname(hostname, sizeof(char)*hostlength);
    std::cout << hostname << std::endl;
    // make array for all hostnames
    char * allhostnames = (char*)calloc(hostlength*_commsize, sizeof(char));
    // copy my name into the big array
    char * host_index = allhostnames + (hostlength * _commrank);
    strncpy(host_index, hostname, hostlength);
    // get all hostnames
    PMPI_Allgather(hostname, hostlength, MPI_CHAR, allhostnames, 
                   hostlength, MPI_CHAR, MPI_COMM_WORLD);
    _daemon_rank = 0;
    // point to the head of the array
    host_index = allhostnames;
    // find the lowest rank with my hostname
    for (i = 0 ; i < _commsize ; i++) {
        //printf("%d:%d comparing '%s' to '%s'\n", _commrank, _commsize, hostname, host_index);
        if (strncmp(hostname, host_index, hostlength) == 0) {
            _daemon_rank = i;
            break;
        }
        host_index = host_index + hostlength;
    }
    // fork the daemon
    if (_commrank == _daemon_rank) {
        _shutdown_daemon = true;
        int pid = vfork();
        if (pid == 0) {
            char* forkCommand = NULL;
            char* ranks_per_node = NULL;
            char* offset = NULL;
            forkCommand = getenv ("SOS_FORK_COMMAND");
            std::cout << "forkCommand " << forkCommand << std::endl;
            ranks_per_node = getenv ("SOS_APP_RANKS_PER_NODE");
            std::cout << "ranks_per_node " << ranks_per_node << std::endl;
            offset = getenv ("SOS_LISTENER_RANK_OFFSET");
            std::cout << "offset " << offset << std::endl;
            if (forkCommand) {
                std::string custom_command(forkCommand);
                size_t index = 0;
                index = custom_command.find("@LISTENER_RANK@", index);
                if (index != std::string::npos) {
                    if (ranks_per_node) {
                        int rpn = atoi(ranks_per_node);
                        int listener_rank = _commrank / rpn;
                        if(offset) {
                            int off = atoi(offset);
                            listener_rank = listener_rank + off;
                        }
                        std::stringstream ss;
                        ss << listener_rank;
                        custom_command.replace(index,15,ss.str());
                    }
                }
                std::cout << "Rank " << _commrank << " spawning SOS daemon(s): " << custom_command << std::endl;
                do_fork(custom_command);
            } else {
                std::cerr << "Please set the SOS_FORK_COMMAND environment variable to spawn SOS in the background." << std::endl;
            }
        }
    }
    //
    // wait until it is running
    //
    //wait(2);
}

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
