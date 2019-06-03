#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <argp.h>
#include <sys/resource.h>
#include "threshold_functions.h"

#define MAX_NAME 1024

char commandline[]= "mrtV2 <-c configfile> [input_filename]";

static struct argp_option options[] = {
  {"configfile", 'c', "configuration_file",    0,   "Specify name of configuration file"},
  { 0 }
};

struct arguments{
  char config_filename[MAX_NAME];
  char *i_filename;
};

static error_t parse_opt(int key, char *arg, struct argp_state *state){
  struct arguments *arguments = state->input;
  switch(key){
  case 'c':
    if (arg[0] != '-')
      strncpy(arguments->config_filename,arg,MAX_NAME);
    else {
      printf("arguments should not begin with a - sign.  Exiting...\n");
      exit(0);
    }
    break;
  case ARGP_KEY_ARG:
    arguments->i_filename = arg;
    break;
  case ARGP_KEY_END:
    if(state->arg_num<1)
      argp_usage(state);
    break;
  default:
    return ARGP_ERR_UNKNOWN;
  }
  return 0;
}

static struct argp argp = {options,parse_opt,"$Id: mrtV2.c,v 1.0 2019-03-02 15:43:15 elyons Exp $"};

int main(int argc, char* argv[])
{
  if(argc == 1) {
    printf("Usage: %s\n", commandline);
    exit(-1);
  }
  
  struct arguments arguments;
  arguments.config_filename[0] = '\0';
  argp_parse(&argp,argc,argv,0,0,&arguments);

  /*locate and verify the config file*/
  if (arguments.config_filename[0] == '\0') {
    /*locate and verify the config file*/
    char *mrt_home = getenv("MRTHOME");
    if (mrt_home == NULL) {
      printf("please declare MRTHOME or specify configfile on command line\n");
      exit(0);
    }
    else 
      sprintf(arguments.config_filename, "%s/mrtV2_config.txt", mrt_home); 
  }

  printf("config_filename: %s\n", arguments.config_filename);

  const rlim_t sSz = 64L * 1024L * 1024L;   // min stack size = 64 Mb
  struct rlimit rlim;
  int stat;

  stat = getrlimit(RLIMIT_STACK, &rlim);
  if (stat == 0) {
    if (rlim.rlim_cur < sSz) {
      rlim.rlim_cur = sSz;
      stat = setrlimit(RLIMIT_STACK, &rlim);
      if (stat != 0) {
	fprintf(stderr, "setrlimit returned result = %d\n", stat);
      }
    }
  }
  int mrtV2( char i_filename[], char config_fn[]  );
  mrtV2(arguments.i_filename, arguments.config_filename);

  return 0;
}
