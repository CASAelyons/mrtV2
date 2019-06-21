/* $Id: mrt.c,v 0.1 2015-03-01 12:02:34 elyons Exp $ */
/* Copyright 2015 University of Massachusetts Amherst all rights reserved*/

/* MRT (merged reflectivity thresholding) is an algorithm used to identify areas of fixed reflectivity levels. It generates various forms of output, json, kml, text, and netcdf.
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <unistd.h>
#include <time.h>
#include <netcdf.h>
#include <math.h>
#include <jansson.h>
#include <libconfig.h>
#include "threshold_functions.h"

int mrtV2(char i_filename[], char config_fn[]) {

  /*locate and verify the config file*/
  config_t mrt_config;
  config_setting_t *gns_params;
  config_setting_t *alg_params;

  config_init(&mrt_config);
  if (!config_read_file(&mrt_config, config_fn)) {
    printf("config file not found or contains errors.  Should be in $MRTHOME directory and called mrt_config.txt. exiting...\n");
    config_destroy(&mrt_config);
    return -1;
  }
  
  /* Check whether this will report to the MCC via sockets */
  int mcc_output = 0;
  struct sockaddr_in mcc_addr;
  int mcc_sock;
  int mcc_port;
  const char *mcc_host;
  int connected = 0;
  
  if (config_lookup_bool(&mrt_config, "mcc_output", &mcc_output))
    printf("MCC output: %d\n", mcc_output);
  else
    printf("No mcc_output definition in mrt_config.  MCC output disabled.\n");

  if (mcc_output == 1) {
    if (!config_lookup_int(&mrt_config, "mcc_port", &mcc_port)) {
      printf("No mcc_port definition in mrt_config. Exiting...\n");
      exit(1);
    }
    
    if (!config_lookup_string(&mrt_config, "mcc_host", &mcc_host)) {
      printf("No mcc_host definition in mrt_config. Exiting...\n");
      exit(1);
    }
  }
  
  /* Check whether this will report to the ACS */
  int gns_output = 0;
  int num_gns_hosts = 0;
  int num_gns_ports = 0;
  const config_setting_t *gns_host_arr;
  const config_setting_t *gns_port_arr;
  
  if (config_lookup_bool(&mrt_config, "gns_output", &gns_output))
    printf("GNS output: %d\n", gns_output);
  else
    printf("No gns_output definition in mrt_config.  GNS output disabled.\n");
  
  if (gns_output == 1) {
    gns_params = config_lookup(&mrt_config, "gns_params");

    if (gns_params != NULL) {
      gns_host_arr = config_setting_lookup(gns_params, "gns_host");
      if (gns_host_arr != NULL) {
	num_gns_hosts = config_setting_length(gns_host_arr);
      }
      else {
	printf("GNS_output is true, but no GNS hosts defined.  Exiting.\n");
	exit(-1);
      }
      gns_port_arr = config_setting_lookup(gns_params, "gns_port");
      if (gns_port_arr != NULL) {
	num_gns_ports = config_setting_length(gns_port_arr);
      }
      else {
	printf("GNS_output is true, but no GNS ports defined.  Exiting.\n");
	exit(-1);
      }
      if (num_gns_ports != num_gns_hosts) {
	printf("Number of GNS ports does not match number of GNS hosts.  Please ensure you have one port entry for each host in the configuration file.  Exiting.\n");
	exit(-1);
      }  
    }
    else {
      printf("No gns_params in config file.  Exiting.\n");
      exit(-1);
    }
  }

  struct sockaddr_in gns_addr[num_gns_hosts];
  int gns_sock[num_gns_hosts];
  long int gns_port[num_gns_ports];
  const char *gns_host[num_gns_hosts];
  int sock_connected[num_gns_hosts];
  int gh;
  for (gh = 0; gh < num_gns_hosts; gh++) {
    sock_connected[gh] = 0;
    gns_port[gh] = config_setting_get_int_elem(gns_port_arr, gh);
    gns_host[gh] = config_setting_get_string_elem(gns_host_arr, gh);
    if ((!gns_port[gh]) || (!gns_host[gh])) {
      printf("GNS port or host not properly defined.  Please check config file.  Exiting...\n");
      exit(-1);
    }
    if ((gns_sock[gh]=socket(AF_INET, SOCK_STREAM, 0))==-1){
      perror("socket");
      exit(-1);
    }
    
    struct timeval timeout = {2, 0};
    
    if (setsockopt(gns_sock[gh], SOL_SOCKET, SO_SNDTIMEO, (struct timeval *)&timeout,sizeof(struct timeval)) < 0) {
      perror("setsockopt failed\n");
      //exit(-1);
    }
    
    memset((char *) &gns_addr[gh],'\0',sizeof(gns_addr[gh]));
    gns_addr[gh].sin_family = AF_INET;
    gns_addr[gh].sin_port = htons(gns_port[gh]);
    gns_addr[gh].sin_addr.s_addr = inet_addr(gns_host[gh]);
    if (connect(gns_sock[gh],(struct sockaddr *) &gns_addr[gh],
                sizeof (struct sockaddr_in))<0){
      perror("connect");
      //exit(-1);
    }
    else {
      sock_connected[gh] = 1;
    }
  }

  /* Check whether netcdf file output will be generated */
  int netcdf_output = 0;

  if (config_lookup_bool(&mrt_config, "netcdf_output", &netcdf_output))
    printf("netcdf output: %d\n", netcdf_output);
  else
    printf("No netcdf_output definition in mrt_config.  Netcdf output disabled.\n");
  
  /* Check whether a json data file will be generated */
  int json_output = 0;
  char *jsonbuf;
  //FILE *jsonfile;
  char jsonname[128];
  const char *jsondir;
  //json_t *json_coords_arr;//  = json_array();
  json_t *jsonroot = json_object();
  json_t *features_json_arr = json_array();
  json_object_set_new( jsonroot, "type", json_string( "FeatureCollection" ) );
  json_object_set_new( jsonroot, "features", features_json_arr );
    
  if (config_lookup_bool(&mrt_config, "json_output", &json_output)){
    printf("JSON output: %d\n", json_output);
  }
  else 
    printf("No json_output definition in mrt_config.  JSON output disabled.\n");

  if (config_lookup_string(&mrt_config, "json_directory", &jsondir)){
    printf("JSON directory: %s\n", jsondir);
  }
  else {
    printf("No json_directory definition in mrt_config.  JSON output disabled.\n");
    json_output = 0;
  }
  
  /*read the alg_params */
  int num_contour_levels = 0;
  alg_params = config_lookup(&mrt_config, "alg_params");
  int min_contour_points;
  int min_contour_points_present;
  const config_setting_t *contour_levels;
  if (alg_params != NULL) {
    contour_levels = config_setting_lookup(alg_params, "contour_levels");
    if (contour_levels != NULL) {
      num_contour_levels = config_setting_length(contour_levels);
    }
    else {
      printf("No contour_levels defined.  Exiting.\n");
      exit(-1);
    }
    min_contour_points_present = config_setting_lookup_int(alg_params, "min_contour_points", &min_contour_points);
    if (min_contour_points_present == CONFIG_FALSE) {
      printf("min_contour_points not defined. No size threshold will be used.");
      min_contour_points = 0;
    }
    else {
      printf("Minimum contour points: %d\n", min_contour_points);
    }
  }
  else {
    printf("No alg_params in config file.  Exiting.\n");
    exit(-1);
  }
  float contours[num_contour_levels];
  int cl;
  for (cl = 0; cl < num_contour_levels; cl++) {
    contours[cl] = config_setting_get_float_elem(contour_levels, cl);
    printf("Contour level %d: %f\n", cl, contours[cl]);
  }

  //Read in input data from netcdf
  float **ref_in;
  float *lat_in;
  float *lon_in;
  
  float n_lat, w_lon, lat_spac, lon_spac;
  
  int status,ncid,latid,lonid,ref_id,h,i,j;
  char varname[128];
  size_t num_lats, num_lons;
  
  status = nc_open(i_filename, NC_NOWRITE, &ncid);
  handle_error(status);
  if (status != NC_NOERR)
    exit(-1);
  
  status = nc_inq_dimid(ncid, "Lat", &latid);    
  handle_error(status);

  status = nc_inq_dimlen(ncid, latid, &num_lats);
  handle_error(status);
  
  status = nc_inq_dimid(ncid, "Lon", &lonid);    
  handle_error(status);

  status = nc_inq_dimlen(ncid, lonid, &num_lons);
  handle_error(status);
  
  lat_in = malloc(num_lats*sizeof(float));
  lon_in = malloc(num_lons*sizeof(float));

  ref_in = malloc(num_lons*sizeof(float *));
  ref_in[0] = malloc(num_lons*num_lats*sizeof(float));
  for(i=0; i < num_lons; ++i)
    ref_in[i] = ref_in[0]+i*num_lats;

  status = nc_inq_varname(ncid, 0, varname);
  handle_error(status);

  status = nc_inq_varid(ncid,varname,&ref_id);
  handle_error(status);

  status = nc_get_var_float(ncid, ref_id, &ref_in[0][0]);
  handle_error(status);
    
  status = nc_get_att_float(ncid, NC_GLOBAL, "Latitude", &n_lat);
  handle_error(status);

  status = nc_get_att_float(ncid, NC_GLOBAL, "Longitude", &w_lon);
  handle_error(status);

  status = nc_get_att_float(ncid, NC_GLOBAL, "LatGridSpacing", &lat_spac);
  handle_error(status);

  status = nc_get_att_float(ncid, NC_GLOBAL, "LonGridSpacing", &lon_spac);
  handle_error(status);

  int v;
  for (v=1; v < num_lons+1; v++) 
    lon_in[v-1] = (w_lon + (lon_spac*(v-1)));
  
  //for (v=1; v < num_lats+1; v++)
  //  lat_in[v-1] = (n_lat - (lat_spac*(v-1))); 
  for (v=0; v < num_lats; v++)
    lat_in[v] = (n_lat - (lat_spac*((num_lats - v)-1)));
  


  size_t filename_len;
  filename_len = strlen(i_filename);
  char starttime[16];
  char yyyy[4];
  char mon[3];
  char dd[3];
  char hh[3];
  char mm[3];
  char ss[3];
  int ymd_start = (int)filename_len - 18;
  memcpy(yyyy, &i_filename[ymd_start], 4);
  yyyy[4] = '\0';
  memcpy(mon, &i_filename[ymd_start+4], 2);
  mon[2] = '\0';
  memcpy(dd, &i_filename[ymd_start+6], 2);
  dd[2] = '\0';
  memcpy(hh, &i_filename[ymd_start+9], 2);
  hh[2] = '\0';
  memcpy(mm, &i_filename[ymd_start+11], 2);
  mm[2] = '\0';
  memcpy(ss, &i_filename[ymd_start+13], 2);
  ss[2] = '\0';
  starttime[15] = '\0';
  sprintf(starttime, "%s%s%s-%s%s%s", yyyy, mon, dd, hh, mm, ss);
  
  char alt_starttime[32];
  char mcctime[16];
  
  sprintf(mcctime,"%s%s%s%s%s%s", yyyy, mon, dd, hh, mm, ss);
  sprintf(alt_starttime, "%s-%s-%sT%s:%s:%sZ", yyyy, mon, dd, hh, mm, ss);
  
  char name_string[32];
  char filemin[8];
  char *s_tok = strtok(i_filename, "_");
  int s_tok_ind = 0;
  while (s_tok != NULL) {
    if (s_tok_ind == 1) {
      size_t minstrsz = strlen(s_tok);
      size_t minsz = minstrsz - 3;
      char minstr[minsz+1];
      memcpy(minstr, &s_tok[0], minsz);
      minstr[minsz] = '\0';
      sprintf(filemin, "%s", minstr);
    }
    s_tok = strtok(NULL, "_");
    s_tok_ind++;
  }
  
  sprintf(name_string, "STORM_CASA_%s", filemin);

  //Only connect to MCC if this is the 1 minute file
  if ((strstr(varname, "PredictedReflectivity_1min") != NULL) && (mcc_output == 1)) {
    if ((mcc_sock=socket(AF_INET, SOCK_STREAM, 0))==-1){
      perror("socket");
    }
    else {
      memset((char *) &mcc_addr,'\0',sizeof(mcc_addr));
      mcc_addr.sin_family = AF_INET;
      mcc_addr.sin_port = htons(mcc_port);
      mcc_addr.sin_addr.s_addr = inet_addr(mcc_host);
      if (connect(mcc_sock,(struct sockaddr *) &mcc_addr, sizeof (struct sockaddr_in))<0) {
	perror("connect");
      }
      else
	connected = 1;
    }
  }
  
  if (json_output == 1) {
    sprintf(jsonname, "%s/mrt_%s_%s.geojson", jsondir, name_string, starttime);
  }


  //prepare for conrec implementation
  
  //declare, allocate memory, set lataxis/lonaxis
  int *lonaxis, *lataxis;
  lataxis = malloc(num_lats*sizeof(int));
  lonaxis = malloc(num_lons*sizeof(int));
  
  for(i=0; i<num_lats; ++i)
    lataxis[i]=i;

  for(i=0; i<num_lons; ++i) 
    lonaxis[i]=i;


  //declare, allocate memory, initialize up/down/left/right arrays
  struct index_pair ***up;
  struct index_pair ***down;
  struct index_pair ***left;
  struct index_pair ***right;
  
  up = (struct index_pair ***)malloc(sizeof(struct index_pair **) * num_contour_levels);
  for (i=0; i < num_contour_levels; i++) {
    up[i] = (struct index_pair **)malloc(sizeof(struct index_pair *) * num_lons);
    for (j=0; j < num_lons; ++j)
      up[i][j] = (struct index_pair *)malloc(sizeof(struct index_pair) * num_lats);
  }

  left = (struct index_pair ***)malloc(sizeof(struct index_pair **) *num_contour_levels);
  for (i=0; i <num_contour_levels; i++) {
    left[i] =(struct index_pair **)malloc(sizeof(struct index_pair *) * num_lons);
    for(j=0; j< num_lons; ++j)
      left[i][j] = (struct index_pair*)malloc(sizeof(struct index_pair) * num_lats);
  }

  down = (struct index_pair ***)malloc(sizeof(struct index_pair **) *num_contour_levels);
  for (i=0; i <num_contour_levels; i++) {
    down[i] =(struct index_pair **)malloc(sizeof(struct index_pair *) * num_lons);
    for(j=0; j< num_lons; ++j)
      down[i][j] = (struct index_pair*)malloc(sizeof(struct index_pair) * num_lats);
  }

  right = (struct index_pair ***)malloc(sizeof(struct index_pair **) *num_contour_levels);
  for (i=0; i <num_contour_levels; i++) {
    right[i] =(struct index_pair **)malloc(sizeof(struct index_pair *) * num_lons);
    for(j=0; j< num_lons; ++j)
      right[i][j] = (struct index_pair*)malloc(sizeof(struct index_pair) * num_lats);
  }
  
  for (h=0; h<num_contour_levels; ++h) {
    for (i=0; i<num_lons; ++i) {
      for (j=0; j<num_lats; ++j) {
      up[h][i][j].xind = -1;
      up[h][i][j].yind = -1;
      down[h][i][j].xind = -1;
      down[h][i][j].yind = -1;
      left[h][i][j].xind = -1;
      left[h][i][j].yind = -1;
      right[h][i][j].xind = -1;
      right[h][i][j].yind = -1;
      }
    }
  }

  //now call the Conrec Contouring function
  
  Contour(ref_in, 0, (num_lats - 1), 0, (num_lons - 1), lataxis, lonaxis, num_contour_levels, contours, up, down, left, right);

  //begin processing Conrec output
  
  //keep track of your number of valid polygons
  int num_valid_contours = 0;

  //create some names in advance to call these polygons
  for (h=0; h<num_contour_levels; h++) {
    printf("processing contour level %d: %f\n", h, contours[h]);

    char hazard_string[128];
    if (strstr(varname, "PredictedReflectivity_15min") != NULL) {
      if(contours[h] > 39.9)
	sprintf(hazard_string, "STORM_INTENSE_CASA_15");
      else
	sprintf(hazard_string, "STORM_CASA_15");
    }
    else if (strstr(varname, "PredictedReflectivity_10min") != NULL) {
      if(contours[h] > 39.9)
	sprintf(hazard_string, "STORM_INTENSE_CASA_10");
      else
	sprintf(hazard_string, "STORM_CASA_10");
    }
    else if (strstr(varname, "PredictedReflectivity_5min") != NULL) {
      if(contours[h] > 39.9)
	sprintf(hazard_string, "STORM_INTENSE_CASA_5");
      else
	sprintf(hazard_string, "STORM_CASA_5");
    }
    else if (strstr(varname, "PredictedReflectivity_1min") != NULL) {
      if(contours[h] > 39.9)
	sprintf(hazard_string, "STORM_INTENSE_CASA_0");
      else
	sprintf(hazard_string, "STORM_CASA_0");
    }
    else {
      sprintf(hazard_string, "STORM_INTENSE_CASA");
    }

    
    int k; //counter
    int numpoints=0; //number of points in this particular polygon
    int numcontours=0; //number of polygons
    int numwaypoints = 0; //number of points in this particular linestring
    int numlinestrings=0; //number of linestrings
    
    struct index_pair reference; //the i&j index reference of the cell that led us to the current cell
    reference.xind = -1; 
    reference.yind = -1;
    
    struct index_pair *contour_ips; //an array of the index points in the current contour
    struct index_pair *reference_ips; //an array of index points that led us to this point
    
    contour_ips = (struct index_pair*)malloc(sizeof(struct index_pair));
    reference_ips = (struct index_pair*)malloc(sizeof(struct index_pair));

    struct index_pair **linestrings = NULL;
    int *waypoints = NULL;
    
    struct index_pair null_ip;
    null_ip.xind = -1;
    null_ip.yind = -1;
    
    for (i=0; i<num_lats;i++) {
      for (j=0; j<num_lons;j=j+0) {

	//Print out the arrays...
	//printf("%d %d %d %d %d %d %d %d %d %d\n", i, j, up[h][i][j].xind, up[h][i][j].yind, down[h][i][j].xind, down[h][i][j].yind, left[h][i][j].xind, left[h][i][j].yind, right[h][i][j].xind, right[h][i][j].yind);
	
	//current_ip index_pair representing this i,j pair                                                                                                             
	struct index_pair current_ip;
	current_ip.xind = i;
	current_ip.yind = j;
	
	//if all the arrays are uninitialized we either go to the reference or move on
	if (((up[h][i][j].xind == -1) && (down[h][i][j].xind == -1) && (left[h][i][j].xind == -1) && (right[h][i][j].xind == -1)) && (ip_is_equal(reference, null_ip))) {
	  //case 1
	  ++j;
	  continue;
	}
	
	//view the arrays
	//printf("%d %d %d %d %d %d %d %d %d %d \n", i, j, up[h][i][j].xind, up[h][i][j].yind, down[h][i][j].xind, down[h][i][j].yind, left[h][i][j].xind, left[h][i][j].yind, right[h][i][j].xind, right[h][i][j].yind);
	//printf("reference %d %d\n", reference.xind, reference.yind);
	
	//now that that's out of the way, let's see if we have a point that referred us here...
	//if not this must be the first point in a new contour
	if(reference.xind == -1) {
	  //case 2
	  int nexti;
	  int nextj;
	  
	  //printf("this is the first point in a new contour\n");
	  //printf("numpoints: %d\n", numpoints);
	  //set the first point in this contour to this index point
	  contour_ips[0] = current_ip;

	  //set the reference of point 0 to null
	  reference_ips[0] = null_ip;

	  //now we start down the path
	  //subcase a
	  if(down[h][i][j].xind != -1) {
	    nexti = down[h][i][j].xind;
	    nextj = down[h][i][j].yind;
	    
	    //wipe out the path so we don't traverse it again
	    down[h][i][j].xind = -1;
	    down[h][i][j].yind = -1;
	  }
	  else if (up[h][i][j].xind != -1){
	    nexti = up[h][i][j].xind;
	    nextj = up[h][i][j].yind;
	  
	    up[h][i][j].xind = -1;
	    up[h][i][j].yind = -1;
	  }
	  else if (right[h][i][j].xind != -1){
	    nexti = right[h][i][j].xind;
	    nextj = right[h][i][j].yind;
	  
	    right[h][i][j].xind = -1;
	    right[h][i][j].yind = -1;
	  }
	  else if (left[h][i][j].xind != -1){
	    nexti = left[h][i][j].xind;
	    nextj = left[h][i][j].yind;
	  
	    left[h][i][j].xind = -1;
	    left[h][i][j].yind = -1;
	  }
	  else {
	    //subcase b
	    printf("should have been handled by case 1a\n");
	  }

	  //set the reference point and move on
	  reference = current_ip;
	  i = nexti;
	  j = nextj;

	  continue;
	}
	
	//ok so if we made it here we have a contour with at least one point already, and a reference
	//if this point is not already contained in the ip array of this contour, we are not closed yet
	//or it could be a multipath and we intentionally didn't add it to the array
	if (!contained(numpoints+1, current_ip, contour_ips)) {
	  //case 3
	  //printf("case 3\n");
	  int nexti;
	  int nextj;

	  //let's wipe out the path to the reference
	  if (ip_is_equal(up[h][i][j], reference)) 
	    up[h][i][j] = null_ip;
	  else if (ip_is_equal(down[h][i][j], reference))
	    down[h][i][j] = null_ip;
	  else if (ip_is_equal(left[h][i][j], reference))
	    left[h][i][j] = null_ip;
	  else if (ip_is_equal(right[h][i][j], reference))
	    right[h][i][j] = null_ip;
	  
	  //if there is a fork in the road, take it
	  //subcase a
	  
	  if (up[h][i][j].xind != -1) {
	    //printf("up %d %d ref %d %d\n", up[h][i][j].xind, up[h][i][j].yind, reference.xind, reference.yind); 
	    nexti = up[h][i][j].xind;
	    nextj = up[h][i][j].yind;
	    up[h][i][j] = null_ip;
	  }
	  else if (left[h][i][j].xind != -1) {
	    //printf("left %d %d ref %d %d\n", left[h][i][j].xind, left[h][i][j].yind, reference.xind, reference.yind);
	    nexti = left[h][i][j].xind;
	    nextj = left[h][i][j].yind;
	    left[h][i][j] = null_ip;
	  }
	  else if (down[h][i][j].xind != -1) {
	    //printf("down %d %d ref %d %d\n", down[h][i][j].xind, down[h][i][j].yind, reference.xind, reference.yind);
	    nexti = down[h][i][j].xind;
	    nextj = down[h][i][j].yind;
	    down[h][i][j] = null_ip;
	  }
	  else if (right[h][i][j].xind != -1) {
	    //printf("right %d %d ref %d %d\n", right[h][i][j].xind, right[h][i][j].yind, reference.xind, reference.yind);
	    nexti = right[h][i][j].xind;
	    nextj = right[h][i][j].yind;
	    right[h][i][j] = null_ip;
	  }
	  else {
	    //if no paths exist, create linestring back to last point with a possible path forward
	    //subcase b
	    //printf("subcase b\n");
	    linestrings = (struct index_pair**)realloc(linestrings, (numlinestrings + 1) * sizeof(struct index_pair*));
	    linestrings[numlinestrings] = (struct index_pair*)malloc((numwaypoints+1)*sizeof(struct index_pair));
	    linestrings[numlinestrings][numwaypoints] = current_ip;
	    //printf("numlinestrings %d numwaypoints %d ip %d %d\n", numlinestrings, numwaypoints, linestrings[numlinestrings][numwaypoints].xind, linestrings[numlinestrings][numwaypoints].yind);
	    //printf("numpoints %d    no further path.  Return to reference %d %d\n", numpoints, reference.xind, reference.yind);
	    numwaypoints++;
	    
	    i = reference.xind;
	    j = reference.yind;
	    continue;
	  }

	  //printf("seems like a legit point\n");
	  
	  numpoints++;
	  
	  reference_ips = (struct index_pair*)realloc(reference_ips, (numpoints+1)*sizeof(struct index_pair));
	  reference_ips[numpoints] = reference;
	  
	  //increase the size of the dynamic contour array
	  contour_ips = (struct index_pair*)realloc(contour_ips, (numpoints+1)*sizeof(struct index_pair));

	  //add this ip to our list of ips in this contour
	  contour_ips[numpoints] = current_ip;
	  //printf("next i %d next j %d\n", nexti, nextj);
	  i = nexti;
	  j = nextj;
	  
	  reference = current_ip;
	  continue;
	}

	else if (contained(numpoints+1, current_ip, contour_ips)) {
	  //case 4
	  //printf("case 4\n");
	  //if it's here we've returned to a point already contained in the contour
	  
	  if (ip_is_equal(reference, current_ip)) {
	    //subcase a
	    //printf("subcase a\n");
	    //if it's here we've been sent back from a previous contour attempt
	    
	    //add this point to linestring
	    linestrings[numlinestrings] = (struct index_pair*)realloc(linestrings[numlinestrings], (numwaypoints+1)*sizeof(struct index_pair));
	    linestrings[numlinestrings][numwaypoints] = current_ip;
	    //printf("numpoints %d numlinestrings %d numwaypoints %d ip %d %d\n", numpoints, numlinestrings, numwaypoints, linestrings[numlinestrings][numwaypoints].xind, linestrings[numlinestrings][numwaypoints].yind);
	    numwaypoints++;
	    //if we have more paths available, or if this is the first point, break off the linestring
	    if (((up[h][i][j].xind > -1) || (down[h][i][j].xind > -1) || (left[h][i][j].xind > -1) || (right[h][i][j].xind > -1)) || (ip_is_equal(reference_ips[numpoints], null_ip))) {
	      //printf("subsubcase 1\n");
	      waypoints = (int*)realloc(waypoints, (numlinestrings+1)*sizeof(int));
	      waypoints[numlinestrings] = numwaypoints;
	      int a;
	      //for (a=0; a<numwaypoints; a++) 
	      //printf("linestring point %d index %d %d\n", a, linestrings[numlinestrings][a].xind, linestrings[numlinestrings][a].yind);
	      ++numlinestrings;
	      numwaypoints = 0;
	     
	    }

	    
	    
	    //if we have more paths, stay right here
	    if((up[h][i][j].xind > -1) || (down[h][i][j].xind > -1) || (left[h][i][j].xind > -1) || (right[h][i][j].xind > -1)) {
	      reference = reference_ips[numpoints];
	    }
	    else if (ip_is_equal(reference_ips[numpoints], null_ip)){
	      reference = null_ip;
	      j = j+1;
	    }
	    else {
	      reference = reference_ips[numpoints];
	      i = reference.xind;
	      j = reference.yind;
	    }
	    
	    //remove the point from the arrays
	    if (numpoints > 0) {
	      numpoints--;
	      reference_ips = (struct index_pair*)realloc(reference_ips, (numpoints+1)*sizeof(struct index_pair));
	      contour_ips = (struct index_pair*)realloc(contour_ips, (numpoints+1)*sizeof(struct index_pair));
	    }
	    continue;
	  }

	  else {
	    //subcase b
	    //if the point is contained and this was not the reference, we have a polygon
	    //but not necessarily the polygon from which we started

	    if (ip_is_equal(current_ip, contour_ips[0])) {
	      //subsubcase 1 
	      //printf("closing original polygon...\n");
	      json_t *json_coords_arr = json_array();
	      //increment the number of contours
	      ++numcontours;
	      //printf("numcontours: %d\n", numcontours);
	      for (k=0; k<numpoints+1; ++k) {
		//print out the contour
		printf("contour %d pair %d i %d j %d\n", numcontours, k, contour_ips[k].xind, contour_ips[k].yind);
	      }
	      
	      //dump contour to output 
	      //double dxbetween;
	      int a;
	      if (numpoints > min_contour_points) {
		++num_valid_contours;
		for (a=0; a < numpoints+1; a++) {
		  if (json_output == 1) {
		    json_t * json_coords = json_pack("[f,f]", lon_in[contour_ips[a].yind], lat_in[contour_ips[a].xind]);
		    json_array_append_new(json_coords_arr, json_coords);
		  }
		}
	      
		/* now append the first point back to the end to conform with geoJSON polygon standard */
		if (json_output == 1) {
		  json_t * json_coords = json_pack("[f,f]", lon_in[contour_ips[0].yind], lat_in[contour_ips[0].xind]);
		  json_array_append_new(json_coords_arr, json_coords);
		}
		
		char contour_level_string[16];
		sprintf(contour_level_string, "%.2fdBZ", contours[h]);
		
		char id_string[128];
		if (strstr(varname, "PredictedReflectivity_15min") != NULL) {
		  sprintf(id_string, "STORM_CASA_15_%s_%s_%d", starttime, contour_level_string, numcontours);
		}
		else if (strstr(varname, "PredictedReflectivity_10min") != NULL) {
		  sprintf(id_string, "STORM_CASA_10_%s_%s_%d", starttime, contour_level_string, numcontours);
		}
		else if (strstr(varname, "PredictedReflectivity_5min") != NULL) {
		  sprintf(id_string, "STORM_CASA_5_%s_%s_%d", starttime, contour_level_string, numcontours);
		}
		else if (strstr(varname, "PredictedReflectivity_1min") != NULL) {
		  sprintf(id_string, "STORM_CASA_0_%s_%s_%d", starttime, contour_level_string, numcontours);
		}
		else {
		  sprintf(id_string, "STORM_CASA_%s_%s_%d", starttime, contour_level_string, numcontours);
		}
		if (json_output == 1) {
		  json_t *polygon = json_pack("{s:s,s:s,s:s,s:{s:s,s:f},s:{s:s,s:[o]}}", "type", "Feature", "id", id_string, "hazardType", hazard_string, "properties", "timestamp", alt_starttime, "ReflectivityLevel", contours[h], "geometry", "type", "Polygon", "coordinates", json_coords_arr);
		  json_array_append_new(features_json_arr, polygon);
		}
		
		//if ((mcc_output == 1) && (numpoints > 8)) {
		//  for (y=0; y < numpoints; ++y) {
		//    char detect[128];
		//    sprintf(detect, "reflectivity :lat %.4f :long %.4f :height 1000.0 :time %s\n", ll_array[y].lat, ll_array[y].lon, mcctime);
		//    if (connected) {
		//      write(mcc_sock,detect,strlen(detect));
		//    }
		// }
		//}
	      }
	      
	      //return to the start of the first contour 
	      i = contour_ips[0].xind;
	      j = contour_ips[0].yind + 1;

	      //clear the path to the reference if it exists
	      if (ip_is_equal(up[h][i][j], reference))
		up[h][i][j] = null_ip;
	      else if (ip_is_equal(down[h][i][j], reference))
		down[h][i][j] = null_ip;
	      else if (ip_is_equal(left[h][i][j], reference))
		left[h][i][j] = null_ip;
	      else if (ip_is_equal(right[h][i][j], reference))
		right[h][i][j] = null_ip;
	      
	      //clear the reference 
	      reference.xind = -1;
	      reference.yind = -1;
	    
	      //reset the numpoints to zero
	      numpoints = 0;
	      
	      //free up and malloc reference_ips, contours_ips, and ll_array to the minimum
	      free(reference_ips);
	      reference_ips = (struct index_pair*)malloc(sizeof(struct index_pair));
	      
	      free(contour_ips);
	      contour_ips = (struct index_pair*)malloc(sizeof(struct index_pair));
	      
	      //free(ll_array);
	      //ll_array = (struct latLon*)malloc(sizeof(struct latLon));
	      
	      continue;
	    }
	    else {
	      //subsubcase 2
	      //printf("closing some incidental polygon we came across on the original path\n");
	      json_t *json_coords_arr = json_array();
	      ++numcontours;
	      //printf("numcontours: %d\n", numcontours);

	      int this_index;
	      for (k=0; k<numpoints+1; ++k) {
		if (ip_is_equal(current_ip, contour_ips[k])) 
		  this_index = k; 
	      }
	      for (k=this_index; k<numpoints+1; ++k) {
		//print out the contour
		printf("contour %d pair %d i %d j %d\n", numcontours, k, contour_ips[k].xind, contour_ips[k].yind);
	      }
	      
	      int a;
	      
	      if (numpoints > min_contour_points) {
		++num_valid_contours;
		for (a=this_index; a < numpoints+1; a++) {
		  if (json_output == 1) {
		    json_t * json_coords = json_pack("[f,f]", lon_in[contour_ips[a].yind], lat_in[contour_ips[a].xind]);
		    json_array_append_new(json_coords_arr, json_coords);
		  }
		}

		/* now append the first point of this contour back to the end to conform with geoJSON polygon standard */
		if (json_output == 1) {
		  json_t * json_coords = json_pack("[f,f]", lon_in[contour_ips[this_index].yind], lat_in[contour_ips[this_index].xind]);
		  json_array_append_new(json_coords_arr, json_coords);
		}

		char contour_level_string[16];
		sprintf(contour_level_string, "%.2fdBZ", contours[h]);

		char id_string[128];
		if (strstr(varname, "PredictedReflectivity_15min") != NULL) {
		  sprintf(id_string, "STORM_CASA_15_%s_%s_%d", starttime, contour_level_string, numcontours);
		}
		else if (strstr(varname, "PredictedReflectivity_10min") != NULL) {
		  sprintf(id_string, "STORM_CASA_10_%s_%s_%d", starttime, contour_level_string, numcontours);
		}
		else if (strstr(varname, "PredictedReflectivity_5min") != NULL) {
		  sprintf(id_string, "STORM_CASA_5_%s_%s_%d", starttime, contour_level_string, numcontours);
		}
		else if (strstr(varname, "PredictedReflectivity_1min") != NULL) {
		  sprintf(id_string, "STORM_CASA_0_%s_%s_%d", starttime, contour_level_string, numcontours);
		}
		else {
		  sprintf(id_string, "STORM_CASA_%s_%s_%d", starttime, contour_level_string, numcontours);
		}
		if (json_output == 1) {
                  json_t *polygon = json_pack("{s:s,s:s,s:s,s:{s:s,s:f},s:{s:s,s:[o]}}", "type", "Feature", "id", id_string, "hazardType", hazard_string, "properties", "timestamp",alt_starttime, "ReflectivityLevel", contours[h], "geometry", "type", "Polygon", "coordinates", json_coords_arr);
                  json_array_append_new(features_json_arr, polygon);
                }

                //if ((mcc_output == 1) && (numpoints > 8)) {
                //  for (y=0; y < numpoints; ++y) {
		//   char detect[128];
                //    sprintf(detect, "reflectivity :lat %.4f :long %.4f :height 1000.0 :time %s\n", ll_array[y].lat, ll_array[y].lon, mcctime);
                //    if (connected) {
                //      write(mcc_sock,detect,strlen(detect));
                //    }
                //  }
                //}
              }
	      
	      //dump the points from this incidental polygon from the original contour list
	      //for (k=this_index+1; k<numpoints+1; ++k) {
	      int this_poly_numpoints = (numpoints+1) - this_index;
	      //printf("this_poly_numpoints %d  numpoints %d this_index %d\n", this_poly_numpoints, numpoints, this_index);
	      numpoints = numpoints - this_poly_numpoints;
	      
	      i = contour_ips[this_index].xind;
	      j = contour_ips[this_index].yind;
	      //set the reference to whatever ip had led us to the first point
	      reference = reference_ips[this_index];

	      //clear out the reference_ips and contour_ips that were in this incidental polygon
	      reference_ips = (struct index_pair*)realloc(reference_ips, (numpoints+1)*sizeof(struct index_pair));
	      contour_ips = (struct index_pair*)realloc(contour_ips, (numpoints+1)*sizeof(struct index_pair));

	      continue;
	      
	    }
	  }
	}
	printf("if you see this message, something unexpected has occurred\n");
      }
    }

    //ok lets arrange our linestrings
    //some counters
    int a; int b; int c;

    //we have to keep track of which linestrings have been combined already
    short *addressed = (short*)malloc(numlinestrings * sizeof(short));
    for(a=0; a<numlinestrings; a++)
      addressed[a] = 0;

    //printf("numlinestrings: %d\n", numlinestrings);
    
    for (a=0; a<numlinestrings; a++) {
      //if we've already addressed this one, just skip it
      if (addressed[a])
	continue;

      //otherwise we'll see if we have a linestring to combine it with
      int match = 0;
      //print out the line segment
      //printf("a: %d  waypoints: %d\n", a, waypoints[a]);
      //for (b=0; b<waypoints[a]; b++) {
      //printf("b: %d  x: %d  y: %d\n", b, linestrings[a][b].xind, linestrings[a][b].yind);
      //}

      
      if ((waypoints[a] > 2) && (((linestrings[a][0].xind == 0) && (linestrings[a][waypoints[a]-1].xind == 0)) || ((linestrings[a][0].yind == 0) && (linestrings[a][waypoints[a]-1].yind == 0)) || ((linestrings[a][0].xind == (int)(num_lats - 1)) && (linestrings[a][waypoints[a]-1].xind == (int)(num_lats - 1))) || ((linestrings[a][0].yind == (int)(num_lons - 1)) && (linestrings[a][waypoints[a]-1].xind == (int)(num_lats - 1))) || ((abs(linestrings[a][0].xind - linestrings[a][waypoints[a]-1].xind) < 2) && (abs(linestrings[a][0].yind - linestrings[a][waypoints[a]-1].yind) < 2)))) {
	//printf("waypoints %d\n", waypoints[a]);
	json_t *json_coords_arr = json_array();
	//increment the number of contours
	++numcontours;
	//printf("linestring %d is a polygon\n", a);
	//printf("numcontours: %d\n", numcontours);
	for (c=0; c<waypoints[a]; c++) {
	  printf("polygon c: %d  x: %d  y: %d\n", c, linestrings[a][c].xind, linestrings[a][c].yind);
	}

	int p;
	if (waypoints[p] > min_contour_points) {
	  ++num_valid_contours;
	  for (p=0; p < waypoints[a]; p++) {
	    if (json_output == 1) {
	      json_t * json_coords = json_pack("[f,f]", lon_in[linestrings[a][p].yind], lat_in[linestrings[a][p].xind]);
	      json_array_append_new(json_coords_arr, json_coords);
	    }
	  }

	  // now append the first point back to the end to conform with geoJSON polygon standard 
	  if (json_output == 1) {
	    json_t * json_coords = json_pack("[f,f]", lon_in[linestrings[a][0].yind], lat_in[linestrings[a][0].xind]);
	    json_array_append_new(json_coords_arr, json_coords);
	  }

	  char contour_level_string[16];
	  sprintf(contour_level_string, "%.2fdBZ", contours[h]);

	  char id_string[128];
	  if (strstr(varname, "PredictedReflectivity_15min") != NULL) {
	    sprintf(id_string, "STORM_CASA_15_%s_%s_%d", starttime, contour_level_string, numcontours);
	  }
	  else if (strstr(varname, "PredictedReflectivity_10min") != NULL) {
	    sprintf(id_string, "STORM_CASA_10_%s_%s_%d", starttime, contour_level_string, numcontours);
	  }
	  else if (strstr(varname, "PredictedReflectivity_5min") != NULL) {
	    sprintf(id_string, "STORM_CASA_5_%s_%s_%d", starttime, contour_level_string, numcontours);
	  }
	  else if (strstr(varname, "PredictedReflectivity_1min") != NULL) {
	    sprintf(id_string, "STORM_CASA_0_%s_%s_%d", starttime, contour_level_string, numcontours);
	  }
	  else {
	    sprintf(id_string, "STORM_CASA_%s_%s_%d", starttime, contour_level_string, numcontours);
	  }
	  if (json_output == 1) {
	    json_t *polygon = json_pack("{s:s,s:s,s:s,s:{s:s,s:f},s:{s:s,s:[o]}}", "type", "Feature", "id", id_string, "hazardType", hazard_string, "properties", "timestamp",alt_starttime, "ReflectivityLevel", contours[h], "geometry", "type", "Polygon", "coordinates", json_coords_arr);
	    json_array_append_new(features_json_arr, polygon);
	  }
	}
	addressed[a] = 1;
	continue;
      }  
      
      //otherwise we'll go through and check the other linestrings
      //once you get a match on one end we'll create a new linestring
      //this means we can increment "a" once we get a match
      for (b=a+1; b<numlinestrings; b++) {
	if (!match) {
	  //for (c=0; c<waypoints[b]; c++) {
	    //printf("c: %d  x: %d  y: %d\n", b, linestrings[b][c].xind, linestrings[b][c].yind);
	  //}
	  if ((abs(linestrings[a][0].xind - linestrings[b][0].xind) < 2) && (abs(linestrings[a][0].yind - linestrings[b][0].yind) < 2)) {
	    //printf("linestrings %d and %d startpoints are near eachother\n", a, b);
	    addressed[a] = 1;
	    addressed[b] = 1;
	    addressed =(short*)realloc(addressed, (numlinestrings+1)*sizeof(short));
	    addressed[numlinestrings] = 0;
	    linestrings = (struct index_pair**)realloc(linestrings, (numlinestrings+1)*sizeof(struct index_pair*));
	    linestrings[numlinestrings] = (struct index_pair *)malloc((waypoints[a]+waypoints[b])*sizeof(struct index_pair));
	    waypoints = (int*)realloc(waypoints, (numlinestrings+1)*sizeof(int));
	    waypoints[numlinestrings] = waypoints[a]+waypoints[b];
	    for (c=0; c<waypoints[b]; c++)
	      linestrings[numlinestrings][c] = linestrings[b][(waypoints[b] - c) - 1];
	    for (c=0; c<waypoints[a]; c++)
	      linestrings[numlinestrings][c+waypoints[b]] = linestrings[a][c];
	    numlinestrings++;
	    match = 1;
	    continue;
	  }

	  if ((abs(linestrings[a][0].xind - linestrings[b][waypoints[b]-1].xind) < 2) && (abs(linestrings[a][0].yind - linestrings[b][waypoints[b]-1].yind) < 2)) {
	    //printf("linestrings %d startpoint and linestring %d endpoint are near eachother\n", a, b);
	    addressed[a] = 1;
	    addressed[b] = 1;
	    addressed =(short*)realloc(addressed, (numlinestrings+1)*sizeof(short));
	    addressed[numlinestrings] = 0;
	    linestrings = (struct index_pair**)realloc(linestrings, (numlinestrings+1)*sizeof(struct index_pair*));
	    linestrings[numlinestrings] = (struct index_pair *)malloc((waypoints[a]+waypoints[b])*sizeof(struct index_pair));
	    waypoints = (int*)realloc(waypoints, (numlinestrings+1)*sizeof(int));
	    waypoints[numlinestrings] = waypoints[a]+waypoints[b];
	    for (c=0; c<waypoints[b]; c++)
	      linestrings[numlinestrings][c] = linestrings[b][c];
	    for (c=0; c<waypoints[a]; c++)
	      linestrings[numlinestrings][c+waypoints[b]] = linestrings[a][c];
	    numlinestrings++;
	    match = 1;
	    continue;
	  }

	  if ((abs(linestrings[a][waypoints[a]-1].xind - linestrings[b][waypoints[b]-1].xind) < 2) && (abs(linestrings[a][waypoints[a]-1].yind - linestrings[b][waypoints[b]-1].yind) < 2)) {
	    //printf("linestrings %d and %d endpoints are near eachother\n", a, b);
	    addressed[a] = 1;
	    addressed[b] = 1;
	    addressed =(short*)realloc(addressed, (numlinestrings+1)*sizeof(short));
	    addressed[numlinestrings] = 0;
	    linestrings = (struct index_pair**)realloc(linestrings, (numlinestrings+1)*sizeof(struct index_pair*));
	    linestrings[numlinestrings] = (struct index_pair *)malloc((waypoints[a]+waypoints[b])*sizeof(struct index_pair));
	    waypoints = (int*)realloc(waypoints, (numlinestrings+1)*sizeof(int));
	    waypoints[numlinestrings] = waypoints[a]+waypoints[b];
	    for (c=0; c<waypoints[b]; c++)
	      linestrings[numlinestrings][c] = linestrings[b][c];
	    for (c=0; c<waypoints[a]; c++)
	      linestrings[numlinestrings][c+waypoints[b]] = linestrings[a][(waypoints[a] - c) - 1];
	    numlinestrings++;
	    match = 1;
	    continue;
	  }
	  if ((abs(linestrings[a][waypoints[a]-1].xind - linestrings[b][0].xind) < 2) && (abs(linestrings[a][waypoints[a]-1].yind - linestrings[b][0].yind) < 2)) {
	    //printf("linestrings %d endpoint and linestring %d startpoint are near eachother\n", a, b);
	    addressed[a] = 1;
	    addressed[b] = 1;
	    addressed =(short*)realloc(addressed, (numlinestrings+1)*sizeof(short));
	    addressed[numlinestrings] = 0;
	    linestrings = (struct index_pair**)realloc(linestrings, (numlinestrings+1)*sizeof(struct index_pair*));
	    linestrings[numlinestrings] = (struct index_pair *)malloc((waypoints[a]+waypoints[b])*sizeof(struct index_pair));
	    waypoints = (int*)realloc(waypoints, (numlinestrings+1)*sizeof(int));
	    waypoints[numlinestrings] = waypoints[a]+waypoints[b];
	    for (c=0; c<waypoints[b]; c++)
	      linestrings[numlinestrings][c] = linestrings[b][(waypoints[b] - c) - 1];
	    for (c=0; c<waypoints[a]; c++)
	      linestrings[numlinestrings][c+waypoints[b]] = linestrings[a][(waypoints[a] - c) - 1];
	    numlinestrings++;
	    match = 1;
	    continue;
	  }
	}
      }
    }
    
    if (json_output == 1) {
      size_t flags = 0;
      flags |= JSON_INDENT(1);
      flags |= JSON_PRESERVE_ORDER;
      flags |= JSON_REAL_PRECISION(6);
      //json_error_t error;
      
      //GNS test
      if (gns_output == 1) {
	jsonbuf = json_dumps( jsonroot, flags );
	char outstr[102800];
	sprintf(outstr, "&&&%d&&&%s", (int)strlen(jsonbuf), jsonbuf);
	size_t strn = strnlen(outstr, 102800);
	printf("strn: %d outstr: %s\n", (int)strn, outstr);
	char tmp[strn];
	int b;
	if (num_valid_contours > 0) {
	  for (b=0; b<num_gns_hosts; b++) {
	    if (sock_connected[b]) {
	      send(gns_sock[b], outstr, sizeof(tmp), 0);
	    }
	  }
	}
      }
      
      json_dump_file(jsonroot, jsonname, flags);
      json_decref( jsonroot );
      json_decref( features_json_arr );
    }

    free(linestrings);
    free(waypoints);
  }
  if (mcc_output == 1) {
    if (connected) 
      close(mcc_sock);
  }
  
  if (gns_output == 1) {
    int b;
    for (b=0; b<num_gns_hosts; b++) {
      if (sock_connected[b]) {
	close(gns_sock[b]);
      }
    }
  }
  
  free(ref_in);
  free(up);
  free(down);
  free(left);
  free(right);
  free(lat_in);
  free(lon_in);
  free(lonaxis);
  free(lataxis);
  //free(ll_array);
  
  return(0);
}
