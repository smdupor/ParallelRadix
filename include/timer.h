#ifndef TIMER_H
#define TIMER_H

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/types.h>

enum POS {EXPERIMENT=0, INIT=1, COUNT=2, XFORM=3, MPI_MSG=4, SORT=5, CRUSH=6 };

struct Timer{
  struct timeval start_time;
  struct timeval elapsed_time;
 };

  void Start(struct Timer* timer);
  void Stop(struct Timer* timer);
  double Seconds(struct Timer* timer);
  double Millisecs(struct Timer* timer);
  double Microsecs(struct Timer* timer);


#endif  // TIMER_H