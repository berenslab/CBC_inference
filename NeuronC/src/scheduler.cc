#include "scheduler.h"

#include <vector>
#include <cmath>

extern "C" {
#include <stdio.h>
}

extern double timinc;

scheduler::scheduler()
{
  elapsedTimesteps = 0;
  idcount = 0;
}

void scheduler::init()
{
  //fprintf(stderr, "# [scheduler] initialized: timinc=%g, \n", timinc);
  reset();
}

void scheduler::reset()
{
  elapsedTimesteps = 0;
}

void scheduler::step()
{
  int tmod, delflag;  

  //fprintf(stderr, "# step(): elapsedTimesteps=%d\n", elapsedTimesteps);
  for (int k = 0; k < tasks.size(); k++) {    
    tmod = elapsedTimesteps % tasks[k]->modint;
    //fprintf(stderr, " k=%d, modint=%d, tmod=%d, elapsedTimesteps=%d\n",
    //        k, tasks[k]->modint, tmod, elapsedTimesteps);

    delflag = 0;
    if ((tasks[k]->maxsteps != -1) && (elapsedTimesteps >= tasks[k]->maxsteps)) delflag = 1;

    if (tmod == 0) {
      //fprintf(stderr, "  running task %d\n", k);
      tasks[k]->execute(); 
    }
  
    if (delflag) remqueue.push_back(tasks[k]->id);
  }
  elapsedTimesteps++;

  if (remqueue.size() > 0) {
    for (int k = 0; k < remqueue.size(); k++) {
      fprintf(stderr, "# [scheduler] removing task %d\n", remqueue[k]);
      removeTask(remqueue[k]);
    }
    remqueue.clear();
  }
}


int scheduler::addTask(double interval, void (*taskfunc)(void))
{
  return addTask(interval, taskfunc, -1);
}

int scheduler::addTask(double interval, void (*taskfunc)(void), double duration)
{
  if (interval < timinc) return -1;

  ctask *t = new ctask(taskfunc);
  t->id = idcount++;
  t->interval = interval;
  t->modint = int(interval / timinc);
  if (duration != -1)   t->maxsteps = int(duration / timinc);
  else                  t->maxsteps = -1;

  tasks.push_back(t);  

  //fprintf(stderr, "# [scheduler] added task, id=%d, interval=%g, maxsteps=%d\n", t->id, t->interval, t->maxsteps);

  return t->id;
}

int scheduler::addTask(double interval, task *t)
{
  return addTask(interval, t, -1);
}

int scheduler::addTask(double interval, task *t, double duration)
{
  if (interval < timinc) return -1;

  t->id = idcount++;
  t->interval = interval;
  t->modint = int(interval / timinc);
  if (duration != -1)  t->maxsteps = int(duration / timinc);
  else                 t->maxsteps = -1;
  
  tasks.push_back(t);

  fprintf(stderr, "# [scheduler] added task, id=%d, interval=%g\n", (int)t->id, t->interval);

  return t->id;
}



void scheduler::removeTask(int id)
{
  vector<task *>::iterator iter = tasks.begin();
  do {
    if ((*iter)->id == id) {
      tasks.erase(iter);      
      fprintf(stderr, "# [scheduler] removed task %d\n", id);
      break;
    }
    iter++;
  }
  while (iter != tasks.end());
}

