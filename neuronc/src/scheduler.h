
#ifndef SCHEDULER_H
#define SCHEDULER_H

#include <vector>
using namespace std;

/* A class that represents a task to be
   run at a specified time interval.
*/
class task
{
 public:
  long id;          //id number of task
  double interval;  //interval to execute task
  int modint;       //the interval modulus timinc
  int maxsteps;     //maximum number of steps to keep task around

  /* Dummy function to be overridden by a subclass */
  inline virtual void execute() {;}
};


/* Subclass of task which represents a global function
   to be executed at a specified time interva.
*/
class ctask : public task
{
 protected:
  void (*cfunc)(void);  //pointer to task function
 public:
  inline ctask(void (*func)(void)) { cfunc = func; }
  inline void execute() { cfunc(); }
};


/* Class that contains a list of tasks to execute at specified
   time intervals.
*/
class scheduler {

 protected:

  long elapsedTimesteps;  //time steps since start
  vector<task *> tasks;   //list of tasks
  long idcount;           //next id to assign
  vector<int> remqueue;    //queue of id's to remove

 public:

  scheduler();

  /* Initialize scheduler */
  void init();
  /* Reset scheduler */
  void reset();
  /* Run scheduler for one step (meant to be run every timinc) */
  void step();
  /* Remove a task by its id */
  void removeTask(int id);
  /* Add a function to be executed as a task */
  int  addTask(double interval, void (*taskfunc)(void));
  /* Add a function to be executed as a task, with a specified duration */
  int addTask(double interval, void (*taskfunc)(void), double duration);
  /* Add a task to be executed */
  int addTask(double interval, task *t);
  /* Add a task to be executed, with a specified duration */
  int addTask(double interval, task *t, double duration);
};

#endif
