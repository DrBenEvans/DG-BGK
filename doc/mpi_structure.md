#MPI Implementation notes

In the original implementation, no MPI was used. A domain decomposition with
MPI was implemented later to tackle memory issues: the goal was *not*
performance, but only be able to solve a larger problem. 

A master-slave configuration was chosen, where the master would do the task
relative to IO and to the domain decomposition logic as, e.g., splitting the
mesh and set up the lookup tables for all the processes - see [main
README](README.md).

## Velocity space partitioning.

In order to do velocity space partitioning, the MPI communicator layout is the
one depicted in this figure: ![communicator_layout.](communicators_layout.svg)
(At the moment on Gitlab the svg is not shown in the text but you can
visualise it separately).

Some remarks:
* Velocity space has been partitioned **(Performance problems coming from this can be easily
  fixed.)**

List of relevant variables:
* **VSPACE_FIRST**: First index of vspace that is in the responsibility of the
  current process.
* **VSPACE_LAST**: Last index of vspace that is in the responsibility of the
  current process.





