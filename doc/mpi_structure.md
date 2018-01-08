
#MPI Implementation notes

In the original implementation, no MPI was used. A domain decomposition with
MPI was implemented later to tackle memory issues: the goal was *not*
performance, but only be able to solve a larger problem. 

A master-slave configuration was chosen, where the master would do the task
relative to IO and to the domain decomposition logic as, e.g., splitting the
mesh and set up the lookup tables for all the processes - see [main
README](README.md).

### Some comments on the old MPI implementation

In the original MPI implementation, in most cases a very inefficient
communication pattern was choosen, like the one descripted in the pseudocode
snippets below.

Here is one pattern used in initialization:
```c
for all points in the mesh 
    if on the master rank:
        find the process IG the point belongs to
    
    MPI_BCAST of IG to all processes

    if on the master process:      MPI_SEND to the IG-th process
    if on the IG-th process:       MPI_RECV from the master process

end for on all points in the mesh 
```

Another variant, for inter process comunication of border data:

```c
for IG in all processes 
    for all points in the partition IG
        if on the IG-th process:
            find the process IG2 to which data must be sent

        MPI_BCAST of IG2 to all other processes

    if on the IG-th process:     MPI_SEND to IG2-th process
    if on the IG2-th process:    MPI_RECV from the IG-th process


    end for on all points in the partition
end for on all processes
```

This is computationally inefficient, since for each step of the time marching 
scheme millions of MPI_BCAST calls are performed.



## Velocity space partitioning.

In order to do velocity space partitioning, the MPI communicator layout is the
one depicted in this figure: ![communicator_layout.](communicators_layout.svg)
(At the moment on Gitlab the svg is not shown in the text but you can
visualise it separately).

Some remarks:
* Velocity space has been partitioned **(Performance problems coming from this 
can be easily fixed.)**

List of relevant variables:
* **VSPACE_FIRST**: First index of vspace that is in the responsibility of the
  current process.
* **VSPACE_LAST**: Last index of vspace that is in the responsibility of the
  current process.

* **VCORD**


### Notes

* The **INICON** and **INICON2** subroutines have been modified for VSPACE
  partitioning, but the highly inefficient communication pattern was kept in
place, since this is not a performance-critical part of the program.


