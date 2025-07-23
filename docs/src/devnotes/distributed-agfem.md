
# Distributed cell aggregation

Given a distributed mesh, we want to associate to each CUT cell an UNCUT cell that we will call the root cell.

We consider the set of paths that connect a CUT cell to all its root candidates. We associate a weight to each path. The simplest weight function would its length (i.e the distance between the cell and its root candidate). We define the optimal root cell as the one associated with the path that has the minimum weight.

Assumptions:

- The root cell is unique for each CUT cell.
- The optimal path is fully contained in one of the subdomains sharing the CUT cell.

The second assumption is crucial, and can always by satisfyed by increasing the number of ghost layers in each subdomain.

If we have the above assumptions, the following holds true:

- Applying the serial algorithm to each subdomain will yield, for each CUT cell, a set of root candidates.
- The optimal root cell is among the root candidates, since its path is fully contained in one of the subdomains.
- We can then communicate the root candidates (along with their weights) to the processor owning the CUT cell, select the optimal root and scatter the information back to the neighbors.
- At the end of the process, each CUT cell will have a unique root cell associated to it. This root cell is the optimal, and partition invariant.

If the second assumption is violated, the algorithm will yield a sub-optimal root cell. This will also depend on the partition.

Another thing that would be nice is to have an indicator that tells the user when the second assumption might have been violated. One possibility, that also makes the next two sections possible, is the following: We require that any owned aggregate has at least a layer of local cells around it. I.e no aggregate touches an interface between two subdomains. This the equivalent condition we normally fullfill with our cell partitions (i.e the fact that we require a layer of ghost cells), but for the aggregates.
If the above is fullfilled, it means that no aggregate wanted to expand any further and therefore they are all optimal. Another way of saying this is that all optimal path sizes are below the maximal path size allowed by the ghost layers, so paths are indeed optimal.

## Distributed constrained FESpaces

Once we have our patches, the question is how we distribute the ownership of the DoFs of the constrained FESpaces.

Some preliminary observations:

- The ConstrainedFESpace eliminates the slave DoFs before assembly. That is during assembly all the cell contributions on a patch get assembled to the master DoFs, i.e the DoFs on the root cell. This means that all slave DoFs of a patch de-facto belong to the same processor as the master DoFs they are tied to.
- We can easily create a patch-conforming cell-partition, where cells are owned by the owner of the patch they belong to. The owner of the patch is taken as the owner of the root cell in the original cell-partition. This new partition does NOT require any communication. It is just a re-partitioning of the local cells within a processor.
- Something I don't know, but should be possible and even required for this to work is that local DoFs can be constrained by DoFs that do not belong to the original space. TODO: Check this.

Given the above, this is how we can create the DoF constraints:

- We also allocate the necessary JaggedArrays to hold the SDoF-to-data info for ALL local slave DoFs.
- Each processor creates constraints for it's owned patches. This is easily done by running the local routines on a Triangulation (as we talked about last week). The Triangulation is given by the owned cells according to the patch-conforming partition described above. If we do the serial code in a general enough way, we can reuse it here.
- We do a single step of communication.

## Optimizing communications during assembly

The assembly routines are general enough that we can use the original cell partition to integrate. However this is not optimal, since contributions on ghost cut cells would always have to be communicated to the processor that owns their root cell. This means that the size of communications grows with the patch size.

There is however a way to obtain optimal communications: Using the patch-conforming cell partition we created before. With this cell partition, the only contributions that have to be communicated are the ones corresponding to DoFs at the interface of two patches. This means that the size of the communication is independent of the patch size.
