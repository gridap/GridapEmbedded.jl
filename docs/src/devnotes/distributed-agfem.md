
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

On dof ownership:

- Original DoFs are split into free and Dirichlet DoFs (positive/negative), with a single all-positive numbering where Dirichlet DoFs come after free DoFs.
- Owned sDOFs are:
  - always local
  - always constrained by owned mDOFs, which are themselves always local.
- Ghost sDOFs are:
  - always local
  - but can be constrained by ghost mDOFs, which can be non-local.
- This means that we might have non-local mDOFs that we need to account for. These do not appear during integration, but have contributions coming from the slave sDOFs they constrain. These contributions have to be communicated to the owner of the mDOF during assembly.
- Moreover, these non-local mDOFs can be free or Dirichlet, which have to be treated differently: While free mDOFs appear in the final FESpace PRange, Dirichlet mDOFs do not. However, we still need to create a PRange for them, so that we can communicate Dirichlet values during interpolation (since they can be non-local, Dirichlet values cannot always be computed locally).

How to get the different components we need to build the distributed space:

- Given the original local spaces, the cell ownership and the aggregates, we can colour each local DOF as sDOF, free mDOF (mfdof) or Dirichlet mDOF (mddof).
  - Generate the sDOF_to_dof mapping locally, as well as the sDOF gids.
  - Generate initial mfdof_to_dof and mddof_to_dof mapping locally, as well as the mfdof and mddof gids. All of these will need to be extended later with non-local mDOFs.
- With the sDOF gids and the local cell_to_dof mapping, we can locally count the number of mDOFs constraining owned sDOFs, then communicate the counts.
- The above allows us to allocate the local sDOF_to_mDOFs table. Like before, we can locally fill in the owned rows (sDOFs), map to mDOF global ids, then communicate the non-local rows.
- This generates a consistent sDOF_to_global_mDOF mapping, that we need to renumber into local mDOF ids. This is the tricky part:
  - For each global mDOF id in the table, it could be either local (belongs to the local-to-global mapping of the previously generated mDOF gids) or non-local.
  - If it is local, the communicated global id will be listed within the local global-to-local map. We can therefore directly map it to a local id.
  - If it is non-local, we extend the local global-to-local map (the gid is known, the lid is appended).
  - When this is done, we can create the new PRanges. This does not require further communication, since we have not added any new global id, only added pre-existing ones.
  - Like explained above, we have to do this while differentiating between free and Dirichlet mDOFs.

This allows us to create the local FESpaceWithLinearConstraints spaces, and together with the mfdof PRange this creates our typical DistributedFESpace. On top of this, we also need to keep both the sDOF PRange and the mddof PRange (for interpolation). These two will be stored within a new metadata type, which will also allow us to dispatch when interpolating.

## Optimizing communications during assembly

The assembly routines are general enough that we can use the original cell partition to integrate. However this is not optimal, since contributions on ghost cut cells would always have to be communicated to the processor that owns their root cell. This means that the size of communications grows with the patch size.

There is however a way to obtain optimal communications: Using the patch-conforming cell partition we created before. With this cell partition, the only contributions that have to be communicated are the ones corresponding to DoFs at the interface of two patches. This means that the size of the communication is independent of the patch size.
