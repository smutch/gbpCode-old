#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees_build.h>
#include <gbpTrees_analysis.h>

void move_treenode_to_snap(tree_node_info **halo,int snap_new){
   int snap_i=(*halo)->snap_tree;
   if(snap_i<snap_new){
     while((*halo)!=NULL && snap_i<snap_new){
        (*halo)=(*halo)->descendant;
        if((*halo)!=NULL)
           snap_i=(*halo)->snap_tree;
     }
   }
   else if(snap_i>snap_new){
     while((*halo)!=NULL && snap_i>snap_new){
        (*halo)=(*halo)->progenitor_first;
        if((*halo)!=NULL)
           snap_i=(*halo)->snap_tree;
     }
   }

   // Return NULL if the halo doesn't exist at the new snap
   if((*halo)!=NULL){
      if((*halo)->snap_tree!=snap_new)
        (*halo)=NULL;
   }
}

