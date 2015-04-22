#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees_build.h>
#include <gbpTrees_analysis.h>
#include <assert.h>

void add_to_treenode_hist(tree_info *trees,treenode_hist_info *hist,tree_node_info *current_halo){
   int i_x=-1;
   int i_y=-1;
   int n_x=-1;
   int n_y=-1;
   for(int i_axis=0;i_axis<2;i_axis++){
      int    *i_d;
      switch(i_axis){
         case 0:
           i_d=&i_x;
           n_x=hist->n_x;
           break;
         case 1:
           i_d=&i_y;
           n_y=hist->n_y;
           break;
      }
      (*i_d)=set_treenode_hist_index(trees,hist,current_halo,i_axis);
   }

   // Populate the array 
   if(i_x>=0 && i_x<n_x && i_y>=0 && i_y<n_y)
      hist->array[i_y*n_x+i_x]++;
}

