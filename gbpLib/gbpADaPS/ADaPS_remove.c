#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <gbpADaPS.h>

void ADaPS_remove(ADaPS      **list, 
                  const char  *name_in,...){
  ADaPS *last   =NULL;
  ADaPS *current=NULL;
  ADaPS *remove =NULL;
  va_list    vargs;
  char       name[ADaPS_NAME_LENGTH];

  // Determine the entry name
  va_start(vargs,name_in);
  vsprintf(name,name_in,vargs);

  // Search the linked list
  current=(*list);
  while(current!=NULL && remove==NULL){
    // If we find an entry with the given name, 
    //   remove it from the list
    if(!strcmp(name,current->name)){
      // If it's the first entry...
      if(last==NULL){
        remove =(*list);
        (*list)=(*list)->next;        
      }
      // ... or else if it is not.
      else{
        remove    =current;
        last->next=current->next;
      }
    }
    last   =current;
    current=last->next;
  }

  // Make sure removed items are deallocated
  if(remove!=NULL)
    ADaPS_deallocate(&remove);

  va_end(vargs);
}
