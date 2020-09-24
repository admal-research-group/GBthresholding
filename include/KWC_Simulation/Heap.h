#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include <stdio.h>
#include <assert.h>
#include <iostream>
#include <vector>
#define k_type double

#ifndef BRANCHING_NUMBER
#define BRANCHING_NUMBER 8
#endif

#ifndef maj_mbo_heap_h
#define maj_mbo_heap_h

typedef struct{
    k_type key;
    int originalIndex;
}i_heap_node;

typedef struct{
    i_heap_node *root;
    int *locations;
    int count;
}i_heap;



static void i_heap_swap_nodes(i_heap *heap,int nodeIndex1, int nodeIndex2){
    int i1=heap->root[nodeIndex1].originalIndex;
    int i2=heap->root[nodeIndex2].originalIndex;
    k_type key1=heap->root[nodeIndex1].key;
    k_type key2=heap->root[nodeIndex2].key;
    heap->root[nodeIndex1].originalIndex=i2;
    heap->root[nodeIndex2].originalIndex=i1;
    heap->root[nodeIndex1].key=key2;
    heap->root[nodeIndex2].key=key1;
    heap->locations[i1]=nodeIndex2;
    heap->locations[i2]=nodeIndex1;
    
}

static void i_heap_bubble_up(i_heap *heap, int nodeIndex){
    
    while (nodeIndex>0) {
        k_type myKey=heap->root[nodeIndex].key;
        int parentIndex=nodeIndex/BRANCHING_NUMBER-((nodeIndex%BRANCHING_NUMBER)==0);
        k_type parentKey=heap->root[parentIndex].key;
        if(myKey<parentKey){
            i_heap_swap_nodes(heap,nodeIndex,parentIndex);
            nodeIndex=parentIndex;
        }else{
            break;
        }
    }
}

static void i_heap_push_down(i_heap *heap, int nodeIndex){
    
    int i;
    int count=heap->count;
    int childIndex=nodeIndex*BRANCHING_NUMBER+1;
    while(childIndex<count){
        int minIndex=childIndex;
        k_type myKey=heap->root[nodeIndex].key;
        k_type min=myKey;
        for(i=0;i<BRANCHING_NUMBER;i++){
            if(childIndex+i<count){
                k_type childKey=heap->root[childIndex+i].key;
                if(childKey<min){
                    min=childKey;
                    minIndex=childIndex+i;
                }
            }
        }
        if(min<myKey){
            i_heap_swap_nodes(heap,nodeIndex,minIndex);
            nodeIndex=minIndex;
            childIndex=nodeIndex*BRANCHING_NUMBER+1;
        }else{
            break;
        }
        
    }
}

static void i_heap_insert_node_with_location(i_heap *heap, int originalIndex, k_type key, int location){
    
    int count=heap->count;
    heap->root[count].originalIndex=originalIndex;
    heap->root[count].key=key;
    heap->locations[location]=count;
    heap->count++;
    i_heap_bubble_up(heap, count);
    
}

static void i_heap_add_node_to_bottom_with_location(i_heap *heap, int originalIndex, k_type key, int location){
    
    int count=heap->count;
    heap->root[count].originalIndex=originalIndex;
    heap->root[count].key=key;
    heap->locations[location]=count;
    heap->count++;
}



static void i_heap_delete_min(i_heap *heap){
    int count=heap->count;
    int index=heap->root[0].originalIndex;
    heap->locations[index]=-1;
    i_heap_swap_nodes(heap,count-1,0);
    heap->count--;
    i_heap_push_down(heap,0);
}


static void i_heap_decrease_key(i_heap *heap, int nodeIndex, k_type newKey){
   
    heap->root[nodeIndex].key=newKey;
    i_heap_bubble_up(heap,nodeIndex);
}

static i_heap i_heap_create_empty_heap(int ncount, int locationCount){
    i_heap heap;
    heap.count=0;
    
    heap.root=(i_heap_node *) calloc(ncount,sizeof(i_heap_node));
    heap.locations=(int *) calloc(locationCount,sizeof(int));
    
    memset(heap.locations,-1,locationCount*sizeof(int));
    return heap;
}


static int i_heap_empty(i_heap *heap){
    return heap->count==0;
}

static void i_heap_destroy_heap(i_heap *heap){
   
    free(heap->root);
    free(heap->locations);
}

static void i_heap_clear_heap(i_heap *heap){
    heap->count=0;
}


#endif






