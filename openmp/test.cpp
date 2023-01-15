#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <omp.h>

void function(int x){
    #pragma omp for schedule(dynamic)
    for(int i=0; i<4; i++){
        printf("%d - %d from thread %d\n", x, i, omp_get_thread_num());
    }
    printf("%d - End function.\n", x);
    
}


int main(int argc, char** argv) {

#pragma omp parallel default(shared) num_threads(8)
{
    
    for(int i=0; i<4; i++){
        function(i);
        printf("Ended iteration %d\n", i);
    }
}
}