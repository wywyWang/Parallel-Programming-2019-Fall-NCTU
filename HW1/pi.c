#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <pthread.h>

typedef unsigned long long unlong;
pthread_mutex_t mutex;
unlong sum = 0;
unlong per_toss = 0;

typedef struct thread_info{
	int pid;
}th_info;

void *tossing(void *t_info){
	th_info t = *(th_info *)t_info;

	unlong partial_sum = 0;
	double x,y;
	unsigned int seed = time(NULL);

	for(int i = 0;i < per_toss;i++){
		x = (double) rand_r(&seed) / RAND_MAX;
		y = (double) rand_r(&seed) / RAND_MAX;

		if (x * x + y * y <= 1.0){
		       partial_sum++;
		}
	}

	pthread_mutex_lock(&mutex);
	sum += partial_sum;
	pthread_mutex_unlock(&mutex);

	pthread_exit(NULL);	
}

int main(int argc, char **argv)
{
 
    int ncpu = atoi(argv[1]);
    unlong ntoss = atoi(argv[2]);
    per_toss = ntoss / ncpu;

    pthread_t *th;
    th = (pthread_t *)malloc(ncpu * sizeof(pthread_t));
    pthread_mutex_init(&mutex, NULL);
    th_info t_info[ncpu];

    for(int i = 0;i < ncpu;i++){
	    t_info[i].pid = i;

	    pthread_create(&th[i], NULL, tossing, (void *)&t_info[i]);
    }
    for(int i = 0;i < ncpu; i++){
	    pthread_join(th[i], NULL);
    }

    printf("%f\n", 4.0 * (double)sum / (double)ntoss);
    return 0;
}
