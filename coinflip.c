//------------------------------------
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <ctime>
using namespace std;

int flips();

int main()
{	int result = 0;
	int i;
	//srand(time(NULL));
	
	for (i = 1; i < 999999; i++)
		result += flips();

	cout << "Expected value = "<< result / i << endl;
        return(0);
}
	

int flips()
{	
        int i, j, counter;
	//i = 0;
	srand(0);
	i = rand() % 2;
        srand(100);
        j = rand() % 2;
	counter = 2;

	//while ((i+j) != 2)
        while((i+j)%2 != 0)
	{
      	        i = j;
      	        //i = rand() % 2;
      	        srand(300);
		j = rand() % 2;
		counter++;
	}
	return counter;
}
//---------------------------
